#!/usr/bin/env python3
"""
Two QC checks on the P3b events.tsv (EDF-clock onsets):
  1. Reaction-time histogram, target vs non-target (correct trials only).
  2. Timing-schematic check: onset-to-onset SOA distribution vs the design
     spec (200 ms letter + 1200-1400 ms fixation ITI = 1400-1600 ms/trial).
Pools all subjects. Non-destructive: reads events.tsv, writes figs+csv to OUT.
"""
import glob, os, ast
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Paths are resolved relative to THIS script's location, so it runs from any cwd.
HERE = os.path.dirname(os.path.abspath(__file__))
BEH = os.path.join(HERE, "bids_data_behavior")
OUT = os.path.join(HERE, "timing_rt_checks")
os.makedirs(OUT, exist_ok=True)

def to_scalar_rt(x):
    """response_rt may be a scalar or a stringified list (invalid multi-press)."""
    if pd.isna(x):
        return np.nan
    s = str(x).strip()
    if s.startswith("["):
        try:
            v = ast.literal_eval(s)
            return float(v[0]) if len(v) else np.nan   # first press
        except Exception:
            return np.nan
    try:
        return float(s)
    except Exception:
        return np.nan

files = sorted(glob.glob(os.path.join(BEH, "*_task-p3b_events.tsv")))
print(f"events.tsv files: {len(files)}")

rt_rows = []      # per-trial RT (correct stim trials)
soa_rows = []     # per-consecutive-stimulus SOA
n_sub_ok = 0

for f in files:
    sub = os.path.basename(f).split("_")[0].replace("sub-", "")
    df = pd.read_csv(f, sep="\t")
    if "trial_type" not in df or "onset" not in df:
        continue
    stim = df[df["trial_type"].isin(["target", "non-target"])].copy()
    if stim.empty:
        continue
    stim = stim.sort_values("onset").reset_index(drop=True)
    n_sub_ok += 1

    # --- RT (correct only, scalar) --- (skip EDF-only files w/o response cols)
    if "response_rt" in stim and "response_correct" in stim:
        stim["rt_s"] = stim["response_rt"].map(to_scalar_rt)
        corr = stim[stim["response_correct"].astype(str).isin(["True", "TRUE", "true"])]
        for _, r in corr.iterrows():
            if np.isfinite(r["rt_s"]):
                rt_rows.append({"sub": sub, "trial_type": r["trial_type"],
                                "rt_ms": r["rt_s"] * 1000.0})

    # --- SOA: onset diff between consecutive stimuli ---
    o = stim["onset"].to_numpy(dtype=float)
    d = np.diff(o) * 1000.0   # ms
    for gap in d:
        soa_rows.append({"sub": sub, "soa_ms": gap})

rt = pd.DataFrame(rt_rows)
soa = pd.DataFrame(soa_rows)
print(f"subjects with stim rows: {n_sub_ok}")
print(f"RT trials (correct): {len(rt)}  |  SOA intervals: {len(soa)}")

# ---------------- 1. RT histogram target vs non-target ----------------
def desc(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return dict(n=len(x), mean=np.mean(x), sd=np.std(x, ddof=1),
                median=np.median(x), p5=np.percentile(x, 5),
                p95=np.percentile(x, 95), min=np.min(x), max=np.max(x))

rt_stats = []
fig, ax = plt.subplots(figsize=(8, 5))
bins = np.arange(0, 1600, 25)
for tt, col in [("target", "#1f77b4"), ("non-target", "#d62728")]:
    v = rt.loc[rt["trial_type"] == tt, "rt_ms"].to_numpy()
    ax.hist(v, bins=bins, alpha=0.55, label=f"{tt} (n={len(v)})",
            color=col, density=True)
    s = desc(v); s["trial_type"] = tt; rt_stats.append(s)
    ax.axvline(s["median"], color=col, ls="--", lw=1)
ax.set_xlabel("Reaction time (ms)"); ax.set_ylabel("Density")
ax.set_title("P3b reaction time — correct trials (pooled subjects)")
ax.legend()
fig.tight_layout(); fig.savefig(os.path.join(OUT, "rt_hist_target_vs_nontarget.png"), dpi=140)
plt.close(fig)

rt_stats_df = pd.DataFrame(rt_stats)[["trial_type","n","mean","sd","median","p5","p95","min","max"]]
rt_stats_df.to_csv(os.path.join(OUT, "rt_summary_target_vs_nontarget.csv"), index=False)
print("\n=== RT (ms), correct trials ===")
print(rt_stats_df.round(1).to_string(index=False))

# ---------------- 2. Timing schematic: SOA distribution ----------------
# Design: 200 ms stim + 1200-1400 ms fixation -> SOA 1400-1600 ms within block.
# Big gaps = between-block breaks; separate them out.
soa_v = soa["soa_ms"].to_numpy()
BLOCK_BREAK = 3000.0   # ms; anything above this is a block boundary, not a trial SOA
within = soa_v[soa_v < BLOCK_BREAK]
breaks = soa_v[soa_v >= BLOCK_BREAK]

sw = desc(within)
spec_lo, spec_hi = 1400.0, 1600.0
in_spec = np.mean((within >= spec_lo) & (within <= spec_hi)) * 100.0

fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(within, bins=np.arange(1300, 1750, 10), color="#2ca02c", alpha=0.8)
ax.axvspan(spec_lo, spec_hi, color="orange", alpha=0.20,
           label=f"design 1400-1600 ms ({in_spec:.0f}% inside)")
ax.axvline(sw["median"], color="k", ls="--", lw=1.2,
           label=f"median {sw['median']:.0f} ms")
ax.set_xlabel("Stimulus onset-to-onset interval / SOA (ms)")
ax.set_ylabel("Count")
ax.set_title("P3b trial timing vs design schematic (within-block SOA)")
ax.legend()
fig.tight_layout(); fig.savefig(os.path.join(OUT, "soa_within_block_vs_design.png"), dpi=140)
plt.close(fig)

soa_out = pd.DataFrame([
    {"metric": "within_block_SOA", **{k: sw[k] for k in ["n","mean","sd","median","p5","p95","min","max"]}},
])
soa_out["pct_in_design_1400_1600"] = [round(in_spec, 1)]
soa_out.to_csv(os.path.join(OUT, "soa_timing_summary.csv"), index=False)

print("\n=== Within-block SOA (ms) vs design 1400-1600 ===")
print(pd.DataFrame([sw]).round(1).to_string(index=False))
print(f"% of within-block SOAs inside design window: {in_spec:.1f}%")
print(f"block-break gaps (>3s): n={len(breaks)}, median={np.median(breaks)/1000:.1f}s"
      if len(breaks) else "no block-break gaps")
print(f"\nfigures + csv -> {OUT}")
