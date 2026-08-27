#!/usr/bin/env python3
"""
Summarise AutoReject channel interpolation from the preprocessing reject log.

IMPORTANT: this counts interpolation performed by *AutoReject* (per-epoch,
status == 'interpolated'). It is NOT the manual bad-channel nulling recorded in
manual_channel_qc_nd26.csv — channels a rater flagged as bad are excluded
*before* AutoReject and never appear in the reject log at all.

Reads : preprocessing_output/p3b_reject_log_autoreject_nd_26_allch.csv
        (long: subject x epoch x channel, status in {good, interpolated, bad})
Writes: preprocessing_output/p3b_autoreject_interpolation_by_channel_nd26.csv
        preprocessing_output/p3b_autoreject_interpolation_by_subject_channel_nd26.csv

Paths resolve relative to this script, so it runs from any working directory.
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "preprocessing_output")
REJECT_LOG = os.path.join(OUT_DIR, "p3b_reject_log_autoreject_nd_26_allch.csv")

CH_OUT   = os.path.join(OUT_DIR, "p3b_autoreject_interpolation_by_channel_nd26.csv")
SUBCH_OUT = os.path.join(OUT_DIR, "p3b_autoreject_interpolation_by_subject_channel_nd26.csv")


def status_counts(df, group_cols):
    """Wide good/interpolated/bad counts + totals + percentages per group."""
    g = df.groupby(group_cols)["status"].value_counts().unstack(fill_value=0)
    for s in ["good", "interpolated", "bad"]:
        if s not in g.columns:
            g[s] = 0
    g["n_epochs"] = g[["good", "interpolated", "bad"]].sum(axis=1)
    g["pct_interpolated"] = 100 * g["interpolated"] / g["n_epochs"]
    g["pct_bad"] = 100 * g["bad"] / g["n_epochs"]
    g = g.rename(columns={
        "good": "n_good",
        "interpolated": "n_interpolated_autoreject",
        "bad": "n_bad_rejected",
    })
    return g[["n_interpolated_autoreject", "pct_interpolated",
              "n_bad_rejected", "pct_bad", "n_good", "n_epochs"]].reset_index()


def main():
    df = pd.read_csv(REJECT_LOG)
    n_sub = df["subject"].nunique()
    print(f"reject log: {len(df):,} rows | {n_sub} subjects | status={sorted(df['status'].unique())}")

    # 1) per-channel (pooled across subjects), ranked by AutoReject interpolation
    by_ch = status_counts(df, ["channel"]).sort_values(
        "n_interpolated_autoreject", ascending=False)
    by_ch.round(2).to_csv(CH_OUT, index=False)

    # 2) per-subject x channel
    by_subch = status_counts(df, ["subject", "channel"]).sort_values(
        ["subject", "n_interpolated_autoreject"], ascending=[True, False])
    by_subch.round(2).to_csv(SUBCH_OUT, index=False)

    print("\n=== AutoReject interpolation by channel (pooled) ===")
    print(by_ch[["channel", "n_interpolated_autoreject", "pct_interpolated",
                 "n_bad_rejected", "pct_bad"]].round(2).to_string(index=False))
    print(f"\nwrote:\n  {CH_OUT}\n  {SUBCH_OUT}")


if __name__ == "__main__":
    main()
