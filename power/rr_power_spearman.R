# =============================================================================
# Stage 1 RR — sample size for the confirmatory test
#
#   Outcome   : P3b difference wave (target - non-target, mean of O1/O2, 300-600 ms)
#   Predictor : sm_composite_hrs (time-weighted daily social-media use)
#   Test      : one-tailed Spearman rank correlation, H0: rho = 0, H1: rho < 0
#
# Method: Fisher z-based sample-size formula for a Spearman coefficient,
#   May & Looney (2020), Journal of Biometrics & Biostatistics 11(2):440,
#   eq. (7) and Table 1.  https://doi.org/10.37421/jbmbs.2020.11.440
#
#       n = b + c^2 * [ (z_alpha + z_power) / (z(rho1) - z(rho0)) ]^2
#
#   Spearman: b = 3 and c^2 = 1 + rho_s0^2 / 2, where rho_s0 is the NULL value
#   (1.06 if |rho_s0| >= 0.95). Here the null is 0, so c^2 = 1.
#   One-tailed: z_alpha replaces z_alpha/2. Non-integer results are rounded up.
#
# Planning effect size: safeguard-power approach (Perugini, Gallucci &
#   Costantini, 2014) — the lower bound of a 60% CI around the pilot estimate.
#   The interval uses the variance appropriate for a Spearman coefficient,
#   var(z) = (1 + rho^2/2)/(n-3)  (Bonett & Wright, 2000), not the Pearson
#   1/(n-3). Note this is the ALTERNATIVE rho, unlike c^2 in the sample-size
#   formula above, which takes the NULL value.
#
# Usage: Rscript rr_power_spearman.R
# =============================================================================

# ---- Inputs -----------------------------------------------------------------
PILOT_RHO <- 0.3264   # |rho| observed in the pilot (ND_p3b_composite_O1O2avg.R)
PILOT_N   <- 86       # pilot analysable N
CI_LEVEL  <- 0.60     # CI level for the safeguard anchor
ALPHA     <- 0.05     # one-tailed
POWER     <- 0.95     # >= 0.95 for frequentist plans
RHO_NULL  <- 0        # null value of the Spearman coefficient
RETENTION <- 86/120   # pilot data-retention rate

# ---- Safeguard anchor: lower bound of a CI_LEVEL CI around the pilot rho ----
safeguard_anchor <- function(rho, n, level = CI_LEVEL) {
  z  <- atanh(rho)
  se <- sqrt((1 + rho^2 / 2) / (n - 3))   # Bonett & Wright (2000) Spearman variance
  zc <- qnorm(1 - (1 - level) / 2)
  c(lower = tanh(z - zc * se), upper = tanh(z + zc * se))
}

# ---- May & Looney (2020), eq. (7) + Table 1 ---------------------------------
n_spearman <- function(rho1, rho0 = RHO_NULL, alpha = ALPHA, power = POWER,
                       one_tailed = TRUE) {
  za <- qnorm(1 - if (one_tailed) alpha else alpha / 2)
  zb <- qnorm(power)
  b  <- 3
  c2 <- if (abs(rho0) < 0.95) 1 + rho0^2 / 2 else 1.06
  ceiling(b + c2 * ((za + zb) / (atanh(rho1) - atanh(rho0)))^2)
}

# ---- Report -----------------------------------------------------------------
ci     <- safeguard_anchor(PILOT_RHO, PILOT_N)
anchor <- round(ci["lower"], 2)
n_req  <- n_spearman(anchor)

cat(sprintf("Pilot estimate        : rho = %.3f (N = %d)\n", PILOT_RHO, PILOT_N))
cat(sprintf("%.0f%% CI                : [%.3f, %.3f]\n", 100*CI_LEVEL, ci["lower"], ci["upper"]))
cat(sprintf("Safeguard anchor      : rho = %.2f\n\n", anchor))
cat(sprintf("One-tailed alpha=%.2f, power=%.2f, null rho=%.2f\n", ALPHA, POWER, RHO_NULL))
cat(sprintf("Required analysable N : %d\n", n_req))
cat(sprintf("Recruit (retention %.1f%%): %d\n\n", 100*RETENTION, ceiling(n_req / RETENTION)))

cat("Sensitivity to the anchor:\n")
cat(sprintf("%8s %14s %10s\n", "rho", "analysable N", "recruit"))
for (r in c(0.20, 0.22, 0.24, 0.26, 0.30)) {
  n <- n_spearman(r)
  cat(sprintf("%8.2f %14d %10d\n", r, n, ceiling(n / RETENTION)))
}
