# =============================================================================
# Stage 1 RR — sample size for the confirmatory test
#
#   Outcome   : P3b difference wave (target - non-target, mean of O1/O2, 300-600 ms)
#   Predictor : sm_composite_hrs (time-weighted daily social-media use)

# Assumptions:
# - Independent participants.
# - Social-media score frequencies follow the N = 86 pilot.
# - P3b follows a normal distribution with the pilot mean and SD.

# 1. SETTINGS
set.seed(123)
N_values <- seq(100, 400, 25)
nsim <- 10000
alpha <- 0.05
bound <- 0.20
target_N <- 350
p3b_mean <- 3.066516
p3b_sd <- 4.510596

# 2. SOCIAL-MEDIA SCORE DISTRIBUTION
# Frequencies of the 25 ordered composite scores in the exploratory sample.
counts <- c(1,1,4,2,1,1,8,4,3,2,7,9,3,3,3,6,1,1,5,4,1,1,7,5,3)
prob <- counts / sum(counts)
cumprob <- c(0, cumsum(prob))
cumprob[length(cumprob)] <- 1
cuts <- qnorm(cumprob)
midranks <- head(cumprob, -1) + prob / 2

# 3. CALIBRATE THE POPULATION SPEARMAN CORRELATION WITH TIES
population_rho <- function(a) {
  integrals <- sapply(seq_along(prob), function(j)
    integrate(function(z) pnorm(a*z/sqrt(2-a^2))*dnorm(z),
              cuts[j], cuts[j+1], rel.tol=1e-8)$value)
  (sum(midranks*integrals)-0.25) /
    sqrt(sum(prob*(midranks-0.5)^2)/12)
}
latent_r <- uniroot(function(a) population_rho(a)-bound,
                    c(0, 0.95), tol=1e-8)$root

# 4. SIMULATE ONE STUDY
simulate <- function(n, scenario) {
  a <- if (scenario == "Directional") -latent_r else 0
  z <- rnorm(n)
  x <- findInterval(z, cuts)
  y <- p3b_mean + p3b_sd*(a*z + sqrt(1-a^2)*rnorm(n))
  
  if (scenario == "Directional") {
    # One-sided Spearman test for a negative association.
    passed <- cor.test(x, y, method="spearman",
                       alternative="less", exact=FALSE)$p.value < alpha
  } else {
    # Equivalence: 90% Fisher-z CI with Bonett-Wright standard error.
    r <- cor(x, y, method="spearman")
    r <- max(-1+1e-12, min(1-1e-12, r))
    se <- sqrt((1+r^2/2)/(n-3))
    ci <- tanh(atanh(r) + c(-1,1)*qnorm(1-alpha)*se)
    passed <- ci[1] > -bound && ci[2] < bound
  }
  
  # Both scenarios also require a positive group-level P3b difference.
  passed && t.test(y, mu=0, alternative="greater")$p.value < alpha
}

# 5. ESTIMATE POWER AND 95% MONTE CARLO UNCERTAINTY INTERVALS
results <- list()
for (n in N_values) {
  message("Estimating power at N = ", n)
  for (scenario in c("Directional", "Equivalence")) {
    successes <- sum(replicate(nsim, simulate(n, scenario)))
    ci <- binom.test(successes, nsim)$conf.int
    results[[length(results)+1]] <- data.frame(
      N=n, test=scenario, power=successes/nsim,
      MC_lower=ci[1], MC_upper=ci[2])
  }
}
power_curve <- do.call(rbind, results)
print(power_curve, row.names=FALSE, digits=4)
subset(power_curve, N == 350)

# Results for the planned analysable sample.
print(subset(power_curve, N == target_N), row.names=FALSE, digits=4)

# 6. PLOT POWER CURVES
colours <- c(Directional="#2369A0", Equivalence="#C66B18")
symbols <- c(Directional=16, Equivalence=17)

power_curve_png <- function() {
  par(mar=c(4.5, 4.5, 1, 1))
  plot(NA, xlim=range(N_values), ylim=c(0,1),
     xlab="Analysable participants (N)", ylab="Estimated power",
     yaxt="n", bty="l")
  axis(2, at=seq(0,1,.1), labels=paste0(seq(0,100,10), "%"), las=1)
  abline(h=.95, lty=2, col="grey35")
  abline(v=target_N, lty=3, col="grey60")
  
  for (scenario in names(colours)) {
    d <- subset(power_curve, test == scenario)
    arrows(d$N, d$MC_lower, d$N, d$MC_upper,
           angle=90, code=3, length=.025, col=colours[scenario])
    lines(d$N, d$power, type="b", pch=symbols[scenario],
          lwd=2, col=colours[scenario])
  }
  
  legend("bottomright", bty="n", cex=.85,
         legend=c("Directional: true Spearman rho = -0.20",
                  "Equivalence (+/-0.20): true rho = 0",
                  "95% power", paste("Planned N =", target_N)),
         col=c("#2369A0", "#C66B18", "grey35", "grey60"),
         pch=c(16, 17, NA, NA),
         lty=c(1, 1, 2, 3),
         lwd=c(2, 2, 1, 1),
         pt.cex=1.2, seg.len=2, x.intersp=1)
}

power_curve_png()
png("power_curve.png", width = 9, height = 5, units = "in", res = 300)
power_curve_png()
dev.off()


