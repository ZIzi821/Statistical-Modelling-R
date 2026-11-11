# A/B testing workflow: sample size, SRM check, conversion test, guardrail,
# CUPED variance reduction, and why "peeking" inflates false positives.
#
# The experiment data are simulated with a fixed seed so that the true effect
# is known; the point is to show the method end to end and check it works.
# Run from this folder:  Rscript ab_test_analysis.R

set.seed(2086)
dir.create("output", showWarnings = FALSE)
sink("output/R_output.txt", split = TRUE)
cat("A/B testing workflow - reproducible output\n\n")

# 1. Sample size --------------------------------------------------------------
# Baseline conversion 10%; we want to detect +1 percentage point (10% -> 11%)
# with alpha = 0.05 (two-sided) and power = 0.80.
plan <- power.prop.test(p1 = 0.10, p2 = 0.11, sig.level = 0.05, power = 0.80)
n_per_group <- ceiling(plan$n)
cat("1. Sample size per group for 10% -> 11%:", n_per_group, "\n\n")

# 2. Simulated experiment -----------------------------------------------------
# Users are randomised 50/50. Treatment truly lifts conversion to 11%.
# Minutes spent per user are correlated with each user's pre-experiment minutes;
# treatment truly adds 0.5 minute on average.
N <- 2 * n_per_group
group <- ifelse(runif(N) < 0.5, "control", "treatment")
pre_minutes <- rgamma(N, shape = 4, scale = 5)                  # pre-period minutes
p_convert <- ifelse(group == "treatment", 0.11, 0.10)
converted <- rbinom(N, 1, p_convert)
minutes <- 0.8 * pre_minutes + rgamma(N, shape = 2, scale = 4) + ifelse(group == "treatment", 0.5, 0)
users <- data.frame(group, converted, minutes, pre_minutes)

# 3. Sample ratio mismatch (SRM) check ---------------------------------------
counts <- table(users$group)
srm <- chisq.test(counts, p = c(0.5, 0.5))
cat("3. Users per group:\n"); print(counts)
cat(sprintf("   SRM chi-square p-value = %.4f (%s)\n\n", srm$p.value,
            ifelse(srm$p.value < 0.001, "MISMATCH - stop and investigate", "no mismatch")))

# 4. Primary metric: conversion ---------------------------------------------
x <- tapply(users$converted, users$group, sum)
n <- tapply(users$converted, users$group, length)
test <- prop.test(x[c("treatment", "control")], n[c("treatment", "control")], correct = FALSE)
lift <- unname(x["treatment"] / n["treatment"] - x["control"] / n["control"])
cat("4. Conversion\n")
cat(sprintf("   control = %.4f, treatment = %.4f\n", x["control"] / n["control"], x["treatment"] / n["treatment"]))
cat(sprintf("   absolute lift = %.4f, 95%% CI [%.4f, %.4f], relative lift = %.1f%%\n",
            lift, test$conf.int[1], test$conf.int[2], 100 * lift / (x["control"] / n["control"])))
cat(sprintf("   two-proportion z-test p-value = %.4f\n\n", test$p.value))

# 5. Minutes per user, with and without CUPED -------------------------------
# CUPED: Y_adj = Y - theta * (X - mean(X)), theta = cov(Y, X) / var(X).
theta <- cov(users$minutes, users$pre_minutes) / var(users$pre_minutes)
users$minutes_cuped <- users$minutes - theta * (users$pre_minutes - mean(users$pre_minutes))
t_raw <- t.test(minutes ~ factor(group, levels = c("treatment", "control")), data = users)
t_cuped <- t.test(minutes_cuped ~ factor(group, levels = c("treatment", "control")), data = users)
width <- function(t) diff(t$conf.int)
cat("5. Minutes per user (treatment - control; true effect = 0.5)\n")
cat(sprintf("   raw:   diff = %.3f, 95%% CI [%.3f, %.3f], p = %.4f\n",
            -diff(t_raw$estimate), t_raw$conf.int[1], t_raw$conf.int[2], t_raw$p.value))
cat(sprintf("   CUPED: diff = %.3f, 95%% CI [%.3f, %.3f], p = %.4f\n",
            -diff(t_cuped$estimate), t_cuped$conf.int[1], t_cuped$conf.int[2], t_cuped$p.value))
cat(sprintf("   CUPED theta = %.3f; CI width reduced by %.1f%%\n\n",
            theta, 100 * (1 - width(t_cuped) / width(t_raw))))

# 6. Peeking: repeated looks at an A/A test ----------------------------------
# No true effect. Analysts check the p-value every day for 14 days and stop at
# the first p < 0.05. How often do they wrongly declare a winner?
n_sims <- 2000; days <- 14; users_per_day <- 1000
peek_hit <- logical(n_sims); fixed_hit <- logical(n_sims)
for (s in seq_len(n_sims)) {
  a <- cumsum(rbinom(days, users_per_day, 0.10))
  b <- cumsum(rbinom(days, users_per_day, 0.10))
  m <- users_per_day * seq_len(days)
  p_pool <- (a + b) / (2 * m)
  z <- (b / m - a / m) / sqrt(p_pool * (1 - p_pool) * 2 / m)
  p <- 2 * pnorm(-abs(z))
  peek_hit[s] <- any(p < 0.05)
  fixed_hit[s] <- p[days] < 0.05
}
cat("6. A/A test false-positive rate over", n_sims, "simulations\n")
cat(sprintf("   check once at day %d:      %.1f%%\n", days, 100 * mean(fixed_hit)))
cat(sprintf("   check every day, stop early: %.1f%%\n\n", 100 * mean(peek_hit)))

png("output/peeking_false_positives.png", width = 1400, height = 800, res = 180)
barplot(c(`Fixed horizon\n(check once)` = mean(fixed_hit), `Peek daily,\nstop at p<0.05` = mean(peek_hit)),
        col = c("#2b6cb0", "#c05621"), ylim = c(0, max(0.3, mean(peek_hit) * 1.2)),
        ylab = "False-positive rate", main = "A/A test: peeking inflates false positives")
abline(h = 0.05, lty = 2)
invisible(dev.off())

sink()
