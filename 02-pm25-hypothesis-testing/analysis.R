# PM2.5 comparison: estimation and hypothesis testing


march_may <- read.csv("pm25.march_may_2024.csv")
june_july <- read.csv("pm25.june_july_2024.csv")

# Question 1
n_m <- length(march_may$pm25)
n_j <- length(june_july$pm25)
mean_m <- mean(march_may$pm25)
mean_j <- mean(june_july$pm25)
var_m <- var(march_may$pm25)
var_j <- var(june_july$pm25)

tcrit <- qt(0.975, df = n_j - 1)
se_j <- sqrt(var_j / n_j)
ci_j <- mean_j + c(-1, 1) * tcrit * se_j

mean_diff <- mean_j - mean_m
se_diff <- sqrt(var_j / n_j + var_m / n_m)
ci_diff <- mean_diff + c(-1, 1) * 1.96 * se_diff
z_diff <- mean_diff / se_diff
p_diff <- 2 * pnorm(-abs(z_diff))

pm25 <- c(march_may$pm25, june_july$pm25)
x <- c(rep(0, n_m), rep(1, n_j))
fit <- lm(pm25 ~ x)

# Question 2.5
set.seed(2086)
v <- 1
n <- 10
B <- 10000
y <- matrix(rexp(B * n, rate = exp(-v)), nrow = B, ncol = n)
v_hat <- log(rowMeans(y))
v_bc <- v_hat + 1 / (2 * n)
bias_mle <- mean(v_hat) - v
bias_bc <- mean(v_bc) - v
mc_se_bias <- sd(v_hat) / sqrt(B)

# Question 3
n_A <- 50
n_B <- 50
m_A <- 44
m_B <- 36
p_A <- m_A / n_A
p_B <- m_B / n_B
se_A <- sqrt(p_A * (1 - p_A) / n_A)
ci_A <- p_A + c(-1, 1) * 1.96 * se_A

p_gap <- p_A - p_B
se_gap <- sqrt(p_A * (1 - p_A) / n_A + p_B * (1 - p_B) / n_B)
ci_gap <- p_gap + c(-1, 1) * 1.96 * se_gap

p_pool <- (m_A + m_B) / (n_A + n_B)
se_pool <- sqrt(p_pool * (1 - p_pool) * (1 / n_A + 1 / n_B))
z_prop <- p_gap / se_pool
p_prop <- 2 * pnorm(-abs(z_prop))

nll0 <- -((m_A + m_B) * log(p_pool) +
          (n_A + n_B - m_A - m_B) * log(1 - p_pool))
nll1 <- -(m_A * log(p_A) + (n_A - m_A) * log(1 - p_A) +
          m_B * log(p_B) + (n_B - m_B) * log(1 - p_B))
aic0 <- nll0 + 1
aic1 <- nll1 + 2
bic0 <- nll0 + 0.5 * log(n_A + n_B)
bic1 <- nll1 + log(n_A + n_B)

dir.create("output", recursive = TRUE,
           showWarnings = FALSE)
sink("output/R_output.txt")
options(digits = 12)
cat("PM2.5 hypothesis testing: reproducible output\n\n")
cat("Question 1\n")
cat("March-May: n =", n_m, "mean =", mean_m, "variance =", var_m, "\n")
cat("June-July: n =", n_j, "mean =", mean_j, "variance =", var_j, "\n")
cat("t critical =", tcrit, "SE June-July =", se_j, "CI =", ci_j, "\n")
cat("difference =", mean_diff, "SE difference =", se_diff, "CI =", ci_diff, "\n")
cat("z =", z_diff, "p-value =", p_diff, "\n")
cat("lm coefficients:\n")
print(coef(fit))

cat("\nQuestion 2.5\n")
cat("seed = 2086, simulations =", B, ", n =", n, ", v =", v, "\n")
cat("MLE estimated bias =", bias_mle, "\n")
cat("bias-corrected estimated bias =", bias_bc, "\n")
cat("Monte Carlo SE for either estimated bias =", mc_se_bias, "\n")

cat("\nQuestion 3\n")
cat("p_A =", p_A, "CI =", ci_A, "\n")
cat("difference =", p_gap, "SE =", se_gap, "CI =", ci_gap, "\n")
cat("pooled p =", p_pool, "z =", z_prop, "p-value =", p_prop, "\n")
cat("M0: NLL =", nll0, "AIC =", aic0, "BIC =", bic0, "\n")
cat("M1: NLL =", nll1, "AIC =", aic1, "BIC =", bic1, "\n")

cat("\n")
print(sessionInfo())
sink()
