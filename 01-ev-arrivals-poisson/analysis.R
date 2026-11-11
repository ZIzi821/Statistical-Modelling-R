# Poisson modelling of EV charging arrivals
# R calculations used in Questions 3 and 4

options(stringsAsFactors = FALSE, scipen = 999)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) == 1L) {
  normalizePath(sub("^--file=", "", script_arg), winslash = "/")
} else {
  normalizePath("analysis.R", winslash = "/", mustWork = FALSE)
}
script_dir <- dirname(script_path)
output_dir <- file.path(script_dir, "output")
figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(output_dir, "R_output.txt")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

cat("EV arrivals: reproducible R output\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

# Question 3.1 -------------------------------------------------------------
x <- seq(0, 1, length.out = 1001)
density_x <- function(x, theta) 1 + theta * (2 * x - 1)

q3_figure <- file.path(figure_dir, "Q3_density_functions.png")
png(q3_figure, width = 1800, height = 1100, res = 180, bg = "white")
par(mar = c(4.5, 4.5, 1.2, 1.0))
plot(
  x, density_x(x, -1), type = "l", lwd = 2, lty = 1,
  xlim = c(0, 1), ylim = c(0, 2),
  xlab = "x", ylab = "p(x | theta)"
)
lines(x, density_x(x, 0), lwd = 2, lty = 2)
lines(x, density_x(x, 1), lwd = 2, lty = 3)
legend(
  "top", horiz = TRUE, bty = "n", lwd = 2,
  lty = c(1, 2, 3),
  legend = c("theta = -1", "theta = 0", "theta = 1")
)
dev.off()
cat("Question 3.1 figure:", normalizePath(q3_figure, winslash = "/"), "\n\n")

# Question 4.1 -------------------------------------------------------------
q4_1 <- ppois(6, lambda = 4.8, lower.tail = FALSE)
cat(sprintf("Question 4.1: P(N >= 7) = %.10f\n", q4_1))

# Question 4.2 -------------------------------------------------------------
q4_2 <- ppois(24, lambda = 4 * 4.8)
cat(sprintf("Question 4.2: P(N_hour <= 24) = %.10f\n", q4_2))

# Question 4.3 -------------------------------------------------------------
prob_at_least_five <- ppois(4, lambda = 4.8, lower.tail = FALSE)
q4_3 <- pbinom(2, size = 4, prob = prob_at_least_five, lower.tail = FALSE)
cat(sprintf("Question 4.3: q = P(N >= 5) = %.10f\n", prob_at_least_five))
cat(sprintf("Question 4.3: P(K >= 3), K ~ Bin(4, q) = %.10f\n\n", q4_3))

# Question 4.4 -------------------------------------------------------------
csv_candidates <- c(
  file.path(script_dir, "ev.arrivals.csv")
)
csv_path <- csv_candidates[file.exists(csv_candidates)][1]
if (is.na(csv_path)) {
  stop("ev.arrivals.csv was not found. Place it beside this script.")
}

arrivals_data <- read.csv(csv_path)
if (!identical(names(arrivals_data), "arrivals")) {
  stop("Expected one column named 'arrivals'.")
}
if (nrow(arrivals_data) != 1200L || anyNA(arrivals_data$arrivals)) {
  stop("The arrivals data must contain 1200 non-missing observations.")
}

observed <- arrivals_data$arrivals
k <- 0:max(observed)
empirical_prob <- as.numeric(prop.table(table(factor(observed, levels = k))))
poisson_prob <- dpois(k, lambda = 4.8)

fit_summary <- data.frame(
  arrivals = k,
  empirical_probability = empirical_prob,
  poisson_4_8_probability = poisson_prob
)
write.csv(
  fit_summary,
  file.path(output_dir, "Q4_4_empirical_vs_poisson.csv"),
  row.names = FALSE
)

q4_figure <- file.path(figure_dir, "Q4_empirical_vs_poisson.png")
png(q4_figure, width = 1800, height = 1100, res = 180, bg = "white")
par(mar = c(4.8, 4.8, 1.2, 1.0))
plot(
  k, empirical_prob, type = "h", lwd = 5, col = "grey60",
  xlim = range(k), ylim = c(0, max(c(empirical_prob, poisson_prob)) * 1.12),
  xaxt = "n", xlab = "Number of arrivals in 15 minutes",
  ylab = "Proportion / probability"
)
axis(1, at = seq(0, max(k), by = 2))
lines(k, poisson_prob, type = "b", lwd = 1.5, lty = 2, pch = 1)
legend(
  "topright", bty = "n", lwd = c(5, 1.5), lty = c(1, 2),
  pch = c(NA, 1), col = c("grey60", "black"),
  legend = c("Empirical", "Poi(4.8)")
)
dev.off()

n_obs <- length(observed)
sample_mean <- mean(observed)
sample_variance <- var(observed)
empirical_zero <- mean(observed == 0)
poisson_zero <- dpois(0, lambda = 4.8)
empirical_tail <- mean(observed >= 10)
poisson_tail <- ppois(9, lambda = 4.8, lower.tail = FALSE)

cat("Question 4.4 data file:", normalizePath(csv_path, winslash = "/"), "\n")
cat(sprintf("n = %d\n", n_obs))
cat(sprintf("sample mean = %.7f\n", sample_mean))
cat(sprintf("sample variance = %.7f\n", sample_variance))
cat(sprintf("variance-to-mean ratio = %.7f\n", sample_variance / sample_mean))
cat(sprintf("empirical P(N = 0) = %.7f\n", empirical_zero))
cat(sprintf("Poi(4.8) P(N = 0) = %.7f\n", poisson_zero))
cat(sprintf("empirical P(N >= 10) = %.7f\n", empirical_tail))
cat(sprintf("Poi(4.8) P(N >= 10) = %.7f\n", poisson_tail))
cat("Question 4.4 figure:", normalizePath(q4_figure, winslash = "/"), "\n\n")

# Question 4.5 -------------------------------------------------------------
y <- c(2, 3, 5, 6, 9, 8, 5, 4, 2, 1)
lambda_A <- rep(4.8, length(y))
lambda_B <- c(2.0, 3.0, 4.5, 6.5, 8.0, 7.5, 5.5, 4.0, 2.5, 1.5)

nll_A <- -sum(dpois(y, lambda = lambda_A, log = TRUE))
nll_B <- -sum(dpois(y, lambda = lambda_B, log = TRUE))

nll_summary <- data.frame(
  model = c("Model A", "Model B"),
  total_negative_log_likelihood = c(nll_A, nll_B)
)
write.csv(
  nll_summary,
  file.path(output_dir, "Q4_5_negative_log_likelihood.csv"),
  row.names = FALSE
)

cat(sprintf("Question 4.5: Model A total NLL = %.10f\n", nll_A))
cat(sprintf("Question 4.5: Model B total NLL = %.10f\n", nll_B))
cat("Question 4.5: preferred model = Model B (smaller total NLL)\n\n")

# Reproducibility checks ---------------------------------------------------
stopifnot(
  abs(q4_1 - 0.2091954158) < 1e-10,
  abs(q4_2 - 0.8841740201) < 1e-10,
  abs(q4_3 - 0.3489296233) < 1e-10,
  n_obs == 1200L,
  abs(sample_mean - 4.8825) < 1e-12,
  abs(sample_variance - 11.96366347) < 1e-8,
  nll_B < nll_A,
  file.exists(q3_figure),
  file.exists(q4_figure)
)

cat("All reproducibility checks passed.\n\n")
sessionInfo()
