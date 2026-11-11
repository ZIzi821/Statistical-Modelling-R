# k-NN signal smoothing and RMSE-based selection of k

rm(list = ls())

library(kknn)
library(boot)

ms.measured <- read.csv("data/ms.measured.2026.csv")
ms.truth <- read.csv("data/ms.truth.2026.csv")

predict.spectrum <- function(k, train = ms.measured, test.mz = ms.truth$MZ) {
  fit <- kknn(
    intensity ~ MZ,
    train = train,
    test = data.frame(MZ = test.mz),
    k = k,
    kernel = "optimal"
  )
  as.numeric(fitted(fit))
}

rmse <- function(estimate, truth) {
  sqrt(mean((estimate - truth)^2))
}

# Q3.1 RMSE for k = 1, ..., 25
k.values <- 1:25
spectrum.estimates <- lapply(k.values, predict.spectrum)
rmse.values <- vapply(
  spectrum.estimates,
  rmse,
  numeric(1),
  truth = ms.truth$intensity
)
rmse.table <- data.frame(k = k.values, RMSE = rmse.values)

cat("Q3.1 RMSE for k = 1, ..., 25:\n")
print(rmse.table, row.names = FALSE)

actual.best.k <- rmse.table$k[which.min(rmse.table$RMSE)]
cat("Q3.1 value of k with the smallest actual RMSE:", actual.best.k, "\n")
cat("Q3.1 minimum actual RMSE:", min(rmse.table$RMSE), "\n")

png(
  "output/Q3_RMSE_by_k.png",
  width = 1500,
  height = 1000,
  res = 180
)
plot(
  rmse.table$k,
  rmse.table$RMSE,
  type = "b",
  pch = 1,
  col = "black",
  xlab = "Number of neighbours (k)",
  ylab = "RMSE against the true spectrum",
  main = "Spectrum estimation error for k-NN smoothing"
)
points(
  actual.best.k,
  min(rmse.table$RMSE),
  pch = 19,
  col = "black",
  cex = 1.3
)
text(
  actual.best.k,
  min(rmse.table$RMSE),
  labels = paste0("  minimum at k = ", actual.best.k),
  pos = 4,
  col = "black"
)
grid(col = "grey85")
dev.off()

# Q3.2 Four fitted spectra
requested.k <- c(2, 6, 12, 25)
png(
  "output/Q3_four_spectra.png",
  width = 1800,
  height = 1250,
  res = 180
)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 1, 0))
for (k in requested.k) {
  estimate <- spectrum.estimates[[k]]
  plot(
    ms.measured$MZ,
    ms.measured$intensity,
    pch = 16,
    cex = 0.30,
    col = adjustcolor("grey45", alpha.f = 0.60),
    xlab = "Mass-to-charge ratio (MZ)",
    ylab = "Relative intensity",
    main = paste0("k = ", k, "   (RMSE = ", format(rmse.values[k], digits = 4), ")"),
    ylim = range(c(ms.measured$intensity, ms.truth$intensity))
  )
  lines(ms.truth$MZ, ms.truth$intensity, col = "black", lwd = 2.1, lty = 1)
  lines(ms.truth$MZ, estimate, col = "black", lwd = 1.7, lty = 2)
  legend(
    "topright",
    legend = c("Measurements", "True spectrum", "k-NN estimate"),
    col = c("grey45", "black", "black"),
    pch = c(16, NA, NA),
    lty = c(NA, 1, 2),
    lwd = c(NA, 2.1, 1.7),
    pt.cex = 0.7,
    cex = 0.72,
    bg = "white"
  )
}
dev.off()

# Q3.3 Numerical comparison
comparison.table <- data.frame(
  k = requested.k,
  RMSE = rmse.values[requested.k],
  maximum.estimated.intensity = vapply(
    spectrum.estimates[requested.k],
    max,
    numeric(1)
  )
)
cat("\nQ3.3 comparison of the four requested smoothers:\n")
print(comparison.table, row.names = FALSE)
cat("True maximum intensity:", max(ms.truth$intensity), "\n")

# Q3.4 Cross-validation
set.seed(2086)
cv.fit <- train.kknn(
  intensity ~ MZ,
  data = ms.measured,
  kmax = 25,
  kernel = "optimal"
)
cv.best.k <- cv.fit$best.parameters$k
cv.spectrum.estimate <- predict.spectrum(cv.best.k)
cv.actual.rmse <- rmse(cv.spectrum.estimate, ms.truth$intensity)

cat("\nQ3.4 cross-validation result:\n")
print(cv.fit)
cat("Q3.4 CV-selected k:", cv.best.k, "\n")
cat("Q3.4 actual-RMSE-minimising k:", actual.best.k, "\n")
cat("Q3.4 actual RMSE at the CV-selected k:", cv.actual.rmse, "\n")

# Q3.5 Noise standard deviation from residuals
fitted.at.measured.mz <- predict.spectrum(
  cv.best.k,
  train = ms.measured,
  test.mz = ms.measured$MZ
)
noise.residuals <- ms.measured$intensity - fitted.at.measured.mz
estimated.noise.sd <- sd(noise.residuals)

cat("\nQ3.5 estimated measurement-noise standard deviation:",
    estimated.noise.sd, "\n")

# Q3.6 Peak heights for the four plotted values of k
true.peak.index <- which.max(ms.truth$intensity)
peak.comparison <- data.frame(
  k = requested.k,
  estimate.at.true.peak = vapply(
    spectrum.estimates[requested.k],
    function(x) x[true.peak.index],
    numeric(1)
  ),
  maximum.estimate = vapply(
    spectrum.estimates[requested.k],
    max,
    numeric(1)
  )
)
cat("\nQ3.6 peak-height comparison for the requested k values:\n")
print(peak.comparison, row.names = FALSE)
cat("True peak MZ:", ms.truth$MZ[true.peak.index], "\n")
cat("True peak intensity:", ms.truth$intensity[true.peak.index], "\n")

# Q3.7 Maximum of the CV-selected curve
estimated.peak.index <- which.max(cv.spectrum.estimate)
estimated.peak.mz <- ms.truth$MZ[estimated.peak.index]
estimated.peak.intensity <- cv.spectrum.estimate[estimated.peak.index]
cat("\nQ3.7 MZ at the maximum CV-selected estimated intensity:",
    estimated.peak.mz, "\n")
cat("Q3.7 maximum estimated intensity:", estimated.peak.intensity, "\n")

# Q3.8 Pairs bootstrap at the peak from Q3.7
bootstrap.knn.intensity <- function(data, indices, target.mz, k) {
  bootstrap.data <- data[indices, , drop = FALSE]
  bootstrap.fit <- kknn(
    intensity ~ MZ,
    train = bootstrap.data,
    test = data.frame(MZ = target.mz),
    k = k,
    kernel = "optimal"
  )
  as.numeric(fitted(bootstrap.fit))
}

bootstrap.k.values <- unique(c(cv.best.k, 3, 20))
bootstrap.summary <- data.frame(
  k = bootstrap.k.values,
  estimate = NA_real_,
  lower = NA_real_,
  upper = NA_real_
)

set.seed(208608)
for (i in seq_along(bootstrap.k.values)) {
  k <- bootstrap.k.values[i]
  result <- boot(
    data = ms.measured,
    statistic = bootstrap.knn.intensity,
    R = 5000,
    target.mz = estimated.peak.mz,
    k = k
  )
  interval <- boot.ci(result, conf = 0.95, type = "perc")
  bootstrap.summary$estimate[i] <- result$t0
  bootstrap.summary$lower[i] <- interval$percent[4]
  bootstrap.summary$upper[i] <- interval$percent[5]
}

cat("\nQ3.8 95% percentile bootstrap confidence intervals:\n")
print(bootstrap.summary, row.names = FALSE)
