# Ames house prices: multiple regression, Bonferroni, BIC selection

rm(list = ls())

ames <- read.csv("data/ames.ass3.2026.csv")

# Q1.1 Full multiple regression
full.model <- lm(sale_price ~ ., data = ames)
full.summary <- summary(full.model)
print(full.summary)

coefficient.table <- full.summary$coefficients
predictor.table <- coefficient.table[rownames(coefficient.table) != "(Intercept)", , drop = FALSE]
predictor.table <- predictor.table[order(predictor.table[, "Pr(>|t|)"]), , drop = FALSE]
cat("\nPredictors ordered by their ordinary two-sided p-values:\n")
print(predictor.table)

# Q1.2 Bonferroni adjustment
bonferroni.alpha <- 0.05 / nrow(predictor.table)
bonferroni.table <- data.frame(
  predictor = rownames(predictor.table),
  raw.p.value = predictor.table[, "Pr(>|t|)"],
  bonferroni.p.value = p.adjust(
    predictor.table[, "Pr(>|t|)"],
    method = "bonferroni",
    n = nrow(predictor.table)
  )
)
cat("\nBonferroni per-test threshold:", bonferroni.alpha, "\n")
print(bonferroni.table)

# Q1.3 Effects of overall quality and living area
effects <- coefficient.table[c("overall_qual", "gr_liv_area"), , drop = FALSE]
effect.cis <- confint(full.model, c("overall_qual", "gr_liv_area"), level = 0.95)
cat("\nCoefficient estimates for Q1.3:\n")
print(effects)
cat("\n95% confidence intervals for the two coefficients:\n")
print(effect.cis)

# Q1.4 BIC stepwise selection
bic.model <- step(
  full.model,
  direction = "both",
  k = log(nrow(ames)),
  trace = 0
)
cat("\nFinal BIC-selected model:\n")
print(summary(bic.model))
cat("\nFinal BIC-selected coefficients:\n")
print(coef(bic.model))
cat("\nFinal BIC-selected formula:\n")
print(formula(bic.model))

# Q1.5 Prediction for the house in Table 2
new.house <- data.frame(
  lot_area = 9500,
  overall_qual = 7,
  year_built = 1998,
  year_remod = 2004,
  bsmt_fin_sf = 650,
  gr_liv_area = 1650,
  full_bath = 2,
  half_bath = 1,
  bedrooms = 3,
  rooms = 7,
  fireplaces = 1,
  garage_cars = 2
)

# Q1.5(a): confidence interval for the mean sale price.
mean.price.interval <- predict(
  bic.model,
  newdata = new.house,
  interval = "confidence",
  level = 0.95
)
cat("\nQ1.5(a) 95% confidence interval for mean sale price ($000):\n")
print(mean.price.interval)

# Q1.5(b): prediction interval for this individual house.
house.price.interval <- predict(
  bic.model,
  newdata = new.house,
  interval = "prediction",
  level = 0.95
)
cat("\nQ1.5(b) 95% prediction interval for this house ($000):\n")
print(house.price.interval)

# Q1.5(c): compare the asking price with the prediction interval.
asking.price <- 275
cat("\nQ1.5(c) Asking price ($000):", asking.price, "\n")
cat(
  "Is the asking price inside the 95% prediction interval? ",
  asking.price >= house.price.interval[1, "lwr"] &&
    asking.price <= house.price.interval[1, "upr"],
  "\n",
  sep = ""
)
