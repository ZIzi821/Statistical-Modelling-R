# Hotel booking cancellation: decision tree vs logistic regression

rm(list = ls())

library(rpart)
library(boot)
library(pROC)

source("wrappers.R")
source("my.prediction.stats.R")

hotel.train <- read.csv(
  "data/hotel.train.ass3.2026.csv",
  stringsAsFactors = TRUE
)
hotel.test <- read.csv(
  "data/hotel.test.ass3.2026.csv",
  stringsAsFactors = TRUE
)

# Match factor levels before predicting on the test data.
factor.variables <- c("Market.Segment", "Customer.Type", "Cancelled")
for (variable in factor.variables) {
  hotel.test[[variable]] <- factor(
    hotel.test[[variable]],
    levels = levels(hotel.train[[variable]])
  )
}

# Q2.1 Cross-validated decision tree
set.seed(2086)
tree.cv <- learn.tree.cv(
  Cancelled ~ .,
  data = hotel.train,
  nfolds = 10,
  m = 1000
)
best.tree <- tree.cv$best.tree

tree.variables <- unique(best.tree$frame$var)
tree.variables <- tree.variables[tree.variables != "<leaf>"]
number.of.leaves <- sum(best.tree$frame$var == "<leaf>")

cat("Q2.1 variables used by the CV-selected tree:\n")
print(tree.variables)
cat("Q2.1 number of leaves:", number.of.leaves, "\n")
cat("Q2.1 selected complexity parameter:", tree.cv$best.cp, "\n")
cat("\nAveraged cross-validation results:\n")
print(data.frame(
  cp = tree.cv$cv.stats$cp,
  leaves = tree.cv$cv.stats$n.leaves,
  mean.cv.error = tree.cv$cv.stats$cv.err
))

# Q2.2 Tree plot and terminal-node rules
# White label backgrounds keep the branches clear of the text.
# Long factor splits are stacked so they do not run into nearby leaves.
format.tree.label <- function(label) {
  if (is.na(label) || !startsWith(label, "Market.Segment=")) {
    return(label)
  }

  segment.levels <- strsplit(
    sub("Market.Segment=", "", label, fixed = TRUE),
    ",",
    fixed = TRUE
  )[[1]]
  paste(c("Market.Segment =", segment.levels), collapse = "\n")
}

draw.tree.label <- function(x, y, labels, cex = 1, adj = 0.5, ...) {
  labels <- as.character(labels)
  market.segment.label <- !is.na(labels) & startsWith(labels, "Market.Segment=")
  labels <- vapply(labels, format.tree.label, character(1))
  y[market.segment.label] <- y[market.segment.label] + 0.035
  label.lines <- strsplit(labels, "\n", fixed = TRUE)
  label.width <- vapply(
    label.lines,
    function(z) max(strwidth(z, units = "user", cex = cex)),
    numeric(1)
  )
  line.height <- strheight("M", units = "user", cex = cex)
  label.height <- vapply(
    label.lines,
    function(z) length(z) * line.height * 1.15,
    numeric(1)
  )
  plot.region <- par("usr")
  horizontal.padding <- diff(plot.region[1:2]) * 0.004
  vertical.padding <- diff(plot.region[3:4]) * 0.004
  rect(
    x - label.width / 2 - horizontal.padding,
    y - label.height / 2 - vertical.padding,
    x + label.width / 2 + horizontal.padding,
    y + label.height / 2 + vertical.padding,
    col = "white",
    border = NA,
    xpd = NA
  )
  text(x, y, labels, cex = cex, adj = c(0.5, 0.5), col = "black", ...)
}

png(
  "output/Q2_tree.png",
  width = 1800,
  height = 1200,
  res = 180
)
par(
  mar = c(2, 2, 3, 2),
  fg = "black",
  col = "black",
  col.axis = "black",
  col.lab = "black",
  col.main = "black"
)
plot(best.tree, uniform = TRUE, branch = 0.5, margin = 0.08, lwd = 1.15)
text(
  best.tree,
  use.n = TRUE,
  all = FALSE,
  pretty = 0,
  cex = 0.96,
  FUN = draw.tree.label
)
title("Decision tree selected by 10-fold cross-validation")
dev.off()

leaf.nodes <- as.numeric(row.names(best.tree$frame))[best.tree$frame$var == "<leaf>"]
leaf.paths <- path.rpart(best.tree, nodes = leaf.nodes, print.it = FALSE)

training.node <- as.numeric(row.names(best.tree$frame))[best.tree$where]
training.probability <- predict(best.tree, hotel.train, type = "prob")[, "Y"]
leaf.probabilities <- aggregate(
  training.probability,
  by = list(node = training.node),
  FUN = unique
)
names(leaf.probabilities)[2] <- "probability.cancelled"
leaf.probabilities$predicted.class <- ifelse(
  leaf.probabilities$probability.cancelled > 0.5,
  "Y",
  "N"
)

cat("\nQ2.2 terminal-node rules and cancellation probabilities:\n")
for (i in seq_along(leaf.nodes)) {
  node <- leaf.nodes[i]
  probability <- leaf.probabilities$probability.cancelled[
    leaf.probabilities$node == node
  ]
  predicted.class <- leaf.probabilities$predicted.class[
    leaf.probabilities$node == node
  ]
  cat("\nNode", node, "- P(Cancelled = Y) =", probability,
      "- predicted class =", predicted.class, "\n")
  cat(paste(leaf.paths[[i]], collapse = " AND "), "\n")
}

# Q2.3 Lowest cancellation probability
lowest.row <- which.min(leaf.probabilities$probability.cancelled)
lowest.node <- leaf.probabilities$node[lowest.row]
lowest.path <- leaf.paths[[match(lowest.node, leaf.nodes)]]
cat("\nQ2.3 lowest-probability node:", lowest.node, "\n")
cat("Q2.3 cancellation probability:",
    leaf.probabilities$probability.cancelled[lowest.row], "\n")
cat("Q2.3 rule:", paste(lowest.path, collapse = " AND "), "\n")

# Q2.4 Logistic regression with KIC selection (k = 3)
full.logistic <- glm(
  Cancelled ~ .,
  data = hotel.train,
  family = binomial
)
kic.logistic <- step(
  full.logistic,
  direction = "both",
  k = 3,
  trace = 0
)

cat("\nQ2.4 final KIC-selected logistic model:\n")
print(summary(kic.logistic))
cat("\nQ2.4 final model terms:\n")
print(attr(terms(kic.logistic), "term.labels"))
cat("\nQ2.4 final model coefficients:\n")
print(coef(kic.logistic))

coefficient.summary <- summary(kic.logistic)$coefficients
coefficient.summary <- coefficient.summary[
  rownames(coefficient.summary) != "(Intercept)",
  ,
  drop = FALSE
]
cat("\nQ2.4 non-intercept coefficients ordered by p-value:\n")
print(coefficient.summary[order(coefficient.summary[, "Pr(>|z|)"]), , drop = FALSE])

# Q2.5 Coefficients for the fitted equation
cat("\nQ2.5 fitted log-odds coefficients:\n")
print(coef(kic.logistic), digits = 12)

# Q2.6 Effect of Special.Requests on the odds
special.requests.coefficient <- coef(kic.logistic)["Special.Requests"]
special.requests.odds.ratio <- exp(special.requests.coefficient)
special.requests.ci <- exp(confint.default(kic.logistic, "Special.Requests"))
cat("\nQ2.6 Special.Requests log-odds coefficient:",
    special.requests.coefficient, "\n")
cat("Q2.6 odds ratio per additional request:",
    special.requests.odds.ratio, "\n")
cat("Q2.6 approximate 95% confidence interval for the odds ratio:\n")
print(special.requests.ci)

# Q2.7 Prediction statistics on the test set
tree.test.probability <- predict(
  best.tree,
  newdata = hotel.test,
  type = "prob"
)[, "Y"]
logistic.test.probability <- predict(
  kic.logistic,
  newdata = hotel.test,
  type = "response"
)

png(
  "output/Q2_ROC_curves.png",
  width = 1800,
  height = 800,
  res = 180
)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
cat("\nQ2.7 prediction statistics for the CV-selected tree:\n")
my.pred.stats(tree.test.probability, hotel.test$Cancelled)
title("CV-selected tree")
cat("\nQ2.7 prediction statistics for the KIC logistic model:\n")
my.pred.stats(logistic.test.probability, hotel.test$Cancelled)
title("KIC logistic regression")
dev.off()

# Q2.8 Odds for booking 10
tree.probability.10 <- tree.test.probability[10]
logistic.probability.10 <- logistic.test.probability[10]
tree.odds.10 <- tree.probability.10 / (1 - tree.probability.10)
logistic.odds.10 <- logistic.probability.10 / (1 - logistic.probability.10)
cat("\nQ2.8 tree probability and odds for test booking 10:\n")
print(c(probability = tree.probability.10, odds = tree.odds.10))
cat("Q2.8 logistic probability and odds for test booking 10:\n")
print(c(probability = logistic.probability.10, odds = logistic.odds.10))

# Q2.9 BCa interval for booking 65 minus booking 66
fixed.formula <- formula(kic.logistic)
booking.65.and.66 <- hotel.test[c(65, 66), , drop = FALSE]

log.odds.difference <- function(data, indices, model.formula, new.bookings) {
  bootstrap.fit <- glm(
    model.formula,
    data = data[indices, , drop = FALSE],
    family = binomial
  )
  predicted.log.odds <- predict(
    bootstrap.fit,
    newdata = new.bookings,
    type = "link"
  )
  unname(predicted.log.odds[1] - predicted.log.odds[2])
}

observed.log.odds <- predict(
  kic.logistic,
  newdata = booking.65.and.66,
  type = "link"
)
observed.difference <- unname(observed.log.odds[1] - observed.log.odds[2])

set.seed(208609)
bootstrap.log.odds <- boot(
  data = hotel.train,
  statistic = log.odds.difference,
  R = 5000,
  model.formula = fixed.formula,
  new.bookings = booking.65.and.66
)
bootstrap.ci <- boot.ci(bootstrap.log.odds, conf = 0.95, type = "bca")

cat("\nQ2.9 observed predicted log-odds (bookings 65 and 66):\n")
print(observed.log.odds)
cat("Q2.9 observed difference, booking 65 minus booking 66:",
    observed.difference, "\n")
cat("Q2.9 BCa 95% confidence interval:\n")
print(bootstrap.ci)
cat("Q2.9 point estimate and BCa interval on the odds-ratio scale:\n")
print(c(
  odds.ratio = exp(observed.difference),
  lower = exp(bootstrap.ci$bca[4]),
  upper = exp(bootstrap.ci$bca[5])
))
