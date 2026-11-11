# Statistical-Modelling-R

用 R 完成的统计建模与预测分析练习合集：概率建模、参数估计与假设检验、回归、分类、非参数平滑。每个子目录可独立运行，输出与图表已一并提交。

A collection of statistical modelling exercises in R — probability models, estimation and hypothesis testing, regression, classification and k-NN smoothing. Every folder runs on its own and its outputs are committed.

作者：Yanyu Zhu

## 目录

| 目录 | 问题 | 方法 | 主要结果 |
|---|---|---|---|
| [`01-ev-arrivals-poisson`](01-ev-arrivals-poisson) | 电动车充电站每小时到达数能否用 Poisson 描述？ | Poisson 概率计算、经验分布 vs 理论分布、负对数似然模型比较 | 1,200 小时样本方差/均值 = 2.45，明显过度离散；零到达与 ≥10 次到达的经验频率远高于 Poisson(4.8) |
| [`02-pm25-hypothesis-testing`](02-pm25-hypothesis-testing) | 6–7 月 PM2.5 是否低于 3–5 月？ | t 区间、两样本 z 检验、MLE 偏差的 Monte Carlo 模拟（10,000 次）、比例检验、AIC/BIC | 均值差 −2.27 μg/m³，95% CI [−4.08, −0.45]，p = 0.014 |
| [`03-predictive-modelling`](03-predictive-modelling) | 房价预测、酒店预订取消分类、质谱信号平滑 | 多元线性回归 + Bonferroni + BIC 逐步选择；决策树（10 折 CV 剪枝）vs Logistic 回归（KIC 选择）、ROC、BCa bootstrap；k-NN 不同 k 的 RMSE 与 CV 选 k | 房价模型 Adj. R² = 0.87；测试集准确率 树 0.736 / Logistic 0.768；k-NN 最优 k = 7（RMSE 1.957），CV 选 k = 6 |

## 运行

R ≥ 4.3。进入对应目录后：

```r
# 01 / 02
source("analysis.R")

# 03（需要 rpart、boot、pROC、kknn）
install.packages(c("boot", "pROC", "kknn"))
source("q1_house_price_regression.R")
source("q2_hotel_cancellation_classification.R")
source("q3_knn_signal_smoothing.R")
```

结果写入各目录下的 `output/`。

## 说明

- `03-predictive-modelling/wrappers.R`、`my.prediction.stats.R` 为分类评估辅助函数（混淆矩阵、准确率、灵敏度、特异度等）。
- 数据均为小型教学数据集（Ames 房价子集、酒店预订子集、模拟质谱信号、EV 到达计数与 PM2.5 观测）。
- 新的独立建模项目会单独开仓库；本仓库只收录练习规模的建模分析。
