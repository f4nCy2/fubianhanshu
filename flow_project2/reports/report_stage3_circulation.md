# 项目二：环量调控与升力优化报告

## 参数
- U = 1.0
- a = 1.0
- rho = 1.225

## 关键结论
- 临界环量：|Gamma_critical| = 12.566371
- 推荐安全环量区间：Gamma in [-8.377580, 8.377580]
- 升力数值验证最大相对误差：6.661338e-16

## 工程解释
- 当 |Gamma| > 4*pi*U*a 时，表面驻点消失，流场进入非安全区。
- 在安全环量区间内，升力与 Gamma 近线性关系满足 Kutta-Joukowski 定理。

## 图像
![Gamma-Lift 曲线](/Users/fancy2/work/code/fubianhanshu/flow_project2/reports/figures/stage3_gamma_lift.png)
