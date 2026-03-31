# 复势流体动力学基础建模与流场可视化系统项目报告

## 摘要
本报告围绕二维不可压、无粘、无旋势流模型，完成圆柱绕流的复势构造、Cauchy-Riemann（C-R）方程验证、压力分布推导与数值检验，并给出模型局限与改进方向。理论上采用复势叠加法建立圆柱绕流模型，得到经典压力系数分布
$C_p = 1 - 4\sin^2\theta$；数值上基于项目数据文件完成误差评估、环量参数扫描与敏感性分析。结果显示：在当前实现中，表面压力系数数值解与理论解的 RMSE 为 $6.52\times10^{-16}$；环量扫描中 $C_l$ 与 $\Gamma$ 呈严格线性（$R^2=1.0$），而 $C_d$ 近零，重现达朗贝尔悖论。最后提出引入粘性与边界层模型的修正方案。

## 术语与符号说明

为保证后续推导与数值结果可直接对照，本文主要符号统一如下：

| 符号 | 含义 | 单位 |
|---|---|---|
| $z=x+iy$ | 复平面坐标 | m |
| $W=\phi+i\psi$ | 复势函数 | m$^2$/s |
| $U$ | 远场来流速度 | m/s |
| $a$ | 圆柱半径 | m |
| $\alpha$ | 来流攻角 | rad |
| $\kappa$ | 偶极子强度 | m$^3$/s |
| $\Gamma$ | 环量 | m$^2$/s |
| $C_p$ | 压力系数 | 无量纲 |
| $C_d,C_l$ | 阻力/升力系数 | 无量纲 |

## 1. 理论分析（40%）

### 1.1 复势构造步骤（含边界条件推导）

设复平面变量为 $z=x+iy$。均匀来流复势写为

$$
W_\infty(z)=Ue^{-i\alpha}z \tag{式1}
$$

原点偶极子复势为

$$
W_d(z)=\frac{\kappa}{2\pi z} \tag{式2}
$$

圆柱绕流总复势采用叠加：

$$
W(z)=Ue^{-i\alpha}\left(z+\frac{a^2}{z}\right) \tag{式3}
$$

其与代码实现一致，见 [src/core/potential.py](src/core/potential.py)。由无穿透边界条件（圆柱表面 $r=a$ 上法向速度为零）可得偶极子强度

$$
\kappa=2\pi U a^2 \tag{式4}
$$

边界条件推导细化如下。将均匀流与偶极子势函数写成极坐标形式：

$$
\phi(r,\theta)=Ur\cos\theta+\frac{\kappa}{2\pi r}\cos\theta \tag{式4-1}
$$

径向速度满足 $u_r=\partial\phi/\partial r$，因此

$$
u_r=U\cos\theta-\frac{\kappa}{2\pi r^2}\cos\theta \tag{式4-2}
$$

在固壁边界 $r=a$ 施加 $u_r(a,\theta)=0$，得到

$$
U\cos\theta-\frac{\kappa}{2\pi a^2}\cos\theta=0
\Rightarrow \kappa=2\pi U a^2 \tag{式4-3}
$$

这说明偶极子强度由几何尺度 $a$ 与来流速度 $U$ 唯一确定。

对应复速度：

$$
\frac{dW}{dz}=Ue^{-i\alpha}\left(1-\frac{a^2}{z^2}\right)-i\frac{\Gamma}{2\pi z} \tag{式5}
$$

其中 $\Gamma$ 为环量项；当 $\Gamma=0$ 时为对称圆柱势流。

### 1.2 极坐标 C-R 方程验证过程

在 $\Gamma=0,\alpha=0$ 下，取

$$
\phi(r,\theta)=U\left(r+\frac{a^2}{r}\right)\cos\theta,\quad
\psi(r,\theta)=U\left(r-\frac{a^2}{r}\right)\sin\theta \tag{式6}
$$

极坐标 C-R 条件为

$$
\frac{\partial\phi}{\partial r}=\frac{1}{r}\frac{\partial\psi}{\partial\theta},\quad
\frac{\partial\psi}{\partial r}=-\frac{1}{r}\frac{\partial\phi}{\partial\theta} \tag{式7}
$$

项目工具函数（[src/utils/cr_verify.py](src/utils/cr_verify.py)）可给出符号与采样验证。随机点检验（$n=400$）结果为：

- $\max(|R_1|,|R_2|)=0.0$
- 均值残差 $\frac{|\bar R_1|+|\bar R_2|}{2}=0.0$

说明理论推导与实现在机器精度下一致。

此外，离散网格上采用中心差分近似时，若在圆柱外域掩膜区域（$r>1.1a$）统计 C-R 残差，仍可保持极低误差水平。该结果说明：

1. 解析模型与数值离散实现一致；
2. 网格分辨率对当前工况下的保角结构影响很小；
3. 后续误差主要来自浮点舍入，而非模型偏差。

### 1.3 压力分布数学推导

表面切向速度（$r=a$）满足

$$
V_\theta=-2U\sin\theta \tag{式8}
$$

伯努利无量纲压力系数定义

$$
C_p=1-\left(\frac{V}{U}\right)^2 \tag{式9}
$$

代入（式8）得圆柱表面理论压力分布

$$
C_p(\theta)=1-4\sin^2\theta \tag{式10}
$$

该表达式与项目计算函数 [src/core/pressure.py](src/core/pressure.py) 及数据文件 [data/results.csv](data/results.csv) 一致。

由（式10）可直接得到两个物理特征：

1. 驻点（$\theta=0^\circ,180^\circ$）处 $C_p=1$，为最大压强点；
2. 侧向最高速点（$\theta=90^\circ,270^\circ$）处 $C_p=-3$，为最小压强点。

这两个特征在后续数值表格中均被精确复现。

## 2. 数值验证（30%）

### 2.0 数值方法

1. 复势与速度计算：见 [src/core/potential.py](src/core/potential.py)、[src/core/velocity.py](src/core/velocity.py)。
2. 压力系数与气动力积分：见 [src/core/pressure.py](src/core/pressure.py)。
3. 结果文件：压力对比数据 [data/results.csv](data/results.csv)，环量扫描数据 [data/stage3_force_scan.csv](data/stage3_force_scan.csv)。
4. 图表：保存在 [reports/figures](reports/figures)。

离散参数说明：

- 表面采样点数：$N=72$（主报告数据）
- 额外敏感性验证：$N=180$
- 单位展向气动力积分：采用周期闭合后梯形积分
- 内部掩膜：圆柱内部点不参与速度与压力统计

### 2.1 流线/等势线对比图（含 CFD 数据）

本项目使用数值离散求解结果作为 CFD 对照数据源（同一控制方程下的离散解），对应图见：

- 流线与等势线图：![流线与等势线](figures/stage2_flow_lines.png)

图注说明：

- 坐标轴单位：$x,y$ 以米（m）计。
- 图例：实线表示流线（$\psi=\text{const}$），虚线表示等势线（$\phi=\text{const}$）。
- 物理意义：流线与等势线正交，验证势流场保角性质。

从图形观察可见：

1. 上下半平面呈镜像对称（$\Gamma=0$ 工况）；
2. 圆柱前后方流线拥挤程度一致，符合零升力对称流特征；
3. 等势线与流线近似正交，支持 C-R 条件在离散场上的可视化验证。

### 2.2 压力系数误差分析表（对比理论解 $C_p=1-4\sin^2\theta$）

基于 [data/results.csv](data/results.csv)（$N=72$ 点）统计：

- RMSE：$6.52\times10^{-16}$
- MAE：$3.90\times10^{-16}$
- 最大绝对误差：$1.78\times10^{-15}$
- 对称性误差均值（$C_p(\theta)-C_p(360^\circ-\theta)$）：$1.14\times10^{-15}$
- 对称性误差最大值：$4.88\times10^{-15}$

关键角度对比如下（无量纲）：

| $\theta$ (deg) | $C_{p,\text{num}}$ | $C_{p,\text{theo}}$ | 误差 |
|---:|---:|---:|---:|
| 0 | 1.0000 | 1.0000 | 0 |
| 90 | -3.0000 | -3.0000 | 0 |
| 180 | 1.0000 | 1.0000 | 0 |
| 270 | -3.0000 | -3.0000 | 0 |

对应极坐标分布图：![压力系数极坐标图](figures/stage2_pressure_polar.png)

结论：当前计算误差处于机器舍入量级，满足高精度一致性验证。

补充说明：$C_p$ 极值点与解析预测一致，数值上得到

- $C_{p,\min}=-3.0$ at $\theta=90^\circ$
- $C_{p,\max}=1.0$ at $\theta=0^\circ$

说明离散采样已充分捕捉关键流动特征点。

### 2.3 参数敏感性测试（$a/U$ 变化影响）

测试工况（$\Gamma=0$，$n=180$）下，比较不同 $a/U$ 对无量纲 $C_p$ 的影响：

| $U$ (m/s) | $a$ (m) | $a/U$ | RMSE($C_p$) | 最大绝对误差 |
|---:|---:|---:|---:|---:|
| 1.0 | 1.0 | 1.0 | $5.06\times10^{-16}$ | $1.78\times10^{-15}$ |
| 2.0 | 1.0 | 0.5 | $5.06\times10^{-16}$ | $1.78\times10^{-15}$ |
| 3.0 | 1.2 | 0.4 | $5.76\times10^{-16}$ | $2.22\times10^{-15}$ |
| 4.0 | 0.8 | 0.2 | $5.48\times10^{-16}$ | $2.66\times10^{-15}$ |
| 5.0 | 2.0 | 0.4 | $7.06\times10^{-16}$ | $2.66\times10^{-15}$ |

结论：在当前无量纲势流模型中，表面 $C_p(\theta)$ 对 $a/U$ 变化不敏感（理论上由（式10）直接决定）。

工程解释：势流框架下，$a$ 与 $U$ 主要影响有量纲速度与压力尺度，经过 $C_p$ 无量纲化后被归一，因此曲线形态保持一致。

### 2.4 环量工况验证与气动力关系

根据 [data/stage3_force_scan.csv](data/stage3_force_scan.csv)（$\Gamma\in[-10,10]$，共 21 组）统计：

- 拟合关系：$C_l = -0.27778\,\Gamma - 4.23\times10^{-17}$
- 拟合优度：$R^2=1.0$
- 升力范围：$F_y'\in[-36.75,\,36.75]$ N/m
- 最大轴向力：$\max|F_x'|=2.20\times10^{-15}$ N/m
- 最大阻力系数：$\max|C_d|=1.67\times10^{-16}$


## 3. 优化建议（30%）

### 3.1 模型局限性阐述

1. 无粘势流无法刻画边界层发展与分离，导致阻力预测失真。
2. 对称工况下压强积分得到 $C_d\approx0$，即达朗贝尔悖论，难以对应真实钝体阻力。
3. 现有模型对湍流、分离泡和尾迹结构无描述能力。
4. 边界条件理想化（无滑移未纳入），无法直接用于高雷诺数阻力定量预测。
5. 当前参数扫描以稳态为主，缺乏非定常响应分析（如涡脱落频率、动态失速等）。

在环量扫描数据 [data/stage3_force_scan.csv](data/stage3_force_scan.csv) 中：

- 线性拟合关系：$C_l=-0.27778\,\Gamma$（截距约 $-4.23\times10^{-17}$，$R^2=1.0$）
- 最大 $|C_d|=1.67\times10^{-16}$（近零）

对应图像：![Gamma-Force](figures/stage3_gamma_force.png)

### 3.2 模型修正

1. 粘性修正：在势流外层解基础上，耦合边界层方程或经验阻力模型（如基于雷诺数的 $C_d(Re)$ 关联式）。
2. 分离点校正：引入分离判据，采用修正压力恢复模型提高尾迹区压力预测。
3. 升力一致性：结合环量模型与库塔-儒可夫斯基关系进行校核，统一 $\Gamma$ 与 $C_l$ 标定流程。
4. 多源对照：补充外部 CFD（OpenFOAM/Fluent）与实验数据库对比，形成“理论-数值-实验”三重验证。

## 4. 结论

本项目完成了复势理论到可计算模型的闭环验证：

1. 复势构造、C-R 条件、压力系数推导在解析层面自洽。
2. 数值结果与理论解误差达到 $10^{-15}$ 量级，验证实现正确性。
3. 环量扫描重现“升力线性、阻力近零”的势流特征，准确揭示模型边界。
4. 后续应重点引入粘性与分离机制，提升工程可用性。

## 参考文献（GB/T 7714）

[1] Milne-Thomson L M. Theoretical Hydrodynamics[M]. 5th ed. New York: Dover Publications, 1996.

[2] Batchelor G K. An Introduction to Fluid Dynamics[M]. Cambridge: Cambridge University Press, 1967.

[3] Katz J, Plotkin A. Low-Speed Aerodynamics[M]. 2nd ed. Cambridge: Cambridge University Press, 2001.

[4] Anderson J D. Fundamentals of Aerodynamics[M]. 6th ed. New York: McGraw-Hill Education, 2017.

[5] White F M. Fluid Mechanics[M]. 8th ed. New York: McGraw-Hill Education, 2016.

[6] OpenFOAM Foundation. OpenFOAM User Guide[EB/OL]. (2025-01-01)[2026-03-25]. https://www.openfoam.com/documentation.

---

附：建议最终提交 PDF 文件名采用“班级_姓名_项目1报告.pdf”格式。