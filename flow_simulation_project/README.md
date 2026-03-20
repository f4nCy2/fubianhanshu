# flow_simulation_project

基于复势理论的二维理想流体仿真项目，包含：
- 复势推导与 C-R 方程验证
- 流场可视化（流线、等势线）
- 压力系数分析
- PyQt 交互式数字孪生界面

## 项目结构

- `notebooks/`: 教学与实验 Notebook
- `src/`: 核心源码
- `data/`: 结果数据
- `reports/`: 图表与报告

## 快速开始

1. 在项目根目录创建并激活 Python 环境
2. 安装依赖：

```bash
python -m pip install -r requirements.txt
```

3. （可选）注册 Notebook 内核：

```bash
python -m ipykernel install --user --name flow-sim --display-name "flow-sim"
```

4. 运行示例：

```python
from src import FlowSimulator
sim = FlowSimulator(U=1.0, a=1.0, n=150)
result = sim.run()
print(result["cp"].shape)
```

## 阶段 2 运行方式

1. 在项目根目录启动 Jupyter：

```bash
python -m jupyter lab
```

2. 打开并运行 Notebook：
- `notebooks/2_流场可视化.ipynb`
- `notebooks/3_压力分析.ipynb`

3. 启动交互式数字孪生界面：

```bash
python -m src.ui.simulator
```

默认控制项包括圆柱半径、来流速度、环量、等势线/驻点/压力分布显示，以及雷诺数实时显示。

## 阶段 3 产物

- 气动力扫描数据：`data/stage3_force_scan.csv`
- 环量-气动力图：`reports/figures/stage3_gamma_force.png`
- 分析报告：`reports/report_stage3_paradox.md`

## 1.4 提交规范对应

- Notebook 分阶段文件：
	- `notebooks/1_复势推导.ipynb`
	- `notebooks/2_流场可视化.ipynb`
	- `notebooks/3_压力分析.ipynb`
- 工具模块：`utils.py`
	- `verify_cr_equations`: C-R 方程验证函数
	- `BernoulliCalculator`: 伯努利计算类
	- `FlowFieldGenerator`: 流场生成器

所有 Notebook 的第 1 个代码单元均输出依赖版本信息，以支持可重复性检查。

## 默认参数

- 圆柱半径 `a = 1.0`
- 来流速度 `U = 1.0`
