# flow_project2

基于项目二任务书的二维理想流体仿真项目，聚焦环量调控、驻点迁移、压力不对称与升力验证。

## 项目结构

- src: 核心源码
- data: 扫描结果数据
- reports: 分析报告与图像
- notebooks: 分阶段 Notebook 入口
- p2.md: 项目二任务书转写稿

## 快速开始

1. 在项目根目录创建并激活 Python 环境
2. 安装依赖：

python -m pip install -r requirements.txt

3. 运行基础仿真：

python - <<'PY'
from src import FlowSimulator
sim = FlowSimulator(U=1.0, a=1.0, gamma=2.0, n=151)
out = sim.run()
print(out['cp'].shape)
PY

4. 启动交互式界面：

python -m src.ui.simulator

项目二界面特性：
- 采用归一化环量控制（Gamma/Gamma_critical）而非固定绝对值滑块。
- 实时显示 Gamma_critical=4*pi*U*a 与工程安全区 Gamma_safe。
- 提供安全区/警戒区/超临界状态提示，直接对应任务 3.1。

5. 运行阶段三参数扫描与报告生成：

python -m src.run_pipeline

## 对应任务映射

- 任务 1.2: 复势构造在 src/core/potential.py
- 任务 2.1: 驻点迁移在 src/core/velocity.py
- 任务 2.2: 压力分布与积分在 src/core/pressure.py
- 任务 2.3: 动态交互可视化在 src/ui/simulator.py
- 任务 3.1/3.2: 安全环量与升力验证在 src/run_pipeline.py

## 产物文件

- data/stage3_force_scan.csv
- reports/figures/stage3_gamma_lift.png
- reports/figures/stage2_flow_snapshot.png
- reports/figures/stage2_stagnation_scan.png
- reports/figures/stage2_surface_cp.png
- reports/figures/stage3_gamma_lift.png
- reports/report_stage3_circulation.md

说明：
- Notebook 中的 matplotlib 图默认也会以内嵌输出保存在 ipynb 文件里。
- 现在相关绘图单元已同时 savefig 到 reports/figures，便于单独提交与报告引用。
