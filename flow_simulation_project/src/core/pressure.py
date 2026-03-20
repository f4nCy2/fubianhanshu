"""Pressure coefficient helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def pressure_coefficient(velocity_mag: np.ndarray, U_inf: float = 1.0) -> np.ndarray:
    """功能：根据速度场计算压力系数 Cp。

    输入：
        velocity_mag: 速度模长数组 |V|。
        U_inf: 远场参考速度，单位 m/s。

    输出：
        Cp 数组，公式 Cp = 1 - (V/U_inf)^2。

    物理参数说明：
        U_inf 是伯努利无量纲化的参考速度尺度。
    """
    if U_inf <= 0:
        raise ValueError("U_inf must be positive.")
    return 1.0 - (velocity_mag / U_inf) ** 2


def cylinder_surface_cp(
    U: float,
    a: float,
    gamma: float = 0.0,
    n_points: int = 72,
) -> dict[str, np.ndarray]:
    """功能：采样圆柱表面并计算数值/理论压力系数及误差。

    输入：
        U: 远场来流速度，单位 m/s。
        a: 圆柱半径，单位 m。
        gamma: 环量，单位 m^2/s。
        n_points: 圆周采样点数。

    输出：
        包含 theta、速度分量、cp_num、cp_theo、error 的字典。

    物理参数说明：
        gamma=0 时理论分布为 Cp = 1 - 4*sin^2(theta)。
    """
    if U <= 0:
        raise ValueError("U must be positive.")
    if a <= 0:
        raise ValueError("a must be positive.")
    if n_points <= 0:
        raise ValueError("n_points must be positive.")

    theta = np.linspace(0.0, 2.0 * np.pi, n_points, endpoint=False)
    z_surf = a * np.exp(1j * theta)
    v_surf = U * (1.0 - (a**2) / (z_surf**2)) - 1j * gamma / (2.0 * np.pi * z_surf)
    v_mag = np.abs(v_surf)

    cp_num = 1.0 - (v_mag / U) ** 2
    cp_theo = 1.0 - 4.0 * (np.sin(theta) ** 2)
    err = cp_num - cp_theo

    return {
        "theta": theta,
        "z_surf": z_surf,
        "u": np.real(v_surf),
        "v": -np.imag(v_surf),
        "v_mag": v_mag,
        "cp_num": cp_num,
        "cp_theo": cp_theo,
        "error": err,
    }


def rmse(a: np.ndarray, b: np.ndarray) -> float:
    """功能：计算两个数组的均方根误差。

    输入：
        a, b: 待比较数组。

    输出：
        RMSE 标量值。

    物理参数说明：
        常用于评估数值结果与理论结果偏差。
    """
    return float(np.sqrt(np.mean((a - b) ** 2)))


def aerodynamic_force_from_surface_cp(
    theta: np.ndarray,
    cp: np.ndarray,
    U: float,
    a: float,
    rho: float = 1.225,
) -> dict[str, float]:
    """功能：由表面压力系数积分得到单位展向气动力及系数。

    输入：
        theta: 圆周角数组，单位 rad。
        cp: 对应压力系数数组。
        U: 远场来流速度，单位 m/s。
        a: 圆柱半径，单位 m。
        rho: 流体密度，单位 kg/m^3。

    输出：
        包含 Fx_per_span、Fy_per_span、Cd、Cl 的字典。

    物理参数说明：
        q_inf = 0.5*rho*U^2 为动压；参考面积取单位展向 2a。
    """
    if U <= 0:
        raise ValueError("U must be positive.")
    if a <= 0:
        raise ValueError("a must be positive.")
    if rho <= 0:
        raise ValueError("rho must be positive.")
    if theta.size != cp.size:
        raise ValueError("theta and cp must have the same length.")

    # Integrate in ascending theta over one period.
    order = np.argsort(theta)
    th = theta[order]
    cp_sorted = cp[order]

    # Ensure periodic closure point for trapezoidal integration.
    th_ext = np.concatenate([th, [th[0] + 2.0 * np.pi]])
    cp_ext = np.concatenate([cp_sorted, [cp_sorted[0]]])

    int_cx = np.trapezoid(cp_ext * np.cos(th_ext), th_ext)
    int_cy = np.trapezoid(cp_ext * np.sin(th_ext), th_ext)

    q_inf = 0.5 * rho * U**2
    fx = -q_inf * a * int_cx
    fy = -q_inf * a * int_cy

    # Reference area per unit span: 2a*1
    cd = fx / (q_inf * 2.0 * a)
    cl = fy / (q_inf * 2.0 * a)

    return {
        "Fx_per_span": float(fx),
        "Fy_per_span": float(fy),
        "Cd": float(cd),
        "Cl": float(cl),
    }


def circulation_force_study(
    U: float,
    a: float,
    gamma_values: np.ndarray,
    n_points: int = 720,
    rho: float = 1.225,
) -> pd.DataFrame:
    """功能：扫描不同环量，统计对应气动力变化。

    输入：
        U: 远场来流速度，单位 m/s。
        a: 圆柱半径，单位 m。
        gamma_values: 环量采样数组，单位 m^2/s。
        n_points: 圆周离散点数。
        rho: 流体密度，单位 kg/m^3。

    输出：
        含 gamma、Fx_per_span、Fy_per_span、Cd、Cl 的 DataFrame。

    物理参数说明：
        可用于观察升力系数 Cl 与 gamma 的关系趋势。
    """
    rows: list[dict[str, float]] = []
    for gamma in gamma_values:
        surf = cylinder_surface_cp(U=U, a=a, gamma=float(gamma), n_points=n_points)
        force = aerodynamic_force_from_surface_cp(surf["theta"], surf["cp_num"], U=U, a=a, rho=rho)
        rows.append(
            {
                "gamma": float(gamma),
                "Fx_per_span": force["Fx_per_span"],
                "Fy_per_span": force["Fy_per_span"],
                "Cd": force["Cd"],
                "Cl": force["Cl"],
            }
        )
    return pd.DataFrame(rows)


def export_surface_cp_csv(output_csv: str | Path, surface_data: dict[str, np.ndarray]) -> pd.DataFrame:
    """功能：将表面压力采样数据导出为 CSV。

    输入：
        output_csv: 输出文件路径。
        surface_data: 表面采样结果字典。

    输出：
        与导出内容一致的 DataFrame。

    物理参数说明：
        导出字段包含 Cp 数值/理论值及速度分量，便于后处理。
    """
    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame(
        {
            "theta_deg": np.degrees(surface_data["theta"]),
            "cp_num": surface_data["cp_num"],
            "cp_theo": surface_data["cp_theo"],
            "error": surface_data["error"],
            "u": surface_data["u"],
            "v": surface_data["v"],
            "v_mag": surface_data["v_mag"],
        }
    )
    df.to_csv(output_csv, index=False)
    return df


def write_pressure_report(
    report_path: str | Path,
    figure_path: str | Path,
    surface_data: dict[str, np.ndarray],
    rmse_value: float,
    U: float,
    a: float,
    gamma: float,
) -> None:
    """功能：生成阶段 2 压力分析 Markdown 报告。

    输入：
        report_path: 报告输出路径。
        figure_path: 图像路径。
        surface_data: 表面采样结果字典。
        rmse_value: 数值与理论 Cp 的 RMSE。
        U, a, gamma: 工况参数。

    输出：
        无返回值；在磁盘写入报告文件。

    物理参数说明：
        U、a、gamma 定义了当前势流与环量工况。
    """
    report_path = Path(report_path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    figure_path = Path(figure_path)

    content = f"""# 阶段 2.2 压力分析报告

## 参数
- U = {U}
- a = {a}
- Gamma = {gamma}
- 采样点数 N = {surface_data['theta'].size}

## 结果
- RMSE = {rmse_value:.6f}
- 驻点理论值: theta=0°,180° 处 Cp=1
- 最小压力理论值: theta=90°,270° 处 Cp=-3

## 图像
![压力分布图]({figure_path.as_posix()})

## 结论
{'满足 RMSE < 0.01 要求。' if rmse_value < 0.01 else '未满足 RMSE < 0.01 要求，请提高网格或检查参数。'}
"""
    report_path.write_text(content, encoding="utf-8")


def write_stage3_report(
    report_path: str | Path,
    force_df: pd.DataFrame,
    figure_path: str | Path,
    U: float,
    a: float,
    rho: float,
) -> None:
    """功能：生成阶段 3 悖论与环量修正分析报告。

    输入：
        report_path: 报告输出路径。
        force_df: 气动力扫描结果表。
        figure_path: 图像路径。
        U, a, rho: 工况参数。

    输出：
        无返回值；在磁盘写入报告文件。

    物理参数说明：
        通过 gamma 扫描可分析达朗贝尔悖论下升力/阻力预测特征。
    """
    report_path = Path(report_path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    figure_path = Path(figure_path)

    zero_gamma = force_df.loc[np.isclose(force_df["gamma"], 0.0)]
    if zero_gamma.empty:
        fx0 = float("nan")
        fy0 = float("nan")
        cd0 = float("nan")
        cl0 = float("nan")
    else:
        row = zero_gamma.iloc[0]
        fx0, fy0, cd0, cl0 = row["Fx_per_span"], row["Fy_per_span"], row["Cd"], row["Cl"]

    content = f"""# 阶段 3.1 气动力悖论与环量修正分析

## 工况参数
- U = {U}
- a = {a}
- rho = {rho}

## 理想模型验证（Gamma=0）
- Fx' = {fx0:.6e} N/m
- Fy' = {fy0:.6e} N/m
- Cd = {cd0:.6e}
- Cl = {cl0:.6e}

结论：理想势流模型预测阻力约为 0，升力约为 0（对称工况），体现达朗贝尔悖论。

## 环量修正预研
- 通过引入 Gamma，升力系数 Cl 随 Gamma 呈近线性变化。
- 阻力系数 Cd 仍接近 0，说明仅靠无粘势流+环量仍无法预测真实阻力。

## 图像
![Gamma-Force]({figure_path.as_posix()})

## 模型局限与后续方向
- 局限：当前模型不含粘性和边界层分离，无法给出实际阻力。
- 修正思路：可引入库塔-儒可夫斯基关系评估升力，并结合粘性模型/经验阻力模型改进阻力预测。
"""
    report_path.write_text(content, encoding="utf-8")
