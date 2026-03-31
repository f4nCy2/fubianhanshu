"""Pressure and aerodynamic-force analysis for Project 2."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .velocity import surface_tangential_velocity


def pressure_coefficient(velocity_mag: np.ndarray, U_inf: float = 1.0) -> np.ndarray:
    """功能：计算压力系数 Cp。

    输入：
        velocity_mag: 局部速度模长数组。
        U_inf: 远场参考速度，单位 m/s。

    输出：
        Cp = 1 - (V/U_inf)^2。

    物理参数说明：
        Cp 来自伯努利关系，反映压力相对变化。
    """
    if U_inf <= 0:
        raise ValueError("U_inf must be positive.")
    return 1.0 - (velocity_mag / U_inf) ** 2


def cylinder_surface_pressure(U: float, a: float, gamma: float, rho: float = 1.225, n_points: int = 720) -> dict[str, np.ndarray]:
    """功能：计算圆柱表面压力与压力系数分布。

    输入：
        U: 来流速度，单位 m/s。
        a: 圆柱半径，单位 m。
        gamma: 环量，单位 m^2/s。
        rho: 密度，单位 kg/m^3。
        n_points: 圆周离散点数。

    输出：
        含 theta、v_theta、cp、p_minus_p_inf 数组的字典。

    物理参数说明：
        环量通过切向速度项影响 Cp 对称性与升力积分结果。
    """
    if U <= 0:
        raise ValueError("U must be positive.")
    if a <= 0:
        raise ValueError("a must be positive.")
    if rho <= 0:
        raise ValueError("rho must be positive.")

    theta = np.linspace(0.0, 2.0 * np.pi, n_points, endpoint=False)
    v_theta = surface_tangential_velocity(theta=theta, U=U, a=a, gamma=gamma)
    cp = 1.0 - (v_theta / U) ** 2
    q_inf = 0.5 * rho * U**2
    p_minus_p_inf = q_inf * cp

    return {
        "theta": theta,
        "v_theta": v_theta,
        "cp": cp,
        "p_minus_p_inf": p_minus_p_inf,
    }


def integrate_lift_from_cp(theta: np.ndarray, cp: np.ndarray, U: float, a: float, rho: float = 1.225) -> float:
    """功能：通过表面 Cp 数值积分计算单位展向升力。

    输入：
        theta: 圆周角数组，单位 rad。
        cp: 压力系数数组。
        U: 来流速度，单位 m/s。
        a: 圆柱半径，单位 m。
        rho: 密度，单位 kg/m^3。

    输出：
        升力 L'（单位 N/m）。

    物理参数说明：
        公式 L' = -q_inf*a*∮Cp*sin(theta)dtheta。
    """
    if theta.size != cp.size:
        raise ValueError("theta and cp must have the same size.")

    order = np.argsort(theta)
    th = theta[order]
    cp_sorted = cp[order]

    th_ext = np.concatenate([th, [th[0] + 2.0 * np.pi]])
    cp_ext = np.concatenate([cp_sorted, [cp_sorted[0]]])

    q_inf = 0.5 * rho * U**2
    lift = -q_inf * a * np.trapezoid(cp_ext * np.sin(th_ext), th_ext)
    return float(lift)


def lift_theory_kj(rho: float, U: float, gamma: float) -> float:
    """功能：计算库塔-茹科夫斯基理论升力。

    输入：
        rho: 密度，单位 kg/m^3。
        U: 来流速度，单位 m/s。
        gamma: 环量，单位 m^2/s。

    输出：
        理论升力 L' = rho*U*Gamma（单位 N/m）。

    物理参数说明：
        该关系是环量与升力的核心定量联系。
    """
    return float(rho * U * gamma)


def safe_gamma_range(U: float, a: float, safety_factor: float = 1.5) -> tuple[float, float, float]:
    """功能：计算临界环量及工程安全环量。

    输入：
        U: 来流速度，单位 m/s。
        a: 圆柱半径，单位 m。
        safety_factor: 安全系数 n。

    输出：
        (gamma_critical, gamma_safe, gamma_safe_negative)。

    物理参数说明：
        |Gamma|>4*pi*U*a 时表面驻点消失；工程上采用 Gamma_safe=Gamma_critical/n。
    """
    if U <= 0 or a <= 0:
        raise ValueError("U and a must be positive.")
    if safety_factor <= 1.0:
        raise ValueError("safety_factor must be greater than 1.0.")

    gamma_critical = 4.0 * np.pi * U * a
    gamma_safe = gamma_critical / safety_factor
    return float(gamma_critical), float(gamma_safe), float(-gamma_safe)


def scan_gamma_lift(
    U: float,
    a: float,
    rho: float,
    gamma_values: np.ndarray,
    n_points: int = 1000,
) -> pd.DataFrame:
    """功能：扫描环量并比较理论/数值升力。

    输入：
        U: 来流速度，单位 m/s。
        a: 圆柱半径，单位 m。
        rho: 密度，单位 kg/m^3。
        gamma_values: 待扫描环量数组。
        n_points: 圆周积分点数。

    输出：
        含 gamma、lift_num、lift_theory、relative_error 的 DataFrame。

    物理参数说明：
        用于验证 L'-Gamma 的线性关系和误差水平。
    """
    rows: list[dict[str, float]] = []
    for gamma in gamma_values:
        surf = cylinder_surface_pressure(U=U, a=a, gamma=float(gamma), rho=rho, n_points=n_points)
        lift_num = integrate_lift_from_cp(theta=surf["theta"], cp=surf["cp"], U=U, a=a, rho=rho)
        lift_theory = lift_theory_kj(rho=rho, U=U, gamma=float(gamma))

        if np.isclose(lift_theory, 0.0):
            relative_error = 0.0 if np.isclose(lift_num, 0.0) else np.nan
        else:
            relative_error = abs(lift_num / lift_theory - 1.0)

        rows.append(
            {
                "gamma": float(gamma),
                "lift_num": float(lift_num),
                "lift_theory": float(lift_theory),
                "relative_error": float(relative_error),
            }
        )

    return pd.DataFrame(rows)


def export_gamma_scan_csv(output_csv: str | Path, df: pd.DataFrame) -> None:
    """功能：导出环量扫描结果 CSV。

    输入：
        output_csv: 输出 CSV 路径。
        df: 扫描结果 DataFrame。

    输出：
        无返回值；将文件写入磁盘。

    物理参数说明：
        结果用于工程参数选型与报告复现。
    """
    output_path = Path(output_csv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)


def write_project2_report(
    report_path: str | Path,
    figure_path: str | Path,
    U: float,
    a: float,
    rho: float,
    gamma_critical: float,
    gamma_safe: float,
    max_relative_error: float,
) -> None:
    """功能：生成项目二阶段三分析报告。

    输入：
        report_path: 报告路径。
        figure_path: 图像路径。
        U, a, rho: 基础物理参数。
        gamma_critical: 临界环量。
        gamma_safe: 推荐安全环量绝对值。
        max_relative_error: 升力验证最大相对误差。

    输出：
        无返回值；写入 Markdown 报告。

    物理参数说明：
        报告综合展示驻点约束与升力验证结果。
    """
    report = Path(report_path)
    report.parent.mkdir(parents=True, exist_ok=True)
    fig = Path(figure_path)

    content = f"""# 项目二：环量调控与升力优化报告

## 参数
- U = {U}
- a = {a}
- rho = {rho}

## 关键结论
- 临界环量：|Gamma_critical| = {gamma_critical:.6f}
- 推荐安全环量区间：Gamma in [{-gamma_safe:.6f}, {gamma_safe:.6f}]
- 升力数值验证最大相对误差：{max_relative_error:.6e}

## 工程解释
- 当 |Gamma| > 4*pi*U*a 时，表面驻点消失，流场进入非安全区。
- 在安全环量区间内，升力与 Gamma 近线性关系满足 Kutta-Joukowski 定理。

## 图像
![Gamma-Lift 曲线]({fig.as_posix()})
"""
    report.write_text(content, encoding="utf-8")
