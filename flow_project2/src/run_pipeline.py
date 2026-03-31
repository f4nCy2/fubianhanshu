"""End-to-end pipeline for Project 2 deliverables."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .core.pressure import (
    export_gamma_scan_csv,
    safe_gamma_range,
    scan_gamma_lift,
    write_project2_report,
)
from .visualization.plot_pressure import plot_gamma_lift_curve


def run_project2_pipeline(
    U: float = 1.0,
    a: float = 1.0,
    rho: float = 1.225,
    safety_factor: float = 1.5,
) -> dict[str, str]:
    """功能：执行项目二阶段三参数扫描并导出图表与报告。

    输入：
        U: 来流速度，单位 m/s。
        a: 圆柱半径，单位 m。
        rho: 流体密度，单位 kg/m^3。
        safety_factor: 安全系数 n。

    输出：
        关键产物路径字典。

    物理参数说明：
        扫描区间基于临界环量 |Gamma|<=4*pi*U*a。
    """
    root = Path(__file__).resolve().parents[1]
    data_dir = root / "data"
    report_dir = root / "reports"
    fig_dir = report_dir / "figures"

    data_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    gamma_critical, gamma_safe, _ = safe_gamma_range(U=U, a=a, safety_factor=safety_factor)

    gamma_values = np.linspace(-gamma_critical, gamma_critical, 81)
    scan_df = scan_gamma_lift(U=U, a=a, rho=rho, gamma_values=gamma_values)

    csv_path = data_dir / "stage3_force_scan.csv"
    export_gamma_scan_csv(csv_path, scan_df)

    fig, _ = plot_gamma_lift_curve(scan_df)
    fig_path = fig_dir / "stage3_gamma_lift.png"
    fig.savefig(fig_path, dpi=180, bbox_inches="tight")

    max_err = float(np.nanmax(scan_df["relative_error"].to_numpy()))
    report_path = report_dir / "report_stage3_circulation.md"
    write_project2_report(
        report_path=report_path,
        figure_path=fig_path,
        U=U,
        a=a,
        rho=rho,
        gamma_critical=gamma_critical,
        gamma_safe=gamma_safe,
        max_relative_error=max_err,
    )

    return {
        "csv": str(csv_path),
        "figure": str(fig_path),
        "report": str(report_path),
    }


if __name__ == "__main__":
    outputs = run_project2_pipeline()
    for key, value in outputs.items():
        print(f"{key}: {value}")
