"""Interactive simulator and high-level orchestration for Project 2."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from matplotlib.patches import Circle

from ..core.potential import total_potential
from ..core.pressure import pressure_coefficient, safe_gamma_range
from ..core.velocity import stagnation_points, velocity_field_numba
from ..utils.grid import create_grid

NU_AIR = 1.5e-5

try:
    from matplotlib.animation import FuncAnimation
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import (
        QApplication,
        QCheckBox,
        QGridLayout,
        QHBoxLayout,
        QLabel,
        QMainWindow,
        QPushButton,
        QSlider,
        QWidget,
    )

    PYQT_AVAILABLE = True
except Exception:  # pragma: no cover
    PYQT_AVAILABLE = False


@dataclass
class FlowSimulator:
    """High-level API for circulation-controlled cylinder flow."""

    U: float = 1.0
    a: float = 1.0
    gamma: float = 0.0
    n: int = 201
    x_min: float = -3.0
    x_max: float = 3.0
    y_min: float = -3.0
    y_max: float = 3.0

    def reynolds_number(self, nu: float = NU_AIR) -> float:
        """Return Reynolds number Re = 2*U*a/nu."""
        if nu <= 0:
            raise ValueError("nu must be positive.")
        return (2.0 * self.U * self.a) / nu

    def run(self) -> dict[str, np.ndarray]:
        """Compute potential, velocity and pressure fields."""
        X, Y, Z = create_grid(
            x_min=self.x_min,
            x_max=self.x_max,
            y_min=self.y_min,
            y_max=self.y_max,
            n=self.n,
        )
        W = total_potential(Z, U=self.U, a=self.a, gamma=self.gamma)

        u, v, v_mag, valid_mask = velocity_field_numba(X, Y, U=self.U, a=self.a, gamma=self.gamma)
        cp = pressure_coefficient(v_mag, U_inf=self.U)

        phi = np.where(valid_mask, np.real(W), np.nan)
        psi = np.where(valid_mask, np.imag(W), np.nan)
        cp = np.where(valid_mask, cp, np.nan)

        return {
            "X": X,
            "Y": Y,
            "Z": Z,
            "W": W,
            "phi": phi,
            "psi": psi,
            "u": u,
            "v": v,
            "V": v_mag,
            "cp": cp,
            "valid_mask": valid_mask,
            "stagnation": stagnation_points(U=self.U, a=self.a, gamma=self.gamma),
        }


if PYQT_AVAILABLE:

    class FlowSimulatorWindow(QMainWindow):
        """PyQt window for dynamic flow visualization."""

        def __init__(self) -> None:
            super().__init__()
            self.setWindowTitle("Project2 Circulation Simulator")
            self.setMinimumSize(1200, 780)

            self.sim = FlowSimulator()
            self.show_potential = True
            self.show_stagnation = True
            self.show_pressure = False
            self.show_quiver = False
            self.safety_factor = 1.5
            self.auto_sweep = False
            self.sweep_direction = 1

            self.figure = Figure(figsize=(8.8, 6.6), dpi=100)
            self.canvas = FigureCanvas(self.figure)
            self.ax = self.figure.add_subplot(111)

            root = QWidget()
            self.setCentralWidget(root)
            layout = QHBoxLayout(root)

            controls = self._build_controls()
            layout.addWidget(controls, 0)
            layout.addWidget(self.canvas, 1)

            self._dirty = True
            self.animation = FuncAnimation(self.figure, self._animate, interval=40, blit=False)
            self._draw_scene()

        def _build_slider(self, min_v: int, max_v: int, value: int, step: int = 1) -> QSlider:
            slider = QSlider(Qt.Orientation.Horizontal)
            slider.setMinimum(min_v)
            slider.setMaximum(max_v)
            slider.setSingleStep(step)
            slider.setValue(value)
            return slider

        def _build_controls(self) -> QWidget:
            panel = QWidget()
            grid = QGridLayout(panel)

            title = QLabel("环量调控参数")
            title.setStyleSheet("font-size: 18px; font-weight: 600;")
            grid.addWidget(title, 0, 0, 1, 2)

            self.lbl_a = QLabel()
            self.sld_a = self._build_slider(5, 20, 10)
            self.sld_a.valueChanged.connect(self._on_slider_change)

            self.lbl_u = QLabel()
            self.sld_u = self._build_slider(2, 20, 2)
            self.sld_u.valueChanged.connect(self._on_slider_change)

            self.lbl_gamma = QLabel()
            self.lbl_gamma_norm = QLabel()
            self.sld_gamma_norm = self._build_slider(-200, 200, 0)
            self.sld_gamma_norm.valueChanged.connect(self._on_slider_change)
            self.lbl_gamma_critical = QLabel()
            self.lbl_gamma_status = QLabel()
            self.lbl_gamma_status.setStyleSheet("font-size: 13px; font-weight: 600;")

            self.chk_phi = QCheckBox("显示等势线")
            self.chk_phi.setChecked(True)
            self.chk_phi.stateChanged.connect(self._on_toggle)

            self.chk_stag = QCheckBox("显示驻点")
            self.chk_stag.setChecked(True)
            self.chk_stag.stateChanged.connect(self._on_toggle)

            self.chk_cp = QCheckBox("显示压力分布")
            self.chk_cp.setChecked(False)
            self.chk_cp.stateChanged.connect(self._on_toggle)

            self.chk_quiver = QCheckBox("显示速度箭头")
            self.chk_quiver.setChecked(False)
            self.chk_quiver.stateChanged.connect(self._on_toggle)

            self.lbl_re = QLabel()
            self.lbl_re.setStyleSheet("font-size: 14px; color: #17496e;")

            self.chk_auto = QCheckBox("自动扫 Gamma")
            self.chk_auto.setChecked(False)
            self.chk_auto.stateChanged.connect(self._on_toggle_auto)

            self.btn_auto_toggle = QPushButton("开始扫频")
            self.btn_auto_toggle.clicked.connect(self._toggle_auto_sweep)

            self.lbl_sweep_speed = QLabel()
            self.sld_sweep_speed = self._build_slider(1, 20, 6)
            self.sld_sweep_speed.valueChanged.connect(self._on_sweep_speed_change)

            btn_refresh = QPushButton("立即刷新")
            btn_refresh.clicked.connect(self._force_redraw)

            row = 1
            grid.addWidget(QLabel("圆柱半径 a (m)"), row, 0)
            grid.addWidget(self.sld_a, row, 1)
            row += 1
            grid.addWidget(self.lbl_a, row, 1)

            row += 1
            grid.addWidget(QLabel("来流速度 U (m/s)"), row, 0)
            grid.addWidget(self.sld_u, row, 1)
            row += 1
            grid.addWidget(self.lbl_u, row, 1)

            row += 1
            grid.addWidget(QLabel("环量 Gamma (m^2/s)"), row, 0)
            grid.addWidget(self.sld_gamma_norm, row, 1)
            row += 1
            grid.addWidget(self.lbl_gamma, row, 1)
            row += 1
            grid.addWidget(self.lbl_gamma_norm, row, 1)
            row += 1
            grid.addWidget(self.lbl_gamma_critical, row, 0, 1, 2)
            row += 1
            grid.addWidget(self.lbl_gamma_status, row, 0, 1, 2)

            row += 1
            grid.addWidget(self.chk_auto, row, 0, 1, 2)
            row += 1
            grid.addWidget(self.btn_auto_toggle, row, 0, 1, 2)
            row += 1
            grid.addWidget(QLabel("扫频速度（步/帧）"), row, 0)
            grid.addWidget(self.sld_sweep_speed, row, 1)
            row += 1
            grid.addWidget(self.lbl_sweep_speed, row, 1)

            row += 1
            grid.addWidget(self.chk_phi, row, 0, 1, 2)
            row += 1
            grid.addWidget(self.chk_stag, row, 0, 1, 2)
            row += 1
            grid.addWidget(self.chk_cp, row, 0, 1, 2)
            row += 1
            grid.addWidget(self.chk_quiver, row, 0, 1, 2)
            row += 1
            grid.addWidget(self.lbl_re, row, 0, 1, 2)
            row += 1
            grid.addWidget(btn_refresh, row, 0, 1, 2)

            self._sync_labels()
            return panel

        def _on_slider_change(self) -> None:
            self.sim.a = self.sld_a.value() / 10.0
            self.sim.U = self.sld_u.value() / 2.0

            gamma_critical = 4.0 * np.pi * self.sim.U * self.sim.a
            gamma_ratio = self.sld_gamma_norm.value() / 100.0
            self.sim.gamma = gamma_ratio * gamma_critical

            self._sync_labels()
            self._dirty = True

        def _on_toggle(self) -> None:
            self.show_potential = self.chk_phi.isChecked()
            self.show_stagnation = self.chk_stag.isChecked()
            self.show_pressure = self.chk_cp.isChecked()
            self.show_quiver = self.chk_quiver.isChecked()
            self._dirty = True

        def _on_toggle_auto(self) -> None:
            self.auto_sweep = self.chk_auto.isChecked()
            self.btn_auto_toggle.setText("停止扫频" if self.auto_sweep else "开始扫频")
            self._dirty = True

        def _toggle_auto_sweep(self) -> None:
            self.auto_sweep = not self.auto_sweep
            self.chk_auto.blockSignals(True)
            self.chk_auto.setChecked(self.auto_sweep)
            self.chk_auto.blockSignals(False)
            self.btn_auto_toggle.setText("停止扫频" if self.auto_sweep else "开始扫频")
            self._dirty = True

        def _on_sweep_speed_change(self) -> None:
            self.lbl_sweep_speed.setText(f"speed = {self.sld_sweep_speed.value()} step/frame")

        def _sync_labels(self) -> None:
            gamma_critical, gamma_safe_pos, _ = safe_gamma_range(
                U=self.sim.U,
                a=self.sim.a,
                safety_factor=self.safety_factor,
            )
            gamma_ratio = self.sim.gamma / gamma_critical if gamma_critical > 0 else 0.0

            self.lbl_a.setText(f"a = {self.sim.a:.1f} m")
            self.lbl_u.setText(f"U = {self.sim.U:.1f} m/s")
            self.lbl_gamma.setText(f"Gamma = {self.sim.gamma:.1f} m^2/s")
            self.lbl_gamma_norm.setText(f"Gamma / Gamma_critical = {gamma_ratio:+.2f}")
            self.lbl_gamma_critical.setText(
                f"Gamma_critical = 4*pi*U*a = {gamma_critical:.3f}; "
                f"Gamma_safe(n={self.safety_factor}) = +/-{gamma_safe_pos:.3f}"
            )

            if abs(self.sim.gamma) > gamma_critical:
                status_text = "状态: 超临界（驻点脱离表面）"
                status_color = "#b00020"
            elif abs(self.sim.gamma) > gamma_safe_pos:
                status_text = "状态: 临界前警戒区"
                status_color = "#b26b00"
            else:
                status_text = "状态: 安全环量区"
                status_color = "#0f7b0f"

            self.lbl_gamma_status.setText(status_text)
            self.lbl_gamma_status.setStyleSheet(f"font-size: 13px; font-weight: 600; color: {status_color};")
            self.lbl_re.setText(f"Re = {self.sim.reynolds_number():.2e} (nu=1.5e-5)")
            self.lbl_sweep_speed.setText(f"speed = {self.sld_sweep_speed.value()} step/frame")

        def _force_redraw(self) -> None:
            self._dirty = True
            self._draw_scene()

        def _animate(self, _frame_idx: int):
            if self.auto_sweep:
                speed = self.sld_sweep_speed.value()
                curr = self.sld_gamma_norm.value()
                next_value = curr + self.sweep_direction * speed

                if next_value >= self.sld_gamma_norm.maximum():
                    next_value = self.sld_gamma_norm.maximum()
                    self.sweep_direction = -1
                elif next_value <= self.sld_gamma_norm.minimum():
                    next_value = self.sld_gamma_norm.minimum()
                    self.sweep_direction = 1

                self.sld_gamma_norm.setValue(next_value)

            if self._dirty:
                self._draw_scene()
            return []

        def _draw_scene(self) -> None:
            out = self.sim.run()
            gamma_critical = 4.0 * np.pi * self.sim.U * self.sim.a
            gamma_ratio = self.sim.gamma / gamma_critical if gamma_critical > 0 else 0.0

            self.ax.clear()
            self.ax.set_facecolor("#0b2340")

            if self.show_pressure:
                self.ax.contourf(out["X"], out["Y"], out["cp"], levels=60, cmap="turbo", alpha=0.72)

            self.ax.contour(
                out["X"],
                out["Y"],
                out["psi"],
                levels=42,
                colors="white",
                linewidths=1.0,
                alpha=0.82,
            )

            if self.show_potential:
                self.ax.contour(
                    out["X"],
                    out["Y"],
                    out["phi"],
                    levels=42,
                    colors="#ffd166",
                    linewidths=0.9,
                    linestyles="--",
                    alpha=0.65,
                )

            if self.show_quiver:
                step = max(1, self.sim.n // 25)
                self.ax.quiver(
                    out["X"][::step, ::step],
                    out["Y"][::step, ::step],
                    out["u"][::step, ::step],
                    out["v"][::step, ::step],
                    color="cyan",
                    scale=45,
                    alpha=0.66,
                )

            circle = Circle((0.0, 0.0), self.sim.a, facecolor="#f94144", edgecolor="black", linewidth=1.1, alpha=0.42)
            self.ax.add_patch(circle)

            if self.show_stagnation and out["stagnation"].size > 0:
                self.ax.plot(np.real(out["stagnation"]), np.imag(out["stagnation"]), "b*", markersize=11, zorder=9)

            self.ax.set_xlim(self.sim.x_min, self.sim.x_max)
            self.ax.set_ylim(self.sim.y_min, self.sim.y_max)
            self.ax.set_aspect("equal", adjustable="box")
            self.ax.set_xlabel("x")
            self.ax.set_ylabel("y")
            self.ax.set_title(
                f"Project2 Circulation Control | Gamma/Gamma_critical={gamma_ratio:+.2f}"
            )

            self.ax.text(
                0.02,
                0.02,
                f"Gamma={self.sim.gamma:.2f}, Gamma_critical={gamma_critical:.2f}",
                transform=self.ax.transAxes,
                color="white",
                fontsize=10,
                bbox={"boxstyle": "round,pad=0.25", "facecolor": "black", "alpha": 0.35, "edgecolor": "none"},
            )

            self.canvas.draw_idle()
            self._dirty = False


def launch_simulator_app() -> None:
    """Launch PyQt interactive simulator."""
    if not PYQT_AVAILABLE:
        raise RuntimeError("PyQt5 or Qt backend is not available in current environment.")

    app = QApplication.instance() or QApplication([])
    win = FlowSimulatorWindow()
    win.show()
    app.exec_()


if __name__ == "__main__":
    launch_simulator_app()
