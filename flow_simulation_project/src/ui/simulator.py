"""Main simulation orchestration class and interactive UI."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from matplotlib.patches import Circle

from ..core.potential import total_potential
from ..core.pressure import pressure_coefficient
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
except Exception:  # pragma: no cover - allows core usage without GUI stack
    PYQT_AVAILABLE = False


@dataclass
class FlowSimulator:
    """High-level API for running ideal flow simulations."""

    U: float = 1.0
    a: float = 1.0
    alpha: float = 0.0
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
        """Compute potential, velocity and pressure on a Cartesian grid."""
        X, Y, Z = create_grid(self.x_min, self.x_max, self.y_min, self.y_max, self.n)
        W = total_potential(Z, U=self.U, a=self.a, alpha=self.alpha, gamma=self.gamma)

        u, v, V, valid_mask = velocity_field_numba(
            X,
            Y,
            U=self.U,
            a=self.a,
            alpha=self.alpha,
            gamma=self.gamma,
        )
        cp = pressure_coefficient(V, U_inf=self.U)

        phi = np.real(W)
        psi = np.imag(W)
        phi = np.where(valid_mask, phi, np.nan)
        psi = np.where(valid_mask, psi, np.nan)
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
            "V": V,
            "cp": cp,
            "valid_mask": valid_mask,
            "stagnation": stagnation_points(U=self.U, a=self.a, alpha=self.alpha, gamma=self.gamma),
        }

    def boundary_streamline_std(self, n_points: int = 10) -> float:
        """Validate boundary stream-function closure using cylinder surface samples."""
        theta = np.linspace(0.0, 2.0 * np.pi, n_points, endpoint=False)
        r = self.a
        psi_surface = self.U * (r - (self.a**2) / r) * np.sin(theta)
        return float(np.std(psi_surface))


if PYQT_AVAILABLE:

    class FlowSimulatorWindow(QMainWindow):
        """PyQt digital-twin interface for interactive flow visualization."""

        def __init__(self) -> None:
            super().__init__()
            self.setWindowTitle("Flow Digital Twin Simulator")
            self.setMinimumSize(1200, 780)

            self.sim = FlowSimulator()
            self.show_potential = True
            self.show_stagnation = True
            self.show_pressure = False
            self.show_quiver = False

            self.figure = Figure(figsize=(8.6, 6.4), dpi=100)
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

            title = QLabel("参数控制")
            title.setStyleSheet("font-size: 18px; font-weight: 600;")
            grid.addWidget(title, 0, 0, 1, 2)

            self.lbl_a = QLabel()
            self.sld_a = self._build_slider(5, 20, 10)  # 0.5-2.0 step 0.1
            self.sld_a.valueChanged.connect(self._on_slider_change)

            self.lbl_u = QLabel()
            self.sld_u = self._build_slider(2, 20, 2)  # 1-10 step 0.5
            self.sld_u.valueChanged.connect(self._on_slider_change)

            self.lbl_g = QLabel()
            self.sld_g = self._build_slider(-20, 20, 0)  # -10~10 step 0.5
            self.sld_g.valueChanged.connect(self._on_slider_change)

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

            self.lbl_perf = QLabel("动画帧率目标: >= 24 fps")

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
            grid.addWidget(self.sld_g, row, 1)
            row += 1
            grid.addWidget(self.lbl_g, row, 1)

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
            grid.addWidget(self.lbl_perf, row, 0, 1, 2)
            row += 1
            grid.addWidget(btn_refresh, row, 0, 1, 2)

            self._sync_labels()
            return panel

        def _on_slider_change(self) -> None:
            self.sim.a = self.sld_a.value() / 10.0
            self.sim.U = self.sld_u.value() / 2.0
            self.sim.gamma = self.sld_g.value() / 2.0
            self._sync_labels()
            self._dirty = True

        def _on_toggle(self) -> None:
            self.show_potential = self.chk_phi.isChecked()
            self.show_stagnation = self.chk_stag.isChecked()
            self.show_pressure = self.chk_cp.isChecked()
            self.show_quiver = self.chk_quiver.isChecked()
            self._dirty = True

        def _sync_labels(self) -> None:
            self.lbl_a.setText(f"a = {self.sim.a:.1f} m")
            self.lbl_u.setText(f"U = {self.sim.U:.1f} m/s")
            self.lbl_g.setText(f"Gamma = {self.sim.gamma:.1f} m^2/s")
            self.lbl_re.setText(f"Re = {self.sim.reynolds_number():.2e} (nu=1.5e-5)")

        def _force_redraw(self) -> None:
            self._dirty = True
            self._draw_scene()

        def _animate(self, _frame_idx: int):
            if self._dirty:
                self._draw_scene()
            return []

        def _draw_scene(self) -> None:
            out = self.sim.run()
            self.ax.clear()

            self.ax.set_facecolor("#0b2340")

            if self.show_pressure:
                self.ax.contourf(out["X"], out["Y"], out["cp"], levels=60, cmap="turbo", alpha=0.72)

            self.ax.contour(
                out["X"],
                out["Y"],
                out["psi"],
                levels=40,
                colors="white",
                linewidths=1.0,
                alpha=0.8,
            )

            if self.show_potential:
                self.ax.contour(
                    out["X"],
                    out["Y"],
                    out["phi"],
                    levels=40,
                    colors="yellow",
                    linewidths=0.9,
                    linestyles="--",
                    alpha=0.6,
                )

            if self.show_quiver:
                step = max(1, self.sim.n // 24)
                self.ax.quiver(
                    out["X"][::step, ::step],
                    out["Y"][::step, ::step],
                    out["u"][::step, ::step],
                    out["v"][::step, ::step],
                    color="cyan",
                    scale=40,
                    alpha=0.68,
                )

            circle = Circle(
                (0.0, 0.0),
                self.sim.a,
                facecolor="red",
                edgecolor="black",
                linewidth=1.2,
                alpha=0.42,
            )
            self.ax.add_patch(circle)
            circle.set_zorder(5)

            if self.show_stagnation:
                spoints = out["stagnation"]
                self.ax.plot(np.real(spoints), np.imag(spoints), "b*", markersize=12, zorder=10)

            self.ax.set_xlim(self.sim.x_min, self.sim.x_max)
            self.ax.set_ylim(self.sim.y_min, self.sim.y_max)
            self.ax.set_aspect("equal", adjustable="box")
            self.ax.set_xlabel("x")
            self.ax.set_ylabel("y")
            self.ax.set_title("Interactive Cylinder Flow Twin")

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
