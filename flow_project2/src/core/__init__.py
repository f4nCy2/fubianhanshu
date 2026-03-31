"""Core flow-physics modules for Project 2."""

from .potential import total_potential
from .pressure import (
    cylinder_surface_pressure,
    export_gamma_scan_csv,
    integrate_lift_from_cp,
    lift_theory_kj,
    pressure_coefficient,
    safe_gamma_range,
    scan_gamma_lift,
    write_project2_report,
)
from .velocity import (
    complex_velocity,
    stagnation_angles,
    stagnation_points,
    surface_tangential_velocity,
    velocity_field_numba,
)

__all__ = [
    "total_potential",
    "pressure_coefficient",
    "cylinder_surface_pressure",
    "integrate_lift_from_cp",
    "lift_theory_kj",
    "safe_gamma_range",
    "scan_gamma_lift",
    "export_gamma_scan_csv",
    "write_project2_report",
    "complex_velocity",
    "velocity_field_numba",
    "surface_tangential_velocity",
    "stagnation_angles",
    "stagnation_points",
]
