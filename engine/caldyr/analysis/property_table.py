"""Property table: stream properties over a (T, P) grid — the HYSYS
"Stream Analysis → Property Table" tool (Hameed, *Chemical Process Simulations
using Aspen HYSYS*, Wiley 2025, §2.1.4), used there to tabulate/plot n-pentane
mass density over T = 500–600 K and P = 12–18 atm.

A pure function over any :class:`~caldyr.thermo.base.PropertyPackage`: give it
a composition and one-or-many values of T and P, get plot-ready 2-D arrays
back. Failed flash points are skipped gracefully (NaN + a logged failure) so a
grid that wanders out of a backend's validity range still returns the rest.
"""
from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

from ..core.components_db import molar_mass

#: Supported dependent properties -> docstring of what each column holds.
PROPERTIES: dict[str, str] = {
    "mass_density": "bulk mass density, kg/m^3",
    "molar_volume": "bulk molar volume, m^3/mol",
    "enthalpy": "molar enthalpy, J/mol (formation-inclusive engine basis)",
    "entropy": "molar entropy, J/mol/K",
    "vapor_fraction": "molar vapor fraction (0 liquid .. 1 vapor)",
}

DEFAULT_PROPS = ("mass_density", "enthalpy", "vapor_fraction")


def _as_grid(name: str, value: float | Sequence[float]) -> np.ndarray:
    arr = np.atleast_1d(np.asarray(value, dtype=float))
    if arr.ndim != 1 or arr.size == 0:
        raise ValueError(f"{name} must be a scalar or a 1-D sequence of values")
    if not np.all(np.isfinite(arr)) or np.any(arr <= 0.0):
        raise ValueError(f"{name} values must be finite and positive (got {arr})")
    return arr


def _point(pp: Any, t: float, p: float, z: dict[str, float], mw: float,
           props: Sequence[str]) -> dict[str, float]:
    res = pp.flash_pt(t, p, z)
    out: dict[str, float] = {}
    for name in props:
        if name == "enthalpy":
            out[name] = res.H
        elif name == "vapor_fraction":
            out[name] = res.vapor_fraction
        elif name == "entropy":
            out[name] = pp.entropy(t, p, z)
        else:                                     # mass_density / molar_volume
            v = pp.volume(t, p, z)
            out[name] = mw / v if name == "mass_density" else v
    return out


def property_table(
    pp: Any,
    z: dict[str, float],
    *,
    T: float | Sequence[float],
    P: float | Sequence[float],
    props: Sequence[str] = DEFAULT_PROPS,
) -> dict[str, Any]:
    """Evaluate ``props`` for composition ``z`` over the grid ``T`` x ``P``.

    Mirrors the HYSYS Property Table (Hameed 2025, §2.1.4): one or both of
    T/P may be an array (a scalar is a 1-point axis), and every property comes
    back as a 2-D array over the full grid — sweep T at fixed P values, or P at
    fixed T values, or both at once.

    Parameters
    ----------
    pp : PropertyPackage for the components appearing in ``z``.
    z : composition ``{component: mole fraction}`` (normalized internally).
    T, P : grid values, K and Pa (scalar or 1-D sequence each).
    props : property names from :data:`PROPERTIES`.

    Returns
    -------
    dict with keys
      * ``"T"`` — 1-D array (n_T,), K
      * ``"P"`` — 1-D array (n_P,), Pa
      * one ``(n_T, n_P)`` array per requested property (NaN where the flash
        failed)
      * ``"failures"`` — list of ``(T, P, message)`` for skipped points.

    Plot-ready: ``plt.plot(out["T"], out["mass_density"][:, j])`` is one
    isobar, exactly the curves of the book's Figure 2.20.
    """
    bad = [name for name in props if name not in PROPERTIES]
    if bad:
        raise ValueError(
            f"unknown property name(s) {bad}; expected a subset of "
            f"{sorted(PROPERTIES)}"
        )
    if not props:
        raise ValueError("props must name at least one property")

    total = sum(z.values())
    if total <= 0.0:
        raise ValueError(f"composition sums to {total}; expected > 0")
    z_norm = {c: v / total for c, v in z.items()}
    mw = sum(frac * molar_mass(c) for c, frac in z_norm.items())   # kg/mol

    t_grid = _as_grid("T", T)
    p_grid = _as_grid("P", P)

    arrays = {name: np.full((t_grid.size, p_grid.size), np.nan) for name in props}
    failures: list[tuple[float, float, str]] = []
    for i, t in enumerate(t_grid):
        for j, p in enumerate(p_grid):
            try:
                values = _point(pp, float(t), float(p), z_norm, mw, props)
            except Exception as exc:  # noqa: BLE001 — any backend failure skips the point
                failures.append((float(t), float(p), f"{type(exc).__name__}: {exc}"))
                continue
            for name, value in values.items():
                arrays[name][i, j] = value

    out: dict[str, Any] = {"T": t_grid, "P": p_grid, "failures": failures}
    out.update(arrays)
    return out
