"""Petroleum assay characterization — distillation curve to pseudo-components.

The gateway to refining workflows (Hameed, *Chemical Process Simulations using
Aspen HYSYS*, Wiley 2025, ch. 10): a crude oil is described not by a component
list but by a laboratory *assay* — a boiling-point distillation curve (TBP or
ASTM D86) plus bulk/curve density and molecular-weight data.
:func:`characterize_assay` slices that curve into boiling-point cuts and turns
each cut into a *pseudo-component* with a full set of estimated constants
(MW, SG, Tc, Pc, omega, ideal-gas Cp), ready to drop into a Caldyr
:class:`~caldyr.core.flowsheet.Flowsheet` under a cubic-EOS property package
(``thermo:PR`` / ``thermo:SRK`` — see :mod:`caldyr.thermo.thermo_pkg`).

Correlations used (each cited at its function):

* ASTM D86 -> TBP: Riazi-Daubert (API) interconversion.
* MW(Tb, SG): Riazi-Daubert two-parameter (extended 1987 / API TDB 2B2.1 form).
* Tc, Pc (Tb, SG): Kesler-Lee 1976 (Riazi-Daubert 1980 also provided).
* omega: Lee-Kesler 1975 (Tbr <= 0.8, via `chemicals.LK_omega`), Kesler-Lee
  1976 Watson-K form above.
* Ideal-gas Cp: Kesler-Lee 1976 quadratic in T from Watson K.

Unit conventions follow the engine: SI in/out (K, Pa, kg/mol); specific
gravity (60F/60F) and API gravity are the customary dimensionless exceptions.
Distillation curves are sequences of ``(vol_pct, T_K)`` with ``vol_pct`` in
percent (0-100).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Sequence

from .core.component import Component

_RHO_WATER_60F = 999.016        # kg/m3 at 60 F — the SG / API reference density

# Standard ("60F/60F") liquid specific gravities of common light ends, used to
# convert light-end liquid-volume fractions to mass/moles. Values from GPA
# Standard 2145 (table reproduced in the API Technical Data Book, ch. 1, and
# Riazi 2005 Table 2.1); methane/ethane/N2/CO2/H2S are the conventional
# pseudo-liquid values used in petroleum LV% accounting.
LIGHT_END_SG: dict[str, float] = {
    "methane": 0.3000,
    "ethane": 0.35643,
    "propane": 0.50699,
    "isobutane": 0.56287,
    "n-butane": 0.58401,
    "isopentane": 0.62470,
    "n-pentane": 0.63112,
    "n-hexane": 0.66383,
    "nitrogen": 0.80687,
    "carbon dioxide": 0.81802,
    "hydrogen sulfide": 0.78934,
    "water": 1.00000,
}


# -- unit / bulk-property helpers ---------------------------------------------
def api_to_sg(api: float) -> float:
    """Specific gravity (60F/60F) from API gravity: SG = 141.5/(131.5 + API)."""
    return 141.5 / (131.5 + api)


def sg_to_api(sg: float) -> float:
    """API gravity from specific gravity (60F/60F)."""
    return 141.5 / sg - 131.5


def watson_k(tb: float, sg: float) -> float:
    """Watson (UOP) characterization factor K = Tb[R]^(1/3) / SG with the
    boiling point in Rankine (Watson & Nelson, Ind. Eng. Chem. 25 (1933) 880).
    ``tb`` in K."""
    return (1.8 * tb) ** (1.0 / 3.0) / sg


# -- ASTM D86 <-> TBP interconversion ------------------------------------------
# Riazi-Daubert: TBP = a * (D86)^b at fixed vol% knots, both temperatures in
# K. Coefficients from Riazi & Daubert, Hydrocarbon Processing 65(8) (1986)
# (the API Technical Data Book procedure 3A1.1), as tabulated in Riazi,
# "Characterization and Properties of Petroleum Fractions", ASTM MNL50 (2005),
# Table 3.5. Stated D86 validity ranges are roughly 20-320 C per knot; AAD
# ~5 C over 1190 points (Riazi 2005 sec. 3.2.1).
_D86_TBP_COEFFS: tuple[tuple[float, float, float], ...] = (
    # vol%   a       b
    (0.0, 0.9177, 1.0019),
    (10.0, 0.5564, 1.0900),
    (30.0, 0.7617, 1.0425),
    (50.0, 0.9013, 1.0176),
    (70.0, 0.8821, 1.0226),
    (90.0, 0.9552, 1.0110),
    (95.0, 0.8177, 1.0355),
)


def astm_d86_to_tbp(curve: Sequence[tuple[float, float]]) -> list[tuple[float, float]]:
    """Convert an ASTM D86 distillation curve to an atmospheric TBP curve with
    the Riazi-Daubert (API) correlation ``TBP = a * D86^b`` (see
    :data:`_D86_TBP_COEFFS` for the source). The input is interpolated
    (monotone PCHIP) to the correlation's vol% knots that fall inside the data
    range, and the conversion is applied knot-by-knot; the result is the TBP
    curve at those knots, in K.
    """
    pts = _validated_curve(curve, "ASTM D86 curve")
    vols = [p[0] for p in pts]
    interp = _interpolator([p[0] for p in pts], [p[1] for p in pts])
    out: list[tuple[float, float]] = []
    for vol, a, b in _D86_TBP_COEFFS:
        if vols[0] - 1e-9 <= vol <= vols[-1] + 1e-9:
            t_d86 = float(interp(min(max(vol, vols[0]), vols[-1])))
            out.append((vol, a * t_d86**b))
    if len(out) < 3:
        raise ValueError(
            f"ASTM D86 curve spans vol% {vols[0]:.1f}-{vols[-1]:.1f}, covering "
            f"fewer than 3 of the conversion knots "
            f"{[v for v, _, _ in _D86_TBP_COEFFS]}; supply a wider curve"
        )
    return out


def sg_from_distillation(t10: float, t50: float, kind: str = "TBP") -> float:
    """Bulk specific gravity estimated from the 10% and 50% distillation
    temperatures (K): ``SG = a * T10^b * T50^c`` (Riazi-Daubert; Riazi 2005
    Table 3.7, SG range 0.67-1.00). Used as a last-resort bulk-gravity
    fallback when neither an SG curve nor a bulk API gravity is given."""
    coeffs = {"TBP": (0.10431, 0.12550, 0.20862),
              "ASTM_D86": (0.08342, 0.10731, 0.26288)}
    try:
        a, b, c = coeffs[kind]
    except KeyError:
        raise ValueError(f"unknown distillation kind {kind!r}; "
                         f"expected one of {sorted(coeffs)}") from None
    return a * t10**b * t50**c


# -- pseudo-component property correlations -----------------------------------
def riazi_daubert_mw(tb: float, sg: float) -> float:
    """Molar mass (kg/mol) of a petroleum fraction from its mean-average
    boiling point ``tb`` (K) and specific gravity — the Riazi-Daubert
    two-parameter correlation in its extended form (Riazi & Daubert,
    Hydrocarbon Processing 66(3) (1987); API Technical Data Book procedure
    2B2.1), valid roughly MW 70-700 g/mol, Tb 300-850 K; ~2-4% AAD in range
    (Riazi 2005 sec. 2.4.2). Outside that span it is applied as a documented
    extrapolation (heavy resid cuts).
    """
    mw_g = (42.965 * math.exp(2.097e-4 * tb - 7.78712 * sg + 2.08476e-3 * tb * sg)
            * tb**1.26007 * sg**4.98308)
    return mw_g / 1000.0


def riazi_daubert_mw_1980(tb: float, sg: float) -> float:
    """Original Riazi-Daubert two-parameter MW correlation (kg/mol):
    ``MW[g/mol] = 1.6607e-4 Tb^2.1962 SG^-1.0164`` (Riazi & Daubert, Ind. Eng.
    Chem. Process Des. Dev. 19 (1980) 289), valid MW 70-300 g/mol. Kept for
    reference/validation; :func:`riazi_daubert_mw` is the working choice."""
    return 1.6607e-4 * tb**2.1962 * sg**-1.0164 / 1000.0


def riazi_daubert_tc(tb: float, sg: float) -> float:
    """Critical temperature (K): ``Tc = 19.0623 Tb^0.58848 SG^0.3596``
    (Riazi & Daubert 1980, ibid.; API TDB; valid Tb ~ 290-850 K)."""
    return 19.0623 * tb**0.58848 * sg**0.3596


def riazi_daubert_pc(tb: float, sg: float) -> float:
    """Critical pressure (Pa): ``Pc[bar] = 5.53027e7 Tb^-2.3125 SG^2.3201``
    (Riazi & Daubert 1980, ibid.)."""
    return 5.53027e7 * tb**-2.3125 * sg**2.3201 * 1e5


def kesler_lee_tc(tb: float, sg: float) -> float:
    """Critical temperature (K) by Kesler & Lee (Hydrocarbon Processing 55(3)
    (1976) 153; API TDB procedure 4D1.1; Riazi 2005 sec. 2.5.1):
    ``Tc = 189.8 + 450.6 SG + (0.4244 + 0.1174 SG) Tb
    + (0.1441 - 1.0069 SG) 1e5 / Tb``."""
    return (189.8 + 450.6 * sg + (0.4244 + 0.1174 * sg) * tb
            + (0.1441 - 1.0069 * sg) * 1e5 / tb)


def kesler_lee_pc(tb: float, sg: float) -> float:
    """Critical pressure (Pa) by Kesler & Lee (1976, ibid.):
    ln Pc[psia] as a cubic in Tb[R] with SG-dependent coefficients."""
    tb_r = 1.8 * tb
    ln_pc = (8.3634 - 0.0566 / sg
             - (0.24244 + 2.2898 / sg + 0.11857 / sg**2) * 1e-3 * tb_r
             + (1.4685 + 3.648 / sg + 0.47227 / sg**2) * 1e-7 * tb_r**2
             - (0.42019 + 1.6977 / sg**2) * 1e-10 * tb_r**3)
    return math.exp(ln_pc) * 6894.757293            # psia -> Pa


def acentric_factor(tb: float, tc: float, pc: float, kw: float) -> float:
    """Acentric factor of a petroleum fraction. For Tbr = Tb/Tc <= 0.8 the
    Lee-Kesler vapor-pressure correlation (Lee & Kesler, AIChE J. 21 (1975)
    510) via `chemicals.acentric.LK_omega`; above 0.8 (heavy fractions, where
    the LK form degrades) the Kesler-Lee 1976 Watson-K expression
    ``omega = -7.904 + 0.1352 Kw - 0.007465 Kw^2 + 8.359 Tbr
    + (1.408 - 0.01063 Kw)/Tbr`` (Riazi 2005 sec. 2.5.3)."""
    from chemicals.acentric import LK_omega

    tbr = tb / tc
    if tbr <= 0.8:
        return float(LK_omega(tb, tc, pc))
    return (-7.904 + 0.1352 * kw - 0.007465 * kw**2 + 8.359 * tbr
            + (1.408 - 0.01063 * kw) / tbr)


def kesler_lee_cp_ig(mw: float, kw: float) -> tuple[float, float, float]:
    """Ideal-gas heat-capacity polynomial coefficients ``(a0, a1, a2)`` such
    that ``Cp_ig [J/mol/K] = a0 + a1*T + a2*T^2`` (T in K), from the Kesler-Lee
    petroleum-fraction correlation ``Cp_ig [kJ/kg/K] = A0 + A1 T + A2 T^2``
    with A0 = -1.41779 + 0.11828 Kw, A1 = -(6.99724 - 8.69326 Kw
    + 0.27715 Kw^2) 1e-4, A2 = -2.2582e-6 (Kesler & Lee 1976, ibid.; Riazi
    2005 sec. 6.5, stated validity ~255-920 K and Kw 10-12.8; the small
    correction factor for Kw outside that band is neglected). ``mw`` in
    kg/mol. Validated against n-decane's databank ideal-gas Cp (~1% at 300 K,
    see tests/test_m14_assay.py)."""
    mw_g = mw * 1000.0          # kJ/kg/K * g/mol == J/mol/K
    a0 = -1.41779 + 0.11828 * kw
    a1 = -(6.99724 - 8.69326 * kw + 0.27715 * kw**2) * 1e-4
    a2 = -2.2582e-6
    return (mw_g * a0, mw_g * a1, mw_g * a2)


# -- data model ----------------------------------------------------------------
@dataclass(frozen=True)
class PseudoComponent:
    """One boiling-point cut of a characterized assay. SI units (K, Pa,
    kg/mol); fractions are of the WHOLE assay (cuts + light ends)."""
    id: str                     # e.g. "NBP_123C"
    Tb: float                   # normal boiling point, K (mid-volume of cut)
    SG: float                   # specific gravity 60F/60F
    MW: float                   # kg/mol
    Tc: float                   # K
    Pc: float                   # Pa
    omega: float
    watson_k: float
    cp_ig: tuple[float, float, float]   # J/mol/K, ascending powers of T
    vol_frac: float
    mass_frac: float
    mole_frac: float
    T_lo: float                 # cut boundary temperatures, K
    T_hi: float
    vol_mid_pct: float = 0.0    # mid-volume location on the assay curve, LV%

    def pseudo_dict(self) -> dict:
        """The constants payload for :class:`caldyr.core.component.Component`."""
        return {
            "MW": self.MW, "Tb": self.Tb, "SG": self.SG,
            "Tc": self.Tc, "Pc": self.Pc, "omega": self.omega,
            "Cp_ig": list(self.cp_ig),
        }

    def to_component(self) -> Component:
        """A flowsheet :class:`Component` carrying the constants (constructing
        it registers them in the process-wide pseudo registry)."""
        return Component(id=self.id, name=self.id, pseudo=self.pseudo_dict())


@dataclass
class AssayResult:
    """Output of :func:`characterize_assay`.

    ``cuts`` are the pseudo-components (ascending boiling point);
    ``light_end_fractions`` maps light-end component id ->
    {"vol_frac", "mass_frac", "mole_frac"} on the whole-assay basis;
    ``bulk`` carries {"MW" (kg/mol), "SG", "API", "watson_k"} computed from
    the characterized blend; ``curves`` holds plot-ready arrays — keys
    "vol_pct"/"tbp_K" (smoothed working curve) and "input_vol_pct"/"input_T_K"
    (the data as supplied, D86 if that was the input kind).
    """
    cuts: list[PseudoComponent]
    light_end_fractions: dict[str, dict[str, float]] = field(default_factory=dict)
    bulk: dict[str, float] = field(default_factory=dict)
    curves: dict[str, list[float]] = field(default_factory=dict)
    kind: str = "TBP"

    def components(self, include_light_ends: bool = True) -> list[Component]:
        """Flowsheet-ready Component list (registers all pseudo cuts)."""
        from .core.components_db import resolve_component

        comps = [cut.to_component() for cut in self.cuts]
        if include_light_ends:
            comps = [resolve_component(cid) for cid in self.light_end_fractions] + comps
        return comps

    def mole_fractions(self, include_light_ends: bool = True) -> dict[str, float]:
        """Whole-assay mole fractions, renormalized over the included set —
        ready to use as a feed composition."""
        z = {cut.id: cut.mole_frac for cut in self.cuts}
        if include_light_ends:
            z.update({cid: f["mole_frac"] for cid, f in self.light_end_fractions.items()})
        total = sum(z.values())
        if total <= 0.0:
            raise ValueError("assay has no material in the requested fraction set")
        return {cid: v / total for cid, v in z.items()}


# -- characterization ----------------------------------------------------------
def _validated_curve(curve: Sequence[tuple[float, float]],
                     label: str) -> list[tuple[float, float]]:
    pts = sorted((float(v), float(t)) for v, t in curve)
    if len(pts) < 5:
        raise ValueError(
            f"{label} needs at least 5 points (got {len(pts)}) — the same "
            f"minimum HYSYS enforces for assay curves (Hameed 2025 fig. 10.24)"
        )
    for (v0, t0), (v1, t1) in zip(pts, pts[1:]):
        if v1 <= v0:
            raise ValueError(f"{label} vol% values must be strictly increasing "
                             f"(got {v0} then {v1})")
        if t1 <= t0:
            raise ValueError(f"{label} temperatures must be strictly increasing "
                             f"with vol% (got {t0:.2f} K at {v0}% then "
                             f"{t1:.2f} K at {v1}%)")
    if pts[0][0] < 0.0 or pts[-1][0] > 100.0:
        raise ValueError(f"{label} vol% must lie in [0, 100]")
    return pts


def _interpolator(x: list[float], y: list[float]):
    """Monotone (PCHIP) interpolant — no Runge overshoot between assay points,
    which matters because the working curve must stay monotone for the
    vol%<->T inversion to be well-defined. Inputs are clamped to the knot span
    (constant end-extrapolation), so round-tripping vol%<->T at the curve
    edges can never fall off the domain by a floating-point epsilon."""
    from scipy.interpolate import PchipInterpolator

    pchip = PchipInterpolator(x, y, extrapolate=False)
    lo, hi = x[0], x[-1]

    def f(value: float) -> float:
        return float(pchip(min(max(float(value), lo), hi)))

    return f


def _light_end_sg(cid: str) -> float:
    """Standard liquid SG of a light-end component: the GPA 2145 table first,
    else the thermo databank liquid density at 60 F."""
    if cid in LIGHT_END_SG:
        return LIGHT_END_SG[cid]
    import thermo as _thermo

    rhol = _thermo.Chemical(cid, T=288.706, P=101325.0).rhol
    if rhol is None:
        raise ValueError(
            f"no standard liquid density available for light-end component "
            f"{cid!r} (not in the GPA 2145 table and the databank has no "
            f"liquid density at 60 F); supply it via a different id or extend "
            f"caldyr.assay.LIGHT_END_SG"
        )
    return float(rhol) / _RHO_WATER_60F


def _unique_name(tb_k: float, taken: set[str]) -> str:
    base = f"NBP_{round(tb_k - 273.15)}C"
    name, n = base, 1
    while name in taken:
        n += 1
        name = f"{base}_{n}"
    taken.add(name)
    return name


def characterize_assay(
    curve: Sequence[tuple[float, float]],
    *,
    kind: str = "TBP",
    api_gravity: float | None = None,
    sg_curve: Sequence[tuple[float, float]] | None = None,
    mw_curve: Sequence[tuple[float, float]] | None = None,
    n_cuts: int = 20,
    cut_points: Sequence[float] | None = None,
    light_ends: dict[str, float] | None = None,
) -> AssayResult:
    """Characterize a petroleum assay into pseudo-components.

    Parameters
    ----------
    curve : sequence of (vol_pct, T_K)
        The distillation curve on a liquid-volume basis (at least 5 points,
        strictly increasing in both axes).
    kind : "TBP" | "ASTM_D86"
        Curve type. An ASTM D86 curve is first converted to TBP with the
        Riazi-Daubert (API) correlation (:func:`astm_d86_to_tbp`).
    api_gravity : float, optional
        Bulk API gravity (60F) of the whole crude. Used (with the constant-
        Watson-K assumption, standard practice when no density curve is
        measured — Riazi 2005 sec. 3.3.4; HYSYS does the same) to distribute
        gravity over the cuts when no ``sg_curve`` is given.
    sg_curve : sequence of (vol_pct, SG), optional
        Specific-gravity distribution on a mid-volume basis (convert an API
        curve point-wise with :func:`api_to_sg`). Takes precedence over the
        Watson-K assumption.
    mw_curve : sequence of (vol_pct, MW_g_per_mol), optional
        Measured molecular-weight distribution (mid-volume basis, the
        customary g/mol). When given, cut MWs are interpolated from it (the
        HYSYS "Molecular Wt. Curve" treatment); otherwise they are estimated
        with the Riazi-Daubert correlation (:func:`riazi_daubert_mw`).
    n_cuts : int
        Number of equal-temperature-range cuts (recommended 10-30) when
        ``cut_points`` is not given.
    cut_points : sequence of float (K), optional
        Explicit cut boundary temperatures (ascending, K) — e.g. HYSYS-style
        ranges "IBP-800F every 25F, 800-1200F every 50F, 1200-1400F every
        100F" (Hameed 2025 sec. 10.2.1). Boundaries are clipped to the curve's
        span; empty slices (< 0.01 vol%) are dropped.
    light_ends : dict id -> LV%, optional
        Defined light-end components (databank ids) with their liquid-volume
        percent of the whole assay. They are taken to occupy the front of the
        curve (the cut region starts after them), the HYSYS "Input
        Composition" treatment.

    Returns
    -------
    AssayResult
        Pseudo-components with constants + whole-assay volume/mass/mole
        fractions, light-end fractions, bulk properties and plot-ready curves.
    """
    if kind not in ("TBP", "ASTM_D86"):
        raise ValueError(f"unknown assay curve kind {kind!r}; "
                         f"expected 'TBP' or 'ASTM_D86'")
    input_pts = _validated_curve(curve, f"{kind} curve")
    tbp_pts = astm_d86_to_tbp(input_pts) if kind == "ASTM_D86" else input_pts

    vols = [p[0] for p in tbp_pts]
    temps = [p[1] for p in tbp_pts]
    t_of_v = _interpolator(vols, temps)
    # T -> vol% must be the *inverse* of the working curve (two independently
    # fitted PCHIPs are only inverses at the knots, which would make the cut
    # volumes inconsistent with the cut NBPs). Invert t_of_v numerically on a
    # fine grid; linear interpolation between grid points keeps it monotone
    # and inverse-consistent to ~1e-4 of the curve span.
    import numpy as np

    _grid_v = np.linspace(vols[0], vols[-1], 4001)
    _grid_t = np.array([t_of_v(v) for v in _grid_v])

    def v_of_t(t: float) -> float:
        return float(np.interp(t, _grid_t, _grid_v))

    light_ends = dict(light_ends or {})
    le_total = sum(light_ends.values())
    if le_total < 0.0 or le_total >= 100.0:
        raise ValueError(f"light ends sum to {le_total} LV%; expected [0, 100)")

    v_start = max(vols[0], le_total)
    v_end = vols[-1]
    if v_start >= v_end:
        raise ValueError(
            f"light ends ({le_total:.2f} LV%) cover the whole supplied curve "
            f"(vol% {vols[0]:.1f}-{v_end:.1f}); nothing left to cut"
        )
    t_start, t_end = float(t_of_v(v_start)), float(t_of_v(v_end))

    # -- cut boundaries (temperature-based slicing, as HYSYS cuts blends) ----
    if cut_points is not None:
        bounds = sorted(float(t) for t in cut_points)
        if len(bounds) < 2:
            raise ValueError("cut_points needs at least 2 boundary temperatures")
        bounds = [max(b, t_start) for b in bounds]
        bounds = [min(b, t_end) for b in bounds]
        # Material outside the explicit boundaries is merged into the first /
        # last cut (extending their ranges) rather than spawning sliver cuts —
        # this matches HYSYS's bookkeeping, where e.g. the book's 28+8+2
        # ranges yield exactly 38 hypocomponents (Hameed 2025 sec. 10.2.1
        # step 20) even though the curve starts below the starting cut point.
        bounds[0] = t_start
        bounds[-1] = t_end
    else:
        if not 2 <= n_cuts <= 60:
            raise ValueError(f"n_cuts={n_cuts} out of range [2, 60] "
                             f"(recommended 10-30)")
        bounds = [t_start + (t_end - t_start) * i / n_cuts for i in range(n_cuts + 1)]

    # -- bulk gravity setup ---------------------------------------------------
    sg_interp = None
    sg_lo = sg_hi = 0.0
    if sg_curve is not None:
        sg_pts = sorted((float(v), float(s)) for v, s in sg_curve)
        if len(sg_pts) < 2:
            raise ValueError("sg_curve needs at least 2 points")
        sg_interp = _interpolator([p[0] for p in sg_pts], [p[1] for p in sg_pts])
        sg_lo, sg_hi = sg_pts[0][0], sg_pts[-1][0]
    if sg_interp is None:
        if api_gravity is not None:
            sg_bulk = api_to_sg(api_gravity)
        else:
            # Last resort: estimate bulk SG from the curve itself (documented
            # fallback; supply api_gravity or sg_curve for real work).
            sg_bulk = sg_from_distillation(float(t_of_v(min(max(10.0, v_start), v_end))),
                                           float(t_of_v(min(50.0, v_end))), "TBP")

    # -- slice ----------------------------------------------------------------
    raw: list[dict] = []
    for t_lo, t_hi in zip(bounds, bounds[1:]):
        if t_hi - t_lo <= 1e-9:
            continue
        v_lo = float(v_of_t(min(max(t_lo, t_start), t_end)))
        v_hi = float(v_of_t(min(max(t_hi, t_start), t_end)))
        dv = v_hi - v_lo
        if dv < 0.01:                                  # empty slice (<0.01 LV%)
            continue
        nbp = float(t_of_v(0.5 * (v_lo + v_hi)))       # mid-volume NBP
        raw.append({"T_lo": t_lo, "T_hi": t_hi, "dv": dv, "v_mid": 0.5 * (v_lo + v_hi),
                    "Tb": nbp})

    if not raw:
        raise ValueError("no non-empty cuts were produced; check cut_points vs "
                         "the curve's temperature span")

    # -- per-cut gravity ------------------------------------------------------
    if sg_interp is not None:
        for c in raw:
            c["SG"] = float(sg_interp(min(max(c["v_mid"], sg_lo), sg_hi)))
    else:
        # Constant Watson K from the bulk gravity and the volume-average
        # boiling point (VABP as the MeABP proxy — the small MeABP correction
        # of Riazi 2005 sec. 3.4.1 is neglected, documented), then a single
        # multiplicative rescale so the blended gravity reproduces the bulk
        # value exactly.
        vabp = sum(c["Tb"] * c["dv"] for c in raw) / sum(c["dv"] for c in raw)
        kw_bulk = watson_k(vabp, sg_bulk)
        for c in raw:
            c["SG"] = (1.8 * c["Tb"]) ** (1.0 / 3.0) / kw_bulk
        blend = sum(c["SG"] * c["dv"] for c in raw) / sum(c["dv"] for c in raw)
        scale = sg_bulk / blend
        for c in raw:
            c["SG"] *= scale

    # -- per-cut molecular weight ----------------------------------------------
    mw_interp = None
    mw_lo = mw_hi = 0.0
    if mw_curve is not None:
        mw_pts = sorted((float(v), float(m)) for v, m in mw_curve)
        if len(mw_pts) < 2:
            raise ValueError("mw_curve needs at least 2 points")
        mw_interp = _interpolator([p[0] for p in mw_pts], [p[1] for p in mw_pts])
        mw_lo, mw_hi = mw_pts[0][0], mw_pts[-1][0]

    # -- light ends -----------------------------------------------------------
    from .core.components_db import molar_mass

    le_rows: list[tuple[str, float, float, float]] = []   # id, LV%, SG, MW
    for cid, lv in light_ends.items():
        if lv < 0.0:
            raise ValueError(f"light end {cid!r} has negative LV% {lv}")
        le_rows.append((cid, lv, _light_end_sg(cid), molar_mass(cid)))

    # -- assemble pseudo-components -------------------------------------------
    # Whole-assay basis: volumes normalized over light ends + cuts (a curve
    # ending below 100% — e.g. the book's 98% — simply characterizes what it
    # reports). mass ∝ LV * SG; moles ∝ mass / MW.
    total_v = le_total + sum(c["dv"] for c in raw)
    le_mass = [(lv / total_v) * sg for _, lv, sg, _ in le_rows]
    cut_sg_mw: list[tuple[float, float]] = []
    for c in raw:
        if mw_interp is not None:
            mw = float(mw_interp(min(max(c["v_mid"], mw_lo), mw_hi))) / 1000.0
        else:
            mw = riazi_daubert_mw(c["Tb"], c["SG"])
        cut_sg_mw.append((c["SG"], mw))
    cut_mass = [(c["dv"] / total_v) * c["SG"] for c in raw]
    total_mass = sum(le_mass) + sum(cut_mass)
    le_mol = [m / mw for m, (_, _, _, mw) in zip(le_mass, le_rows)]
    cut_mol = [m / mw for m, (_, mw) in zip(cut_mass, cut_sg_mw)]
    total_mol = sum(le_mol) + sum(cut_mol)

    taken: set[str] = set()
    cuts: list[PseudoComponent] = []
    for c, (sg, mw), m, n in zip(raw, cut_sg_mw, cut_mass, cut_mol):
        tb = c["Tb"]
        kw = watson_k(tb, sg)
        tc = kesler_lee_tc(tb, sg)
        pc = kesler_lee_pc(tb, sg)
        cuts.append(PseudoComponent(
            id=_unique_name(tb, taken),
            Tb=tb, SG=sg, MW=mw, Tc=tc, Pc=pc,
            omega=acentric_factor(tb, tc, pc, kw),
            watson_k=kw,
            cp_ig=kesler_lee_cp_ig(mw, kw),
            vol_frac=c["dv"] / total_v,
            mass_frac=m / total_mass,
            mole_frac=n / total_mol,
            T_lo=c["T_lo"], T_hi=c["T_hi"],
            vol_mid_pct=c["v_mid"],
        ))

    le_fracs = {
        cid: {"vol_frac": lv / total_v, "mass_frac": m / total_mass,
              "mole_frac": n / total_mol}
        for (cid, lv, _, _), m, n in zip(le_rows, le_mass, le_mol)
    }

    # -- bulk + plot curves -----------------------------------------------------
    # mass_i carries units of (vol frac)*(SG) and mol_i = mass_i / MW[kg/mol],
    # so total_mass/total_mol is the blend molar mass in kg/mol directly; and
    # because the volume fractions sum to 1, total_mass IS the volume-weighted
    # blend specific gravity.
    mw_bulk = total_mass / total_mol
    sg_blend = total_mass
    vabp_all = sum(c.Tb * c.vol_frac for c in cuts) / sum(c.vol_frac for c in cuts)
    bulk = {
        "MW": mw_bulk,
        "SG": sg_blend,
        "API": sg_to_api(sg_blend),
        "watson_k": watson_k(vabp_all, sg_blend),
    }

    grid_n = 101
    grid_v = [v_start + (v_end - v_start) * i / (grid_n - 1) for i in range(grid_n)]
    curves = {
        "vol_pct": grid_v,
        "tbp_K": [float(t_of_v(v)) for v in grid_v],
        "input_vol_pct": [p[0] for p in input_pts],
        "input_T_K": [p[1] for p in input_pts],
    }

    return AssayResult(cuts=cuts, light_end_fractions=le_fracs, bulk=bulk,
                       curves=curves, kind=kind)
