"""Reactive acid-gas property package for amine gas sweetening.

Wraps the modified Kent-Eisenberg acid-gas solubility model
(:mod:`caldyr.thermo.amine`) behind the :class:`~caldyr.thermo.base.PropertyPackage`
protocol, so a :class:`~caldyr.unitops.absorber.Absorber` (or stripper /
regenerator) can sweeten a sour gas with aqueous DEA/MDEA. This is the open
analogue of HYSYS's "Acid Gas - Chemical Solvents" package (Hameed 2025,
*Chemical Process Simulations using Aspen HYSYS*, §15.3).

**Why a dedicated package.** A cubic EOS or activity model cannot represent the
*chemical* uptake of CO2/H2S by an amine — most of the absorbed acid gas sits in
the liquid as ionic species (carbamate, bicarbonate, bisulfide, protonated
amine) with no vapour pressure of its own. The signature is a K-value that
*collapses at low loading* (very favourable absorption) and rises steeply as the
solvent saturates. Here the acid-gas K-values come from the reactive model:

    K_i = p_i(T, loading_i, M) / (x_i * P)            for CO2, H2S

where ``p_i`` is the equilibrium partial pressure the amine model predicts at the
stage's *apparent* liquid loading (loading_i = mol acid gas / mol amine, read off
the apparent liquid mole fractions) and bulk amine molarity ``M`` (computed from
the apparent composition and a solution molar volume). CO2 and H2S are inverted
*jointly* (:func:`caldyr.thermo.amine.acid_gas_partial_pressures`) so they
compete for the amine through the shared charge balance.

The remaining species are conventional:
  * **water** — Raoult's law, ``K = Psat(T) / P`` (Antoine; activity ~ 1, the
    bulk solvent);
  * **amine** — essentially non-volatile, ``K = Psat_amine(T) / P`` from a
    Clausius-Clapeyron vapour pressure (tiny);
  * **light gases** (CH4, C2H6, C3H8, N2, CO, H2 ...) — *physically* dissolved
    only, Henry's law ``K = H_i(T) / P`` (large, sparingly soluble; van't Hoff
    Henry constants from the Sander 2015 compilation).

**Enthalpy** is on a per-component reference basis (each pure ideal gas at
298.15 K = 0): vapour is sensible only; the liquid carries a per-component
reference offset — heat of vaporization for the condensables (water, amine) and
the **heat of absorption** for the reactive acid gases (CO2 ≈ -55..-70, H2S ≈
-38 kJ/mol). Transferring an acid gas from vapour to the absorbed liquid state
therefore releases ``-ΔH_abs`` — which is exactly the absorber temperature
bulge. Per-component references cancel in conserved-composition balances, so the
energy balance closes while the heat of absorption is carried correctly.

This is a steady-state VLE + energy package: it implements the column contract
(``k_values``, ``enthalpy_liquid/vapor``, ``enthalpy``, the PT/PH flashes,
bubble point). Entropy-based operations (``entropy``, ``flash_ps``) and
three-phase (VLLE) flashes raise :class:`NotImplementedError` with a clear
message — out of scope for a gas-sweetening column.

Validation: ``tests/test_m16_amine_package.py`` (reactive K-value collapse at
low loading, energy-balanced flashes, and a §15.3-style DEA absorber that
removes the acid gas, closes the mass balance, and lands a rich loading near the
solubility model's value).
"""
from __future__ import annotations

import math
from functools import lru_cache

from scipy.optimize import brentq

from . import amine as ke
from .base import PhaseResult, ThreePhaseResult

_TREF = 298.15            # K, enthalpy reference temperature
_R = 8.314462618          # J/mol/K
_X_FLOOR = 1e-15          # floor for apparent mole fractions
_KPA_PER_PA = 1.0e-3      # kPa per Pa (amine model speaks kPa)
_PA_PER_KPA = 1.0e3
_T_LO, _T_HI = 200.0, 600.0   # K, flash/bubble temperature search window

# CAS numbers used to classify a flowsheet component into its physical role.
_CAS_WATER = "7732-18-5"
_CAS_CO2 = "124-38-9"
_CAS_H2S = "7783-06-4"
_CAS_AMINE = {"DEA": "111-42-2", "MDEA": "105-59-9"}


# -- pure-component property data ---------------------------------------------
# Aqueous-solvent species: molar mass (kg/mol), liquid molar volume (m^3/mol,
# from density at ~298 K), liquid heat capacity (J/mol/K), heat of vaporization
# (J/mol), ideal-gas heat capacity (J/mol/K), and normal boiling point (K, for
# the Clausius-Clapeyron amine vapour pressure). Values from standard references
# (CRC / DIPPR / Yaws); the package is validated by the absorber behaviour, not
# by these constants to high precision.
_WATER = dict(MW=18.015e-3, Vm=18.07e-6, CpL=75.3, Hvap=44.0e3, Cpig=33.6, Tb=373.15)
_AMINE_DATA = {
    "DEA": dict(MW=105.136e-3, Vm=95.8e-6, CpL=290.0, Hvap=60.9e3, Cpig=200.0, Tb=541.5),
    "MDEA": dict(MW=119.163e-3, Vm=114.8e-6, CpL=350.0, Hvap=62.0e3, Cpig=230.0, Tb=520.2),
}

# Reactive acid gases: ideal-gas Cp (J/mol/K) and heat of absorption (J/mol,
# negative, exothermic). CO2 heat of absorption is amine-specific; H2S's is the
# weaker (instantaneous protonation) reaction, taken constant. Representative
# literature values (Kohl & Nielsen, *Gas Purification* 5e; the absorber bulge
# depends on them but the acid-gas removal does not).
_ACID_CPIG = {"CO2": 37.1, "H2S": 34.0}
_ACID_HABS = {
    "CO2": {"DEA": -70.0e3, "MDEA": -54.0e3},
    "H2S": {"DEA": -38.0e3, "MDEA": -38.0e3},
}

# Light gases — physical (Henry's law) solubility in water. Per species:
# Sander (2015) intrinsic Henry solubility H^cp (mol/(m^3.Pa)) at 298.15 K and
# its temperature coefficient d(ln H^cp)/d(1/T) (K). The mole-fraction Henry
# constant is K_x = 1 / (Vm_water * H^cp) [Pa]; solubility falls with T, so K_x
# rises. CAS-keyed so any recognised id maps in.
_HENRY: dict[str, tuple[float, float]] = {
    "74-82-8": (1.4e-5, 1900.0),    # methane
    "74-84-0": (1.9e-5, 2400.0),    # ethane
    "74-98-6": (1.5e-5, 2700.0),    # propane
    "7727-37-9": (6.4e-6, 1600.0),  # nitrogen
    "630-08-0": (9.7e-6, 1300.0),   # carbon monoxide
    "1333-74-0": (7.8e-6, 500.0),   # hydrogen
}
_HENRY_DEFAULT = (5.0e-6, 1500.0)   # unknown light gas: sparingly soluble
# Physical (non-reactive) Henry solubility of the acid gases themselves —
# the fallback when no aqueous amine is present (e.g. the dry sour-gas feed
# flash). CO2/H2S are far more soluble than the light gases, but with no amine
# there is no chemical uptake. Sander (2015).
_ACID_HENRY = {"CO2": (3.3e-4, 2400.0), "H2S": (1.0e-3, 2100.0)}


def _antoine_water_pa(T: float) -> float:
    """Saturation pressure of water (Pa) from the NIST Antoine equation
    (Stull 1947; mmHg, deg C), with the >99 deg C coefficient set above 372 K."""
    tc = T - 273.15
    if T <= 372.15:
        a, b, c = 8.07131, 1730.63, 233.426
    else:
        a, b, c = 8.14019, 1810.94, 244.485
    return 10.0 ** (a - b / (tc + c)) * 133.322


def _classify(cas: str | None, cid: str, amine_name: str) -> str:
    """Map a component to its role: 'amine', 'water', 'CO2', 'H2S' or 'gas'."""
    if cas == _CAS_WATER:
        return "water"
    if cas == _CAS_CO2:
        return "CO2"
    if cas == _CAS_H2S:
        return "H2S"
    if cas == _CAS_AMINE.get(amine_name) or cid.upper() == amine_name:
        return "amine"
    return "gas"


@lru_cache(maxsize=64)
def _roles(components: tuple[str, ...], amine_name: str) -> dict[str, str]:
    """Per-component role map (cached — the chemicals lookup is the slow part)."""
    from ..core.components_db import resolve_component

    roles: dict[str, str] = {}
    for cid in components:
        cas = None
        try:
            cas = resolve_component(cid).cas
        except Exception:
            cas = None
        roles[cid] = _classify(cas, cid, amine_name)
    return roles


@lru_cache(maxsize=64)
def _henry_params(components: tuple[str, ...]) -> dict[str, tuple[float, float]]:
    """Henry (H^cp, dlnH/d(1/T)) per light-gas component, keyed by CAS."""
    from ..core.components_db import resolve_component

    out: dict[str, tuple[float, float]] = {}
    for cid in components:
        try:
            cas = resolve_component(cid).cas
        except Exception:
            cas = None
        out[cid] = _HENRY.get(cas or "", _HENRY_DEFAULT)
    return out


class AmineAcidGasPackage:
    """Reactive acid-gas (Kent-Eisenberg) property package for amine sweetening.

    Built for a fixed, ordered component list that must contain the selected
    amine and water plus at least one acid gas (CO2 and/or H2S); any further
    components are treated as physically-dissolved light gases. See the module
    docstring for the model.
    """

    SUPPORTED = ("DEA", "MDEA")

    def __init__(self, components: list[str], amine: str = "DEA") -> None:
        amine = amine.upper()
        if amine not in self.SUPPORTED:
            raise ValueError(
                f"unsupported amine {amine!r}; the acid-gas package supports "
                f"{self.SUPPORTED}"
            )
        if not components:
            raise ValueError("AmineAcidGasPackage requires at least one component")
        self.amine = amine
        self.components = list(components)
        self._sys = ke.amine_system(amine)
        self._roles = _roles(tuple(components), amine)
        self._henry = _henry_params(tuple(components))
        amine_id = self._role_id("amine", required=True)
        water_id = self._role_id("water", required=True)
        assert amine_id is not None and water_id is not None  # required=True
        self._amine_id: str = amine_id
        self._water_id: str = water_id
        self._co2_id = self._role_id("CO2")
        self._h2s_id = self._role_id("H2S")
        if self._co2_id is None and self._h2s_id is None:
            raise ValueError(
                "AmineAcidGasPackage needs at least one acid gas (CO2 or H2S) "
                f"in the component list; got {components}"
            )
        self._amine_data = _AMINE_DATA[amine]

    def _role_id(self, role: str, required: bool = False) -> str | None:
        ids = [c for c, r in self._roles.items() if r == role]
        if len(ids) > 1:
            raise ValueError(
                f"AmineAcidGasPackage: more than one {role!r} component "
                f"({ids}); exactly one is expected"
            )
        if not ids:
            if required:
                raise ValueError(
                    f"AmineAcidGasPackage: no {role!r} component found in "
                    f"{self.components}; the amine package needs aqueous "
                    f"{self.amine} (amine + water) present"
                )
            return None
        return ids[0]

    @classmethod
    def from_spec(cls, spec: str, components: list[str]) -> "AmineAcidGasPackage":
        """Build from a flowsheet ``property_package`` string like
        ``"amine:DEA"``."""
        backend, _, amine = spec.partition(":")
        if backend != "amine":
            raise ValueError(
                f"AmineAcidGasPackage cannot build backend {backend!r} "
                f"(got {spec!r})"
            )
        return cls(components, amine or "DEA")

    # -- composition helpers ----------------------------------------------
    def _zs(self, z: dict[str, float]) -> dict[str, float]:
        unknown = set(z) - set(self.components)
        if unknown:
            raise ValueError(
                f"composition has components not in package: {sorted(unknown)}")
        raw = {c: max(float(z.get(c, 0.0)), 0.0) for c in self.components}
        total = sum(raw.values())
        if total <= 0.0:
            raise ValueError(f"composition sums to {total}; expected > 0")
        return {c: v / total for c, v in raw.items()}

    def _amine_molarity(self, x: dict[str, float]) -> float:
        """Bulk amine molarity (mol/L) from the apparent liquid composition,
        using an ideal-mixing solvent molar volume (water + amine)."""
        xa = max(x[self._amine_id], _X_FLOOR)
        xw = max(x[self._water_id], 0.0)
        vm = xw * _WATER["Vm"] + xa * self._amine_data["Vm"]   # m^3 per mol soln
        if vm <= 0.0:
            return 0.0
        return xa / (vm * 1000.0)        # mol amine / L solution

    # -- K-values (the reactive core) -------------------------------------
    def k_values(self, T: float, P: float, x: dict[str, float],
                 y: dict[str, float]) -> dict[str, float]:
        """Phase-equilibrium K = y/x per component: reactive (amine model) for
        CO2/H2S, Raoult for water/amine, Henry for light gases. The vapour
        composition ``y`` is not used (the vapour is treated as ideal)."""
        xs = {c: max(v, _X_FLOOR) for c, v in self._zs(x).items()}
        xa = xs[self._amine_id]
        M = self._amine_molarity(xs)
        K: dict[str, float] = {}

        # Reactive acid gases, inverted jointly so they compete for the amine.
        co2_id, h2s_id = self._co2_id, self._h2s_id
        a_co2 = xs[co2_id] / xa if co2_id else 0.0
        a_h2s = xs[h2s_id] / xa if h2s_id else 0.0
        reactive = M > 1e-6 and xs[self._water_id] > 1e-6
        if reactive:
            try:
                p_co2, p_h2s = ke.acid_gas_partial_pressures(
                    T, a_co2, a_h2s, M, self._sys)
            except (ValueError, RuntimeError, OverflowError):
                # The Kent-Eisenberg charge balance is fitted to ~290-400 K; far
                # outside it (e.g. the hot end of a flash_ph temperature search)
                # the ionic equilibrium can fail to bracket. There the solution
                # is well into the vapour region anyway, so fall back to the
                # acid gases' physical (Henry) volatility.
                reactive = False
                p_co2 = p_h2s = 0.0
        else:
            p_co2 = p_h2s = 0.0
        for cid, p, alpha, role in ((co2_id, p_co2, a_co2, "CO2"),
                                    (h2s_id, p_h2s, a_h2s, "H2S")):
            if cid is None:
                continue
            if reactive and alpha > 0.0:
                K[cid] = p * _PA_PER_KPA / (xs[cid] * P)
            else:
                # No aqueous amine present (e.g. the dry gas feed flash): the
                # acid gas is a normal physical gas -> Henry's law.
                K[cid] = self._henry_k(cid, T, P, params=_ACID_HENRY[role])

        # Condensable solvents: Raoult's law.
        K[self._water_id] = _antoine_water_pa(T) / P
        K[self._amine_id] = self._amine_psat(T) / P

        # Light gases: Henry's law.
        for cid, role in self._roles.items():
            if role == "gas":
                K[cid] = self._henry_k(cid, T, P)
        return K

    def _henry_k(self, cid: str, T: float, P: float,
                 params: tuple[float, float] | None = None) -> float:
        hcp298, slope = params or self._henry.get(cid, _HENRY_DEFAULT)
        hcp = hcp298 * math.exp(slope * (1.0 / T - 1.0 / _TREF))
        kx = 1.0 / (_WATER["Vm"] * hcp)        # mole-fraction Henry constant, Pa
        return kx / P

    def _amine_psat(self, T: float) -> float:
        """Amine vapour pressure (Pa) by Clausius-Clapeyron from its normal
        boiling point and heat of vaporization — small, so the amine is
        effectively non-volatile."""
        d = self._amine_data
        return 101325.0 * math.exp(-d["Hvap"] / _R * (1.0 / T - 1.0 / d["Tb"]))

    # -- enthalpy ---------------------------------------------------------
    def _hL_ref(self, cid: str) -> float:
        """Pure-liquid enthalpy of component ``cid`` at the reference T, relative
        to its ideal gas (J/mol)."""
        role = self._roles[cid]
        if role == "water":
            return -_WATER["Hvap"]
        if role == "amine":
            return -self._amine_data["Hvap"]
        if role in ("CO2", "H2S"):
            return _ACID_HABS[role][self.amine]
        return 0.0       # dissolved light gas ~ ideal-gas enthalpy

    def _cpL(self, cid: str) -> float:
        role = self._roles[cid]
        if role == "water":
            return _WATER["CpL"]
        if role == "amine":
            return self._amine_data["CpL"]
        if role in ("CO2", "H2S"):
            return _ACID_CPIG[role]
        return self._cpig(cid)

    def _cpig(self, cid: str) -> float:
        role = self._roles[cid]
        if role == "water":
            return _WATER["Cpig"]
        if role == "amine":
            return self._amine_data["Cpig"]
        if role in ("CO2", "H2S"):
            return _ACID_CPIG[role]
        return _cpig_gas(cid)

    def enthalpy_liquid(self, T: float, P: float, x: dict[str, float]) -> float:
        """Molar enthalpy (J/mol) of ``x`` as a liquid: sensible + per-component
        liquid reference (heat of vaporization for solvents, heat of absorption
        for the reactive acid gases)."""
        xs = self._zs(x)
        return sum(v * (self._hL_ref(c) + self._cpL(c) * (T - _TREF))
                   for c, v in xs.items())

    def enthalpy_vapor(self, T: float, P: float, y: dict[str, float]) -> float:
        """Molar enthalpy (J/mol) of ``y`` as an ideal-gas vapour: sensible
        only (pure ideal gas at the reference T is the zero)."""
        ys = self._zs(y)
        return sum(v * self._cpig(c) * (T - _TREF) for c, v in ys.items())

    def enthalpy(self, T: float, P: float, z: dict[str, float]) -> float:
        """Bulk molar enthalpy at (T, P): a PT flash, phase-weighted."""
        return self.flash_pt(T, P, z).H

    # -- flashes ----------------------------------------------------------
    def _flash_beta(self, T: float, P: float, zs: dict[str, float]
                    ) -> tuple[float, dict[str, float], dict[str, float]]:
        """Isothermal flash: vapour fraction beta and phase compositions, with
        the reactive K-values re-evaluated as the liquid composition settles
        (outer fixed-point because K depends on x)."""
        comps = self.components
        x = dict(zs)
        y = dict(zs)
        beta = 0.5
        for _ in range(60):
            K = self.k_values(T, P, x, y)
            g0 = sum(zs[c] * (K[c] - 1.0) for c in comps)
            g1 = sum(zs[c] * (K[c] - 1.0) / K[c] for c in comps)
            if g0 <= 0.0:                       # below bubble: all liquid
                beta = 0.0
            elif g1 >= 0.0:                     # above dew: all vapour
                beta = 1.0
            else:
                beta = brentq(
                    lambda b: sum(zs[c] * (K[c] - 1.0) / (1.0 + b * (K[c] - 1.0))
                                  for c in comps), 0.0, 1.0, xtol=1e-12)
            x_new = {c: zs[c] / (1.0 + beta * (K[c] - 1.0)) for c in comps}
            s = sum(x_new.values()) or 1.0
            x_new = {c: v / s for c, v in x_new.items()}
            y_new = {c: K[c] * x_new[c] for c in comps}
            sy = sum(y_new.values()) or 1.0
            y_new = {c: v / sy for c, v in y_new.items()}
            if max(abs(x_new[c] - x[c]) for c in comps) < 1e-12:
                x, y = x_new, y_new
                break
            x, y = x_new, y_new
        return beta, x, y

    def flash_pt(self, T: float, P: float, z: dict[str, float]) -> PhaseResult:
        zs = self._zs(z)
        beta, x, y = self._flash_beta(T, P, zs)
        hL: float | None = None
        hV: float | None = None
        if beta <= 0.0:
            phase = "liquid"
            H = hL = self.enthalpy_liquid(T, P, zs)
        elif beta >= 1.0:
            phase = "vapor"
            H = hV = self.enthalpy_vapor(T, P, zs)
        else:
            phase = "VLE"
            hL = self.enthalpy_liquid(T, P, x)
            hV = self.enthalpy_vapor(T, P, y)
            H = beta * hV + (1.0 - beta) * hL
        return PhaseResult(
            T=float(T), P=float(P), H=float(H), phase=phase,
            vapor_fraction=float(beta),
            x=dict(x), y=dict(y), H_liquid=hL, H_vapor=hV)

    def flash_ph(self, P: float, H: float, z: dict[str, float]) -> PhaseResult:
        """Resolve the temperature whose PT-flash enthalpy equals ``H`` (Brent),
        then return that flash. Used to recover the rich-solvent product T from
        the column's exact energy balance."""
        zs = self._zs(z)

        def resid(T: float) -> float:
            return self.flash_pt(T, P, zs).H - H

        lo, hi = _T_LO, _T_HI
        f_lo, f_hi = resid(lo), resid(hi)
        if f_lo > 0.0:                 # colder than the search floor
            return self.flash_pt(lo, P, zs)
        if f_hi < 0.0:                 # hotter than the search ceiling
            return self.flash_pt(hi, P, zs)
        T = float(brentq(resid, lo, hi, xtol=1e-6))
        return self.flash_pt(T, P, zs)

    def flash_ps(self, P: float, S: float, z: dict[str, float]) -> PhaseResult:
        raise NotImplementedError(
            "AmineAcidGasPackage does not implement entropy-based (PS) flashes; "
            "it is a VLE + energy package for amine absorbers/regenerators")

    # -- bubble / dew -----------------------------------------------------
    def bubble_point(self, P: float, x: dict[str, float]) -> PhaseResult:
        """Saturated-liquid state of ``x`` at ``P``: the temperature where
        sum_i K_i x_i = 1, with the incipient vapour y = K x."""
        xs = self._zs(x)

        def resid(T: float) -> float:
            K = self.k_values(T, P, xs, xs)
            return sum(K[c] * xs[c] for c in self.components) - 1.0

        T = self._bracket_solve(resid)
        K = self.k_values(T, P, xs, xs)
        y = {c: K[c] * xs[c] for c in self.components}
        sy = sum(y.values()) or 1.0
        y = {c: v / sy for c, v in y.items()}
        return PhaseResult(
            T=T, P=float(P), H=self.enthalpy_liquid(T, P, xs),
            phase="liquid", vapor_fraction=0.0, x=dict(xs), y=y,
            H_liquid=self.enthalpy_liquid(T, P, xs),
            H_vapor=self.enthalpy_vapor(T, P, y))

    def bubble_dew(self, P: float, z: dict[str, float]) -> tuple[float, float]:
        zs = self._zs(z)

        def bub(T: float) -> float:
            K = self.k_values(T, P, zs, zs)
            return sum(K[c] * zs[c] for c in self.components) - 1.0

        def dew(T: float) -> float:
            K = self.k_values(T, P, zs, zs)
            return sum(zs[c] / K[c] for c in self.components) - 1.0

        t_bub = self._bracket_solve(bub)
        try:
            t_dew = self._bracket_solve(dew)
        except ValueError:
            t_dew = _T_HI
        return min(t_bub, t_dew), max(t_bub, t_dew)

    def _bracket_solve(self, resid) -> float:
        lo, hi = _T_LO, _T_HI
        f_lo, f_hi = resid(lo), resid(hi)
        if f_lo == 0.0:
            return lo
        if f_hi == 0.0:
            return hi
        if f_lo * f_hi > 0.0:
            # No sign change across the window: clamp to the nearer edge so an
            # all-volatile / all-condensable mixture degrades gracefully.
            return lo if abs(f_lo) < abs(f_hi) else hi
        return float(brentq(resid, lo, hi, xtol=1e-6))

    # -- entropy / volume -------------------------------------------------
    def entropy(self, T: float, P: float, z: dict[str, float]) -> float:
        raise NotImplementedError(
            "AmineAcidGasPackage does not implement entropy; it is a VLE + "
            "energy package for amine absorbers/regenerators")

    def volume(self, T: float, P: float, z: dict[str, float]) -> float:
        res = self.flash_pt(T, P, z)
        beta = res.vapor_fraction
        vL = self.volume_liquid(T, P, res.x or z) if beta < 1.0 else 0.0
        vV = self.volume_vapor(T, P, res.y or z) if beta > 0.0 else 0.0
        return beta * vV + (1.0 - beta) * vL

    def volume_liquid(self, T: float, P: float, x: dict[str, float]) -> float:
        """Liquid molar volume (m^3/mol): ideal mixing of solvent volumes
        (water, amine); dissolved gases contribute a water-like partial volume."""
        xs = self._zs(x)
        vm = 0.0
        for c, v in xs.items():
            role = self._roles[c]
            if role == "water":
                vm += v * _WATER["Vm"]
            elif role == "amine":
                vm += v * self._amine_data["Vm"]
            else:
                vm += v * _WATER["Vm"]      # dissolved species ~ water volume
        return vm

    def volume_vapor(self, T: float, P: float, y: dict[str, float]) -> float:
        """Ideal-gas molar volume, m^3/mol."""
        return _R * T / P

    # -- three-phase (not supported) --------------------------------------
    def flash_pt_3p(self, T: float, P: float, z: dict[str, float]) -> ThreePhaseResult:
        raise NotImplementedError(
            "AmineAcidGasPackage does not support three-phase (VLLE) flashes")

    def flash_ph_3p(self, P: float, H: float, z: dict[str, float]) -> ThreePhaseResult:
        raise NotImplementedError(
            "AmineAcidGasPackage does not support three-phase (VLLE) flashes")


@lru_cache(maxsize=256)
def _cpig_gas(cid: str) -> float:
    """Ideal-gas Cp (J/mol/K, near ambient) of a light gas, from a small table
    keyed by CAS with a sane default."""
    table = {
        "74-82-8": 35.7,    # methane
        "74-84-0": 52.5,    # ethane
        "74-98-6": 73.6,    # propane
        "7727-37-9": 29.1,  # nitrogen
        "630-08-0": 29.1,   # carbon monoxide
        "1333-74-0": 28.8,  # hydrogen
    }
    try:
        from ..core.components_db import resolve_component
        cas = resolve_component(cid).cas
    except Exception:
        cas = None
    return table.get(cas or "", 33.0)
