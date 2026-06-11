"""Steam-tables property package: pure water on CoolProp (IAPWS-95).

Selected as ``"coolprop:Water"``. This is the HYSYS "NBS Steam" / "ASME Steam"
equivalent (cf. Hameed, *Chemical Process Simulations using Aspen HYSYS*, Wiley
2025, §2.2 "Steam Table"): a reference-quality pure-water equation of state
(CoolProp's default for ``Water`` is the IAPWS-95 Helmholtz formulation, the
basis of the modern steam tables), far more accurate for water/steam than a
cubic EOS. **Single-component water only** — anything else raises a typed
error pointing at the thermo backends.

Enthalpy basis
--------------
The engine's enthalpy basis is *formation-inclusive* (see
:mod:`caldyr.thermo._flasher`): the `thermo` packages report
``H = H_sensible(T, P; ref ideal gas at 298.15 K) + Hf_ig(298.15 K)``.
CoolProp's water enthalpy uses the IAPWS reference (zero internal energy and
entropy for saturated liquid at the triple point), so we shift it by a single
constant chosen to make the *ideal-gas state at 298.15 K* equal to water's
ideal-gas formation enthalpy — exactly the anchor the thermo packages use:

    offset = Hf_ig(water) - Hmolar_CoolProp(298.15 K, P -> 0)

``Hf_ig`` comes from the same `chemicals` database `thermo` reads, so streams
flashed by this package and by ``thermo:PR`` sit on the same absolute basis
(they differ only by the two models' physical accuracy — PR's water enthalpy
departures are a few percent off the IAPWS values).

Entropy keeps CoolProp's native reference: it is only used internally
(``flash_ps`` round-trips against this package's own ``entropy``), so no
cross-package offset is needed.
"""
from __future__ import annotations

from functools import lru_cache

from .base import PhaseResult, ThreePhaseResult

_WATER_CAS = "7732-18-5"
_T_REF = 298.15        # K, formation-enthalpy anchor (matches _flasher.py)
_P_IDEAL = 1.0         # Pa — low enough that water is ideal-gas to < 0.01 J/mol
_VF_EPS = 1e-9


@lru_cache(maxsize=1)
def _water_constants() -> tuple[float, float, float, float, float, float]:
    """(MW kg/mol, Hf_ig J/mol, enthalpy offset J/mol, Tc K, Pc Pa,
    rho_c mol/m^3).

    Cached once: the offset anchors CoolProp's IAPWS reference state to the
    engine's formation-inclusive basis (see module docstring). ``Hf_ig`` is
    read from `chemicals` — the same source `thermo`'s constants use — so the
    two backends share the basis exactly.
    """
    from chemicals.reaction import Hfg
    from CoolProp.CoolProp import PropsSI

    hf = float(Hfg(_WATER_CAS))                      # -241,822 J/mol
    h_ig_ref = float(PropsSI("Hmolar", "T", _T_REF, "P", _P_IDEAL, "Water"))
    return (
        float(PropsSI("M", "Water")),                # kg/mol
        hf,
        hf - h_ig_ref,                               # J/mol offset
        float(PropsSI("Tcrit", "Water")),
        float(PropsSI("pcrit", "Water")),
        float(PropsSI("rhomolar_critical", "Water")),
    )


class CoolPropWaterPackage:
    """Pure-water :class:`~caldyr.thermo.base.PropertyPackage` on CoolProp
    (IAPWS-95 — steam-table quality). Built by ``make_package("coolprop:Water",
    ["water"])``."""

    def __init__(self, components: list[str]) -> None:
        from ..core.components_db import resolve_component

        if len(components) != 1:
            raise ValueError(
                f"coolprop:Water is a pure-water steam-tables package and "
                f"supports exactly one component; got {len(components)}: "
                f"{list(components)}. Use a thermo:* package for mixtures."
            )
        comp = resolve_component(components[0])
        if comp.cas != _WATER_CAS:
            raise ValueError(
                f"coolprop:Water supports only water (CAS {_WATER_CAS}); component "
                f"{components[0]!r} resolved to {comp.name!r} (CAS {comp.cas}). "
                f"Use a thermo:* package for other species."
            )
        self.components = list(components)
        self._cid = components[0]
        (self.mw, self._hf, self._h_offset,
         self._t_crit, self._p_crit, self._rho_crit) = _water_constants()

    # -- internals -----------------------------------------------------------
    def _props(self, output: str, n1: str, v1: float, n2: str, v2: float) -> float:
        from CoolProp.CoolProp import PropsSI

        return float(PropsSI(output, n1, v1, n2, v2, "Water"))

    def _check_z(self, z: dict[str, float]) -> dict[str, float]:
        unknown = set(z) - {self._cid}
        if unknown:
            raise ValueError(
                f"composition has components not in package: {sorted(unknown)} "
                f"(coolprop:Water is pure water)"
            )
        if sum(z.values()) <= 0.0:
            raise ValueError(f"composition sums to {sum(z.values())}; expected > 0")
        return {self._cid: 1.0}

    def _phase_label(self, T: float, P: float) -> str:
        """"vapor" | "liquid" for a single-phase (T, P) point. Supercritical
        states are classified by density relative to the critical density (the
        same density-based convention the cubic packages use)."""
        from CoolProp.CoolProp import PhaseSI

        ph = PhaseSI("T", T, "P", P, "Water")
        if ph in ("liquid", "supercritical_liquid"):
            return "liquid"
        if ph in ("gas", "supercritical_gas"):
            return "vapor"
        rho = self._props("Dmolar", "T", T, "P", P)
        return "liquid" if rho > self._rho_crit else "vapor"

    def _single_phase_result(self, T: float, P: float, H: float) -> PhaseResult:
        comp = {self._cid: 1.0}
        phase = self._phase_label(T, P)
        if phase == "vapor":
            return PhaseResult(T=T, P=P, H=H, phase="vapor", vapor_fraction=1.0,
                               x=None, y=comp, H_liquid=None, H_vapor=H)
        return PhaseResult(T=T, P=P, H=H, phase="liquid", vapor_fraction=0.0,
                           x=comp, y=None, H_liquid=H, H_vapor=None)

    def _flash_px(self, P: float, name: str, value: float) -> PhaseResult:
        """Flash at (P, Hmolar) or (P, Smolar): the workhorse for PH/PS. A pure
        fluid at fixed P is two-phase only on the saturation line, where
        CoolProp reports the quality Q directly."""
        q = self._props("Q", "P", P, name, value)
        T = self._props("T", "P", P, name, value)
        h = self._props("Hmolar", "P", P, name, value) + self._h_offset
        if 0.0 <= q <= 1.0:                       # on the saturation line
            comp = {self._cid: 1.0}
            hl = self._props("Hmolar", "P", P, "Q", 0.0) + self._h_offset
            hv = self._props("Hmolar", "P", P, "Q", 1.0) + self._h_offset
            if q <= _VF_EPS:
                phase = "liquid"
            elif q >= 1.0 - _VF_EPS:
                phase = "vapor"
            else:
                phase = "VLE"
            return PhaseResult(T=T, P=P, H=h, phase=phase, vapor_fraction=q,
                               x=comp, y=comp, H_liquid=hl, H_vapor=hv)
        return self._single_phase_result(T, P, h)

    # -- PropertyPackage protocol ---------------------------------------------
    def enthalpy(self, T: float, P: float, z: dict[str, float]) -> float:
        """Molar enthalpy, J/mol, on the engine's formation-inclusive basis."""
        self._check_z(z)
        return self._props("Hmolar", "T", T, "P", P) + self._h_offset

    def entropy(self, T: float, P: float, z: dict[str, float]) -> float:
        """Molar entropy, J/mol/K (CoolProp's IAPWS reference; internally
        consistent with this package's flash_ps)."""
        self._check_z(z)
        return self._props("Smolar", "T", T, "P", P)

    def volume(self, T: float, P: float, z: dict[str, float]) -> float:
        """Molar volume, m^3/mol."""
        self._check_z(z)
        return 1.0 / self._props("Dmolar", "T", T, "P", P)

    def flash_pt(self, T: float, P: float, z: dict[str, float]) -> PhaseResult:
        self._check_z(z)
        h = self._props("Hmolar", "T", T, "P", P) + self._h_offset
        return self._single_phase_result(T, P, h)

    def flash_ph(self, P: float, H: float, z: dict[str, float]) -> PhaseResult:
        # H arrives on the absolute (formation-inclusive) basis; shift back to
        # CoolProp's reference before flashing — mirror of _flasher.flash_ph.
        self._check_z(z)
        return self._flash_px(P, "Hmolar", H - self._h_offset)

    def flash_ps(self, P: float, S: float, z: dict[str, float]) -> PhaseResult:
        self._check_z(z)
        return self._flash_px(P, "Smolar", S)

    def bubble_dew(self, P: float, z: dict[str, float]) -> tuple[float, float]:
        """For a pure fluid the bubble and dew temperatures coincide at the
        saturation temperature Tsat(P) (Hameed 2025, §2.3: "for a single
        component, bubble point = dew point")."""
        self._check_z(z)
        if P >= self._p_crit:
            raise ValueError(
                f"bubble_dew: P={P:.4g} Pa is at/above water's critical pressure "
                f"({self._p_crit:.4g} Pa); no saturation temperature exists"
            )
        t_sat = self._props("T", "P", P, "Q", 0.0)
        return t_sat, t_sat

    def bubble_point(self, P: float, x: dict[str, float]) -> PhaseResult:
        """Saturated liquid at P: T = Tsat(P), with both saturated phase
        enthalpies populated (so h_fg = H_vapor - H_liquid)."""
        self._check_z(x)
        t_sat, _ = self.bubble_dew(P, x)
        comp = {self._cid: 1.0}
        hl = self._props("Hmolar", "P", P, "Q", 0.0) + self._h_offset
        hv = self._props("Hmolar", "P", P, "Q", 1.0) + self._h_offset
        return PhaseResult(T=t_sat, P=P, H=hl, phase="liquid", vapor_fraction=0.0,
                           x=comp, y=comp, H_liquid=hl, H_vapor=hv)

    # -- forced-phase properties (per-phase protocol extensions) ----------------
    # For a pure fluid only one phase is stable at any (T, P) off the saturation
    # line; a "forced" property of the other phase is metastable and outside
    # IAPWS-95's domain via CoolProp's PT interface. We return the stable-phase
    # value where the requested phase IS stable, and clamp to the saturated
    # phase at P otherwise — continuous, and exact everywhere a flash actually
    # places that phase.
    def _forced(self, T: float, P: float, output: str, want_vapor: bool) -> float:
        if P < self._p_crit:
            t_sat = self._props("T", "P", P, "Q", 0.0)
            stable_ok = (T >= t_sat) if want_vapor else (T <= t_sat)
            if not stable_ok:
                return self._props(output, "P", P, "Q", 1.0 if want_vapor else 0.0)
        return self._props(output, "T", T, "P", P)

    def k_values(self, T: float, P: float, x: dict[str, float],
                 y: dict[str, float]) -> dict[str, float]:
        """Pure-component K-value: K = Psat(T) / P (exact at saturation, the
        standard pure-fluid limit elsewhere)."""
        self._check_z(x)
        self._check_z(y)
        if T >= self._t_crit:
            raise ValueError(
                f"k_values: T={T:.2f} K is at/above water's critical temperature "
                f"({self._t_crit:.2f} K); no vapor pressure exists"
            )
        p_sat = self._props("P", "T", T, "Q", 0.0)
        return {self._cid: p_sat / P}

    def enthalpy_liquid(self, T: float, P: float, x: dict[str, float]) -> float:
        self._check_z(x)
        return self._forced(T, P, "Hmolar", want_vapor=False) + self._h_offset

    def enthalpy_vapor(self, T: float, P: float, y: dict[str, float]) -> float:
        self._check_z(y)
        return self._forced(T, P, "Hmolar", want_vapor=True) + self._h_offset

    def volume_liquid(self, T: float, P: float, x: dict[str, float]) -> float:
        self._check_z(x)
        return 1.0 / self._forced(T, P, "Dmolar", want_vapor=False)

    def volume_vapor(self, T: float, P: float, y: dict[str, float]) -> float:
        self._check_z(y)
        return 1.0 / self._forced(T, P, "Dmolar", want_vapor=True)

    # -- not supported ---------------------------------------------------------
    def lnKeq(self, stoich: dict[str, float], T: float) -> float:
        raise NotImplementedError(
            "coolprop:Water is a pure-water steam-tables package; reaction "
            "equilibrium (lnKeq) is not available — select a thermo:* property "
            "package (thermo:PR / thermo:SRK / thermo:NRTL) for reactors."
        )

    def flash_pt_3p(self, T: float, P: float, z: dict[str, float]) -> ThreePhaseResult:
        raise NotImplementedError(
            "coolprop:Water does not support three-phase (VLLE) flashes — a pure "
            "fluid cannot split into two liquids; select 'thermo:PR' or "
            "'thermo:SRK' for three-phase separators."
        )

    def flash_ph_3p(self, P: float, H: float, z: dict[str, float]) -> ThreePhaseResult:
        raise NotImplementedError(
            "coolprop:Water does not support three-phase (VLLE) flashes — a pure "
            "fluid cannot split into two liquids; select 'thermo:PR' or "
            "'thermo:SRK' for three-phase separators."
        )

    @classmethod
    def from_spec(cls, spec: str, components: list[str]) -> "CoolPropWaterPackage":
        """Build from a flowsheet ``property_package`` string ``"coolprop:Water"``."""
        backend, _, method = spec.partition(":")
        if backend != "coolprop" or method.lower() != "water":
            raise ValueError(
                f"CoolPropWaterPackage builds only 'coolprop:Water' (got {spec!r})"
            )
        return cls(components)
