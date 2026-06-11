"""Pipe segment: single-phase pressure drop along a pipe (Darcy-Weisbach).

The HYSYS "Pipe Segment" equivalent (Hameed, *Chemical Process Simulations
using Aspen HYSYS*, Wiley 2025, ch. 7 / Eq. 7.1): the pressure gradient is the
sum of a friction term, a gravity (elevation) term and minor (fitting) losses,

    dP = f (L/D) rho v^2 / 2  +  rho g dz  +  K rho v^2 / 2

with the Darcy friction factor ``f`` from **Churchill's all-regime correlation**
(S.W. Churchill, "Friction-factor equation spans all fluid-flow regimes",
*Chemical Engineering* 84(24), 91-92, 1977), which reproduces f = 64/Re in the
laminar limit and tracks Colebrook within ~2% in the turbulent regime — one
smooth expression, no regime branching. The pipe is marched in ``segments``
increments (HYSYS does the same): density and velocity are re-evaluated from
the property package at each increment's local pressure, so gas lines pick up
the rising velocity (and friction gradient) as the gas expands.

Viscosity is NOT part of the engine's PropertyPackage protocol, so this unit
sources it from the `thermo`/`chemicals` correlation stack directly:
pure-component temperature-dependent viscosities from
``thermo.ChemicalConstantsPackage.from_IDs`` (the same DIPPR-quality databank
the thermo property packages read), mixed with the Arrhenius mole-fraction
log-mean rule for liquids (Grunberg-Nissan with zero interaction; Poling,
Prausnitz & O'Connell, *The Properties of Gases and Liquids* 5e, ch. 9) and
Herning-Zipperer for gases (ibid.). The pressure dependence of viscosity is
neglected (small for liquids and for gases below ~0.9 Tc / moderate P).

Modeling scope (v1, honest limitations):

* **Single-phase only.** A two-phase (VLE) inlet — or a liquid that flashes as
  the pressure falls along the pipe — raises :class:`PipeFlowError`; two-phase
  pressure drop (Beggs & Brill et al.) is future work.
* **Isothermal.** The outlet is at the inlet temperature; heat loss to the
  surroundings (the HYSYS heat-transfer page) is out of scope. The enthalpy is
  re-flashed at the outlet pressure, so downstream energy balances stay exact.
* The kinetic/acceleration term of Eq. 7.1 is neglected (negligible for
  liquids and low-Mach gas flow; Perry's 8e sec. 6).
"""
from __future__ import annotations

import math
from functools import lru_cache
from typing import Any

from ..core import Port, Stream, UnitOp
from ..core.components_db import molar_mass
from ..core.unitop import PortStream
from .base import register

_G = 9.80665                  # m/s^2, standard gravity
_VF_EPS = 1e-9


class PipeFlowError(ValueError):
    """A PipeSegment cannot model the requested flow (two-phase fluid, or the
    pressure drop exhausts the inlet pressure)."""


def churchill_friction_factor(re: float, rel_roughness: float) -> float:
    """Darcy friction factor, Churchill (1977) all-regime equation.

    f = 8 [ (8/Re)^12 + (A + B)^-1.5 ]^(1/12)
    A = (2.457 ln(1 / ((7/Re)^0.9 + 0.27 eps/D)))^16,  B = (37530/Re)^16

    Valid for all Re and relative roughness: reduces to f = 64/Re for laminar
    flow and follows Colebrook in the turbulent regime (within ~2%).
    """
    if re <= 0.0:
        raise ValueError(f"Reynolds number must be positive, got {re}")
    a = (2.457 * math.log(1.0 / ((7.0 / re) ** 0.9 + 0.27 * rel_roughness))) ** 16
    b = (37530.0 / re) ** 16
    return 8.0 * ((8.0 / re) ** 12 + 1.0 / (a + b) ** 1.5) ** (1.0 / 12.0)


@lru_cache(maxsize=64)
def _viscosity_correlations(components: tuple[str, ...]):
    """Pure-component viscosity correlation objects (and molar masses, g/mol)
    for an ordered component tuple, from the thermo/chemicals databank. Cached:
    ``from_IDs`` hits a database and is slow."""
    from thermo import ChemicalConstantsPackage

    constants, props = ChemicalConstantsPackage.from_IDs(list(components))
    return constants.MWs, props.ViscosityLiquids, props.ViscosityGases


def mixture_viscosity(components: tuple[str, ...], zs: list[float], T: float,
                      phase: str) -> float:
    """Dynamic viscosity (Pa*s) of a single-phase mixture at temperature ``T``.

    Pure-component values come from the `thermo` correlation databank; mixing
    is Arrhenius (mole-fraction log-mean) for liquids and Herning-Zipperer for
    gases (Poling et al. 5e ch. 9). Pressure dependence is neglected.
    """
    mws, visc_liq, visc_gas = _viscosity_correlations(components)
    models = visc_liq if phase == "liquid" else visc_gas
    mus: list[float] = []
    for comp, model in zip(components, models):
        mu = model.T_dependent_property(T)
        if mu is None or not math.isfinite(mu) or mu <= 0.0:
            raise PipeFlowError(
                f"no {phase} viscosity correlation value for {comp!r} at "
                f"{T:.1f} K (thermo databank returned {mu!r})"
            )
        mus.append(float(mu))
    if phase == "liquid":
        return math.exp(sum(z * math.log(mu) for z, mu in zip(zs, mus)))
    num = sum(z * mu * math.sqrt(mw) for z, mu, mw in zip(zs, mus, mws))
    den = sum(z * math.sqrt(mw) for z, mw in zip(zs, mws))
    return num / den


def _flow_regime(re: float) -> str:
    if re < 2300.0:
        return "laminar"
    if re < 4000.0:
        return "transitional"
    return "turbulent"


@register("PipeSegment")
class PipeSegment(UnitOp):
    """Single-phase pipe segment (see module docstring for the model).

    Parameters
    ----------
    length : m, pipe length (required, > 0)
    diameter : m, inner diameter (required, > 0)
    roughness : m, absolute roughness (default 4.5e-5, commercial steel —
        Crane TP-410 / Perry's 8e Table 6-1 give 0.045 mm)
    elevation_change : m, outlet height minus inlet height (default 0)
    fittings_K : extra velocity heads of minor losses (elbows, valves...; sum
        of Crane TP-410 K factors; default 0)
    segments : number of marching increments for the pressure profile
        (default 10)

    After a solve, ``unit.design`` carries dP_friction / dP_elevation /
    dP_fittings / dP_total (Pa), inlet velocity (m/s), Re, friction_factor,
    flow_regime, density (kg/m^3), viscosity (Pa*s) and the plot-ready
    pressure-vs-length profile (``L_profile`` m, ``P_profile`` Pa).
    """

    design: dict[str, Any] | None = None

    def define_ports(self) -> list[Port]:
        return [Port("in1", "inlet"), Port("out", "outlet")]

    def _param(self, name: str) -> float:
        value = self.params.get(name)
        if value is None:
            raise ValueError(f"PipeSegment {self.id!r}: parameter {name!r} is required")
        return float(value)

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        inlet = inlets.get("in1")
        if inlet is None or not inlet.molar_flow:
            raise ValueError(f"PipeSegment {self.id!r}: missing or empty inlet on 'in1'")

        T_in, P_in, n = inlet.require_state()
        z = inlet.normalized_z()

        length = self._param("length")
        diameter = self._param("diameter")
        roughness = float(self.params.get("roughness", 4.5e-5))
        dz_total = float(self.params.get("elevation_change", 0.0))
        fittings_k = float(self.params.get("fittings_K", 0.0))
        nseg = int(self.params.get("segments", 10))
        if length <= 0.0 or diameter <= 0.0:
            raise ValueError(
                f"PipeSegment {self.id!r}: length and diameter must be positive "
                f"(got length={length}, diameter={diameter})"
            )
        if roughness < 0.0 or fittings_k < 0.0:
            raise ValueError(
                f"PipeSegment {self.id!r}: roughness and fittings_K must be >= 0 "
                f"(got roughness={roughness}, fittings_K={fittings_k})"
            )
        if nseg < 1:
            raise ValueError(f"PipeSegment {self.id!r}: segments must be >= 1 (got {nseg})")

        res_in = pp.flash_pt(T_in, P_in, z)
        if _VF_EPS < res_in.vapor_fraction < 1.0 - _VF_EPS:
            raise PipeFlowError(
                f"PipeSegment {self.id!r}: inlet is two-phase (vapor fraction "
                f"{res_in.vapor_fraction:.3f} at {T_in:.1f} K / {P_in:.3g} Pa). "
                f"This v1 models single-phase flow only; two-phase pressure "
                f"drop (Beggs-Brill) is not implemented yet."
            )
        phase = "liquid" if res_in.vapor_fraction <= _VF_EPS else "vapor"

        comps = tuple(c for c in inlet.components if z.get(c, 0.0) > 0.0)
        zs = [z[c] for c in comps]
        mu = mixture_viscosity(comps, zs, T_in, phase)
        mw = sum(zi * molar_mass(c) for c, zi in zip(comps, zs))   # kg/mol
        area = math.pi * diameter * diameter / 4.0

        # March down the pipe re-evaluating density/velocity at local pressure.
        d_len = length / nseg
        d_z = dz_total / nseg
        p = P_in
        l_profile = [0.0]
        p_profile = [P_in]
        dp_fric = dp_elev = dp_fit = 0.0
        v_in = re_in = f_in = rho_in = 0.0
        for i in range(nseg):
            vm = pp.volume(T_in, p, z)                  # m^3/mol at local P
            rho = mw / vm
            vel = n * vm / area
            re = rho * vel * diameter / mu
            f = churchill_friction_factor(re, roughness / diameter)
            if i == 0:
                v_in, re_in, f_in, rho_in = vel, re, f, rho
            dyn = rho * vel * vel / 2.0                 # dynamic pressure, Pa
            seg_fric = f * (d_len / diameter) * dyn
            seg_elev = rho * _G * d_z
            seg_fit = (fittings_k / nseg) * dyn
            p_next = p - seg_fric - seg_elev - seg_fit
            if p_next <= 0.0 or P_in - p_next >= P_in:
                raise PipeFlowError(
                    f"PipeSegment {self.id!r}: pressure drop exceeds the inlet "
                    f"pressure ({P_in:.4g} Pa) {((i + 1) * d_len):.1f} m down a "
                    f"{length:.1f} m pipe — the specified flow cannot be driven "
                    f"through this line (D={diameter:.4g} m, "
                    f"friction so far {dp_fric + seg_fric:.4g} Pa)"
                )
            dp_fric += seg_fric
            dp_elev += seg_elev
            dp_fit += seg_fit
            p = p_next
            l_profile.append((i + 1) * d_len)
            p_profile.append(p)

        res_out = pp.flash_pt(T_in, p, z)
        if _VF_EPS < res_out.vapor_fraction < 1.0 - _VF_EPS:
            raise PipeFlowError(
                f"PipeSegment {self.id!r}: the fluid flashes inside the pipe — "
                f"the outlet at {T_in:.1f} K / {p:.4g} Pa is two-phase (vapor "
                f"fraction {res_out.vapor_fraction:.3f}). Two-phase flow "
                f"(Beggs-Brill) is not implemented yet."
            )
        phase_out = "liquid" if res_out.vapor_fraction <= _VF_EPS else "vapor"
        if phase_out != phase:
            # A pure(-ish) component can jump straight across the saturation
            # line (VF 0 -> 1) without ever flashing to an intermediate VF.
            raise PipeFlowError(
                f"PipeSegment {self.id!r}: the fluid changes phase inside the "
                f"pipe ({phase} in, {phase_out} out at {T_in:.1f} K / "
                f"{p:.4g} Pa) — the pressure falls below saturation along the "
                f"line. Two-phase/flashing flow is not implemented yet."
            )

        vm_out = pp.volume(T_in, p, z)
        self.design = {
            "dP_total": P_in - p,
            "dP_friction": dp_fric,
            "dP_elevation": dp_elev,
            "dP_fittings": dp_fit,
            "velocity": v_in,
            "velocity_out": n * vm_out / area,
            "Re": re_in,
            "friction_factor": f_in,
            "flow_regime": _flow_regime(re_in),
            "phase": phase,
            "density": rho_in,
            "viscosity": mu,
            "length": length,
            "diameter": diameter,
            "L_profile": l_profile,
            "P_profile": p_profile,
        }

        out = Stream(
            id=f"{self.id}.out", components=list(inlet.components),
            T=res_out.T, P=res_out.P, molar_flow=n, z=z,
            H=res_out.H, phase=res_out.phase, vapor_fraction=res_out.vapor_fraction,
        )
        return {"out": out}
