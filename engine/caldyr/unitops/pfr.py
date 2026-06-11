"""Plug-flow reactor (PFR) with power-law kinetics.

Integrates the steady-state molar balances along the reactor volume

    d n_i / dV = Σ_j ν_ij · r_j(C, T)        [mol/(s·m^3)]

with ``scipy.integrate.solve_ivp`` (LSODA — power-law networks are routinely
stiff). Rates and concentrations are exactly the CSTR's: power-law
:class:`~caldyr.unitops.reaction.KineticReaction` with C_i = z_i / v(T,P,z)
from the property package's bulk molar volume (one documented basis for vapor
and liquid; see :func:`~caldyr.unitops.reaction.concentrations`).

Energy: isothermal if ``params['T']`` is given (the total duty to hold it is
reported on the energy port); otherwise adiabatic. The adiabatic temperature
profile is **solved from enthalpy conservation at each volume node** — total
enthalpy flow Σ n_i·h is constant down an adiabatic PFR, so at every
right-hand-side evaluation T is recovered from the current composition with a
PH flash (rather than integrating a separate dT/dV equation, which would need
heat capacities and could drift off the enthalpy surface). Formation-inclusive
enthalpies carry the heat of reaction automatically.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register
from .cstr import KineticSolveError
from .reaction import KineticReaction, concentrations, reactor_outlet

_DEFAULT_N_STEPS = 50


@register("PFR")
class PFR(UnitOp):
    """Power-law kinetic plug-flow reactor.

    Params:
      * ``V`` (required) — reactor volume, m^3.
      * ``reactions`` (required) — list of kinetic reaction dicts
        (see :class:`KineticReaction`).
      * ``T`` (optional) — isothermal operating temperature, K. Absent →
        adiabatic (T from enthalpy conservation at each volume node).
      * ``n_steps`` (optional, default 50) — caps the integrator step at
        ``V / n_steps`` (resolution control; accuracy itself is governed by
        the integrator's rtol of 1e-9).
      * ``dP`` (optional) — pressure drop, Pa, applied as a single outlet
        drop; the kinetics are evaluated at the outlet pressure.
    """

    def define_ports(self) -> list[Port]:
        return [Port("in1", "inlet"), Port("out", "outlet"),
                Port("duty", "outlet", "energy")]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        inlet = inlets.get("in1")
        if inlet is None or not inlet.molar_flow:
            raise ValueError(f"PFR {self.id!r}: missing/empty inlet on 'in1'")
        rxn_dicts = self.params.get("reactions")
        if not rxn_dicts:
            raise ValueError(f"PFR {self.id!r}: params['reactions'] is required "
                             f"(a list of kinetic reaction dicts)")
        rxns = [KineticReaction.from_param(d) for d in rxn_dicts]
        volume = float(self.params["V"])
        if volume <= 0.0:
            raise ValueError(f"PFR {self.id!r}: params['V'] must be > 0 m^3")
        n_steps = int(self.params.get("n_steps", _DEFAULT_N_STEPS))

        T_in, P_in, n_in = inlet.require_state()
        P_out = P_in - float(self.params.get("dP", 0.0))
        z_in = inlet.normalized_z()
        comps = list(inlet.components)
        t_spec_raw = self.params.get("T")
        t_spec: float | None = None if t_spec_raw is None else float(t_spec_raw)
        h_in = inlet.H if inlet.H is not None else pp.enthalpy(T_in, P_in, z_in)
        h_total_in = n_in * h_in

        # stoichiometry matrix: nu[j][i] for reaction j, component i
        nu = np.array([[rxn.stoich.get(c, 0.0) for c in comps] for rxn in rxns])
        y0 = np.array([n_in * z_in.get(c, 0.0) for c in comps])

        def temperature(moles: dict[str, float]) -> float:
            if t_spec is not None:
                return t_spec
            n_tot = sum(moles.values())
            return pp.flash_ph(P_out, h_total_in / n_tot, moles).T

        def rhs(_v: float, y: np.ndarray) -> np.ndarray:
            moles = {c: max(float(val), 0.0) for c, val in zip(comps, y)}
            temp = temperature(moles)
            conc = concentrations(pp, temp, P_out, moles)
            rates = np.array([rxn.rate(conc, temp) for rxn in rxns])
            return nu.T @ rates

        sol = solve_ivp(rhs, (0.0, volume), y0, method="LSODA",
                        rtol=1e-9, atol=1e-12 * n_in, max_step=volume / n_steps)
        if not sol.success:
            raise KineticSolveError(
                f"PFR {self.id!r}: volume integration failed ({sol.message})"
            )
        moles_out = {c: max(float(v), 0.0) for c, v in zip(comps, sol.y[:, -1])}

        out, duty = reactor_outlet(self.id, inlet, pp, moles_out, P_out, t_spec)
        return {"out": out, "duty": EnergyStream(id=f"{self.id}.duty", duty=duty)}
