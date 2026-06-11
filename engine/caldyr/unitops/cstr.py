"""Continuous stirred-tank reactor (CSTR) with power-law kinetics.

The classic mixing-cell balance: the whole volume sits at the *outlet*
composition and temperature, so each reaction's molar extent is

    ξ_j = V · r_j(C_out, T_out)          [mol/s]

with the power-law rate from :class:`~caldyr.unitops.reaction.KineticReaction`
and concentrations C_i = z_i / v(T,P,z) from the property package's bulk molar
volume (see :func:`~caldyr.unitops.reaction.concentrations` — one documented
concentration basis for vapor and liquid alike). The resulting algebraic
system in the extents (plus T, if adiabatic) is solved with ``scipy`` root.

Energy: isothermal if ``params['T']`` is given (duty to hold it reported on the
energy port); otherwise adiabatic — total enthalpy flow is conserved, and
because stream enthalpies are formation-inclusive the heat of reaction is
carried automatically (same mechanism as :class:`ConversionReactor`). The
adiabatic case is solved as a *nested* problem for robustness: the extents are
solved at fixed T (the well-behaved isothermal system), and a 1-D outer search
brackets and solves the energy balance in T (Brent). An adiabatic CSTR can
have multiple steady states (the classic ignition/extinction multiplicity);
the outward search from the feed temperature finds the steady state nearest
the feed condition. The outlet state then comes from a PH flash on the same
enthalpy surface, so the energy balance closes exactly.
"""
from __future__ import annotations

from scipy.optimize import brentq, root

from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register
from .reaction import (
    KineticReaction,
    apply_extents,
    concentrations,
    reactor_outlet,
)


class KineticSolveError(ValueError):
    """The kinetic reactor's balance equations failed to converge."""


@register("CSTR")
class CSTR(UnitOp):
    """Power-law kinetic CSTR.

    Params:
      * ``V`` (required) — reactor volume, m^3.
      * ``reactions`` (required) — list of kinetic reaction dicts
        (``{"stoich": ..., "key": ..., "k0": ..., "Ea": ..., "orders": ...}``;
        see :class:`KineticReaction`).
      * ``T`` (optional) — isothermal operating temperature, K. Absent →
        adiabatic.
      * ``dP`` (optional) — pressure drop, Pa.
    """

    def define_ports(self) -> list[Port]:
        return [Port("in1", "inlet"), Port("out", "outlet"),
                Port("duty", "outlet", "energy")]

    def _reactions(self) -> list[KineticReaction]:
        rxn_dicts = self.params.get("reactions")
        if not rxn_dicts:
            raise ValueError(f"{type(self).__name__} {self.id!r}: params['reactions'] "
                             f"is required (a list of kinetic reaction dicts)")
        return [KineticReaction.from_param(d) for d in rxn_dicts]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        inlet = inlets.get("in1")
        if inlet is None or not inlet.molar_flow:
            raise ValueError(f"CSTR {self.id!r}: missing/empty inlet on 'in1'")
        rxns = self._reactions()
        volume = float(self.params["V"])
        if volume <= 0.0:
            raise ValueError(f"CSTR {self.id!r}: params['V'] must be > 0 m^3")

        T_in, P_in, n_in = inlet.require_state()
        P_out = P_in - float(self.params.get("dP", 0.0))
        z_in = inlet.normalized_z()
        n0 = {c: n_in * z_in.get(c, 0.0) for c in inlet.components}
        t_spec_raw = self.params.get("T")
        t_spec: float | None = None if t_spec_raw is None else float(t_spec_raw)
        h_in = inlet.H if inlet.H is not None else pp.enthalpy(T_in, P_in, z_in)
        h_total_in = n_in * h_in
        m = len(rxns)

        def solve_extents(temp: float, guess: list[float]) -> list[float]:
            """Extents [mol/s] of the isothermal mixing-cell balance at temp."""
            def residuals(x) -> list[float]:
                xis = [xi * n_in for xi in x]
                moles = apply_extents(n0, rxns, xis)
                conc = concentrations(pp, temp, P_out, moles)
                return [(xis[j] - volume * rxns[j].rate(conc, temp)) / n_in
                        for j in range(m)]

            sol = root(residuals, [g / n_in for g in guess], method="hybr", tol=1e-12)
            # Judge convergence by the scaled residuals, not sol.success alone:
            # hybr reports "not making progress" when it cannot improve *below*
            # machine precision, which is success, not failure.
            if max(abs(r) for r in sol.fun) > 1e-9:
                raise KineticSolveError(
                    f"CSTR {self.id!r}: extent balance did not converge at "
                    f"T={temp:.2f} K ({sol.message}; residuals {list(sol.fun)})"
                )
            return [xi * n_in for xi in sol.x]

        if t_spec is not None:
            xis = solve_extents(t_spec, [0.0] * m)
        else:
            xis = self._solve_adiabatic(solve_extents, n0, rxns, pp, P_out,
                                        T_in, h_total_in)
        moles_out = {c: max(v, 0.0)
                     for c, v in apply_extents(n0, rxns, xis).items()}

        out, duty = reactor_outlet(self.id, inlet, pp, moles_out, P_out, t_spec)
        return {"out": out, "duty": EnergyStream(id=f"{self.id}.duty", duty=duty)}

    def _solve_adiabatic(self, solve_extents, n0, rxns, pp, P_out,
                         T_in: float, h_total_in: float) -> list[float]:
        """Outer 1-D energy balance: bracket the outlet temperature walking
        outward from the feed temperature, then solve with Brent. The extent
        solve is warm-started along the walk (continuation)."""
        warm = [0.0] * len(rxns)

        def energy_residual(temp: float) -> float:
            nonlocal warm
            warm = solve_extents(temp, warm)
            moles = {c: max(v, 0.0)
                     for c, v in apply_extents(n0, rxns, warm).items()}
            n_tot = sum(moles.values())
            return n_tot * pp.enthalpy(temp, P_out, moles) - h_total_in

        g_prev = energy_residual(T_in)
        if g_prev == 0.0:
            return warm
        direction = 1.0 if g_prev < 0 else -1.0     # exothermic -> look hotter
        t_prev, step, bracket = T_in, max(5.0, 0.01 * T_in), None
        for _ in range(80):
            t_next = max(t_prev + direction * step, 50.0)
            g_next = energy_residual(t_next)
            if g_prev * g_next <= 0.0:
                bracket = (min(t_prev, t_next), max(t_prev, t_next))
                break
            if t_next <= 50.0:
                break
            t_prev, g_prev, step = t_next, g_next, step * 1.5
        if bracket is None:
            raise KineticSolveError(
                f"CSTR {self.id!r}: could not bracket the adiabatic outlet "
                f"temperature (searched {'up' if direction > 0 else 'down'} "
                f"from {T_in:.1f} K)"
            )
        t_out = float(brentq(energy_residual, *bracket, xtol=1e-9))
        return solve_extents(t_out, warm)
