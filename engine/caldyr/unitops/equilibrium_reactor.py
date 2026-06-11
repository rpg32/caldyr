import math

from scipy.optimize import brentq

from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register
from .reaction import Reaction, apply_extents, reactor_outlet

_P_REF = 1e5     # Pa, ideal-gas standard state for the equilibrium constant


@register("EquilibriumReactor")
class EquilibriumReactor(UnitOp):
    """Single-reaction gas-phase equilibrium reactor, isothermal at ``params['T']``.

    Solves the reaction extent ξ so the ideal-gas mass-action expression equals
    K(T) (from the property package's :meth:`lnKeq`):

        K = Π y_i^ν_i · (P / P_ref)^Σν

    Activities are ideal-gas partial pressures referenced to 1 bar — consistent
    with K built from ideal-gas formation Gibbs energies. Duty to hold T is
    reported on the energy port.
    """

    def define_ports(self) -> list[Port]:
        return [Port("in1", "inlet"), Port("out", "outlet"),
                Port("duty", "outlet", "energy")]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        inlet = inlets.get("in1")
        if inlet is None or not inlet.molar_flow:
            raise ValueError(f"EquilibriumReactor {self.id!r}: missing/empty inlet on 'in1'")

        rxn = Reaction.from_param(self.params.get("reaction") or self.params["reactions"][0])
        T = float(self.params["T"])
        _, P_in, n_in = inlet.require_state()
        P_out = P_in - float(self.params.get("dP", 0.0))
        z_in = inlet.normalized_z()
        n0 = {c: n_in * z_in.get(c, 0.0) for c in inlet.components}

        xi = self._solve_extent(rxn, n0, pp.lnKeq(rxn.stoich, T), P_out)
        moles_out = apply_extents(n0, [rxn], [xi])

        out, duty = reactor_outlet(self.id, inlet, pp, moles_out, P_out, T)
        return {"out": out, "duty": EnergyStream(id=f"{self.id}.duty", duty=duty)}

    @staticmethod
    def _solve_extent(rxn: Reaction, n0: dict[str, float], lnK: float, P: float) -> float:
        stoich = rxn.stoich
        dn = rxn.dn

        # Extent window: reactants cannot go below zero (upper bound); products
        # cannot reverse below zero (lower bound).
        xi_hi = min((n0.get(c, 0.0) / -nu for c, nu in stoich.items() if nu < 0), default=0.0)
        xi_lo = max((-n0.get(c, 0.0) / nu for c, nu in stoich.items() if nu > 0), default=0.0)
        span = xi_hi - xi_lo
        if span <= 1e-12:
            return 0.0
        eps = 1e-9 * max(span, 1.0)

        def residual(xi: float) -> float:
            n = {c: n0.get(c, 0.0) + stoich.get(c, 0.0) * xi for c in n0}
            n_tot = sum(n.values())
            log_q = sum(nu * math.log(max(n[c] / n_tot, 1e-300)) for c, nu in stoich.items())
            log_q += dn * math.log(P / _P_REF)
            return log_q - lnK            # solve in log space for numerical range

        lo, hi = xi_lo + eps, xi_hi - eps
        f_lo, f_hi = residual(lo), residual(hi)
        if f_lo * f_hi > 0:               # no sign change: equilibrium at a bound
            return lo if abs(f_lo) < abs(f_hi) else hi
        return float(brentq(residual, lo, hi, xtol=1e-12, rtol=1e-12))
