import math

from ..core import Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register


@register("Mixer")
class Mixer(UnitOp):
    """Adiabatic mixer: combine inlets, conserve mass and energy, flash to the
    outlet phase at the lowest inlet pressure minus ``params['dP']``.

    Energy balance: H_out (per mole) = Σ(n_i·H_i) / Σ n_i, with each inlet's
    molar enthalpy taken from the property package at its own (T, P, z). The
    outlet is then a PH flash, so its T/phase are consistent with that enthalpy.
    """

    def define_ports(self) -> list[Port]:
        return [Port("in1", "inlet"), Port("in2", "inlet"), Port("out", "outlet")]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        streams = [s for s in inlets.values() if s.molar_flow]
        if not streams:
            raise ValueError(f"Mixer {self.id!r}: no inlet streams with flow")

        components = streams[0].components
        n_total = 0.0
        moles: dict[str, float] = {c: 0.0 for c in components}
        H_flow = 0.0      # W (J/s): Σ n_i · H_i
        P_min = math.inf  # outlet sits at the lowest inlet pressure (minus dP)
        for s in streams:
            T_i, P_i, n_i = s.require_state()
            zi = s.normalized_z()
            n_total += n_i
            for c in components:
                moles[c] += n_i * zi.get(c, 0.0)
            H_i = s.H if s.H is not None else pp.enthalpy(T_i, P_i, zi)
            H_flow += n_i * H_i
            P_min = min(P_min, P_i)

        z_out = {c: moles[c] / n_total for c in components}
        H_out = H_flow / n_total                      # J/mol
        dP = float(self.params.get("dP", 0.0))
        P_out = P_min - dP

        res = pp.flash_ph(P_out, H_out, z_out)
        out = Stream(
            id=f"{self.id}.out",
            components=list(components),
            T=res.T, P=res.P, molar_flow=n_total, z=z_out,
            H=res.H, phase=res.phase, vapor_fraction=res.vapor_fraction,
        )
        return {"out": out}
