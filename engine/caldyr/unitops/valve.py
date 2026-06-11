from ..core import Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register


@register("Valve")
class Valve(UnitOp):
    """Throttling valve: an isenthalpic (Joule-Thomson) pressure drop. The
    outlet pressure is ``params['P_out']`` if given, else ``P_in - params['dP']``.
    Enthalpy is conserved, so the outlet phase/temperature follow from a PH flash
    (a real fluid generally cools, and may partially vaporize).
    """

    def define_ports(self) -> list[Port]:
        return [Port("in1", "inlet"), Port("out", "outlet")]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        inlet = inlets.get("in1")
        if inlet is None or not inlet.molar_flow:
            raise ValueError(f"Valve {self.id!r}: missing or empty inlet on 'in1'")

        T_in, P_in, n = inlet.require_state()
        z = inlet.normalized_z()
        if self.params.get("P_out") is not None:
            P_out = float(self.params["P_out"])
        else:
            P_out = P_in - float(self.params.get("dP", 0.0))
        if P_out >= P_in:
            raise ValueError(
                f"Valve {self.id!r}: outlet P {P_out} must be below inlet P {P_in}"
            )

        H_in = inlet.H if inlet.H is not None else pp.enthalpy(T_in, P_in, z)
        res = pp.flash_ph(P_out, H_in, z)
        out = Stream(
            id=f"{self.id}.out", components=list(inlet.components),
            T=res.T, P=res.P, molar_flow=n, z=z,
            H=res.H, phase=res.phase, vapor_fraction=res.vapor_fraction,
        )
        return {"out": out}
