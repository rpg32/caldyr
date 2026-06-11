from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register


@register("Heater")
class Heater(UnitOp):
    """Heater/Cooler: bring a stream to ``params['T_out']`` (or apply a fixed
    duty ``params['Q']``), dropping pressure by ``params['dP']``. Duty is
    Q = n·(H_out − H_in); a positive Q heats, negative cools. The duty is
    reported on the energy ``duty`` port as an :class:`EnergyStream`.
    """

    def define_ports(self) -> list[Port]:
        return [Port("in1", "inlet"), Port("out", "outlet"), Port("duty", "outlet", "energy")]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        inlet = inlets.get("in1")
        if inlet is None or not inlet.molar_flow:
            raise ValueError(f"Heater {self.id!r}: missing or empty inlet on 'in1'")

        T_in, P_in, n = inlet.require_state()
        z = inlet.normalized_z()
        dP = float(self.params.get("dP", 0.0))
        P_out = P_in - dP
        H_in = inlet.H if inlet.H is not None else pp.enthalpy(T_in, P_in, z)

        has_T = self.params.get("T_out") is not None
        has_Q = self.params.get("Q") is not None
        if has_T == has_Q:
            raise ValueError(
                f"Heater {self.id!r}: specify exactly one of 'T_out' or 'Q' "
                f"(got T_out={self.params.get('T_out')}, Q={self.params.get('Q')})"
            )

        if has_T:
            res = pp.flash_pt(float(self.params["T_out"]), P_out, z)
            Q = n * (res.H - H_in)
        else:
            Q = float(self.params["Q"])
            H_out = H_in + Q / n
            res = pp.flash_ph(P_out, H_out, z)

        out = Stream(
            id=f"{self.id}.out",
            components=list(inlet.components),
            T=res.T, P=res.P, molar_flow=n, z=z,
            H=res.H, phase=res.phase, vapor_fraction=res.vapor_fraction,
        )
        return {"out": out, "duty": EnergyStream(id=f"{self.id}.duty", duty=Q)}
