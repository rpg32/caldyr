from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register


@register("Compressor")
class Compressor(UnitOp):
    """Adiabatic gas compressor. Raises pressure to ``params['P_out']`` with
    isentropic efficiency ``params['eta']`` (default 0.75): the ideal path is
    constant-entropy to P_out, and the real enthalpy rise is the ideal rise
    divided by η. Shaft work is reported on the energy ``work`` port.
    """

    def define_ports(self) -> list[Port]:
        return [Port("in1", "inlet"), Port("out", "outlet"),
                Port("work", "outlet", "energy")]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        inlet = inlets.get("in1")
        if inlet is None or not inlet.molar_flow:
            raise ValueError(f"Compressor {self.id!r}: missing or empty inlet on 'in1'")

        T_in, P_in, n = inlet.require_state()
        z = inlet.normalized_z()
        P_out = float(self.params["P_out"])
        eta = float(self.params.get("eta", 0.75))
        if P_out <= P_in:
            raise ValueError(
                f"Compressor {self.id!r}: outlet P {P_out} must exceed inlet P {P_in}"
            )
        if not 0.0 < eta <= 1.0:
            raise ValueError(f"Compressor {self.id!r}: eta={eta} must be in (0, 1]")

        H_in = inlet.H if inlet.H is not None else pp.enthalpy(T_in, P_in, z)
        s_in = pp.entropy(T_in, P_in, z)
        h_ideal = pp.flash_ps(P_out, s_in, z).H        # isentropic outlet enthalpy
        w_molar = (h_ideal - H_in) / eta               # actual shaft work, J/mol
        res = pp.flash_ph(P_out, H_in + w_molar, z)

        out = Stream(
            id=f"{self.id}.out", components=list(inlet.components),
            T=res.T, P=res.P, molar_flow=n, z=z,
            H=res.H, phase=res.phase, vapor_fraction=res.vapor_fraction,
        )
        return {"out": out, "work": EnergyStream(id=f"{self.id}.work", duty=n * w_molar)}
