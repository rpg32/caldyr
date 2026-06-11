import math

from ..core import Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register

_DT = 0.5   # K, half-step for the constant-Cp estimate used by effectiveness-NTU


@register("HeatExchanger")
class HeatExchanger(UnitOp):
    """Two-stream heat exchanger (no heat lost to surroundings: the duty the hot
    side loses is exactly what the cold side gains).

    Exactly one operating spec in ``params``:
      * ``duty`` — heat transferred hot→cold, W;
      * ``T_hot_out`` or ``T_cold_out`` — a target outlet temperature;
      * ``UA`` — the conductance, rated by the effectiveness-NTU method with
        ``arrangement`` ``"counterflow"`` (default) or ``"cocurrent"``.

    Outlet states always come from a rigorous PH flash, so phase changes are
    captured; only the effectiveness-NTU duty assumes constant heat capacities
    (estimated at the inlets). ``dP_hot`` / ``dP_cold`` set pressure drops.
    """

    _SPECS = ("duty", "T_hot_out", "T_cold_out", "UA")

    def define_ports(self) -> list[Port]:
        return [
            Port("hot_in", "inlet"), Port("cold_in", "inlet"),
            Port("hot_out", "outlet"), Port("cold_out", "outlet"),
        ]

    @staticmethod
    def lmtd(dT1: float, dT2: float) -> float:
        """Log-mean temperature difference for the two end approaches."""
        if dT1 <= 0 or dT2 <= 0:
            raise ValueError(f"non-positive approach temperature(s): {dT1}, {dT2}")
        if abs(dT1 - dT2) < 1e-9:
            return dT1
        return (dT1 - dT2) / math.log(dT1 / dT2)

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        hot, cold = inlets.get("hot_in"), inlets.get("cold_in")
        if hot is None or not hot.molar_flow:
            raise ValueError(f"HeatExchanger {self.id!r}: missing/empty 'hot_in'")
        if cold is None or not cold.molar_flow:
            raise ValueError(f"HeatExchanger {self.id!r}: missing/empty 'cold_in'")

        present = [s for s in self._SPECS if self.params.get(s) is not None]
        if len(present) != 1:
            raise ValueError(
                f"HeatExchanger {self.id!r}: give exactly one of {self._SPECS} (got {present})"
            )
        spec = present[0]

        Th, Ph, nh = hot.require_state()
        Tc, Pc, nc = cold.require_state()
        zh, zc = hot.normalized_z(), cold.normalized_z()
        Ph_out = Ph - float(self.params.get("dP_hot", 0.0))
        Pc_out = Pc - float(self.params.get("dP_cold", 0.0))
        h_hot_in = hot.H if hot.H is not None else pp.enthalpy(Th, Ph, zh)
        h_cold_in = cold.H if cold.H is not None else pp.enthalpy(Tc, Pc, zc)

        if spec == "duty":
            Q = float(self.params["duty"])
        elif spec == "T_hot_out":
            Q = nh * (h_hot_in - pp.enthalpy(float(self.params["T_hot_out"]), Ph_out, zh))
        elif spec == "T_cold_out":
            Q = nc * (pp.enthalpy(float(self.params["T_cold_out"]), Pc_out, zc) - h_cold_in)
        else:  # UA via effectiveness-NTU
            Q = self._duty_from_ua(pp, Th, Ph, zh, nh, Tc, Pc, zc, nc)

        res_hot = pp.flash_ph(Ph_out, h_hot_in - Q / nh, zh)
        res_cold = pp.flash_ph(Pc_out, h_cold_in + Q / nc, zc)
        return {
            "hot_out": Stream(
                id=f"{self.id}.hot_out", components=list(hot.components),
                T=res_hot.T, P=res_hot.P, molar_flow=nh, z=zh,
                H=res_hot.H, phase=res_hot.phase, vapor_fraction=res_hot.vapor_fraction,
            ),
            "cold_out": Stream(
                id=f"{self.id}.cold_out", components=list(cold.components),
                T=res_cold.T, P=res_cold.P, molar_flow=nc, z=zc,
                H=res_cold.H, phase=res_cold.phase, vapor_fraction=res_cold.vapor_fraction,
            ),
        }

    def _duty_from_ua(self, pp, Th, Ph, zh, nh, Tc, Pc, zc, nc) -> float:
        ua = float(self.params["UA"])
        c_hot = nh * self._cp(pp, Th, Ph, zh)        # heat capacity rates, W/K
        c_cold = nc * self._cp(pp, Tc, Pc, zc)
        c_min, c_max = sorted((c_hot, c_cold))
        cr = c_min / c_max
        ntu = ua / c_min

        arrangement = self.params.get("arrangement", "counterflow")
        if arrangement == "cocurrent":
            eps = (1.0 - math.exp(-ntu * (1.0 + cr))) / (1.0 + cr)
        elif arrangement == "counterflow":
            if abs(cr - 1.0) < 1e-9:
                eps = ntu / (1.0 + ntu)
            else:
                e = math.exp(-ntu * (1.0 - cr))
                eps = (1.0 - e) / (1.0 - cr * e)
        else:
            raise ValueError(f"HeatExchanger {self.id!r}: unknown arrangement {arrangement!r}")

        return eps * c_min * (Th - Tc)

    @staticmethod
    def _cp(pp, T, P, z) -> float:
        """Constant-pressure molar heat capacity, J/mol/K, by central difference."""
        return (pp.enthalpy(T + _DT, P, z) - pp.enthalpy(T - _DT, P, z)) / (2.0 * _DT)
