from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register


@register("Flash")
class FlashDrum(UnitOp):
    """Two-phase (vapor-liquid) flash drum / separator.

    Operating spec:
      * isothermal — both ``params['T']`` and ``params['P']`` given: a PT flash,
        with the heat needed to hold T reported on the ``duty`` energy port.
      * adiabatic — only ``params['P']`` given: a PH flash at the inlet enthalpy
        (duty = 0); the drum finds its own temperature.

    Outlets ``vapor`` and ``liquid`` carry the equilibrium phase compositions and
    flows; a single-phase feed yields a zero-flow stream on the absent phase.
    """

    def define_ports(self) -> list[Port]:
        return [
            Port("in1", "inlet"),
            Port("vapor", "outlet"),
            Port("liquid", "outlet"),
            Port("duty", "outlet", "energy"),
        ]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        inlet = inlets.get("in1")
        if inlet is None or not inlet.molar_flow:
            raise ValueError(f"Flash {self.id!r}: missing or empty inlet on 'in1'")

        T_in, P_in, n = inlet.require_state()
        z = inlet.normalized_z()
        P = float(self.params.get("P", P_in))
        H_in = inlet.H if inlet.H is not None else pp.enthalpy(T_in, P_in, z)

        if self.params.get("T") is not None:                 # isothermal
            res = pp.flash_pt(float(self.params["T"]), P, z)
            duty = n * (res.H - H_in)
        else:                                                # adiabatic
            res = pp.flash_ph(P, H_in, z)
            duty = 0.0

        vf = res.vapor_fraction
        y = res.y if res.y is not None else z
        x = res.x if res.x is not None else z
        h_vap = res.H_vapor if res.H_vapor is not None else res.H
        h_liq = res.H_liquid if res.H_liquid is not None else res.H

        vapor = Stream(
            id=f"{self.id}.vapor", components=list(inlet.components),
            T=res.T, P=P, molar_flow=n * vf, z=dict(y),
            H=h_vap, phase="vapor", vapor_fraction=1.0,
        )
        liquid = Stream(
            id=f"{self.id}.liquid", components=list(inlet.components),
            T=res.T, P=P, molar_flow=n * (1.0 - vf), z=dict(x),
            H=h_liq, phase="liquid", vapor_fraction=0.0,
        )
        return {
            "vapor": vapor,
            "liquid": liquid,
            "duty": EnergyStream(id=f"{self.id}.duty", duty=duty),
        }
