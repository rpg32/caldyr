from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register


@register("FiredHeater")
class FiredHeater(UnitOp):
    """Direct-fired heater (furnace): bring a stream to ``params['T_out']`` (or
    apply a fixed process duty ``params['Q']``), dropping pressure by
    ``params['dP']`` — the same spec contract as :class:`Heater`, but heat is
    supplied by burning fuel at a fired efficiency ``params['efficiency']``
    (default 0.85, a typical value for a modern process furnace; Turton 4e
    Ch. 8 uses 0.80-0.90 for fired heaters).

    The energy ``duty`` port reports the PROCESS duty Q = n·(H_out − H_in), so
    flowsheet energy balances close exactly like a Heater's. The *fuel* duty —
    what the burners must release, Q_fuel = Q / efficiency — is published on
    ``unit.design`` after each solve::

        unit.design = {"process_duty": Q, "fuel_duty": Q / efficiency,
                       "efficiency": efficiency}

    The economics layer sizes the heater on the process duty (the Turton
    fired-heater correlation capacity) and books fuel cost on the fuel duty.
    A fired heater only heats: a spec that implies cooling raises.
    """

    design: dict[str, float] | None = None

    def define_ports(self) -> list[Port]:
        return [Port("in1", "inlet"), Port("out", "outlet"), Port("duty", "outlet", "energy")]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        inlet = inlets.get("in1")
        if inlet is None or not inlet.molar_flow:
            raise ValueError(f"FiredHeater {self.id!r}: missing or empty inlet on 'in1'")

        eta = float(self.params.get("efficiency", 0.85))
        if not 0.0 < eta <= 1.0:
            raise ValueError(
                f"FiredHeater {self.id!r}: efficiency={eta} must be in (0, 1]"
            )

        T_in, P_in, n = inlet.require_state()
        z = inlet.normalized_z()
        dP = float(self.params.get("dP", 0.0))
        P_out = P_in - dP
        H_in = inlet.H if inlet.H is not None else pp.enthalpy(T_in, P_in, z)

        has_T = self.params.get("T_out") is not None
        has_Q = self.params.get("Q") is not None
        if has_T == has_Q:
            raise ValueError(
                f"FiredHeater {self.id!r}: specify exactly one of 'T_out' or 'Q' "
                f"(got T_out={self.params.get('T_out')}, Q={self.params.get('Q')})"
            )

        if has_T:
            res = pp.flash_pt(float(self.params["T_out"]), P_out, z)
            Q = n * (res.H - H_in)
        else:
            Q = float(self.params["Q"])
            H_out = H_in + Q / n
            res = pp.flash_ph(P_out, H_out, z)

        if Q < 0.0:
            raise ValueError(
                f"FiredHeater {self.id!r}: process duty {Q:.3g} W is negative — a "
                f"fired heater only heats (T_out below the inlet temperature?). "
                f"Use a Heater or AirCooler for cooling service."
            )

        self.design = {"process_duty": Q, "fuel_duty": Q / eta, "efficiency": eta}

        out = Stream(
            id=f"{self.id}.out",
            components=list(inlet.components),
            T=res.T, P=res.P, molar_flow=n, z=z,
            H=res.H, phase=res.phase, vapor_fraction=res.vapor_fraction,
        )
        return {"out": out, "duty": EnergyStream(id=f"{self.id}.duty", duty=Q)}
