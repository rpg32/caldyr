from typing import Any

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

    **Radiant/convective design split (opt-in).** Set ``params['design_split']``
    truthy (or supply any of ``fuel``/``excess_air``/``bridgewall_T``) and the
    unit additionally runs the firebox combustion + radiant/convective design of
    :mod:`caldyr.economics.fired_heater_design` (Hameed §4.3): it predicts the
    fuel and air molar flows, the flue-gas composition and temperature, the
    adiabatic flame and bridgewall/stack temperatures, the radiant vs convective
    duty split, and the radiant/convective tube areas — all published as nested
    dicts on ``unit.design`` (``combustion`` and ``firing``). Design knobs:
    ``fuel`` (composition dict, default pure methane), ``excess_air`` (fraction,
    default 0.15), ``fuel_T``/``air_T`` (K), ``bridgewall_T`` (K),
    ``loss_fraction``, ``radiant_flux`` (W/m^2), ``convective_U`` (W/m^2K).
    """

    design: dict[str, Any] | None = None

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

        if self._wants_design_split():
            self._add_design_split(Q, eta, T_in, res.T)

        out = Stream(
            id=f"{self.id}.out",
            components=list(inlet.components),
            T=res.T, P=res.P, molar_flow=n, z=z,
            H=res.H, phase=res.phase, vapor_fraction=res.vapor_fraction,
        )
        return {"out": out, "duty": EnergyStream(id=f"{self.id}.duty", duty=Q)}

    def _wants_design_split(self) -> bool:
        triggers = ("design_split", "fuel", "excess_air", "bridgewall_T")
        return any(self.params.get(k) is not None for k in triggers)

    def _add_design_split(self, Q: float, eta: float,
                          T_in: float, T_out: float) -> None:
        """Run the radiant/convective firebox design and fold the results onto
        ``unit.design`` (assumes self.design already holds the basic duties)."""
        from ..economics.fired_heater_design import design_fired_heater

        p = self.params
        kw: dict[str, object] = {}
        for key in ("excess_air", "fuel_T", "air_T", "bridgewall_T",
                    "loss_fraction", "radiant_flux", "convective_U"):
            if p.get(key) is not None:
                kw[key] = float(p[key])
        fuel = p.get("fuel")
        if isinstance(fuel, str):
            fuel = {fuel: 1.0}
        if fuel:
            kw["fuel"] = fuel

        d = design_fired_heater(Q, eta, T_in, T_out, **kw)  # type: ignore[arg-type]
        c = d.combustion
        assert self.design is not None
        self.design["combustion"] = {
            "lhv_mix": c.lhv_mix,
            "fuel_flow": c.fuel_flow,
            "fuel_flows": c.fuel_flows,
            "air_flow": c.air_flow,
            "o2_stoich": c.o2_stoich,
            "excess_air": c.excess_air,
            "flue_flow": c.flue_flow,
            "flue_composition": c.flue_composition,
            "flue_flows": c.flue_flows,
        }
        self.design["firing"] = {
            "fired_duty": d.fired_duty,
            "heat_available": d.heat_available,
            "efficiency_gross": d.efficiency_gross,
            "flame_temperature": d.flame_temperature,
            "bridgewall_temperature": d.bridgewall_temperature,
            "stack_temperature": d.stack_temperature,
            "radiant_duty": d.radiant_duty,
            "convective_duty": d.convective_duty,
            "radiant_fraction": d.radiant_fraction,
            "casing_loss": d.casing_loss,
            "stack_loss": d.stack_loss,
            "radiant_area": d.radiant_area,
            "convective_area": d.convective_area,
            "radiant_flux": d.radiant_flux,
            "convective_U": d.convective_U,
            "convective_lmtd": d.convective_lmtd,
            "notes": d.notes,
        }
