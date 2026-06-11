from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register


class AirCoolerApproachError(ValueError):
    """The specified outlet temperature is closer to the ambient-air inlet than
    the minimum approach allows — an air cooler cannot reach it."""


@register("AirCooler")
class AirCooler(UnitOp):
    """Air-cooled exchanger ("fin-fan"): cool a stream to ``params['T_out']``
    (required), rejecting the heat to ambient air instead of cooling water.

    Parameters:
      * ``T_out`` — outlet temperature, K (required).
      * ``t_air_in`` — ambient (design) air inlet temperature, K; default
        308.15 K (35 C, a common hot-day design ambient; GPSA Engineering Data
        Book, air-cooled exchangers section).
      * ``approach`` — minimum process-outlet-to-air-inlet approach, K; default
        10 K (typical economic minimum for air coolers; Towler & Sinnott,
        *Chemical Engineering Design*, 2e, Ch. 12). ``T_out`` below
        ``t_air_in + approach`` raises :class:`AirCoolerApproachError`.
      * ``dP`` — tube-side pressure drop, Pa; default 0.
      * ``fan_power_frac`` — fan electricity per unit heat rejected, kW/kW;
        default 0.02. (Rule of thumb: ACHE fans draw on the order of 1-3% of
        the rejected duty; cf. GPSA Engineering Data Book, Section 10 fan-power
        examples. A sweepable assumption, used by the economics layer.)

    Duty Q = n·(H_out − H_in) is negative (heat leaves the process) and is
    reported on the energy ``duty`` port, exactly like a Heater in cooling
    service. No cooling-water utility is drawn — the operating cost is the fan
    electricity, booked during sizing from ``fan_power_frac``.
    """

    def define_ports(self) -> list[Port]:
        return [Port("in1", "inlet"), Port("out", "outlet"), Port("duty", "outlet", "energy")]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        inlet = inlets.get("in1")
        if inlet is None or not inlet.molar_flow:
            raise ValueError(f"AirCooler {self.id!r}: missing or empty inlet on 'in1'")

        if self.params.get("T_out") is None:
            raise ValueError(f"AirCooler {self.id!r}: 'T_out' is required")
        T_out = float(self.params["T_out"])
        t_air_in = float(self.params.get("t_air_in", 308.15))
        approach = float(self.params.get("approach", 10.0))
        if T_out < t_air_in + approach:
            raise AirCoolerApproachError(
                f"AirCooler {self.id!r}: T_out={T_out:.2f} K is below the air inlet "
                f"({t_air_in:.2f} K) plus the minimum approach ({approach:.1f} K) = "
                f"{t_air_in + approach:.2f} K — ambient air cannot cool that far. "
                f"Raise T_out, lower t_air_in/approach, or use a refrigerated cooler."
            )

        T_in, P_in, n = inlet.require_state()
        if T_out > T_in:
            raise ValueError(
                f"AirCooler {self.id!r}: T_out={T_out:.2f} K is above the inlet "
                f"temperature {T_in:.2f} K — an air cooler only cools."
            )
        z = inlet.normalized_z()
        P_out = P_in - float(self.params.get("dP", 0.0))
        H_in = inlet.H if inlet.H is not None else pp.enthalpy(T_in, P_in, z)

        res = pp.flash_pt(T_out, P_out, z)
        Q = n * (res.H - H_in)                    # <= 0: heat rejected to air

        out = Stream(
            id=f"{self.id}.out",
            components=list(inlet.components),
            T=res.T, P=res.P, molar_flow=n, z=z,
            H=res.H, phase=res.phase, vapor_fraction=res.vapor_fraction,
        )
        return {"out": out, "duty": EnergyStream(id=f"{self.id}.duty", duty=Q)}
