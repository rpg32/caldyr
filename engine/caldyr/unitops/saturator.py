"""Stream saturator: bring a gas to a target humidity with a condensable liquid.

The open analogue of HYSYS's **Stream Saturator** (Hameed 2025, *Chemical
Process Simulations using Aspen HYSYS*, sec. 10.4, where it humidifies the acid
gas and combustion-air feeds of a sulfur-recovery unit — "used to quickly
calculate saturation of acid gas and air streams"). A dry (or partly humid) gas
is contacted with a condensable component — water by default — and leaves
carrying that component up to its saturation partial pressure at the operating
temperature and pressure (relative humidity 100 %, or a chosen ``relative_
humidity``).

**Ports**: ``gas_in`` (the gas to saturate), an optional ``water_in`` (the
liquid saturant supply — if absent the unit supplies exactly the amount needed,
reported in ``design['saturant_added']``), ``gas_out`` (the humidified gas),
``liquid_out`` (any unevaporated / condensed liquid saturant) and ``duty`` (the
heat to hold the operating temperature — saturation is endothermic, so an
isothermal saturator absorbs the latent heat; reported honestly like a
ComponentSplitter's duty).

**Model**. The saturation vapor mole fraction ``y*`` of the saturant at
``(T, P)`` is read from a rigorous flash of a saturant-rich probe mixture (the
vapor in equilibrium with liquid saturant). The gas is then loaded with saturant
vapour to ``relative_humidity · y*``:

    n_sat_vap = y_target · n_dry / (1 - y_target),   y_target = RH · y*

where ``n_dry`` is the gas's non-saturant moles (conserved). The saturant added
is ``n_sat_vap - (saturant already in the gas)``; if a ``water_in`` supply is
short of that, the gas leaves **sub-saturated** (all the supply evaporates) and
the achieved humidity is reported. Any surplus supply — or saturant condensed
out of an already-supersaturated feed — leaves as ``liquid_out``.

**Scope (v1)**: isothermal at the gas inlet temperature (or ``T``), with the
saturation duty reported on the ``duty`` port. The saturant is a single
condensable component (``saturant``, default ``"water"``). An *adiabatic*
saturator (solve the outlet temperature for zero duty — the adiabatic-saturation
/ wet-bulb temperature) is the natural next form and is left as a documented
follow-up, not silently dropped.
"""
from __future__ import annotations

from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register


class SaturatorError(ValueError):
    """A saturator could not be configured or solved (missing saturant, a
    saturant that does not condense at the operating conditions, or a bad
    relative-humidity spec). Raised with diagnostics — never a silent answer."""


@register("Saturator")
class Saturator(UnitOp):
    """Gas saturator / humidifier. See the module docstring."""

    def define_ports(self) -> list[Port]:
        return [
            Port("gas_in", "inlet"),
            Port("water_in", "inlet"),          # optional saturant-liquid supply
            Port("gas_out", "outlet"),
            Port("liquid_out", "outlet"),
            Port("duty", "outlet", "energy"),
        ]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        gas = inlets.get("gas_in")
        if gas is None or not gas.molar_flow:
            raise SaturatorError(f"Saturator {self.id!r}: missing/empty 'gas_in'")

        saturant = str(self.params.get("saturant", "water"))
        if saturant not in gas.components:
            raise SaturatorError(
                f"Saturator {self.id!r}: saturant {saturant!r} is not in the "
                f"component list {list(gas.components)}")
        rh = float(self.params.get("relative_humidity", 1.0))
        if not 0.0 < rh <= 1.0:
            raise SaturatorError(
                f"Saturator {self.id!r}: relative_humidity must be in (0, 1], "
                f"got {rh}")

        Tg, Pg, ng = gas.require_state()
        zg = gas.normalized_z()
        comps = list(gas.components)
        T_out = float(self.params.get("T", Tg))
        P_out = float(self.params.get("P", Pg)) - float(self.params.get("dP", 0.0))

        y_sat = self._saturation_y(pp, comps, saturant, T_out, P_out)
        y_target = rh * y_sat
        if y_target >= 1.0:
            raise SaturatorError(
                f"Saturator {self.id!r}: target saturant fraction {y_target:.4g} "
                f">= 1 at {T_out:.1f} K, {P_out:.0f} Pa — the gas cannot hold it")

        # Dry (non-saturant) gas is conserved; size the vapour saturant load.
        n_dry_i = {c: ng * zg.get(c, 0.0) for c in comps if c != saturant}
        n_dry = sum(n_dry_i.values())
        n_sat_in = ng * zg.get(saturant, 0.0)
        n_sat_vap_target = y_target * n_dry / (1.0 - y_target)
        needed = n_sat_vap_target - n_sat_in        # >0 add, <0 must condense

        z_sat_pure = {c: (1.0 if c == saturant else 0.0) for c in comps}

        # Saturant-liquid supply (optional). v1: a pure-saturant liquid stream.
        # When no supply is wired the unit delivers exactly what is needed as a
        # saturated liquid at the operating temperature (its enthalpy is part of
        # the energy balance — saturation absorbs the latent heat).
        water = inlets.get("water_in")
        h_supply = 0.0
        if water is not None and water.molar_flow:
            Tw, Pw, nw = water.require_state()
            zw = water.normalized_z()
            if abs(zw.get(saturant, 0.0) - 1.0) > 1e-6:
                raise SaturatorError(
                    f"Saturator {self.id!r}: 'water_in' must be pure {saturant} "
                    f"(v1), got {zw}")
            available = nw
            h_supply = nw * (water.H if water.H is not None
                             else pp.enthalpy_liquid(Tw, Pw, zw))
        else:
            available = max(0.0, needed)
            h_supply = available * pp.enthalpy_liquid(T_out, P_out, z_sat_pure)
        sub_saturated = False
        if needed > 0.0 and available < needed - 1e-12:
            # Not enough saturant: everything supplied evaporates, gas leaves
            # below the target humidity.
            n_sat_vap = n_sat_in + available
            liquid_out = 0.0
            sub_saturated = True
        else:
            n_sat_vap = n_sat_vap_target
            # surplus supply + any condensed-out saturant leave as liquid
            liquid_out = max(0.0, available - needed)

        # -- assemble products ------------------------------------------------
        n_gas_out = n_dry + n_sat_vap
        z_gas_out = {**{c: n_dry_i[c] / n_gas_out for c in n_dry_i},
                     saturant: n_sat_vap / n_gas_out}
        h_gas_out = pp.enthalpy_vapor(T_out, P_out, z_gas_out)
        gas_out = Stream(
            id=f"{self.id}.gas_out", components=comps, T=T_out, P=P_out,
            molar_flow=n_gas_out, z=z_gas_out, H=h_gas_out,
            phase="vapor", vapor_fraction=1.0)

        h_liq = pp.enthalpy_liquid(T_out, P_out, z_sat_pure) if liquid_out > 0 else 0.0
        liquid = Stream(
            id=f"{self.id}.liquid_out", components=comps, T=T_out, P=P_out,
            molar_flow=liquid_out, z=z_sat_pure, H=h_liq,
            phase="liquid", vapor_fraction=0.0)

        # Duty closes the energy balance (Heater sign: + heats the process). The
        # saturant supply enthalpy (wired or auto) is part of the inlet side.
        h_gas_in = ng * (gas.H if gas.H is not None else pp.enthalpy(Tg, Pg, zg))
        duty = (n_gas_out * h_gas_out + liquid_out * h_liq) - (h_gas_in + h_supply)

        achieved_rh = (n_sat_vap / n_gas_out) / y_sat if y_sat > 0 else 0.0
        self.design = {
            "saturant": saturant,
            "y_saturation": y_sat,
            "relative_humidity_target": rh,
            "relative_humidity_achieved": achieved_rh,
            "saturant_added": max(0.0, needed) if water is None else min(available, max(0.0, needed)),
            "saturant_condensed": max(0.0, -needed),
            "liquid_out": liquid_out,
            "sub_saturated": sub_saturated,
            "duty": duty,
            "T": T_out, "P": P_out,
        }
        return {
            "gas_out": gas_out,
            "liquid_out": liquid,
            "duty": EnergyStream(id=f"{self.id}.duty", duty=duty),
        }

    def _saturation_y(self, pp, comps, saturant, T, P) -> float:
        """Saturant vapour mole fraction in equilibrium with liquid saturant at
        ``(T, P)`` — read from a flash of a saturant-rich probe mixture so the
        result folds in mixture non-ideality, not just the pure vapour pressure.
        """
        # The saturant can only saturate the gas if it can form a liquid at
        # (T, P): pure saturant must NOT be all-vapour there (P > its vapour
        # pressure). This pre-check also distinguishes the saturant condensing
        # from some *other* component condensing in the probe below.
        z_pure = {c: (1.0 if c == saturant else 0.0) for c in comps}
        try:
            pure = pp.flash_pt(T, P, z_pure)
        except Exception as exc:                       # noqa: BLE001 - re-typed
            raise SaturatorError(
                f"Saturator {self.id!r}: could not flash pure {saturant!r} at "
                f"{T:.1f} K, {P:.0f} Pa ({type(exc).__name__})") from exc
        if pure.vapor_fraction >= 1.0 - 1e-9:
            raise SaturatorError(
                f"Saturator {self.id!r}: {saturant!r} does not condense at "
                f"{T:.1f} K, {P:.0f} Pa (its vapour pressure exceeds the operating "
                f"pressure), so the gas cannot be saturated with it here")

        others = [c for c in comps if c != saturant]
        # 70 % saturant guarantees a saturant-rich liquid phase; the remaining
        # 30 % spreads over the other gases so y* folds in mixture non-ideality.
        z = {c: 0.0 for c in comps}
        z[saturant] = 0.7
        for c in others:
            z[c] = 0.3 / len(others)
        if not others:
            z[saturant] = 1.0
        try:
            res = pp.flash_pt(T, P, z)
        except Exception as exc:                       # noqa: BLE001 - re-typed
            raise SaturatorError(
                f"Saturator {self.id!r}: could not flash the saturation probe for "
                f"{saturant!r} at {T:.1f} K, {P:.0f} Pa ({type(exc).__name__})"
            ) from exc
        if res.vapor_fraction >= 1.0 - 1e-9 or res.y is None:
            # No liquid even with 70 % saturant -> treat pure-saturant Psat/P as
            # the saturation fraction (Raoult limit).
            return min(pure.P / P if pure.P else 0.0, 1.0)
        return res.y.get(saturant, 0.0)
