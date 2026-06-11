"""Size equipment from a *solved* flowsheet.

Each sizer reads the resolved streams and duties around a unit and returns an
:class:`EquipmentSize` carrying the capacity attribute the costing correlations
need (area, power, volume, or duty) plus the utility it draws (for operating
cost).

Sizers are looked up in a **registry** keyed by the unit-type name (the unit's
class name): decorate a function with ``@register_sizer("MyUnit")`` to make
``size_flowsheet`` handle it. Every sizer has the same signature
``(unit, ctx: SizerContext) -> list[EquipmentSize]`` — an empty list means the
unit contributes no significant equipment cost (mixers, splitters).

Design heuristics (overall U for heater/cooler service, vessel residence times,
L/D) are explicit options with documented defaults — they are engineering
assumptions, not physics, so they live here in the economics layer rather than in
the solve units.
"""
from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

from ..unitops.heat_exchanger import HeatExchanger
from . import data

_MIN_APPROACH = 10.0   # K, minimum temperature approach when picking a utility


@dataclass
class EquipmentSize:
    unit_id: str
    equipment_type: str          # key into economics.data tables
    attribute: float             # capacity attribute (area_m2 / power_kW / volume_m3)
    attribute_name: str
    pressure_barg: float = 0.0
    material: str = "CS"
    diameter_m: float = 1.0      # for vessel pressure factor
    utility: str | None = None   # utility drawn (for opex), e.g. "fired_heat"
    utility_duty_W: float = 0.0  # magnitude of the utility duty, W
    quantity: int = 1            # stacked identical items (column trays)
    notes: list[str] = field(default_factory=list)


@dataclass
class SizingOptions:
    material: str = "CS"
    hx_overall_U: float = 500.0          # W/m^2/K, process-process exchanger
    vessel_residence_s: float = 300.0    # separator liquid+vapor holdup
    reactor_residence_s: float = 15.0    # gas-phase reactor space-time
    ld_ratio: float = 3.0                # vessel length/diameter
    # Air-cooled exchanger overall U on a bare-tube-equivalent area basis.
    # Gas/light-duty process service runs ~30-60 W/m^2/K bare-tube (GPSA
    # Engineering Data Book, 13e, Section 10, air-cooled exchangers; Towler &
    # Sinnott, Chemical Engineering Design, 2e, Ch. 12 gives similar ranges).
    air_cooler_U: float = 40.0           # W/m^2/K, bare-tube equivalent
    # Design air-side temperature rise across the tube bundle (typical ACHE
    # design rise of ~10-15 K; GPSA Section 10).
    air_cooler_air_rise: float = 10.0    # K


def _pa_to_barg(p_pa: float) -> float:
    return p_pa / 1e5 - 1.01325


def _streams_around(fs, uid):
    ins, outs = {}, {}
    for conn in fs.connections:
        if conn.to_unit == uid and conn.stream_id in fs.streams:
            ins[conn.to_port] = fs.streams[conn.stream_id]
        if conn.from_unit == uid and conn.stream_id in fs.streams:
            outs[conn.from_port] = fs.streams[conn.stream_id]
    return ins, outs


def energy_duty(fs, report, unit) -> float:
    """Signed duty (W) summed over the unit's energy ports. Works whether or not
    the duty port is wired to a boundary: a connected port uses its stream id, an
    unconnected one falls back to the engine's default id (``unit.port``)."""
    total = 0.0
    for p in unit.ports:
        if p.kind != "energy":
            continue
        conn = next((c for c in fs.connections
                     if c.from_unit == unit.id and c.from_port == p.name), None)
        sid = conn.stream_id if conn else f"{unit.id}.{p.name}"
        total += report.duties.get(sid, 0.0)
    return total


# Backwards-compatible private alias (pre-registry name).
_energy_duty = energy_duty


def _select_utility(kind: str, t_target: float) -> str:
    """Cheapest utility able to heat/cool a process stream to ``t_target``."""
    candidates = []
    for name, u in data.UTILITIES.items():
        if u.kind != kind:
            continue
        if kind == "heat" and u.T_supply > t_target + _MIN_APPROACH:
            candidates.append((u.price_per_GJ, name))
        elif kind == "cool" and u.T_supply < t_target - _MIN_APPROACH:
            candidates.append((u.price_per_GJ, name))
    if not candidates:
        raise ValueError(f"no {kind} utility can reach {t_target:.1f} K")
    return min(candidates)[1]


# -- the sizer registry ------------------------------------------------------
@dataclass
class SizerContext:
    """Everything a sizer may need around one unit of a solved flowsheet."""
    fs: Any
    report: Any
    pp: Any
    opts: SizingOptions
    ins: dict[str, Any]      # inlet material streams by port name
    outs: dict[str, Any]     # outlet material streams by port name

    def duty(self, unit) -> float:
        """Signed duty (W) over the unit's energy ports (wired or not)."""
        return energy_duty(self.fs, self.report, unit)


Sizer = Callable[[Any, SizerContext], list[EquipmentSize]]
SIZER_REGISTRY: dict[str, Sizer] = {}


def register_sizer(*unit_types: str) -> Callable[[Sizer], Sizer]:
    """Register a sizer for one or more unit-type names (unit class names)."""
    def deco(fn: Sizer) -> Sizer:
        for name in unit_types:
            SIZER_REGISTRY[name] = fn
        return fn
    return deco


# -- sizers ------------------------------------------------------------------
@register_sizer("Mixer", "Splitter")
def _negligible_size(unit, ctx: SizerContext) -> list[EquipmentSize]:
    """Mixers and splitters contribute no significant equipment cost (just
    piping/headers)."""
    return []


@register_sizer("Heater")
def _heater_size(unit, ctx: SizerContext) -> list[EquipmentSize]:
    inlet = ctx.ins["in1"]
    outlet = ctx.outs["out"]
    duty = ctx.duty(unit)
    q = abs(duty)
    kind = "heat" if duty > 0 else "cool"
    target = outlet.T
    util = _select_utility(kind, target)
    u = data.UTILITIES[util]

    # LMTD between the process stream and a (often isothermal) utility.
    if kind == "heat":
        dt1, dt2 = u.T_supply - outlet.T, u.T_supply - inlet.T
    else:
        dt1, dt2 = inlet.T - u.T_supply, outlet.T - u.T_supply
    lmtd = HeatExchanger.lmtd(max(dt1, 1.0), max(dt2, 1.0))
    area = q / (u.U * lmtd)
    return [EquipmentSize(
        unit_id=unit.id, equipment_type="heat_exchanger",
        attribute=area, attribute_name="area_m2",
        pressure_barg=_pa_to_barg(outlet.P), material=ctx.opts.material,
        utility=util, utility_duty_W=q,
        notes=[f"{kind} with {util}; LMTD={lmtd:.1f} K, U={u.U} W/m^2/K"],
    )]


@register_sizer("FiredHeater")
def _fired_heater_size(unit, ctx: SizerContext) -> list[EquipmentSize]:
    """Fired heater: the Turton correlation capacity is the absorbed (process)
    duty in kW; the *fuel* duty (process / fired efficiency, from the unit's
    ``design`` results) is what the fired-heat utility bills for."""
    outlet = ctx.outs["out"]
    q = ctx.duty(unit)                      # process duty, W (> 0)
    design = getattr(unit, "design", None) or {}
    eta = float(design.get("efficiency", unit.params.get("efficiency", 0.85)))
    fuel = float(design.get("fuel_duty", q / eta))
    return [EquipmentSize(
        unit_id=unit.id, equipment_type="fired_heater",
        attribute=q / 1e3, attribute_name="duty_kW",
        pressure_barg=_pa_to_barg(outlet.P), material=ctx.opts.material,
        utility="fired_heat", utility_duty_W=fuel,
        notes=[f"fired heater: process duty {q / 1e3:.0f} kW at efficiency "
               f"{eta:.2f} -> fuel duty {fuel / 1e3:.0f} kW"],
    )]


@register_sizer("AirCooler")
def _air_cooler_size(unit, ctx: SizerContext) -> list[EquipmentSize]:
    """Air-cooled exchanger: bare-tube-equivalent area from duty/(U·LMTD) vs
    ambient air with a fixed air-side rise; the only utility drawn is fan
    electricity at ``fan_power_frac`` kW per kW rejected (see AirCooler docs)."""
    inlet = ctx.ins["in1"]
    outlet = ctx.outs["out"]
    q = abs(ctx.duty(unit))
    t_air_in = float(unit.params.get("t_air_in", 308.15))
    t_air_out = t_air_in + ctx.opts.air_cooler_air_rise
    lmtd = HeatExchanger.lmtd(max(inlet.T - t_air_out, 1.0),
                              max(outlet.T - t_air_in, 1.0))
    area = q / (ctx.opts.air_cooler_U * lmtd)
    fan_frac = float(unit.params.get("fan_power_frac", 0.02))
    return [EquipmentSize(
        unit_id=unit.id, equipment_type="air_cooler",
        attribute=area, attribute_name="area_m2",
        pressure_barg=_pa_to_barg(outlet.P), material=ctx.opts.material,
        utility="electricity", utility_duty_W=fan_frac * q,
        notes=[f"air cooler: LMTD={lmtd:.1f} K vs air at {t_air_in:.1f} K "
               f"(+{ctx.opts.air_cooler_air_rise:.0f} K rise), "
               f"U={ctx.opts.air_cooler_U} W/m^2/K bare-tube; "
               f"fan power {fan_frac * 100:.1f}% of duty"],
    )]


@register_sizer("HeatExchanger")
def _hx_size(unit, ctx: SizerContext) -> list[EquipmentSize]:
    hot_in, cold_out = ctx.ins["hot_in"], ctx.outs["cold_out"]
    hot_out, cold_in = ctx.outs["hot_out"], ctx.ins["cold_in"]
    q = hot_in.molar_flow * (hot_in.H - hot_out.H)
    lmtd = HeatExchanger.lmtd(max(hot_in.T - cold_out.T, 1.0),
                              max(hot_out.T - cold_in.T, 1.0))
    area = abs(q) / (ctx.opts.hx_overall_U * lmtd)
    return [EquipmentSize(
        unit_id=unit.id, equipment_type="heat_exchanger",
        attribute=area, attribute_name="area_m2",
        pressure_barg=_pa_to_barg(min(hot_in.P, cold_in.P)), material=ctx.opts.material,
        notes=[f"process-process; LMTD={lmtd:.1f} K, U={ctx.opts.hx_overall_U} W/m^2/K"],
    )]


def _rotating_size(unit, ctx: SizerContext, eq_type: str) -> list[EquipmentSize]:
    duty = ctx.duty(unit)
    power_kW = abs(duty) / 1000.0
    return [EquipmentSize(
        unit_id=unit.id, equipment_type=eq_type,
        attribute=power_kW, attribute_name="power_kW",
        pressure_barg=_pa_to_barg(ctx.outs["out"].P), material=ctx.opts.material,
        utility="electricity", utility_duty_W=abs(duty),
    )]


@register_sizer("Pump")
def _pump_size(unit, ctx: SizerContext) -> list[EquipmentSize]:
    return _rotating_size(unit, ctx, "pump_centrifugal")


@register_sizer("Compressor")
def _compressor_size(unit, ctx: SizerContext) -> list[EquipmentSize]:
    return _rotating_size(unit, ctx, "compressor_centrifugal")


@register_sizer("Expander")
def _turbine_size(unit, ctx: SizerContext) -> list[EquipmentSize]:
    """Expander/turbine: sized on the magnitude of the (negative, extracted)
    shaft work. No utility is *drawn* — the power produced is a potential
    credit, which the opex model does not auto-book (flagged in the notes)."""
    power_kW = abs(ctx.duty(unit)) / 1000.0
    return [EquipmentSize(
        unit_id=unit.id, equipment_type="turbine_axial",
        attribute=power_kW, attribute_name="power_kW",
        pressure_barg=_pa_to_barg(ctx.outs["out"].P), material=ctx.opts.material,
        notes=[f"expander: {power_kW:.1f} kW shaft power extracted "
               f"(power credit not auto-booked in opex)"],
    )]


def _vessel_size(unit, ctx: SizerContext, residence_s: float,
                 orientation: str = "vessel_vertical") -> list[EquipmentSize]:
    inlet = next(iter(ctx.ins.values()))
    vol_flow = inlet.molar_flow * ctx.pp.volume(inlet.T, inlet.P, inlet.normalized_z())
    volume = vol_flow * residence_s
    diameter = (4.0 * volume / (math.pi * ctx.opts.ld_ratio)) ** (1.0 / 3.0)
    return [EquipmentSize(
        unit_id=unit.id, equipment_type=orientation,
        attribute=volume, attribute_name="volume_m3",
        pressure_barg=_pa_to_barg(inlet.P), material=ctx.opts.material,
        diameter_m=diameter,
        notes=[f"V = {vol_flow:.3g} m^3/s x {residence_s:.0f} s; D={diameter:.2f} m"],
    )]


@register_sizer("FlashDrum", "ComponentSplitter")
def _separator_size(unit, ctx: SizerContext) -> list[EquipmentSize]:
    """Flash drums, and black-box component splitters costed as a vertical
    vessel by residence time, like a flash drum (a deliberate placeholder
    heuristic for the latter)."""
    return _vessel_size(unit, ctx, ctx.opts.vessel_residence_s)


@register_sizer("ThreePhaseSeparator")
def _three_phase_size(unit, ctx: SizerContext) -> list[EquipmentSize]:
    """Three-phase separators are horizontal drums (liquid-liquid disengagement
    needs settling length, not height)."""
    return _vessel_size(unit, ctx, ctx.opts.vessel_residence_s,
                        orientation="vessel_horizontal")


@register_sizer("ConversionReactor", "EquilibriumReactor", "GibbsReactor")
def _reactor_size(unit, ctx: SizerContext) -> list[EquipmentSize]:
    return _vessel_size(unit, ctx, ctx.opts.reactor_residence_s)


@register_sizer("CSTR", "PFR")
def _kinetic_reactor_size(unit, ctx: SizerContext) -> list[EquipmentSize]:
    """CSTR/PFR: the reaction volume is a design parameter (``params['V']``),
    so size the vessel on it directly (Turton vertical process vessel — a PFR
    shell is a vessel too, for now), not a residence-time heuristic."""
    inlet = next(iter(ctx.ins.values()))
    volume = float(unit.params["V"])
    diameter = (4.0 * volume / (math.pi * ctx.opts.ld_ratio)) ** (1.0 / 3.0)
    return [EquipmentSize(
        unit_id=unit.id, equipment_type="vessel_vertical",
        attribute=volume, attribute_name="volume_m3",
        pressure_barg=_pa_to_barg(inlet.P), material=ctx.opts.material,
        diameter_m=diameter,
        notes=[f"kinetic reactor: V={volume:.3g} m^3 from params; D={diameter:.2f} m"],
    )]


# -- shortcut distillation column --------------------------------------------
# Design heuristics for the FUG column. The tower shell is costed with the
# Turton vertical-vessel ("towers: tray and packed") correlation, the sieve
# trays separately per tray area, and the condenser/reboiler as utility heat
# exchangers — so both duties also show up as utility opex, exactly like a
# Heater. This mirrors the Turton flow: tower + trays + auxiliary exchangers.
TRAY_EFFICIENCY = 0.7        # overall column efficiency (typical aromatics ~0.7)
TRAY_SPACING_M = 0.61        # 24 in tray spacing
TOWER_EXTRA_HEIGHT_M = 3.0   # sump + disengagement allowance
# Superficial-velocity F-factor u*sqrt(rho_V) ~ 1.2 Pa^0.5 — a typical design
# vapor loading (~80% of flood) for sieve trays at 0.6 m spacing; Kister,
# "Distillation Design" (1992); Towler & Sinnott, "Chemical Engineering
# Design" 2e, Ch. 17. (An F-factor is the Souders-Brown form with the liquid
# density folded into the constant.)
TOWER_F_FACTOR = 1.2         # Pa^0.5


@register_sizer("ShortcutColumn", "RigorousColumn")
def _column_sizes(unit, ctx: SizerContext) -> list[EquipmentSize]:
    """Tower shell + sieve trays + condenser + reboiler for a distillation
    column. Reads the design results from ``unit.design`` (set by solve) —
    both ShortcutColumn (FUG) and RigorousColumn (MESH) publish the same keys
    (``N`` is the equilibrium-stage count excluding a total condenser), so one
    sizer covers both."""
    pp, opts = ctx.pp, ctx.opts
    design = getattr(unit, "design", None)
    if not design:
        raise ValueError(
            f"{type(unit).__name__} {unit.id!r} has no design results; solve "
            f"the flowsheet before sizing"
        )

    n_real = math.ceil(design["N"] / TRAY_EFFICIENCY)
    height = n_real * TRAY_SPACING_M + TOWER_EXTRA_HEIGHT_M
    p = design["P"]
    x_d = design["x_D"]

    # Vapor density of the saturated overhead at the top of the column.
    # (data.molar_mass falls back to the chemicals database for components
    # missing from the local table, so any resolvable component sizes cleanly.)
    mw = sum(frac * data.molar_mass(comp)
             for comp, frac in x_d.items() if frac > 0.0)
    v_vap = pp.volume(design["T_top_dew"], p, x_d)        # m^3/mol, sat. vapor
    rho_v = mw / v_vap                                     # kg/m^3
    u_super = TOWER_F_FACTOR / math.sqrt(rho_v)            # m/s
    q_vap = design["V_top"] * v_vap                        # m^3/s overhead vapor
    area = q_vap / u_super
    diameter = math.sqrt(4.0 * area / math.pi)
    volume = area * height
    pbarg = _pa_to_barg(p)

    tower = EquipmentSize(
        unit_id=unit.id, equipment_type="vessel_vertical",
        attribute=volume, attribute_name="volume_m3",
        pressure_barg=pbarg, material=opts.material, diameter_m=diameter,
        notes=[f"tower: N={design['N']:.1f} ideal stages / {TRAY_EFFICIENCY} eff "
               f"= {n_real} trays; H={height:.1f} m, D={diameter:.2f} m "
               f"(F-factor {TOWER_F_FACTOR} Pa^0.5, rho_V={rho_v:.2f} kg/m^3)"],
    )
    trays = EquipmentSize(
        unit_id=f"{unit.id}.trays", equipment_type="tray_sieve",
        attribute=area, attribute_name="area_m2",
        pressure_barg=pbarg, material=opts.material, diameter_m=diameter,
        quantity=n_real,
        notes=[f"{n_real} sieve trays of {area:.2f} m^2"],
    )

    # Condenser: condensing overhead (T_dew -> T_top) against a cooling utility.
    q_cond = abs(design["Q_condenser"])
    util_c = _select_utility("cool", design["T_top"])
    uc = data.UTILITIES[util_c]
    lmtd_c = HeatExchanger.lmtd(max(design["T_top_dew"] - uc.T_return, 1.0),
                                max(design["T_top"] - uc.T_supply, 1.0))
    condenser = EquipmentSize(
        unit_id=f"{unit.id}.condenser", equipment_type="heat_exchanger",
        attribute=q_cond / (uc.U * lmtd_c), attribute_name="area_m2",
        pressure_barg=pbarg, material=opts.material,
        utility=util_c, utility_duty_W=q_cond,
        notes=[f"condenser: cool with {util_c}; LMTD={lmtd_c:.1f} K, U={uc.U} W/m^2/K"],
    )

    # Reboiler: boiling the bottoms against a (roughly isothermal) hot utility.
    q_reb = abs(design["Q_reboiler"])
    util_r = _select_utility("heat", design["T_bottom"])
    ur = data.UTILITIES[util_r]
    lmtd_r = max(ur.T_supply - design["T_bottom"], 1.0)
    reboiler = EquipmentSize(
        unit_id=f"{unit.id}.reboiler", equipment_type="heat_exchanger",
        attribute=q_reb / (ur.U * lmtd_r), attribute_name="area_m2",
        pressure_barg=pbarg, material=opts.material,
        utility=util_r, utility_duty_W=q_reb,
        notes=[f"reboiler: heat with {util_r}; LMTD={lmtd_r:.1f} K, U={ur.U} W/m^2/K"],
    )
    return [tower, trays, condenser, reboiler]


def size_flowsheet(fs, report, pp, options: SizingOptions | None = None) -> list[EquipmentSize]:
    """Size every cost-bearing unit in a solved flowsheet."""
    opts = options or SizingOptions()
    sizes: list[EquipmentSize] = []
    for uid, unit in fs.units.items():
        kind = type(unit).__name__
        sizer = SIZER_REGISTRY.get(kind)
        if sizer is None:
            raise ValueError(f"no sizer for unit type {kind!r} ({uid})")
        ins, outs = _streams_around(fs, uid)
        ctx = SizerContext(fs=fs, report=report, pp=pp, opts=opts, ins=ins, outs=outs)
        sizes.extend(sizer(unit, ctx))
    return sizes
