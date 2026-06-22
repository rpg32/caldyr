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
from . import data, tray_sizing

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
    # Tray-column hydraulic design point: fraction of the Fair flooding
    # velocity (80% is the classic sieve-tray design point; Seader 3e
    # sec. 6.6) and the downcomer share of the tower cross-section.
    tray_flood_fraction: float = 0.8
    tray_downcomer_frac: float = 0.10
    # Liquid-liquid extraction column capacity: combined (both phases)
    # volumetric throughput per unit of tower cross-section. Sieve-tray
    # extractors run ~20-30 m^3/(m^2 h) total (Seader 3e ch. 8, sieve-tray
    # extraction columns; same order in Towler & Sinnott 2e ch. 11) — there
    # is no vapor, so Fair flooding does not apply.
    extractor_capacity_m3_m2_h: float = 25.0


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
@register_sizer("Mixer", "Splitter", "Valve", "Source", "Makeup")
def _negligible_size(unit, ctx: SizerContext) -> list[EquipmentSize]:
    """Mixers, splitters, throttling valves, boundary sources and make-up
    controllers contribute no significant equipment cost (piping/header
    components)."""
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


@register_sizer("ThreePhaseSeparator", "Decanter")
def _three_phase_size(unit, ctx: SizerContext) -> list[EquipmentSize]:
    """Three-phase separators and decanters are horizontal drums (liquid-liquid
    disengagement needs settling length, not height)."""
    return _vessel_size(unit, ctx, ctx.opts.vessel_residence_s,
                        orientation="vessel_horizontal")


@register_sizer("Balance")
def _balance_size(unit, ctx: SizerContext) -> list[EquipmentSize]:
    """Balance is a logical (bookkeeping) operation — no physical equipment."""
    return []


@register_sizer("PipeSegment")
def _pipe_size(unit, ctx: SizerContext) -> list[EquipmentSize]:
    """Straight pipe, costed per installed length with a diameter power law:
    the ``pipe`` correlation in :mod:`caldyr.economics.data` (Sinnott, C&R
    Vol. 6, 4e, sec. 5.5) takes the composite attribute A = L * D^0.74 in
    metres, and the Sinnott figure is already an *installed* cost (Fbm = 1).
    Length and diameter are the pipe's own design parameters — no stream data
    needed beyond the operating pressure."""
    length = float(unit.params["length"])
    diameter = float(unit.params["diameter"])
    outlet = ctx.outs.get("out")
    pbarg = _pa_to_barg(outlet.P) if outlet is not None else 0.0
    return [EquipmentSize(
        unit_id=unit.id, equipment_type="pipe",
        attribute=length * diameter ** 0.74, attribute_name="L_m_x_D_m^0.74",
        pressure_barg=pbarg, material=ctx.opts.material, diameter_m=diameter,
        notes=[f"pipe: L={length:.1f} m, D={diameter * 1000:.1f} mm "
               f"(installed-cost correlation, Fbm=1)"],
    )]


@register_sizer("Evaporator")
def _evaporator_size(unit, ctx: SizerContext) -> list[EquipmentSize]:
    """Evaporator: a vertical vessel (residence-time heuristic, like a flash
    drum) whose heating duty draws the cheapest hot utility able to reach the
    boiling temperature — typically steam, mirroring the steam-heated
    evaporator of Hameed (2025) §5.2. The heating surface (calandria) is not
    itemized separately; capital = vessel, opex = the heating utility."""
    sizes = _vessel_size(unit, ctx, ctx.opts.vessel_residence_s)
    duty = ctx.duty(unit)
    if duty > 0.0:
        temps = [s.T for s in ctx.outs.values()
                 if getattr(s, "T", None) is not None]
        t_boil = max(temps) if temps else next(iter(ctx.ins.values())).T
        util = _select_utility("heat", t_boil)
        vessel = sizes[0]
        vessel.utility = util
        vessel.utility_duty_W = duty
        vessel.notes.append(
            f"evaporator duty {duty / 1e3:.0f} kW heated with {util} "
            f"(boiling at {t_boil:.1f} K)"
        )
    return sizes


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


# -- tray columns (distillation, absorption, stripping) -----------------------
# Design heuristics for tray columns. The tower shell is costed with the
# Turton vertical-vessel ("towers: tray and packed") correlation, the sieve
# trays separately per tray area, and the condenser/reboiler as utility heat
# exchangers — so both duties also show up as utility opex, exactly like a
# Heater. This mirrors the Turton flow: tower + trays + auxiliary exchangers.
# The diameter comes from Fair's flooding correlation at the governing tray
# (see economics.tray_sizing) at opts.tray_flood_fraction (default 80%) of
# flooding, replacing the earlier fixed F-factor (1.2 Pa^0.5) heuristic —
# the two agree within ~10-15% on the benzene/toluene benchmark, but the
# flooding form responds correctly to liquid load and pressure.
TRAY_EFFICIENCY = 0.7        # overall column efficiency (typical aromatics ~0.7)
TRAY_SPACING_M = 0.61        # 24 in tray spacing
TOWER_EXTRA_HEIGHT_M = 3.0   # sump + disengagement allowance
# Stage efficiencies in absorption service are much poorer than in
# distillation: "the efficiency in the absorption is low (<0.4) but that in
# the stripping (desorption) is moderate" (Hameed 2025, sec. 9.1.3 — same
# message as O'Connell's absorber correlation).
ABSORBER_TRAY_EFFICIENCY = 0.40
STRIPPER_TRAY_EFFICIENCY = 0.50
# Liquid-liquid sieve-tray extractors are poorer still: overall stage
# efficiencies of 10-30% are typical (Seader 3e ch. 8.4).
EXTRACTOR_TRAY_EFFICIENCY = 0.25


def _column_loads(design) -> list[tray_sizing.StageLoad]:
    """Candidate (governing) tray loads for a column: the top-most and
    bottom-most *trays* from the converged stage profiles when the model
    publishes them (RigorousColumn, Absorber, ReboiledAbsorber — the
    condenser/reboiler end stages are skipped where present), else the
    constant-molal-overflow section loads of the FUG design (ShortcutColumn;
    the bottom vapor is approximated as saturated at the bottoms
    composition)."""
    if "V_profile" in design:
        n = design["n_stages"]
        top = 1 if design.get("Q_condenser") is not None else 0
        bottom = n - 2 if design.get("Q_reboiler") is not None else n - 1
        candidates = sorted({max(top, 0), max(bottom, top, 0)})
        return tray_sizing.loads_from_profiles(design, candidates)
    d_rate, b_rate = design["D"], design["B"]
    f = d_rate + b_rate
    v_top = design["V_top"]
    v_bot = max(v_top - (1.0 - design["q"]) * f, 0.05 * f)
    p = design["P"]
    return [
        tray_sizing.StageLoad(stage=2, V=v_top, L=design["R"] * d_rate,
                              T=design["T_top_dew"], P=p,
                              x=design["x_D"], y=design["x_D"]),
        tray_sizing.StageLoad(stage=max(round(design["N"]), 2),
                              V=v_bot, L=v_bot + b_rate,
                              T=design["T_bottom"], P=p,
                              x=design["x_B"], y=design["x_B"]),
    ]


def _tower_and_trays(unit_id: str, design, pp, opts: SizingOptions,
                     efficiency: float) -> list[EquipmentSize]:
    """Tower shell + sieve trays, diameter from Fair flooding at the
    governing tray. Shared by every tray column (distillation, absorber,
    reboiled absorber)."""
    n_real = math.ceil(design["N"] / efficiency)
    height = n_real * TRAY_SPACING_M + TOWER_EXTRA_HEIGHT_M
    hyd = tray_sizing.governing_tray(
        pp, _column_loads(design),
        tray_spacing_m=TRAY_SPACING_M,
        flood_fraction=opts.tray_flood_fraction,
        downcomer_frac=opts.tray_downcomer_frac,
    )
    diameter, area = hyd.diameter_m, hyd.area_m2
    volume = area * height
    pbarg = _pa_to_barg(design["P"])

    tower = EquipmentSize(
        unit_id=unit_id, equipment_type="vessel_vertical",
        attribute=volume, attribute_name="volume_m3",
        pressure_barg=pbarg, material=opts.material, diameter_m=diameter,
        notes=[f"tower: N={design['N']:.1f} ideal stages / {efficiency} eff "
               f"= {n_real} trays; H={height:.1f} m, D={diameter:.2f} m",
               *hyd.notes],
    )
    trays = EquipmentSize(
        unit_id=f"{unit_id}.trays", equipment_type="tray_sieve",
        attribute=area, attribute_name="area_m2",
        pressure_barg=pbarg, material=opts.material, diameter_m=diameter,
        quantity=n_real,
        notes=[f"{n_real} sieve trays of {area:.2f} m^2"],
    )
    return [tower, trays]


def _reboiler_size(unit_id: str, design, opts: SizingOptions,
                   pbarg: float) -> EquipmentSize:
    """Reboiler as a utility heat exchanger boiling the bottoms against a
    (roughly isothermal) hot utility."""
    q_reb = abs(design["Q_reboiler"])
    util_r = _select_utility("heat", design["T_bottom"])
    ur = data.UTILITIES[util_r]
    lmtd_r = max(ur.T_supply - design["T_bottom"], 1.0)
    return EquipmentSize(
        unit_id=f"{unit_id}.reboiler", equipment_type="heat_exchanger",
        attribute=q_reb / (ur.U * lmtd_r), attribute_name="area_m2",
        pressure_barg=pbarg, material=opts.material,
        utility=util_r, utility_duty_W=q_reb,
        notes=[f"reboiler: heat with {util_r}; LMTD={lmtd_r:.1f} K, "
               f"U={ur.U} W/m^2/K"],
    )


def _require_design(unit) -> dict:
    design = getattr(unit, "design", None)
    if not design:
        raise ValueError(
            f"{type(unit).__name__} {unit.id!r} has no design results; solve "
            f"the flowsheet before sizing"
        )
    return design


@register_sizer("ShortcutColumn", "RigorousColumn")
def _column_sizes(unit, ctx: SizerContext) -> list[EquipmentSize]:
    """Tower shell + sieve trays + condenser + reboiler for a distillation
    column. Reads the design results from ``unit.design`` (set by solve) —
    both ShortcutColumn (FUG) and RigorousColumn (MESH) publish the same keys
    (``N`` is the equilibrium-stage count excluding a total condenser), so one
    sizer covers both."""
    pp, opts = ctx.pp, ctx.opts
    design = _require_design(unit)
    pbarg = _pa_to_barg(design["P"])
    tower, trays = _tower_and_trays(unit.id, design, pp, opts, TRAY_EFFICIENCY)

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
    reboiler = _reboiler_size(unit.id, design, opts, pbarg)
    return [tower, trays, condenser, reboiler]


@register_sizer("Absorber")
def _absorber_sizes(unit, ctx: SizerContext) -> list[EquipmentSize]:
    """Gas absorber / stripper: tower shell + sieve trays only (no condenser
    or reboiler). The diameter comes from Fair flooding at the governing tray
    of the converged stage profiles; the tray count uses the (low) absorption
    stage efficiency."""
    design = _require_design(unit)
    return _tower_and_trays(unit.id, design, ctx.pp, ctx.opts,
                            ABSORBER_TRAY_EFFICIENCY)


@register_sizer("ReboiledAbsorber")
def _reboiled_absorber_sizes(unit, ctx: SizerContext) -> list[EquipmentSize]:
    """Reboiled absorber (stripping tower): tower shell + sieve trays + the
    reboiler exchanger (no condenser). ``design['N']`` already excludes the
    reboiler stage; stripping-service tray efficiency."""
    design = _require_design(unit)
    sizes = _tower_and_trays(unit.id, design, ctx.pp, ctx.opts,
                             STRIPPER_TRAY_EFFICIENCY)
    sizes.append(_reboiler_size(unit.id, design, ctx.opts,
                                _pa_to_barg(design["P"])))
    return sizes


@register_sizer("ExtractionColumn")
def _extraction_column_sizes(unit, ctx: SizerContext) -> list[EquipmentSize]:
    """Liquid-liquid extraction column: tower shell + sieve trays (no
    condenser, reboiler or vapor traffic). There is no vapor, so the Fair
    flooding correlation does not apply; the diameter comes from the combined
    liquid throughput of both phases at the design capacity
    ``opts.extractor_capacity_m3_m2_h`` (sieve-tray extractors run ~20-30
    m^3/(m^2 h) total; Seader 3e ch. 8). Tray count from the (low)
    liquid-liquid stage efficiency ``EXTRACTOR_TRAY_EFFICIENCY``."""
    opts = ctx.opts
    design = _require_design(unit)
    t, p = design["T"], design["P"]
    vol_flow = (design["extract_total"]
                * ctx.pp.volume_liquid(t, p, design["x_extract"])
                + design["raffinate_total"]
                * ctx.pp.volume_liquid(t, p, design["x_raffinate"]))
    area = vol_flow / (opts.extractor_capacity_m3_m2_h / 3600.0)
    diameter = math.sqrt(4.0 * area / math.pi)
    n_real = math.ceil(design["N"] / EXTRACTOR_TRAY_EFFICIENCY)
    height = n_real * TRAY_SPACING_M + TOWER_EXTRA_HEIGHT_M
    pbarg = _pa_to_barg(p)
    tower = EquipmentSize(
        unit_id=unit.id, equipment_type="vessel_vertical",
        attribute=area * height, attribute_name="volume_m3",
        pressure_barg=pbarg, material=opts.material, diameter_m=diameter,
        notes=[f"extraction tower: N={design['N']:.1f} ideal stages / "
               f"{EXTRACTOR_TRAY_EFFICIENCY} eff = {n_real} trays; "
               f"H={height:.1f} m, D={diameter:.2f} m from "
               f"{vol_flow * 3600.0:.2f} m^3/h combined liquid at "
               f"{opts.extractor_capacity_m3_m2_h} m^3/(m^2 h)"],
    )
    trays = EquipmentSize(
        unit_id=f"{unit.id}.trays", equipment_type="tray_sieve",
        attribute=area, attribute_name="area_m2",
        pressure_barg=pbarg, material=opts.material, diameter_m=diameter,
        quantity=n_real,
        notes=[f"{n_real} sieve trays of {area:.2f} m^2"],
    )
    return [tower, trays]


# -- solids operations (Hameed 2025 ch. 12) -----------------------------------
@register_sizer("Cyclone")
def _cyclone_size(unit, ctx: SizerContext) -> list[EquipmentSize]:
    """Gas cyclone: costed per cyclone on the actual gas volumetric flow each
    parallel unit handles (the Couper/Walas heavy-duty correlation in
    economics.data — order-of-magnitude, see the source note there). The
    parallel count rides on ``quantity``."""
    design = _require_design(unit)
    inlet = ctx.ins.get("gas_in")
    pbarg = _pa_to_barg(inlet.P) if inlet is not None else 0.0
    n = int(design["n_cyclones"])
    return [EquipmentSize(
        unit_id=unit.id, equipment_type="cyclone",
        attribute=design["Q_per_cyclone_m3_s"], attribute_name="gas_m3_s",
        pressure_barg=pbarg, material=ctx.opts.material,
        diameter_m=design["body_diameter_m"], quantity=n,
        notes=[f"{n} x {design['geometry']} cyclone(s), "
               f"D={design['body_diameter_m']:.2f} m, "
               f"{design['Q_per_cyclone_m3_s']:.2f} m^3/s each, "
               f"dP={design['dP_Pa'] / 1e3:.2f} kPa "
               f"(order-of-magnitude cost correlation)"],
    )]


@register_sizer("RotaryVacuumFilter")
def _rotary_vacuum_filter_size(unit, ctx: SizerContext) -> list[EquipmentSize]:
    """Rotary vacuum drum filter: costed on the filtration area the unit's
    cake-filtration design computed (order-of-magnitude correlation in
    economics.data — see the source note there)."""
    design = _require_design(unit)
    inlet = ctx.ins.get("slurry_in")
    pbarg = _pa_to_barg(inlet.P) if inlet is not None else 0.0
    return [EquipmentSize(
        unit_id=unit.id, equipment_type="rotary_vacuum_filter",
        attribute=design["area_m2"], attribute_name="area_m2",
        pressure_barg=pbarg, material=ctx.opts.material,
        diameter_m=2.0 * design["drum_radius_m"],
        notes=[f"rotary drum filter: A={design['area_m2']:.2f} m^2 "
               f"(R={design['drum_radius_m']:.2f} m x "
               f"W={design['drum_width_m']:.2f} m), "
               f"dP={design['pressure_drop_Pa'] / 1e3:.1f} kPa "
               f"(order-of-magnitude cost correlation)"],
    )]


@register_sizer("BaghouseFilter")
def _baghouse_size(unit, ctx: SizerContext) -> list[EquipmentSize]:
    """Baghouse: costed on the gross cloth area (A = Q / face velocity; EPA
    Cost Manual correlation in economics.data)."""
    design = _require_design(unit)
    inlet = ctx.ins.get("gas_in")
    pbarg = _pa_to_barg(inlet.P) if inlet is not None else 0.0
    return [EquipmentSize(
        unit_id=unit.id, equipment_type="baghouse",
        attribute=design["cloth_area_m2"], attribute_name="area_m2",
        pressure_barg=pbarg, material=ctx.opts.material,
        notes=[f"baghouse: {design['cloth_area_m2']:.0f} m^2 cloth at "
               f"{design['face_velocity_m_s'] * 1e2:.1f} cm/s air-to-cloth, "
               f"~{design['n_bags']} bags"],
    )]


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
