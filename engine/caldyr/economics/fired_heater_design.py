"""Fired-heater (process furnace) thermal design: the radiant / convective split.

The :class:`~caldyr.unitops.fired_heater.FiredHeater` unit op solves the *process*
side as a black-box heater (bring the feed to ``T_out`` at a fired efficiency
``eta`` — fuel duty ``Q/eta``). This module adds the **fired side**: the actual
combustion (fuel and air molar flows, flue-gas composition and temperature) and
the division of the absorbed duty into a **radiant** and a **convective** zone —
the design split the book (Hameed 2025, *Chemical Process Simulations using Aspen
HYSYS*, §4.3) treats with Aspen EDR, and the classic API/Lobo-Evans furnace
design (Wimpress; Perry 8e §27; Towler & Sinnott 2e §19.16).

A direct-fired heater has, gas-side, three zones (Hameed Fig. 4.13): the
**radiant** firebox (flame radiates to the tubes, Stefan-Boltzmann, Hameed
eq. 4.8), the **convective** bank (hot flue gas convects to the tubes, ``Q = U A
ΔT_lm``, eq. 4.9), and an optional economizer. Process fluid flows
counter-current to the flue gas: it enters cold at the top of the convective
bank, is preheated there, then takes the bulk of its duty in the radiant section;
the flue gas leaves the firebox at the **bridgewall temperature** and is cooled
further in the convective bank before going up the stack.

Energy bookkeeping (all enthalpies formation-inclusive, on the NASA ideal-gas
basis of :class:`~caldyr.thermo.nasa_pkg.NasaGasPackage`, so the heat of
combustion falls straight out of a reactant-vs-product enthalpy difference):

* The user fixes the *process* absorbed duty ``Q_abs`` and the LHV-basis fired
  efficiency ``eta`` (HYSYS's definition), so the firing rate on a lower-heating-
  value basis is ``Q_fired = Q_abs / eta`` and the fuel molar flow is
  ``Q_fired / LHV_mix``.
* Stoichiometric O2 and the excess-air specification fix the air flow; the flue
  gas is ``{CO2, H2O, O2_excess, N2}``.
* ``Q_avail`` is the *absolute* heat liberated, i.e. the flue-gas sensible
  enthalpy at the adiabatic flame temperature above the reference temperature;
  it equals ``n_fuel·LHV`` plus any reactant pre-heat (fuel/air entering above
  the reference). The adiabatic flame temperature ``T_ad`` solves
  ``h_sens(T_ad) = Q_avail``.
* A small fixed casing/radiation loss (``loss_fraction·Q_avail``, default 1.5 %,
  Towler & Sinnott) is removed in the hot radiant zone.
* The **bridgewall temperature** ``T_bw`` (the flue temperature leaving the
  firebox) is the design knob that sets the split: the radiant tubes absorb
  ``Q_rad = Q_avail - h_sens(T_bw) - casing_loss``; the convective bank takes the
  remainder ``Q_conv = Q_abs - Q_rad``; the **stack temperature** then solves
  ``h_sens(T_stack) = h_sens(T_bw) - Q_conv``. By construction the balance
  ``Q_avail = Q_abs + casing_loss + stack_loss`` closes exactly.

Areas: the radiant section is sized from an average **radiant heat flux** (design
fluxes are 25-40 kW/m^2 of *cold-plane* tube surface; default 32 kW/m^2, Perry /
Towler), ``A_rad = Q_rad / flux``. The convective bank is sized from an overall
``U`` and the counter-current log-mean ΔT between the flue gas (``T_bw ->
T_stack``) and the process fluid over the fraction of its temperature rise taken
in that bank, ``A_conv = Q_conv / (U_conv·ΔT_lm)`` (gas-side convective ``U`` is
low, ~20-60 W/m^2K; default 40).

**Validation.** The Hameed §4.3.4 worked problem — a 1250 kmol/h C1-C5 feed
heated 50->300 °C with methane fuel at 100 % excess air and 60 % efficiency —
gives a fuel rate of 62.01 kmol/h (book Fig. 4.19). This model returns ~66
kmol/h: the ~6 % gap is the process-duty difference between our PR enthalpy and
HYSYS's (our absorbed duty is 8.82 MW), not the combustion model — fed the
book's implied duty the fuel rate matches. See ``tests/test_m23_fired_heater_
design.py``.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

# Reference temperature for heating values and "loss" bookkeeping (25 C).
_T_REF = 298.15
# Air composition (dry), mole fraction.
_O2_IN_AIR = 0.21
_N2_IN_AIR = 0.79


@dataclass
class FuelCombustion:
    """Combustion mass balance for the fuel burnt in the firebox."""

    lhv_mix: float                 # lower heating value of the fuel, J/mol-fuel
    fuel_flow: float               # total fuel molar flow, mol/s
    fuel_flows: dict[str, float]   # per-component fuel molar flow, mol/s
    air_flow: float                # combustion-air molar flow, mol/s
    o2_stoich: float               # stoichiometric O2 demand, mol/s
    excess_air: float              # excess-air fraction (0.15 = 15 %)
    flue_flow: float               # total flue-gas molar flow, mol/s
    flue_composition: dict[str, float]   # flue mole fractions {CO2,H2O,O2,N2}
    flue_flows: dict[str, float]   # flue component molar flows, mol/s


@dataclass
class FiredHeaterDesign:
    """Radiant/convective thermal design of a fired heater."""

    combustion: FuelCombustion
    process_duty: float            # absorbed (process) duty, W
    fired_duty: float              # LHV-basis firing rate Q/eta, W
    heat_available: float          # absolute heat liberated to the flue, W
    efficiency_lhv: float          # the input LHV-basis efficiency
    efficiency_gross: float        # Q_abs / Q_available (absolute basis)
    flame_temperature: float       # adiabatic flame temperature, K
    bridgewall_temperature: float  # flue temp leaving the firebox, K
    stack_temperature: float       # flue temp up the stack, K
    radiant_duty: float            # heat to the radiant tubes, W
    convective_duty: float         # heat to the convective bank, W
    radiant_fraction: float        # Q_rad / Q_abs
    casing_loss: float             # casing/radiation loss, W
    stack_loss: float              # sensible heat lost up the stack, W
    radiant_area: float            # radiant cold-plane area, m^2
    convective_area: float         # convective bare-tube area, m^2
    radiant_flux: float            # design radiant flux used, W/m^2
    convective_U: float            # convective overall U used, W/m^2K
    convective_lmtd: float         # convective log-mean ΔT, K
    notes: list[str] = field(default_factory=list)


def _atoms(component_id: str) -> dict[str, int]:
    """Element-count of a component from its chemical formula."""
    from chemicals.elements import nested_formula_parser

    from ..core.components_db import resolve_component

    formula = resolve_component(component_id).formula
    if not formula:
        raise ValueError(
            f"no chemical formula for {component_id!r}; cannot compute its "
            f"combustion stoichiometry for fired-heater design"
        )
    return {k: int(v) for k, v in nested_formula_parser(formula).items()}


def _hfg(component_id: str) -> float:
    """Standard gas-phase enthalpy of formation, J/mol (298.15 K)."""
    from chemicals import Hfg
    from chemicals.identifiers import CAS_from_any

    cas = CAS_from_any(component_id)
    h = Hfg(cas)
    if h is None:
        raise ValueError(
            f"no gas-phase formation enthalpy in `chemicals` for {component_id!r}; "
            f"cannot compute its lower heating value for fired-heater combustion"
        )
    return float(h)


# Standard gas-phase formation enthalpies of the combustion products (J/mol,
# 298.15 K) — pulled once from `chemicals`, so the LHV shares the engine basis.
def _combustion_yield(component_id: str) -> tuple[float, float, float, float]:
    """For one fuel component return (O2 demand, CO2 yield, H2O yield, LHV) per
    mole, from its formula and gas-phase formation enthalpies. Burns C->CO2,
    H->H2O(g); fuel-bound O2 reduces the air demand."""
    a = _atoms(component_id)
    nC = a.get("C", 0)
    nH = a.get("H", 0)
    nO = a.get("O", 0)
    if a.get("S") or a.get("N"):
        raise ValueError(
            f"fired-heater combustion supports C/H/O fuels only; {component_id!r} "
            f"contains S or N ({a})"
        )
    if nC == 0 and nH == 0:
        raise ValueError(f"{component_id!r} is not a fuel (no C or H atoms)")
    o2 = nC + nH / 4.0 - nO / 2.0          # stoichiometric O2 per mol fuel
    co2 = float(nC)
    h2o = nH / 2.0
    from chemicals import Hfg
    from chemicals.identifiers import CAS_from_any

    hf_co2 = float(Hfg(CAS_from_any("carbon dioxide")))
    hf_h2o = float(Hfg(CAS_from_any("water")))     # gas-phase -> LHV (not HHV)
    # LHV = Hf(fuel) - nC*Hf(CO2) - (nH/2)*Hf(H2O,g)
    lhv = _hfg(component_id) - co2 * hf_co2 - h2o * hf_h2o
    return o2, co2, h2o, lhv


def combust_fuel(
    fuel: dict[str, float],
    process_duty: float,
    efficiency: float,
    excess_air: float,
) -> FuelCombustion:
    """Combustion mass balance: fuel and air molar flows and the flue-gas
    composition needed to deliver ``process_duty`` (W) at LHV-basis fired
    ``efficiency`` with ``excess_air`` (fraction) over stoichiometric.

    ``fuel`` is a composition (mole fractions, normalised here)."""
    if process_duty <= 0.0:
        raise ValueError(f"process_duty must be > 0 (got {process_duty})")
    if not 0.0 < efficiency <= 1.0:
        raise ValueError(f"efficiency must be in (0, 1] (got {efficiency})")
    if excess_air < 0.0:
        raise ValueError(f"excess_air must be >= 0 (got {excess_air})")

    total = sum(max(v, 0.0) for v in fuel.values())
    if total <= 0.0:
        raise ValueError("fuel composition sums to <= 0")
    y = {k: max(v, 0.0) / total for k, v in fuel.items()}

    # Per-mole-of-fuel yields, mixture-averaged.
    o2_per, co2_per, h2o_per, lhv_mix = 0.0, 0.0, 0.0, 0.0
    for cid, yi in y.items():
        if yi <= 0.0:
            continue
        o2, co2, h2o, lhv = _combustion_yield(cid)
        o2_per += yi * o2
        co2_per += yi * co2
        h2o_per += yi * h2o
        lhv_mix += yi * lhv
    if lhv_mix <= 0.0:
        raise ValueError("fuel mixture has a non-positive heating value")

    fired_duty = process_duty / efficiency
    n_fuel = fired_duty / lhv_mix
    fuel_flows = {cid: n_fuel * yi for cid, yi in y.items() if yi > 0.0}

    o2_stoich = n_fuel * o2_per
    o2_fed = o2_stoich * (1.0 + excess_air)
    air = o2_fed / _O2_IN_AIR
    n2 = air * _N2_IN_AIR

    flue_flows = {
        "carbon dioxide": n_fuel * co2_per,
        "water": n_fuel * h2o_per,
        "oxygen": o2_fed - o2_stoich,
        "nitrogen": n2,
    }
    flue_flows = {k: v for k, v in flue_flows.items() if v > 0.0}
    n_flue = sum(flue_flows.values())
    flue_comp = {k: v / n_flue for k, v in flue_flows.items()}

    return FuelCombustion(
        lhv_mix=lhv_mix, fuel_flow=n_fuel, fuel_flows=fuel_flows,
        air_flow=air, o2_stoich=o2_stoich, excess_air=excess_air,
        flue_flow=n_flue, flue_composition=flue_comp, flue_flows=flue_flows,
    )


def design_fired_heater(
    process_duty: float,
    efficiency: float,
    process_T_in: float,
    process_T_out: float,
    *,
    fuel: dict[str, float] | None = None,
    excess_air: float = 0.15,
    fuel_T: float = 298.15,
    air_T: float = 298.15,
    bridgewall_T: float | None = None,
    loss_fraction: float = 0.015,
    radiant_flux: float = 32_000.0,
    convective_U: float = 40.0,
) -> FiredHeaterDesign:
    """Radiant/convective design split of a fired heater.

    ``process_duty`` (W, absorbed) and LHV-basis ``efficiency`` set the firing
    rate; ``process_T_in``/``process_T_out`` (K) are the feed inlet/outlet used
    for the convective log-mean ΔT. ``fuel`` defaults to pure methane.
    ``bridgewall_T`` (K) is the firebox-exit flue temperature that fixes the
    radiant/convective split; if omitted a typical 1150 K (~880 °C) is used,
    clamped below the adiabatic flame temperature.
    """
    from ..thermo.nasa_pkg import NasaGasPackage

    if process_T_out <= process_T_in:
        raise ValueError(
            f"process outlet T ({process_T_out} K) must exceed inlet "
            f"({process_T_in} K) — a fired heater only heats"
        )
    fuel = dict(fuel) if fuel else {"methane": 1.0}
    comb = combust_fuel(fuel, process_duty, efficiency, excess_air)
    fired_duty = process_duty / efficiency

    # Flue-gas thermo on the NASA ideal-gas (formation-inclusive) basis.
    flue_pkg = NasaGasPackage(sorted(comb.flue_composition))
    P = 101_325.0
    zf = comb.flue_composition
    h_ref = flue_pkg.enthalpy(_T_REF, P, zf)

    def h_sens(T: float) -> float:
        """Flue-gas sensible enthalpy above the reference, W."""
        return comb.flue_flow * (flue_pkg.enthalpy(T, P, zf) - h_ref)

    # Absolute heat liberated = n_fuel*LHV + reactant pre-heat above the ref.
    air_pkg = NasaGasPackage(["oxygen", "nitrogen"])
    z_air = {"oxygen": _O2_IN_AIR, "nitrogen": _N2_IN_AIR}
    air_preheat = comb.air_flow * (
        air_pkg.enthalpy(air_T, P, z_air) - air_pkg.enthalpy(_T_REF, P, z_air)
    )
    fuel_preheat = _fuel_preheat(fuel, comb.fuel_flow, fuel_T)
    q_available = comb.fuel_flow * comb.lhv_mix + air_preheat + fuel_preheat

    # Adiabatic flame temperature: all liberated heat in the flue sensible.
    T_ad = _solve_T(h_sens, q_available, lo=_T_REF, hi=2800.0)

    # Bridgewall temperature sets the radiant/convective split.
    notes: list[str] = []
    if bridgewall_T is None:
        bridgewall_T = min(1150.0, 0.5 * (T_ad + process_T_out))
        notes.append(f"bridgewall T defaulted to {bridgewall_T - 273.15:.0f} C")
    if bridgewall_T >= T_ad:
        raise ValueError(
            f"bridgewall_T ({bridgewall_T} K) must be below the adiabatic flame "
            f"temperature ({T_ad:.0f} K); lower it or reduce excess_air"
        )
    if bridgewall_T <= process_T_out:
        raise ValueError(
            f"bridgewall_T ({bridgewall_T} K) must exceed the process outlet "
            f"temperature ({process_T_out} K) for a positive radiant driving force"
        )

    casing_loss = loss_fraction * q_available
    radiant_duty = q_available - h_sens(bridgewall_T) - casing_loss
    convective_duty = process_duty - radiant_duty
    if radiant_duty <= 0.0 or convective_duty <= 0.0:
        raise ValueError(
            f"the bridgewall temperature gives a non-physical split "
            f"(Q_rad={radiant_duty:.3g} W, Q_conv={convective_duty:.3g} W); "
            f"adjust bridgewall_T"
        )

    # Stack temperature: convective bank cools the flue from T_bw. The stack
    # loss is the analytic sensible heat leaving (= h_sens(T_bw) - Q_conv), so
    # the balance Q_avail = Q_abs + casing + stack closes exactly; stack_T is the
    # temperature that carries it (brentq, to ~1e-3 K).
    stack_loss = h_sens(bridgewall_T) - convective_duty
    stack_T = _solve_T(h_sens, stack_loss, lo=_T_REF, hi=bridgewall_T)

    radiant_fraction = radiant_duty / process_duty

    # Areas. Radiant: design cold-plane flux. Convective: counter-current LMTD
    # between the flue gas (T_bw -> T_stack) and the process fluid over the
    # fraction of its rise taken in the convective bank.
    radiant_area = radiant_duty / radiant_flux
    proc_T_mid = process_T_in + (convective_duty / process_duty) * (
        process_T_out - process_T_in
    )
    dt_hot = bridgewall_T - proc_T_mid       # hot flue in vs warm process out
    dt_cold = stack_T - process_T_in         # cool flue out vs cold process in
    lmtd = _lmtd(dt_hot, dt_cold)
    convective_area = convective_duty / (convective_U * lmtd)

    return FiredHeaterDesign(
        combustion=comb,
        process_duty=process_duty,
        fired_duty=fired_duty,
        heat_available=q_available,
        efficiency_lhv=efficiency,
        efficiency_gross=process_duty / q_available,
        flame_temperature=T_ad,
        bridgewall_temperature=bridgewall_T,
        stack_temperature=stack_T,
        radiant_duty=radiant_duty,
        convective_duty=convective_duty,
        radiant_fraction=radiant_fraction,
        casing_loss=casing_loss,
        stack_loss=stack_loss,
        radiant_area=radiant_area,
        convective_area=convective_area,
        radiant_flux=radiant_flux,
        convective_U=convective_U,
        convective_lmtd=lmtd,
        notes=notes,
    )


def _fuel_preheat(fuel: dict[str, float], n_fuel: float, fuel_T: float) -> float:
    """Sensible heat carried in by fuel entering above the reference (W). Small;
    estimated from the ideal-gas heat capacity in `chemicals` (falls back to 0
    if unavailable — the term is < 1 % of the firing rate near ambient)."""
    if abs(fuel_T - _T_REF) < 1e-6:
        return 0.0
    total = sum(max(v, 0.0) for v in fuel.values())
    if total <= 0.0:
        return 0.0
    try:
        from chemicals import Cpgm  # ideal-gas molar heat capacity, J/mol/K
        from chemicals.identifiers import CAS_from_any
    except Exception:
        return 0.0
    cp = 0.0
    Tm = 0.5 * (fuel_T + _T_REF)
    for cid, v in fuel.items():
        yi = max(v, 0.0) / total
        if yi <= 0.0:
            continue
        try:
            cpi = Cpgm(CAS_from_any(cid), Tm)
        except Exception:
            cpi = None
        if cpi is None:
            return 0.0
        cp += yi * float(cpi)
    return n_fuel * cp * (fuel_T - _T_REF)


def _solve_T(h_sens, target: float, lo: float, hi: float) -> float:
    """Find T in [lo, hi] with h_sens(T) == target (h_sens monotonic increasing)."""
    from scipy.optimize import brentq

    flo = h_sens(lo) - target
    fhi = h_sens(hi) - target
    if flo > 0.0:               # target below the reference floor
        return lo
    if fhi < 0.0:
        raise ValueError(
            f"flue sensible enthalpy cannot reach the target {target:.3g} W below "
            f"{hi} K — combustion heat balance out of range"
        )
    return float(brentq(lambda T: h_sens(T) - target, lo, hi, xtol=1e-3))


def _lmtd(dt1: float, dt2: float) -> float:
    """Log-mean temperature difference of two end approaches (K)."""
    if dt1 <= 0.0 or dt2 <= 0.0:
        raise ValueError(
            f"non-positive convective ΔT ({dt1:.1f}, {dt2:.1f} K) — temperature "
            f"cross in the convective bank"
        )
    if abs(dt1 - dt2) < 1e-9:
        return dt1
    return (dt1 - dt2) / math.log(dt1 / dt2)
