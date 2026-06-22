"""M11 tests: the FiredHeater unit op + Turton fired-heater costing.

Physics reference: a fired heater is a Heater whose heat comes from burning
fuel at a fired efficiency — the fuel (released) duty exceeds the process
(absorbed) duty by exactly 1/efficiency. Turton 4e Ch. 8 uses fired
efficiencies of 0.80-0.90; the default here is 0.85.

Costing reference: **Turton et al., 4e, Appendix A** — non-reactive fired
heater, capacity = absorbed duty Q in kW (Table A.1: K1=7.3488, K2=-1.1666,
K3=0.2028, valid 1,000-100,000 kW). Hand-checked point: at Q = 10,000 kW,
log10 Cp0 = 7.3488 - 1.1666*4 + 0.2028*16 = 5.9272 -> Cp0 ~ $845,800 (CEPCI
397 basis). Bare module Cbm = Cp0 * Fbm * Fp with Fbm = 2.13 (Table A.7, CS)
and the Table A.2 pressure factor. P13 (2026-06-22) verified these against the
authoritative CAPCOST spreadsheet (Turton CD): the non-reactive *process heater*
pressure-factor row is (C1=0.1347, C2=-0.2368, C3=0.1021, 10-200 barg) — the
earlier (0.1017, -0.1957, 0.09403) was the *pyrolysis furnace* (reactive) row.
At 50 barg, log10 Fp = 0.1347 - 0.2368*log10(50) + 0.1021*(log10 50)^2 = 0.0271
-> Fp ~ 1.064.
"""
import pytest

from caldyr.core import Component, Flowsheet
from caldyr.economics.costing import cost_equipment
from caldyr.economics.sizing import EquipmentSize, size_flowsheet
from caldyr.thermo import make_package
from caldyr.unitops import FiredHeater

T_IN = 300.0     # K
T_OUT = 650.0    # K
P = 5.0e5        # Pa
N = 150.0        # mol/s  (-> ~1.5 MW process duty, inside the 1-100 MW range)


def n2_fired_heater(**param_overrides) -> Flowsheet:
    params: dict = {"T_out": T_OUT}
    params.update(param_overrides)
    fs = Flowsheet(components=[Component("nitrogen")], property_package="thermo:PR")
    fs.add(FiredHeater("FH", params))
    fs.feed("FEED", "FH:in1", T=T_IN, P=P, molar_flow=N, z={"nitrogen": 1.0})
    fs.connect("OUT", "FH:out", None)
    fs.connect("Q", "FH:duty", None)
    return fs


def test_fuel_duty_is_process_duty_over_efficiency_exactly():
    fs = n2_fired_heater()                       # default efficiency = 0.85
    rep = fs.solve()
    assert rep.converged
    q = rep.duties["Q"]
    assert q > 0
    unit = fs.units["FH"]
    assert unit.design is not None
    assert unit.design["efficiency"] == pytest.approx(0.85)
    assert unit.design["process_duty"] == pytest.approx(q, rel=1e-12)
    # The deliverable: fuel duty / process duty = 1 / efficiency, exact.
    assert unit.design["fuel_duty"] / q == pytest.approx(1.0 / 0.85, rel=1e-12)


def test_custom_efficiency_scales_fuel_duty():
    fs = n2_fired_heater(efficiency=0.92)
    rep = fs.solve()
    assert fs.units["FH"].design["fuel_duty"] == \
        pytest.approx(rep.duties["Q"] / 0.92, rel=1e-12)


def test_energy_balance_closes_on_the_process_duty():
    """The duty port carries the PROCESS duty (the stream's enthalpy change),
    not the fuel duty — flowsheet energy balances close like a Heater's."""
    fs = n2_fired_heater()
    rep = fs.solve()
    feed, out = fs.streams["FEED"], fs.streams["OUT"]
    assert out.T == pytest.approx(T_OUT)
    assert out.molar_flow * out.H - feed.molar_flow * feed.H == \
        pytest.approx(rep.duties["Q"], rel=1e-9)


def test_fixed_q_spec_and_dp():
    fs = n2_fired_heater(T_out=None, Q=1.5e6, dP=0.5e5)
    rep = fs.solve()
    assert rep.duties["Q"] == pytest.approx(1.5e6)
    assert fs.streams["OUT"].P == pytest.approx(P - 0.5e5)
    assert fs.units["FH"].design["fuel_duty"] == pytest.approx(1.5e6 / 0.85, rel=1e-12)


def test_bad_specs_raise():
    with pytest.raises(ValueError, match="efficiency"):
        n2_fired_heater(efficiency=0.0).solve()
    with pytest.raises(ValueError, match="efficiency"):
        n2_fired_heater(efficiency=1.2).solve()
    with pytest.raises(ValueError, match="exactly one of"):
        n2_fired_heater(Q=1e6).solve()              # both T_out and Q
    with pytest.raises(ValueError, match="only heats"):
        n2_fired_heater(T_out=250.0).solve()        # below the inlet T


# -- economics ---------------------------------------------------------------
def test_sizing_books_fuel_duty_on_the_fired_heat_utility():
    fs = n2_fired_heater()
    rep = fs.solve()
    pp = make_package(fs.property_package, fs.component_ids)
    sizes = size_flowsheet(fs, rep, pp)
    assert len(sizes) == 1
    size = sizes[0]
    assert size.equipment_type == "fired_heater"
    assert size.attribute_name == "duty_kW"
    # Sized on the process duty (the Turton correlation capacity)...
    assert size.attribute == pytest.approx(rep.duties["Q"] / 1e3, rel=1e-12)
    assert 1000.0 < size.attribute < 100_000.0
    # ...but the utility (opex) is billed on the FUEL duty.
    assert size.utility == "fired_heat"
    assert size.utility_duty_W == pytest.approx(rep.duties["Q"] / 0.85, rel=1e-12)

    cost = cost_equipment(size)
    assert cost.bare_module > cost.purchased > 0


def test_purchased_cost_matches_turton_hand_point():
    """Q = 10,000 kW -> log10 Cp0 = 5.9272 -> Cp0 ~ $845,800 (CEPCI 397)."""
    fh = EquipmentSize("F1", "fired_heater", 10_000.0, "duty_kW", pressure_barg=1.0)
    c = cost_equipment(fh, year=2001)              # no escalation at base year
    assert c.purchased == pytest.approx(845_800.0, rel=0.01)
    assert c.factors["Fp"] == pytest.approx(1.0)   # below 10 barg -> Fp = 1
    # Cbm = Cp0 * Fbm * Fp with the direct bare-module factor 2.13 (CS).
    assert c.bare_module == pytest.approx(2.13 * c.purchased, rel=1e-9)
    assert not c.warnings


def test_pressure_factor_applies_above_10_barg():
    """Turton 4e Table A.2 / CAPCOST non-reactive process heater: at 50 barg,
    Fp ~ 1.064 (P13-corrected from the pyrolysis-furnace row's ~1.098)."""
    lo = cost_equipment(EquipmentSize("F", "fired_heater", 10_000.0, "duty_kW",
                                      pressure_barg=5.0), year=2001)
    hi = cost_equipment(EquipmentSize("F", "fired_heater", 10_000.0, "duty_kW",
                                      pressure_barg=50.0), year=2001)
    assert lo.factors["Fp"] == pytest.approx(1.0)
    assert hi.factors["Fp"] == pytest.approx(1.064, abs=0.01)
    assert hi.bare_module > lo.bare_module


def test_fired_heater_factors_are_the_nonreactive_process_heater_row():
    """P13: lock the fired-heater costing factors to the CAPCOST/Table A.2
    *non-reactive process heater* row (the purchased K-triple and the Fp
    C-triple must come from the SAME row, not a mix with the pyrolysis furnace)."""
    from caldyr.economics import data

    # Purchased cost: non-reactive Process Heater K-triple (CAPCOST).
    k = data.PURCHASED["fired_heater"]
    assert (k.K1, k.K2, k.K3) == (7.3488, -1.1666, 0.2028)
    # Pressure factor: the MATCHING Process Heater C-triple, NOT the pyrolysis
    # furnace's (0.1017, -0.1957, 0.09403).
    p = data.PRESSURE["fired_heater"]
    assert (p.C1, p.C2, p.C3) == (0.1347, -0.2368, 0.1021)
