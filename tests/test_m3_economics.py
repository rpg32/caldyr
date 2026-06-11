"""M3 acceptance tests: the techno-economic pipeline.

Costing reference: **Turton et al., 4th ed., Appendix A** worked values — a
floating-head shell-and-tube heat exchanger of area 100 m^2 has purchased cost
Cp0 ~ $23,500 (CEPCI 397, carbon steel, ambient) and bare-module factor
Fbm ~ 3.29 at low pressure. Financial relations (CRF, NPV, payback) are checked
against closed-form hand calculations.
"""
import pytest

from caldyr.core import Component, Flowsheet
from caldyr.economics import (
    TEAConfig,
    analyze,
    capital_recovery_factor,
    cost_equipment,
    estimate_capital,
    monte_carlo,
    profitability,
    purchased_cost,
    six_tenths,
    tornado,
)
from caldyr.economics.costing import pressure_factor, vessel_pressure_factor
from caldyr.economics.sizing import EquipmentSize
from caldyr.unitops import EquilibriumReactor, FlashDrum, Heater, Mixer, Splitter

AMMONIA = {"stoich": {"nitrogen": -1, "hydrogen": -3, "ammonia": 2}, "key": "nitrogen"}


# -- costing correlations vs Turton ----------------------------------------
def test_purchased_cost_matches_turton_heat_exchanger():
    cp0 = purchased_cost(4.3247, -0.3030, 0.1634, 100.0)
    assert cp0 == pytest.approx(23500.0, rel=0.01)        # Turton 4e, A=100 m^2


def test_bare_module_matches_turton_at_low_pressure():
    hx = EquipmentSize(unit_id="E1", equipment_type="heat_exchanger",
                       attribute=100.0, attribute_name="area_m2",
                       pressure_barg=1.0, material="CS")
    c = cost_equipment(hx, year=2001)                     # no escalation at base year
    assert c.factors["Fp"] == pytest.approx(1.0)          # below 5 barg -> Fp = 1
    assert c.factors["Fbm"] == pytest.approx(3.29, abs=0.01)
    assert c.bare_module == pytest.approx(23500.0 * 3.29, rel=0.02)


def test_cepci_escalation():
    base = cost_equipment(EquipmentSize("E", "heat_exchanger", 100.0, "area_m2",
                                        pressure_barg=1.0), year=2001)
    now = cost_equipment(EquipmentSize("E", "heat_exchanger", 100.0, "area_m2",
                                       pressure_barg=1.0), year=2023)
    assert now.bare_module / base.bare_module == pytest.approx(797.9 / 397.0, rel=1e-6)


def test_pressure_factor_clamped_below_minimum():
    assert pressure_factor("heat_exchanger", 2.0) == 1.0          # below pmin=5
    assert pressure_factor("heat_exchanger", 50.0) > 1.0


def test_vessel_pressure_factor_rises_with_pressure():
    lo = vessel_pressure_factor(10.0, 1.5)
    hi = vessel_pressure_factor(200.0, 1.5)
    assert hi > lo >= 1.0


def test_six_tenths_scaling():
    assert six_tenths(1.0e6, 100.0, 200.0) == pytest.approx(1.0e6 * 2 ** 0.6, rel=1e-9)


# -- capital roll-up -------------------------------------------------------
def test_capital_rollup_follows_turton_factors():
    sizes = [EquipmentSize("E1", "heat_exchanger", 100.0, "area_m2", pressure_barg=1.0),
             EquipmentSize("V1", "vessel_vertical", 5.0, "volume_m3",
                           pressure_barg=10.0, diameter_m=1.2)]
    costs = [cost_equipment(s) for s in sizes]
    cap = estimate_capital(costs)
    isbl = sum(c.bare_module for c in costs)
    base = sum(c.bare_module_base for c in costs)
    assert cap.isbl == pytest.approx(isbl)
    assert cap.grassroots == pytest.approx(1.18 * isbl + 0.50 * base)
    assert cap.tci == pytest.approx(cap.grassroots * 1.15)


# -- financial relations ---------------------------------------------------
def test_capital_recovery_factor_closed_form():
    assert capital_recovery_factor(0.10, 20) == pytest.approx(0.117459, abs=1e-6)
    assert capital_recovery_factor(0.0, 10) == pytest.approx(0.1)     # zero-rate limit


def test_profitability_npv_payback_hand_calc():
    pr = profitability(tci=1000.0, annual_opex=0.0, annual_revenue=200.0,
                       annual_production=100.0, discount_rate=0.10, project_years=10)
    # NPV = -1000 + 200 * annuity(10%,10), annuity = 6.14457
    assert pr.npv == pytest.approx(-1000.0 + 200.0 * 6.144567, abs=1.0)
    assert pr.payback_years == pytest.approx(5.0)
    assert pr.irr == pytest.approx(0.15098, abs=1e-3)    # 200/yr for 10 yr on 1000


def test_lcop_is_levelized_break_even():
    """By construction LCOP * production == annualized capital + opex."""
    pr = profitability(tci=1.0e6, annual_opex=2.0e5, annual_revenue=0.0,
                       annual_production=1.0e6, discount_rate=0.10, project_years=20)
    crf = capital_recovery_factor(0.10, 20)
    assert pr.lcop * 1.0e6 == pytest.approx(crf * 1.0e6 + 2.0e5, rel=1e-9)


# -- the M3 deliverable: cost the ammonia loop -----------------------------
def ammonia_loop() -> Flowsheet:
    fs = Flowsheet(
        components=[Component(c) for c in ("nitrogen", "hydrogen", "ammonia", "argon")],
        property_package="thermo:PR")
    fs.add(Mixer("MIX", {"dP": 0.0}))
    fs.add(Heater("PREHEAT", {"T_out": 673.15}))
    fs.add(EquilibriumReactor("RXN", {"reaction": AMMONIA, "T": 673.15}))
    fs.add(Heater("COOL", {"T_out": 250.0}))
    fs.add(FlashDrum("SEP", {"T": 250.0, "P": 2e7}))
    fs.add(Splitter("SPLIT", {"split": 0.90}))
    fs.feed("MAKEUP", "MIX:in1", T=300.0, P=2e7, molar_flow=100.0,
            z={"nitrogen": 0.2475, "hydrogen": 0.7425, "ammonia": 0.0, "argon": 0.01})
    fs.connect("S1", "MIX:out", "PREHEAT:in1")
    fs.connect("S2", "PREHEAT:out", "RXN:in1")
    fs.connect("S3", "RXN:out", "COOL:in1")
    fs.connect("S4", "COOL:out", "SEP:in1")
    fs.connect("PRODUCT", "SEP:liquid", None)
    fs.connect("VAP", "SEP:vapor", "SPLIT:in1")
    fs.connect("RECYCLE", "SPLIT:out1", "MIX:in2")
    fs.connect("PURGE", "SPLIT:out2", None)
    for u in ("PREHEAT", "RXN", "COOL", "SEP"):
        fs.connect(f"Q_{u}", f"{u}:duty", None)
    return fs


def test_analyze_ammonia_loop_end_to_end():
    fs = ammonia_loop()
    report = fs.solve(tol=1e-7, max_iter=400)
    res = analyze(fs, report, TEAConfig())

    # A unit was sized and costed for each cost-bearing block (mixer/splitter free).
    assert {s.unit_id for s in res.sizes} == {"PREHEAT", "RXN", "COOL", "SEP"}
    assert res.capital.tci > res.capital.isbl > 0
    assert res.opex.total > res.opex.raw_materials > 0
    assert res.annual_production_kg > 0
    assert res.profitability.lcop > 0
    # Feed (H2) cost dominates ammonia economics -> top tornado bar.
    bars = tornado(fs, res.sizes, res.config)
    assert bars == sorted(bars, key=lambda b: b.swing, reverse=True)
    assert bars[0].variable.startswith("feed price")


def test_monte_carlo_bands_are_ordered_and_reproducible():
    fs = ammonia_loop()
    report = fs.solve(tol=1e-7, max_iter=400)
    res = analyze(fs, report, TEAConfig())
    mc1 = monte_carlo(fs, res.sizes, res.config, n=500, seed=7)
    mc2 = monte_carlo(fs, res.sizes, res.config, n=500, seed=7)
    assert mc1.lcop["p10"] < mc1.lcop["p50"] < mc1.lcop["p90"]
    assert mc1.lcop["p50"] == mc2.lcop["p50"]      # seeded -> reproducible
