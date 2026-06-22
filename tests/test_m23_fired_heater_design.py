"""M23 / P11 tests: the fired-heater radiant/convective design split.

Reference: Hameed (2025), *Chemical Process Simulations using Aspen HYSYS*,
§4.3 — a direct-fired heater divided into a radiant firebox (Stefan-Boltzmann,
eq. 4.8) and a convective bank (``Q = U A ΔT_lm``, eq. 4.9), plus the classic
API/Lobo-Evans furnace heat balance (Perry 8e §27; Towler & Sinnott 2e §19.16).

The book's §4.3.4 worked problem heats a 1250 kmol/h C1-C5 feed from 50 to
300 °C with methane fuel at 100 % excess air and 60 % efficiency, and reports a
fuel rate of 62.01 kmol/h (Fig. 4.19). Combustion is LHV-based: feeding the
book's implied absorbed duty reproduces 62.01 kmol/h to two decimals; on our PR
process duty (8.82 MW, ~6 % above HYSYS) the fuel comes out ~66 kmol/h — the gap
is the PR-vs-HYSYS enthalpy difference, not the combustion model.
"""
import pytest

from caldyr.core import Component, Flowsheet
from caldyr.economics.fired_heater_design import (
    combust_fuel,
    design_fired_heater,
)
from caldyr.unitops import FiredHeater

# Book §4.3.4: fuel 62.01 kmol/h, methane LHV ~802.6 kJ/mol, eta 0.6.
_LHV_CH4 = 802.58e3
_BOOK_FUEL_KMOLH = 62.01
_BOOK_IMPLIED_DUTY = _BOOK_FUEL_KMOLH / 3.6 * _LHV_CH4 * 0.60   # W


# -- combustion mass balance -------------------------------------------------
def test_methane_lhv_matches_standard():
    c = combust_fuel({"methane": 1.0}, 8.0e6, 0.60, 0.15)
    # Lower heating value of methane is ~802.3 kJ/mol (gaseous water).
    assert c.lhv_mix == pytest.approx(802.6e3, rel=0.005)


def test_book_implied_duty_reproduces_fuel_rate():
    """Fed the book's implied absorbed duty, the model returns 62.01 kmol/h."""
    c = combust_fuel({"methane": 1.0}, _BOOK_IMPLIED_DUTY, 0.60, excess_air=1.0)
    assert c.fuel_flow * 3.6 == pytest.approx(_BOOK_FUEL_KMOLH, abs=0.05)


def test_stoichiometric_air_and_flue_composition():
    """100 % excess air on methane: O2 fed is twice stoichiometric, and the flue
    gas carries CO2 : H2O : O2 : N2 in the expected combustion ratios."""
    c = combust_fuel({"methane": 1.0}, 8.0e6, 0.60, excess_air=1.0)
    # CH4 + 2 O2 -> CO2 + 2 H2O ; per mol fuel stoich O2 = 2.
    assert c.o2_stoich == pytest.approx(2.0 * c.fuel_flow, rel=1e-9)
    # O2 fed = 2 x stoich (100 % excess) -> air = 2*stoich/0.21.
    assert c.air_flow == pytest.approx(2.0 * c.o2_stoich / 0.21, rel=1e-9)
    f = c.flue_flows
    assert f["carbon dioxide"] == pytest.approx(c.fuel_flow, rel=1e-9)
    assert f["water"] == pytest.approx(2.0 * c.fuel_flow, rel=1e-9)
    assert f["oxygen"] == pytest.approx(c.o2_stoich, rel=1e-9)   # half the fed O2 left
    assert sum(c.flue_composition.values()) == pytest.approx(1.0, rel=1e-12)


def test_excess_air_increases_air_and_dilutes_flue():
    lo = combust_fuel({"methane": 1.0}, 8.0e6, 0.60, excess_air=0.15)
    hi = combust_fuel({"methane": 1.0}, 8.0e6, 0.60, excess_air=1.0)
    assert hi.air_flow > lo.air_flow
    # Same fuel -> same fired duty -> same fuel flow; more air dilutes O2 frac up.
    assert hi.fuel_flow == pytest.approx(lo.fuel_flow, rel=1e-9)
    assert hi.flue_composition["oxygen"] > lo.flue_composition["oxygen"]


def test_heavier_fuel_burns_via_formula():
    """A C2-C3 fuel mixture is handled from its formula/formation enthalpy
    (heavier hydrocarbons are not in the NASA flue package, only the products)."""
    c = combust_fuel({"ethane": 0.5, "propane": 0.5}, 8.0e6, 0.80, excess_air=0.2)
    # Heavier hydrocarbons have higher per-mole LHV than methane.
    assert c.lhv_mix > 1.2e6
    assert c.flue_flows["carbon dioxide"] > c.fuel_flow   # >1 C per fuel molecule


# -- radiant / convective design split ---------------------------------------
def book_design(**kw):
    base = dict(
        process_duty=8.8167e6, efficiency=0.60,
        process_T_in=50 + 273.15, process_T_out=300 + 273.15,
        fuel={"methane": 1.0}, excess_air=1.0,
        fuel_T=40 + 273.15, air_T=40 + 273.15,
    )
    base.update(kw)
    return design_fired_heater(**base)


def test_energy_balance_closes_exactly():
    """Q_available = Q_absorbed + casing loss + stack loss, to machine precision
    (the radiant/convective split is constructed to conserve energy)."""
    d = book_design()
    assert d.heat_available == pytest.approx(
        d.process_duty + d.casing_loss + d.stack_loss, rel=1e-9
    )
    # Radiant + convective = the absorbed process duty.
    assert d.radiant_duty + d.convective_duty == pytest.approx(d.process_duty, rel=1e-12)


def test_temperatures_ordered_and_physical():
    d = book_design()
    # T_ad > T_bridgewall > T_stack > reference; bridgewall > process outlet.
    assert d.flame_temperature > d.bridgewall_temperature > d.stack_temperature
    assert d.bridgewall_temperature > 300 + 273.15
    assert d.stack_temperature > 298.15
    # 100 % excess air on methane: adiabatic flame around ~1200 C, not stoich-hot.
    assert 1100 < d.flame_temperature - 273.15 < 1350


def test_radiant_fraction_in_typical_range():
    """A process furnace takes the bulk of its absorbed duty in the radiant
    section: the radiant fraction is ~0.5-0.75 (Towler & Sinnott)."""
    d = book_design()
    assert 0.45 < d.radiant_fraction < 0.80
    # gross (absolute-basis) efficiency is close to the LHV-basis 0.60.
    assert d.efficiency_gross == pytest.approx(0.60, abs=0.02)


def test_areas_positive_and_flux_consistent():
    d = book_design()
    assert d.radiant_area == pytest.approx(d.radiant_duty / d.radiant_flux, rel=1e-12)
    assert d.convective_area == pytest.approx(
        d.convective_duty / (d.convective_U * d.convective_lmtd), rel=1e-9
    )
    assert d.radiant_area > 0 and d.convective_area > 0


def test_higher_bridgewall_shifts_duty_to_convective():
    """Raising the bridgewall temperature leaves more sensible heat in the flue
    entering the convective bank -> less radiant, more convective duty. The stack
    temperature is set by the overall efficiency (casing + stack loss), so the
    split does not move it."""
    lo = book_design(bridgewall_T=950.0)
    hi = book_design(bridgewall_T=1150.0)
    assert hi.radiant_duty < lo.radiant_duty
    assert hi.convective_duty > lo.convective_duty
    assert hi.stack_temperature == pytest.approx(lo.stack_temperature, abs=0.01)


def test_bridgewall_validation():
    # Above the adiabatic flame temperature -> impossible.
    with pytest.raises(ValueError, match="adiabatic flame"):
        book_design(bridgewall_T=2000.0)
    # Below the process outlet -> no radiant driving force.
    with pytest.raises(ValueError, match="process outlet"):
        book_design(bridgewall_T=400.0)


def test_combust_fuel_bad_inputs_raise():
    with pytest.raises(ValueError, match="process_duty"):
        combust_fuel({"methane": 1.0}, -1.0, 0.6, 0.15)
    with pytest.raises(ValueError, match="efficiency"):
        combust_fuel({"methane": 1.0}, 8e6, 1.5, 0.15)
    with pytest.raises(ValueError, match="excess_air"):
        combust_fuel({"methane": 1.0}, 8e6, 0.6, -0.1)
    with pytest.raises(ValueError, match="S or N"):
        combust_fuel({"hydrogen sulfide": 1.0}, 8e6, 0.6, 0.15)


# -- wired into the FiredHeater unit op --------------------------------------
def fired_heater_fs(**params) -> Flowsheet:
    comps = ["methane", "ethane", "propane", "isobutane",
             "n-butane", "isopentane", "n-pentane"]
    flows = {"methane": 300, "ethane": 250, "propane": 200, "isobutane": 150,
             "n-butane": 150, "isopentane": 100, "n-pentane": 100}
    tot = sum(flows.values())
    z = {c: flows[c] / tot for c in comps}
    fs = Flowsheet(components=[Component(c) for c in comps],
                   property_package="thermo:PR")
    p = {"T_out": 300 + 273.15}
    p.update(params)
    fs.add(FiredHeater("FH", p))
    fs.feed("FEED", "FH:in1", T=50 + 273.15, P=5 * 101325.0,
            molar_flow=tot * 1000 / 3600.0, z=z)
    fs.connect("OUT", "FH:out", None)
    fs.connect("Q", "FH:duty", None)
    return fs


def test_unit_op_basic_path_has_no_design_split():
    """Without a design trigger the unit keeps the simple efficiency contract."""
    fs = fired_heater_fs(efficiency=0.60)
    fs.solve()
    d = fs.units["FH"].design
    assert "combustion" not in d and "firing" not in d
    assert d["fuel_duty"] == pytest.approx(d["process_duty"] / 0.60, rel=1e-12)


def test_unit_op_design_split_book_case():
    """The book §4.3 case through the unit op: fuel ~66 kmol/h (PR duty), full
    radiant/convective design published on unit.design."""
    fs = fired_heater_fs(efficiency=0.60, excess_air=1.0,
                         fuel="methane", fuel_T=40 + 273.15, air_T=40 + 273.15)
    fs.solve()
    d = fs.units["FH"].design
    assert "combustion" in d and "firing" in d
    # PR process duty ~8.8 MW -> fuel ~66 kmol/h (book 62.01; ~6 % thermo gap).
    assert d["combustion"]["fuel_flow"] * 3.6 == pytest.approx(66.0, abs=2.0)
    fr = d["firing"]
    assert 0.45 < fr["radiant_fraction"] < 0.80
    assert fr["flame_temperature"] > fr["bridgewall_temperature"] > fr["stack_temperature"]
    # energy balance closes on the published numbers
    assert fr["heat_available"] == pytest.approx(
        d["process_duty"] + fr["casing_loss"] + fr["stack_loss"], rel=1e-9
    )
