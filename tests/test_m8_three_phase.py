"""M8 tests: the ThreePhaseSeparator (VLLE) and the thermo layer's three-phase
flash (`flash_pt_3p`, built on thermo's FlashVLN with two trial cubic-EOS
liquids).

Validation reference: water + n-hexane partial miscibility at ~298 K, 1 atm.
The experimental mutual solubilities are tiny — x(hexane in water) ~ 2e-6 and
x(water in hexane) ~ 5e-4 at 25 C (Tsonopoulos, *Fluid Phase Equilibria* 156
(1999) 21-33; IUPAC-NIST Solubility Data Series Vol. 38, hydrocarbons with
water) — i.e. an essentially pure aqueous phase under an organic phase. PR
with kij = 0 reproduces the *structure* (aqueous phase >99% water, organic
phase hexane-rich, aqueous denser) but overpredicts water-in-hexane (~2 mol%
vs 0.05%), so the assertions check compositional ordering and phase presence
rather than exact x values, per the known limits of a cubic EOS for water.

Adding nitrogen as a non-condensable gives a genuine three-phase V/L/L state
at ambient conditions (a binary at fixed P has VLLE only at a single T, by the
phase rule).
"""
import pytest

from caldyr.core import Component, Flowsheet
from caldyr.economics.costing import cost_equipment
from caldyr.economics.sizing import size_flowsheet
from caldyr.io import from_dict, to_dict
from caldyr.thermo import ThermoPackage, make_package
from caldyr.unitops import ThreePhaseSeparator

P_ATM = 101325.0


def vlle_separator(**param_overrides) -> Flowsheet:
    """N2 + water + n-hexane into a three-phase separator at 298.15 K, 1 atm."""
    params: dict = {"T": 298.15, "P": P_ATM}
    params.update(param_overrides)
    fs = Flowsheet(
        components=[Component("nitrogen"), Component("water"), Component("n-hexane")],
        property_package="thermo:PR")
    fs.add(ThreePhaseSeparator("SEP", params))
    fs.feed("FEED", "SEP:in1", T=300.0, P=P_ATM, molar_flow=10.0,
            z={"nitrogen": 0.2, "water": 0.4, "n-hexane": 0.4})
    fs.connect("VAP", "SEP:vapor", None)
    fs.connect("ORG", "SEP:liquid_light", None)
    fs.connect("AQ", "SEP:liquid_heavy", None)
    fs.connect("Q", "SEP:duty", None)
    return fs


def lle_separator() -> Flowsheet:
    """Binary water/n-hexane, all-liquid at 298 K (no vapor: the sum of pure
    vapor pressures, ~23 kPa, is far below 1 atm)."""
    fs = Flowsheet(components=[Component("water"), Component("n-hexane")],
                   property_package="thermo:PR")
    fs.add(ThreePhaseSeparator("SEP", {"T": 298.15, "P": P_ATM}))
    fs.feed("FEED", "SEP:in1", T=298.15, P=P_ATM, molar_flow=10.0,
            z={"water": 0.5, "n-hexane": 0.5})
    fs.connect("VAP", "SEP:vapor", None)
    fs.connect("ORG", "SEP:liquid_light", None)
    fs.connect("AQ", "SEP:liquid_heavy", None)
    fs.connect("Q", "SEP:duty", None)
    return fs


# -- the water/hexane split ----------------------------------------------------
def test_vlle_produces_three_phases_with_the_right_chemistry():
    fs = vlle_separator()
    assert fs.solve().converged
    vap, org, aq = fs.streams["VAP"], fs.streams["ORG"], fs.streams["AQ"]
    assert vap.molar_flow > 0 and org.molar_flow > 0 and aq.molar_flow > 0
    # Aqueous (heavy) phase: essentially pure water (expt: >99.99% water).
    assert aq.z["water"] > 0.99
    # Organic (light) phase: hexane-rich, nearly water/N2-free.
    assert org.z["n-hexane"] > 0.90
    assert org.z["water"] < 0.10
    # Vapor: carries the nitrogen.
    assert vap.z["nitrogen"] > org.z["nitrogen"] and vap.z["nitrogen"] > aq.z["nitrogen"]


def test_liquids_are_ordered_by_density():
    fs = vlle_separator()
    fs.solve()
    result = fs.units["SEP"].result
    assert result is not None
    # Heavy liquid (aqueous) is denser than light liquid (organic) — that is
    # the *definition* of the port assignment.
    assert result["rho_heavy"] > result["rho_light"]
    # And the heavy one is the water phase.
    assert fs.streams["AQ"].z["water"] > fs.streams["ORG"].z["water"]


def test_mass_balance_closes_per_component():
    fs = vlle_separator()
    fs.solve()
    feed = fs.streams["FEED"]
    outs = [fs.streams[s] for s in ("VAP", "ORG", "AQ")]
    for c in feed.components:
        n_in = feed.molar_flow * feed.z[c]
        n_out = sum(s.molar_flow * s.z.get(c, 0.0) for s in outs)
        assert n_out == pytest.approx(n_in, rel=1e-8, abs=1e-8)


def test_energy_balance_closes_with_the_duty():
    fs = vlle_separator(T=290.0)                   # force a nonzero duty
    rep = fs.solve()
    feed = fs.streams["FEED"]
    h_in = feed.molar_flow * feed.H + rep.duties["Q"]
    h_out = sum(fs.streams[s].molar_flow * fs.streams[s].H
                for s in ("VAP", "ORG", "AQ"))
    assert h_out == pytest.approx(h_in, rel=1e-9)
    assert rep.duties["Q"] != 0.0


# -- degraded phase counts (graceful, not crashes) -------------------------------
def test_lle_without_vapor_gives_an_empty_vapor_stream():
    fs = lle_separator()
    assert fs.solve().converged
    assert fs.streams["VAP"].molar_flow == pytest.approx(0.0, abs=1e-12)
    assert fs.streams["ORG"].molar_flow > 0 and fs.streams["AQ"].molar_flow > 0
    assert fs.streams["AQ"].z["water"] > 0.99
    assert fs.streams["ORG"].z["n-hexane"] > 0.90


def test_fully_miscible_system_gives_an_empty_heavy_liquid():
    """Benzene/toluene at 300 K is one liquid: the separator must degrade to a
    two-phase result (empty heavy-liquid stream), NOT crash."""
    fs = Flowsheet(components=[Component("benzene"), Component("toluene")],
                   property_package="thermo:PR")
    fs.add(ThreePhaseSeparator("SEP", {"T": 300.0, "P": P_ATM}))
    fs.feed("FEED", "SEP:in1", T=300.0, P=P_ATM, molar_flow=10.0,
            z={"benzene": 0.5, "toluene": 0.5})
    fs.connect("VAP", "SEP:vapor", None)
    fs.connect("L1", "SEP:liquid_light", None)
    fs.connect("L2", "SEP:liquid_heavy", None)
    fs.connect("Q", "SEP:duty", None)
    assert fs.solve().converged
    assert fs.streams["L2"].molar_flow == pytest.approx(0.0, abs=1e-12)
    assert fs.streams["VAP"].molar_flow == pytest.approx(0.0, abs=1e-12)
    assert fs.streams["L1"].molar_flow == pytest.approx(10.0, rel=1e-12)


def test_all_vapor_feed_gives_empty_liquids():
    fs = vlle_separator(T=500.0)
    fs.solve()
    assert fs.streams["VAP"].molar_flow == pytest.approx(10.0, rel=1e-12)
    assert fs.streams["ORG"].molar_flow == pytest.approx(0.0, abs=1e-12)
    assert fs.streams["AQ"].molar_flow == pytest.approx(0.0, abs=1e-12)


# -- typed, documented error paths ------------------------------------------------
def test_missing_t_raises_a_clear_error():
    """No adiabatic (PH) mode: upstream enthalpies live on the two-phase flash
    surface, so PH-three-phase would mix enthalpy surfaces (documented)."""
    fs = vlle_separator()
    del fs.units["SEP"].params["T"]
    with pytest.raises(ValueError, match=r"params\['T'\] is required"):
        fs.solve()


def test_activity_package_raises_not_implemented():
    """Three-phase flashes are PR/SRK-only for now; the NRTL activity package
    must fail loudly with a clear message, not mislabel phases."""
    fs = Flowsheet(components=[Component("water"), Component("ethanol")],
                   property_package="thermo:NRTL")
    fs.add(ThreePhaseSeparator("SEP", {"T": 298.15, "P": P_ATM}))
    fs.feed("FEED", "SEP:in1", T=298.15, P=P_ATM, molar_flow=10.0,
            z={"water": 0.5, "ethanol": 0.5})
    fs.connect("VAP", "SEP:vapor", None)
    fs.connect("L1", "SEP:liquid_light", None)
    fs.connect("L2", "SEP:liquid_heavy", None)
    fs.connect("Q", "SEP:duty", None)
    with pytest.raises(NotImplementedError, match="three-phase"):
        fs.solve()


def test_pure_component_three_phase_flash_raises():
    pp = ThermoPackage(["water"])
    with pytest.raises(ValueError, match="at least two components"):
        pp.flash_pt_3p(298.15, P_ATM, {"water": 1.0})


# -- package-level PH consistency ---------------------------------------------------
def test_flash_ph_3p_round_trips_its_own_surface():
    """flash_ph_3p is exposed for callers that hold a VLN-surface enthalpy:
    feeding back the H from a PT three-phase flash recovers the temperature."""
    pp = ThermoPackage(["nitrogen", "water", "n-hexane"])
    z = {"nitrogen": 0.2, "water": 0.4, "n-hexane": 0.4}
    res = pp.flash_pt_3p(298.15, P_ATM, z)
    back = pp.flash_ph_3p(P_ATM, res.H, z)
    assert back.T == pytest.approx(298.15, abs=1e-4)


# -- `.flow` round-trip and economics ----------------------------------------------
def test_flow_round_trip_is_exact():
    fs = vlle_separator()
    fs.solve()
    doc = to_dict(fs)
    assert doc == to_dict(from_dict(doc))


def test_economics_sizes_a_horizontal_vessel():
    """Three-phase separators are horizontal drums (Turton 4e Table A.1
    horizontal process vessel; settling length governs, not height)."""
    fs = lle_separator()
    rep = fs.solve()
    pp = make_package(fs.property_package, fs.component_ids)
    sizes = size_flowsheet(fs, rep, pp)
    assert len(sizes) == 1
    assert sizes[0].equipment_type == "vessel_horizontal"
    assert sizes[0].attribute_name == "volume_m3"
    cost = cost_equipment(sizes[0])
    assert cost.bare_module > 0
