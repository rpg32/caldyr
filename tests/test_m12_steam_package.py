"""M12 tests: the ``coolprop:Water`` steam-tables property package.

References:
* Hameed, *Chemical Process Simulations using Aspen HYSYS* (Wiley 2025), §2.2
  "Steam Table" — a saturated-steam table over 40–400 °C (Psat, mass enthalpy,
  entropy, heat of vaporization); the package reproduces its saturation line.
* Cengel & Boles, *Thermodynamics: An Engineering Approach* (8e), Tables A-4 /
  A-5 (the classic steam tables): Tsat(101.325 kPa) = 99.97 C = 373.12 K with
  h_fg = 2256.5 kJ/kg; Psat(40 C) = 7.3851 kPa with h_fg = 2406.0 kJ/kg.

The enthalpy-basis cross-check (formation-inclusive absolute basis shared with
the `thermo` backends) is the critical test: a Heater duty on a water stream
must agree between ``thermo:PR`` and ``coolprop:Water`` to within PR's physical
accuracy for water. Measured deltas (asserted with headroom):
  * vapor-phase heating 400 -> 500 K at 1 bar: PR is 2.20% below IAPWS-95;
  * liquid heating 300 -> 350 K at 5 bar: PR is 7.9% above IAPWS-95 (a cubic
    EOS is a poor liquid-water Cp model — that *is* why this package exists).
"""
import pytest

from caldyr.core import Component, Flowsheet
from caldyr.thermo import CoolPropWaterPackage, PropertyPackage, make_package
from caldyr.unitops import Heater

Z = {"water": 1.0}
P_ATM = 101_325.0


@pytest.fixture(scope="module")
def pp() -> CoolPropWaterPackage:
    return make_package("coolprop:Water", ["water"])  # type: ignore[return-value]


# -- construction & typed errors ----------------------------------------------
def test_make_package_builds_and_satisfies_protocol(pp):
    assert isinstance(pp, CoolPropWaterPackage)
    assert isinstance(pp, PropertyPackage)


def test_multicomponent_is_a_typed_error():
    with pytest.raises(ValueError, match="exactly one component"):
        make_package("coolprop:Water", ["water", "ethanol"])


def test_non_water_component_is_a_typed_error():
    with pytest.raises(ValueError, match="only water"):
        make_package("coolprop:Water", ["methane"])


def test_unknown_coolprop_method_is_a_typed_error():
    with pytest.raises(ValueError, match="coolprop:Water"):
        make_package("coolprop:R134a", ["water"])


def test_unsupported_calls_raise_not_implemented(pp):
    with pytest.raises(NotImplementedError, match="steam-tables"):
        pp.lnKeq({"water": 1.0}, 400.0)
    with pytest.raises(NotImplementedError, match="three-phase"):
        pp.flash_pt_3p(300.0, 1e5, Z)
    with pytest.raises(NotImplementedError, match="three-phase"):
        pp.flash_ph_3p(1e5, -280000.0, Z)


def test_foreign_composition_is_a_typed_error(pp):
    with pytest.raises(ValueError, match="not in package"):
        pp.enthalpy(300.0, 1e5, {"ethanol": 1.0})


# -- steam-table validation (book §2.2 + Cengel & Boles A-4/A-5) ---------------
def test_saturation_at_one_atm_matches_steam_tables(pp):
    """Tsat(1 atm) = 373.12 K, h_fg = 2256.5 kJ/kg (Cengel & Boles 8e, A-5)."""
    t_bub, t_dew = pp.bubble_dew(P_ATM, Z)
    assert t_bub == t_dew                       # pure fluid: bubble = dew (book §2.3)
    assert t_bub == pytest.approx(373.12, abs=0.01)

    sat = pp.bubble_point(P_ATM, Z)
    h_fg = (sat.H_vapor - sat.H_liquid) / pp.mw / 1e3      # kJ/kg
    assert h_fg == pytest.approx(2256.5, rel=1e-3)


def test_saturation_at_40C_matches_steam_tables(pp):
    """Book §2.2's table starts at 40 C: Psat = 7.3851 kPa, h_fg = 2406.0 kJ/kg
    (Cengel & Boles 8e, A-4)."""
    t_bub, _ = pp.bubble_dew(7385.1, Z)
    assert t_bub == pytest.approx(313.15, abs=0.05)
    sat = pp.bubble_point(7385.1, Z)
    h_fg = (sat.H_vapor - sat.H_liquid) / pp.mw / 1e3
    assert h_fg == pytest.approx(2406.0, rel=1e-3)


def test_above_critical_pressure_has_no_saturation(pp):
    with pytest.raises(ValueError, match="critical pressure"):
        pp.bubble_dew(25e6, Z)


# -- flashes -------------------------------------------------------------------
def test_flash_pt_phase_identification(pp):
    liq = pp.flash_pt(300.0, 1e5, Z)
    assert liq.phase == "liquid" and liq.vapor_fraction == 0.0
    assert liq.x == {"water": 1.0} and liq.y is None
    vap = pp.flash_pt(400.0, 1e5, Z)
    assert vap.phase == "vapor" and vap.vapor_fraction == 1.0
    assert vap.y == {"water": 1.0} and vap.x is None
    # Supercritical states classify by density without crashing.
    sc = pp.flash_pt(700.0, 30e6, Z)
    assert sc.phase in ("vapor", "liquid")


def test_flash_ph_roundtrip_and_two_phase(pp):
    h = pp.enthalpy(400.0, 1e5, Z)
    r = pp.flash_ph(1e5, h, Z)
    assert r.T == pytest.approx(400.0, abs=1e-6)
    assert r.phase == "vapor"

    sat = pp.bubble_point(P_ATM, Z)
    r2 = pp.flash_ph(P_ATM, 0.5 * (sat.H_liquid + sat.H_vapor), Z)
    assert r2.phase == "VLE"
    assert r2.vapor_fraction == pytest.approx(0.5, abs=1e-9)
    assert r2.T == pytest.approx(373.12, abs=0.01)
    assert r2.H_liquid == pytest.approx(sat.H_liquid, rel=1e-12)
    assert r2.H_vapor == pytest.approx(sat.H_vapor, rel=1e-12)


def test_flash_ps_roundtrip(pp):
    s = pp.entropy(400.0, 1e5, Z)
    assert pp.flash_ps(1e5, s, Z).T == pytest.approx(400.0, abs=1e-6)


# -- formation-inclusive enthalpy basis ----------------------------------------
def test_ideal_gas_anchor_is_the_formation_enthalpy(pp):
    """At 298.15 K and P -> 0 the absolute enthalpy is water's ideal-gas
    formation enthalpy, -241,822 J/mol (`chemicals` Hfg — the same constant the
    thermo backends fold in; see _flasher.py)."""
    assert pp.enthalpy(298.15, 100.0, Z) == pytest.approx(-241_822.0, abs=5.0)


def test_absolute_basis_agrees_with_thermo_pr_for_steam(pp):
    """Steam at 400 K / 1 bar: PR and IAPWS-95 absolute enthalpies agree to
    ~100 J/mol (0.04%) once both sit on the formation-inclusive basis."""
    pr = make_package("thermo:PR", ["water"])
    assert pp.enthalpy(400.0, 1e5, Z) == pytest.approx(
        pr.enthalpy(400.0, 1e5, Z), abs=300.0
    )


def _heater_duty(package: str, t_in: float, t_out: float, p: float) -> float:
    fs = Flowsheet(components=[Component("water")], property_package=package)
    fs.add(Heater("H", {"T_out": t_out}))
    fs.feed("FEED", "H:in1", T=t_in, P=p, molar_flow=1.0, z=dict(Z))
    fs.connect("OUT", "H:out", None)
    fs.connect("Q", "H:duty", None)
    rep = fs.solve()
    assert rep.converged
    return rep.duties["Q"]


def test_heater_duty_cross_check_vapor_phase():
    """The cross-check: the same Heater solved under thermo:PR and under
    coolprop:Water. Steam 400 -> 500 K at 1 bar — measured delta 2.20%
    (PR's vapor Cp/departure error vs IAPWS-95); asserted < 3%."""
    q_pr = _heater_duty("thermo:PR", 400.0, 500.0, 1e5)
    q_cp = _heater_duty("coolprop:Water", 400.0, 500.0, 1e5)
    assert q_cp > 0.0
    assert q_pr == pytest.approx(q_cp, rel=0.03)


def test_heater_duty_cross_check_liquid_phase():
    """Liquid water 300 -> 350 K at 5 bar — measured delta 7.9% (a cubic EOS
    is a poor liquid-water heat-capacity model; that is precisely the gap the
    steam-tables package closes); asserted < 10%."""
    q_pr = _heater_duty("thermo:PR", 300.0, 350.0, 5e5)
    q_cp = _heater_duty("coolprop:Water", 300.0, 350.0, 5e5)
    assert q_cp > 0.0
    assert q_pr == pytest.approx(q_cp, rel=0.10)


def test_forced_phase_properties_and_k_values(pp):
    """Per-phase protocol extensions: stable-phase values where the phase is
    stable, clamped to saturation otherwise; K = Psat/P for a pure fluid."""
    sat = pp.bubble_point(P_ATM, Z)
    # in the stable region the forced value equals the bulk value
    assert pp.enthalpy_liquid(300.0, 1e5, Z) == pytest.approx(
        pp.enthalpy(300.0, 1e5, Z), rel=1e-12)
    assert pp.enthalpy_vapor(400.0, 1e5, Z) == pytest.approx(
        pp.enthalpy(400.0, 1e5, Z), rel=1e-12)
    # outside it, clamp to the saturated phase at P
    assert pp.enthalpy_vapor(300.0, P_ATM, Z) == pytest.approx(sat.H_vapor, rel=1e-12)
    assert pp.enthalpy_liquid(450.0, P_ATM, Z) == pytest.approx(sat.H_liquid, rel=1e-12)
    assert pp.volume_vapor(400.0, 1e5, Z) == pytest.approx(
        pp.volume(400.0, 1e5, Z), rel=1e-12)
    assert pp.volume_liquid(300.0, 1e5, Z) > 0
    # K(Tsat(P), P) = 1 for a pure fluid
    t_sat, _ = pp.bubble_dew(P_ATM, Z)
    assert pp.k_values(t_sat, P_ATM, Z, Z)["water"] == pytest.approx(1.0, rel=1e-6)


def test_volume_matches_liquid_water_density(pp):
    """rho(300 K, 1 bar) = 996.56 kg/m^3 (IAPWS; steam tables ~996.5)."""
    rho = pp.mw / pp.volume(300.0, 1e5, Z)
    assert rho == pytest.approx(996.56, rel=1e-3)
