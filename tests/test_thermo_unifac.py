"""Predictive UNIFAC gamma-phi package tests.

The headline: Modified UNIFAC (Dortmund) is a *predictive* simultaneous VLE+LLE
liquid model — it reproduces the ethanol/water minimum-boiling azeotrope (VLE)
AND the liquid-liquid miscibility gap that lets a heterogeneous azeotrope decant
(LLE), from one consistent set of group-contribution parameters with no
per-system fitting. This is the model that unblocks integrated heteroazeotropic
(3-phase) distillation, where a single NRTL set structurally cannot fit both.

References:
  * ethanol/water azeotrope ~89.4 mol% ethanol at 78.2 C / 351.3 K, 1 atm
    (Gmehling, DECHEMA VLE Data Collection; CRC).
  * ethanol/water/cyclohexane ternary decant — the entrainer of anhydrous
    ethanol dehydration (Hameed 2025 sec. 9.5.6).
  * water/n-butanol LLE — mutual solubilities, organic ~50 mol% water, aqueous
    ~98 mol% water at 25 C (Sorensen & Arlt, DECHEMA LLE Data Collection).
"""
import math

import numpy as np
import pytest

from caldyr.core import Component, Flowsheet
from caldyr.thermo import make_package

P_ATM = 101325.0
AZ_T = 351.3      # K, ethanol/water azeotrope temperature at 1 atm


def _min_bubble(spec: str, a: str = "ethanol", b: str = "water"):
    """Return (x_a, T) at the minimum bubble-point temperature of the a/b pair."""
    pp = make_package(spec, [a, b])
    xs = np.linspace(0.02, 0.98, 49)
    temps = [pp.bubble_dew(P_ATM, {a: float(x), b: 1 - float(x)})[0] for x in xs]
    i = int(np.argmin(temps))
    return float(xs[i]), float(temps[i])


def test_unifac_captures_ethanol_water_azeotrope():
    """Modified UNIFAC predicts the minimum-boiling azeotrope without any fitted
    binary data — and lands closer to experiment than ChemSep NRTL."""
    pp = make_package("thermo:UNIFAC", ["ethanol", "water"])
    x_az, t_az = _min_bubble("thermo:UNIFAC")

    assert 0.85 < x_az < 0.93        # experimental ~0.894 mol ethanol
    assert t_az == pytest.approx(AZ_T, abs=1.0)

    # It must boil below both pure components (the defining azeotrope property).
    bp_eth = pp.bubble_dew(P_ATM, {"ethanol": 1.0})[0]
    bp_wat = pp.bubble_dew(P_ATM, {"water": 1.0})[0]
    assert t_az < bp_eth < bp_wat


def test_unifac_ternary_heteroazeotrope_decants():
    """The ethanol/water/cyclohexane overhead splits into a cyclohexane-rich
    organic layer (the recycled entrainer) and an ethanol/water aqueous layer —
    the boundary crossing that azeotropic dehydration relies on."""
    pp = make_package("thermo:UNIFAC", ["ethanol", "water", "cyclohexane"])
    z = {"ethanol": 0.3, "water": 0.1, "cyclohexane": 0.6}
    r = pp.flash_pt_3p(320.0, P_ATM, z)

    # Two liquids, no vapour at 320 K / 1 atm.
    assert r.beta_vapor == pytest.approx(0.0, abs=1e-6)
    assert r.beta_light > 0.1 and r.beta_heavy > 0.1
    assert r.beta_light + r.beta_heavy == pytest.approx(1.0, abs=1e-9)

    # One layer is the entrainer-rich organic; the other carries the ethanol.
    cx = [r.x_light["cyclohexane"], r.x_heavy["cyclohexane"]]
    organic, aqueous = (r.x_light, r.x_heavy) if cx[0] > cx[1] else (r.x_heavy, r.x_light)
    assert organic["cyclohexane"] > 0.8        # entrainer concentrates here
    assert aqueous["ethanol"] > organic["ethanol"]   # ethanol favours the other

    # Overall component balance closes over the two liquids (no vapour).
    for c in z:
        recombined = r.beta_light * r.x_light[c] + r.beta_heavy * r.x_heavy[c]
        assert recombined == pytest.approx(z[c], abs=1e-6)


def test_unifac_ternary_split_is_isoactivity():
    """The split is a true liquid-liquid equilibrium: a_i = gamma_i x_i is equal
    in both layers for every component (the LLE criterion)."""
    pp = make_package("thermo:UNIFAC", ["ethanol", "water", "cyclohexane"])
    z = {"ethanol": 0.3, "water": 0.1, "cyclohexane": 0.6}
    r = pp.flash_pt_3p(320.0, P_ATM, z)

    xI = [r.x_light[c] for c in pp.components]
    xII = [r.x_heavy[c] for c in pp.components]
    gI = pp._gammas(320.0, xI)
    gII = pp._gammas(320.0, xII)
    for i in range(len(pp.components)):
        assert gI[i] * xI[i] == pytest.approx(gII[i] * xII[i], rel=1e-6, abs=1e-9)


def test_unifac_opens_water_butanol_gap_that_nrtl_cannot():
    """Water/n-butanol partially miscible under UNIFAC (predictive LLE) but a
    single mixture, no gap, under the VLE-fit ChemSep NRTL — the structural
    limitation UNIFAC removes."""
    z = {"water": 0.6, "1-butanol": 0.4}
    pp = make_package("thermo:UNIFAC", ["water", "1-butanol"])
    r = pp.flash_pt_3p(298.15, P_ATM, z)
    assert r.beta_light > 0.05 and r.beta_heavy > 0.05      # two liquids
    aqueous = r.x_heavy if r.x_heavy["water"] > r.x_light["water"] else r.x_light
    assert aqueous["water"] > 0.9                           # nearly pure water layer

    # ChemSep NRTL: the same flash returns a single liquid (no miscibility gap).
    pp_nrtl = make_package("thermo:NRTL", ["water", "1-butanol"])
    rn = pp_nrtl.flash_pt_3p(298.15, P_ATM, z)
    assert rn.beta_heavy == pytest.approx(0.0, abs=1e-9)


def test_unifac_lle_selector_also_splits():
    """The secondary 'thermo:UNIFAC-LLE' selector (Magnussen LLE table, built via
    the original->LLE subgroup remap) likewise opens the water/butanol gap."""
    pp = make_package("thermo:UNIFAC-LLE", ["water", "1-butanol"])
    r = pp.flash_pt_3p(298.15, P_ATM, {"water": 0.6, "1-butanol": 0.4})
    assert r.beta_light > 0.05 and r.beta_heavy > 0.05
    aqueous = r.x_heavy if r.x_heavy["water"] > r.x_light["water"] else r.x_light
    assert aqueous["water"] > 0.9


def test_unifac_fully_miscible_degrades_to_one_liquid():
    """No spurious gap: an ideal-ish, fully miscible pair stays one liquid."""
    pp = make_package("thermo:UNIFAC", ["benzene", "toluene"])
    r = pp.flash_pt_3p(360.0, P_ATM, {"benzene": 0.5, "toluene": 0.5})
    assert r.beta_heavy == pytest.approx(0.0, abs=1e-9)


def test_unifac_single_component_boiling_points():
    """Pure-component path (FlashPureVLS) against textbook normal BPs."""
    pp_w = make_package("thermo:UNIFAC", ["water"])
    pp_c = make_package("thermo:UNIFAC", ["cyclohexane"])
    assert pp_w.bubble_dew(P_ATM, {"water": 1.0})[0] == pytest.approx(373.15, abs=1.5)
    # cyclohexane normal BP 80.7 C = 353.9 K
    assert pp_c.bubble_dew(P_ATM, {"cyclohexane": 1.0})[0] == pytest.approx(353.9, abs=1.5)


def test_unifac_unfragmentable_component_raises():
    """A species the DDBST database cannot fragment (a monatomic/permanent gas)
    fails loudly at construction — no silent gamma = 1."""
    with pytest.raises(ValueError, match="group assignment"):
        make_package("thermo:UNIFAC", ["argon", "water"])


def test_unifac_energy_balance_closes_in_a_flowsheet():
    """The package stays consistent with the solver and unit ops."""
    fs = Flowsheet(components=[Component("ethanol"), Component("water"),
                              Component("cyclohexane")],
                   property_package="thermo:UNIFAC")
    from caldyr.unitops import Heater, Mixer
    fs.add(Mixer("MIX", {"dP": 0.0}))
    fs.add(Heater("H", {"T_out": 340.0, "dP": 0.0}))
    fs.feed("S1", "MIX:in1", T=298.15, P=P_ATM, molar_flow=10.0,
            z={"ethanol": 0.5, "water": 0.3, "cyclohexane": 0.2})
    fs.feed("S2", "MIX:in2", T=320.0, P=P_ATM, molar_flow=5.0,
            z={"cyclohexane": 1.0})
    fs.connect("S3", "MIX:out", "H:in1")
    fs.connect("S4", "H:out", None)
    fs.connect("Q", "H:duty", None)
    report = fs.solve()

    h_in = sum(fs.streams[s].molar_flow * fs.streams[s].H for s in ("S1", "S2"))
    h_out = fs.streams["S4"].molar_flow * fs.streams["S4"].H
    assert math.isclose(h_in + report.duties["Q"], h_out, rel_tol=1e-9, abs_tol=1e-3)
