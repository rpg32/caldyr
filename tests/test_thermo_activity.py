"""Activity-coefficient package tests.

The headline: a cubic EOS (PR) cannot represent the ethanol/water minimum-boiling
azeotrope, while the NRTL gamma-phi package does. Reference: the ethanol/water
azeotrope is ~89 mol% ethanol at 78.2 C / 351.3 K, 1 atm (Gmehling, DECHEMA VLE
Data Collection; CRC). NRTL binary parameters come from ChemSep (bundled with
`thermo`).
"""
import math

import numpy as np
import pytest

from caldyr.core import Component, Flowsheet
from caldyr.thermo import make_package
from caldyr.thermo.activity_pkg import _has_offdiagonal

P_ATM = 101325.0
AZ_T = 351.3      # K, ethanol/water azeotrope temperature at 1 atm


def _min_bubble(spec: str):
    """Return (x_ethanol, T) at the minimum bubble-point temperature."""
    pp = make_package(spec, ["ethanol", "water"])
    xs = np.linspace(0.02, 0.98, 49)
    temps = [pp.bubble_dew(P_ATM, {"ethanol": float(x), "water": 1 - float(x)})[0] for x in xs]
    i = int(np.argmin(temps))
    return float(xs[i]), float(temps[i])


def test_nrtl_captures_ethanol_water_azeotrope():
    pp = make_package("thermo:NRTL", ["ethanol", "water"])
    x_az, t_az = _min_bubble("thermo:NRTL")

    # An interior minimum in bubble temperature *is* a minimum-boiling azeotrope.
    assert 0.80 < x_az < 0.95
    assert t_az == pytest.approx(AZ_T, abs=1.5)

    # It must boil below both pure components.
    bp_eth = pp.bubble_dew(P_ATM, {"ethanol": 1.0})[0]
    bp_wat = pp.bubble_dew(P_ATM, {"water": 1.0})[0]
    assert t_az < bp_eth < bp_wat


def test_pr_misses_the_azeotrope():
    """Contrast: the cubic EOS has no interior bubble-T minimum -> no azeotrope."""
    x_az, _ = _min_bubble("thermo:PR")
    assert not (0.05 < x_az < 0.95)   # minimum sits at a (near-)pure endpoint


def test_nrtl_vapor_enriches_then_crosses_at_azeotrope():
    """Below the azeotrope vapor is richer in ethanol; above it, leaner."""
    pp = make_package("thermo:NRTL", ["ethanol", "water"])

    def y_minus_x(x):
        b, d = pp.bubble_dew(P_ATM, {"ethanol": x, "water": 1 - x})
        r = pp.flash_pt((b + d) / 2 if d > b else b + 0.02, P_ATM,
                        {"ethanol": x, "water": 1 - x})
        return r.y["ethanol"] - x

    assert y_minus_x(0.30) > 0.0     # vapor enriched in ethanol
    assert y_minus_x(0.97) < 0.0     # past the azeotrope, vapor leaner


def test_nrtl_energy_balance_closes_in_a_flowsheet():
    """The new package stays consistent with the solver and unit ops."""
    fs = Flowsheet(components=[Component("water"), Component("ethanol")],
                   property_package="thermo:NRTL")
    from caldyr.unitops import Heater, Mixer
    fs.add(Mixer("MIX", {"dP": 0.0}))
    fs.add(Heater("H", {"T_out": 350.0, "dP": 0.0}))
    fs.feed("S1", "MIX:in1", T=298.15, P=P_ATM, molar_flow=10.0,
            z={"water": 0.6, "ethanol": 0.4})
    fs.feed("S2", "MIX:in2", T=320.0, P=P_ATM, molar_flow=5.0, z={"water": 1.0})
    fs.connect("S3", "MIX:out", "H:in1")
    fs.connect("S4", "H:out", None)
    fs.connect("Q", "H:duty", None)
    report = fs.solve()

    h_in = sum(fs.streams[s].molar_flow * fs.streams[s].H for s in ("S1", "S2"))
    h_out = fs.streams["S4"].molar_flow * fs.streams["S4"].H
    assert math.isclose(h_in + report.duties["Q"], h_out, rel_tol=1e-9, abs_tol=1e-3)


def test_nrtl_single_component_boiling_points():
    """Pure-component path (FlashPureVLS) against textbook normal BPs."""
    pp_w = make_package("thermo:NRTL", ["water"])
    pp_e = make_package("thermo:NRTL", ["ethanol"])
    assert pp_w.bubble_dew(P_ATM, {"water": 1.0})[0] == pytest.approx(373.15, abs=1.5)
    assert pp_e.bubble_dew(P_ATM, {"ethanol": 1.0})[0] == pytest.approx(351.44, abs=1.5)


def test_has_offdiagonal_detects_ideal_fallback():
    """The zero-parameter detector that triggers the ideal-solution warning."""
    assert not _has_offdiagonal([[0.0, 0.0], [0.0, 0.0]])
    assert _has_offdiagonal([[0.0, 12.3], [-4.5, 0.0]])
    assert not _has_offdiagonal([[5.0]])   # single component: no off-diagonal

