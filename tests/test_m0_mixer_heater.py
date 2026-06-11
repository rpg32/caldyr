"""M0 acceptance tests.

Property reference: Peng-Robinson EOS via `thermo` (Caleb Bell), with ideal-gas
heat capacities from the same package — enthalpies share one internally
consistent reference state, so flowsheet energy balances close to machine
precision. Spot-checks below are validated against textbook normal boiling
points (CRC/NIST): water 373.15 K, ethanol 351.44 K at 101.325 kPa.
"""
import math

import pytest

from caldyr.core import Component, Flowsheet
from caldyr.io import from_dict, to_dict
from caldyr.thermo import make_package
from caldyr.unitops import Heater, Mixer

P_ATM = 101325.0


def build_flowsheet() -> Flowsheet:
    fs = Flowsheet(
        components=[Component("water"), Component("ethanol")],
        property_package="thermo:PR",
    )
    fs.add(Mixer("MIX1", {"dP": 0.0}))
    fs.add(Heater("H1", {"T_out": 350.0, "dP": 0.0}))
    fs.feed("S1", "MIX1:in1", T=298.15, P=P_ATM, molar_flow=10.0,
            z={"water": 0.6, "ethanol": 0.4})
    fs.feed("S2", "MIX1:in2", T=320.0, P=P_ATM, molar_flow=5.0, z={"water": 1.0})
    fs.connect("S3", "MIX1:out", "H1:in1")
    fs.connect("S4", "H1:out", None)
    fs.connect("Q1", "H1:duty", None)
    return fs


def test_mixer_heater_energy_balance_closes():
    """Σ(inlet enthalpy flow) + heater duty == outlet enthalpy flow."""
    fs = build_flowsheet()
    report = fs.solve()
    assert report.converged
    assert report.order == ["MIX1", "H1"]

    h_in = sum(fs.streams[s].molar_flow * fs.streams[s].H for s in ("S1", "S2"))
    q = report.duties["Q1"]
    h_out = fs.streams["S4"].molar_flow * fs.streams["S4"].H
    # Closes to machine precision (shared enthalpy reference).
    assert math.isclose(h_in + q, h_out, rel_tol=1e-9, abs_tol=1e-3)


def test_mass_balance_closes():
    fs = build_flowsheet()
    fs.solve()
    n_in = fs.streams["S1"].molar_flow + fs.streams["S2"].molar_flow
    assert math.isclose(fs.streams["S4"].molar_flow, n_in, rel_tol=1e-12)
    # per-component water balance: 10*0.6 + 5*1.0 = 11 mol/s water out of 15
    z4 = fs.streams["S4"].normalized_z()
    assert math.isclose(z4["water"] * 15.0, 11.0, rel_tol=1e-9)


def test_heater_duty_matches_independent_enthalpy_difference():
    """Heater duty equals n*(H_out - H_in) recomputed straight from the
    property package — guards the solver's stream propagation."""
    fs = build_flowsheet()
    report = fs.solve()
    pp = make_package("thermo:PR", fs.component_ids)
    s3, s4 = fs.streams["S3"], fs.streams["S4"]
    h_in = pp.enthalpy(s3.T, s3.P, s3.normalized_z())
    h_out = pp.enthalpy(s4.T, s4.P, s4.normalized_z())
    q_independent = s4.molar_flow * (h_out - h_in)
    assert math.isclose(report.duties["Q1"], q_independent, rel_tol=1e-9)


def test_normal_boiling_points_against_textbook():
    """Validate the property package against CRC/NIST normal boiling points."""
    pp = make_package("thermo:PR", ["water", "ethanol"])
    water_bp, _ = pp.bubble_dew(P_ATM, {"water": 1.0})
    ethanol_bp, _ = pp.bubble_dew(P_ATM, {"ethanol": 1.0})
    assert water_bp == pytest.approx(373.15, abs=2.0)      # 100 C
    assert ethanol_bp == pytest.approx(351.44, abs=2.0)    # 78.3 C


def test_heater_cooling_with_fixed_duty():
    """A negative Q cools; T_out/Q specs are interchangeable and consistent."""
    fs = Flowsheet(components=[Component("water")], property_package="thermo:PR")
    fs.add(Heater("H", {"Q": -5.0e3, "dP": 0.0}))  # ~33 K of cooling on 2 mol/s
    fs.feed("F", "H:in1", T=350.0, P=P_ATM, molar_flow=2.0, z={"water": 1.0})
    fs.connect("OUT", "H:out", None)
    fs.connect("D", "H:duty", None)
    report = fs.solve()
    assert report.duties["D"] == pytest.approx(-5.0e3)
    assert fs.streams["OUT"].T < 350.0  # cooled


def test_flow_round_trip_is_exact():
    fs = build_flowsheet()
    fs.solve()
    doc = to_dict(fs)
    again = to_dict(from_dict(doc))
    assert doc == again


def test_acyclic_solve_reports_single_sweep():
    """An acyclic flowsheet takes the fast path: one sweep, no tear streams."""
    fs = build_flowsheet()
    report = fs.solve()
    assert report.iterations == 1
    assert report.tear_streams == []
