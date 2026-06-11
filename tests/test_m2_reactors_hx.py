"""M2 acceptance tests: reactors and the heat exchanger.

Property/thermo references:
  * Heat of reaction N2 + 3H2 -> 2NH3 = -91.8 kJ per 2 mol NH3 at 298 K (CRC).
  * Ammonia equilibrium conversion: ~0.34-0.40 mole-fraction NH3 at 400 C,
    200 bar, stoichiometric feed (classic Haber-Bosch equilibrium data; e.g.
    Smith, Van Ness & Abbott, *Intro. to Chemical Engineering Thermodynamics*).
  * Heat exchanger: effectiveness-NTU and LMTD methods must agree when heat
    capacities are constant (Incropera, *Fundamentals of Heat & Mass Transfer*).
"""
import math

import pytest

from caldyr.core import Component, Flowsheet
from caldyr.thermo import make_package
from caldyr.unitops import HeatExchanger

P_ATM = 101325.0
AMMONIA = {"stoich": {"nitrogen": -1, "hydrogen": -3, "ammonia": 2}, "key": "nitrogen"}


def syngas_flowsheet(reactor) -> Flowsheet:
    fs = Flowsheet(
        components=[Component("nitrogen"), Component("hydrogen"), Component("ammonia")],
        property_package="thermo:PR",
    )
    fs.add(reactor)
    fs.feed("F", "R:in1", T=673.15, P=2e7, molar_flow=4.0,
            z={"nitrogen": 0.25, "hydrogen": 0.75, "ammonia": 0.0})
    fs.connect("O", "R:out", None)
    fs.connect("Q", "R:duty", None)
    return fs


def _atoms_close(stream):
    n = stream.molar_flow
    n_atoms = 2 * n * stream.z["nitrogen"] + n * stream.z["ammonia"]
    h_atoms = 2 * n * stream.z["hydrogen"] + 3 * n * stream.z["ammonia"]
    return n_atoms, h_atoms


# -- conversion reactor ----------------------------------------------------
def test_conversion_reactor_extent_and_atom_balance():
    from caldyr.unitops import ConversionReactor
    fs = syngas_flowsheet(ConversionReactor("R", {"reaction": AMMONIA,
                                                  "conversion": 0.25, "T_out": 673.15}))
    fs.solve()
    out = fs.streams["O"]
    # extent xi = 0.25 * n_N2_in / 1 = 0.25; total moles drop by 2*xi = 0.5
    assert math.isclose(out.molar_flow, 3.5, rel_tol=1e-9)
    assert math.isclose(out.z["ammonia"] * out.molar_flow, 0.5, rel_tol=1e-6)  # 2*xi
    n_atoms, h_atoms = _atoms_close(out)
    assert math.isclose(n_atoms, 2.0, rel_tol=1e-9)   # feed N atoms: 2*1
    assert math.isclose(h_atoms, 6.0, rel_tol=1e-9)   # feed H atoms: 2*3


def test_conversion_reactor_adiabatic_heats_up():
    from caldyr.unitops import ConversionReactor
    fs = syngas_flowsheet(ConversionReactor("R", {"reaction": AMMONIA, "conversion": 0.25}))
    report = fs.solve()
    assert report.duties["Q"] == pytest.approx(0.0, abs=1e-6)   # adiabatic
    assert fs.streams["O"].T > 673.15 + 50                      # exothermic rise


def test_conversion_reactor_isothermal_duty_equals_heat_of_reaction():
    from caldyr.unitops import ConversionReactor
    fs = syngas_flowsheet(ConversionReactor("R", {"reaction": AMMONIA,
                                                  "conversion": 0.25, "T_out": 673.15}))
    report = fs.solve()
    pp = make_package("thermo:PR", ["nitrogen", "hydrogen", "ammonia"])
    T, P, xi = 673.15, 2e7, 0.25
    dH_rxn = (2 * pp.enthalpy(T, P, {"ammonia": 1.0})
              - pp.enthalpy(T, P, {"nitrogen": 1.0})
              - 3 * pp.enthalpy(T, P, {"hydrogen": 1.0}))
    assert report.duties["Q"] < 0                               # exothermic: remove heat
    # ~3% gap vs the pure-component heat of reaction is real-gas departure /
    # enthalpy of mixing at 200 bar — the duty is the rigorous mixture value.
    assert report.duties["Q"] == pytest.approx(xi * dH_rxn, rel=0.05)


# -- equilibrium reactor ---------------------------------------------------
def test_equilibrium_reactor_matches_haber_conversion():
    from caldyr.unitops import EquilibriumReactor
    fs = syngas_flowsheet(EquilibriumReactor("R", {"reaction": AMMONIA, "T": 673.15}))
    fs.solve()
    y_nh3 = fs.streams["O"].z["ammonia"]
    assert y_nh3 == pytest.approx(0.37, abs=0.06)    # ~0.34-0.40 at 400 C, 200 bar
    n_atoms, h_atoms = _atoms_close(fs.streams["O"])
    assert math.isclose(n_atoms, 2.0, rel_tol=1e-9)
    assert math.isclose(h_atoms, 6.0, rel_tol=1e-9)


def test_equilibrium_pressure_and_temperature_trends():
    from caldyr.unitops import EquilibriumReactor

    def y_nh3(T, P):
        fs = Flowsheet(
            components=[Component("nitrogen"), Component("hydrogen"), Component("ammonia")],
            property_package="thermo:PR")
        fs.add(EquilibriumReactor("R", {"reaction": AMMONIA, "T": T}))
        fs.feed("F", "R:in1", T=T, P=P, molar_flow=4.0,
                z={"nitrogen": 0.25, "hydrogen": 0.75, "ammonia": 0.0})
        fs.connect("O", "R:out", None)
        fs.connect("Q", "R:duty", None)
        fs.solve()
        return fs.streams["O"].z["ammonia"]

    # Le Chatelier: more NH3 at higher pressure (dn<0) and lower temperature (exo).
    assert y_nh3(673.15, 3e7) > y_nh3(673.15, 2e7) > y_nh3(673.15, 1e7)
    assert y_nh3(623.15, 2e7) > y_nh3(723.15, 2e7)


def test_equilibrium_constant_falls_with_temperature():
    pp = make_package("thermo:PR", ["nitrogen", "hydrogen", "ammonia"])
    stoich = AMMONIA["stoich"]
    assert pp.lnKeq(stoich, 400.0) > pp.lnKeq(stoich, 700.0)   # exothermic


# -- heat exchanger --------------------------------------------------------
def water_hx(params) -> Flowsheet:
    fs = Flowsheet(components=[Component("water")], property_package="thermo:PR")
    fs.add(HeatExchanger("HX", params))
    fs.feed("H", "HX:hot_in", T=370.0, P=5e5, molar_flow=10.0, z={"water": 1.0})
    fs.feed("C", "HX:cold_in", T=300.0, P=5e5, molar_flow=8.0, z={"water": 1.0})
    fs.connect("HO", "HX:hot_out", None)
    fs.connect("CO", "HX:cold_out", None)
    return fs


def test_hx_duty_mode_conserves_energy():
    fs = water_hx({"duty": 4.0e4})
    fs.solve()
    h, c, ho, co = (fs.streams[k] for k in ("H", "C", "HO", "CO"))
    q_hot = h.molar_flow * (h.H - ho.H)
    q_cold = co.molar_flow * (co.H - c.H)
    assert q_hot == pytest.approx(4.0e4, rel=1e-9)
    assert q_cold == pytest.approx(4.0e4, rel=1e-9)   # all hot duty reaches cold


def test_hx_outlet_temperature_spec():
    fs = water_hx({"T_hot_out": 340.0})
    fs.solve()
    assert fs.streams["HO"].T == pytest.approx(340.0, abs=1e-3)
    h, c, ho, co = (fs.streams[k] for k in ("H", "C", "HO", "CO"))
    assert math.isclose(h.molar_flow * (h.H - ho.H),
                        co.molar_flow * (co.H - c.H), rel_tol=1e-9)


def test_hx_eps_ntu_agrees_with_lmtd_for_constant_cp():
    """Gas streams (no phase change): the two HX methods must agree."""
    fs = Flowsheet(components=[Component("nitrogen")], property_package="thermo:PR")
    fs.add(HeatExchanger("HX", {"UA": 600.0, "arrangement": "counterflow"}))
    fs.feed("H", "HX:hot_in", T=500.0, P=5e5, molar_flow=10.0, z={"nitrogen": 1.0})
    fs.feed("C", "HX:cold_in", T=300.0, P=5e5, molar_flow=12.0, z={"nitrogen": 1.0})
    fs.connect("HO", "HX:hot_out", None)
    fs.connect("CO", "HX:cold_out", None)
    fs.solve()
    h, c, ho, co = (fs.streams[k] for k in ("H", "C", "HO", "CO"))
    q = h.molar_flow * (h.H - ho.H)
    lmtd = HeatExchanger.lmtd(h.T - co.T, ho.T - c.T)   # counterflow approaches
    assert q == pytest.approx(600.0 * lmtd, rel=0.02)   # within Cp(T) drift
    assert co.T <= h.T and ho.T >= c.T                  # second-law limits


def test_lmtd_formula_and_equal_approach_limit():
    assert HeatExchanger.lmtd(50.0, 30.0) == pytest.approx(39.15, abs=0.05)
    assert HeatExchanger.lmtd(20.0, 20.0) == pytest.approx(20.0)   # equal-approach limit
    with pytest.raises(ValueError):
        HeatExchanger.lmtd(10.0, -5.0)


# -- the M2 deliverable: a converging ammonia synthesis loop ----------------
def ammonia_loop() -> Flowsheet:
    from caldyr.unitops import EquilibriumReactor, FlashDrum, Heater, Mixer, Splitter
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


def test_ammonia_loop_converges_and_closes_atom_balances():
    fs = ammonia_loop()
    report = fs.solve(tol=1e-7, max_iter=400)
    assert report.converged
    assert report.tear_streams == ["RECYCLE"]

    s = fs.streams

    def atoms(coeff):
        def total(sid):
            return sum(s[sid].molar_flow * s[sid].z.get(c, 0.0) * k for c, k in coeff.items())
        return total("MAKEUP"), total("PRODUCT") + total("PURGE")

    for coeff in ({"nitrogen": 2, "ammonia": 1},     # N atoms
                  {"hydrogen": 2, "ammonia": 3},     # H atoms
                  {"argon": 1}):                     # inert: in == out
        inp, out = atoms(coeff)
        assert math.isclose(inp, out, rel_tol=1e-5)

    # The separator condenses a high-purity ammonia product.
    assert s["PRODUCT"].z["ammonia"] > 0.95
    # Inert argon is enriched in the loop relative to the makeup (it accumulates
    # until the purge removes it as fast as it enters).
    assert s["RECYCLE"].z["argon"] > s["MAKEUP"].z["argon"]
