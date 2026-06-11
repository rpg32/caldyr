"""M12 tests: the Balance logical unit op (HYSYS Balance, registered "Balance").

Reference: Hameed, *Chemical Process Simulations using Aspen HYSYS* (Wiley
2025), §6.3 "Balance" — the five balance types of the Parameters tab
(Component Mole Flow / Mass Flow / Heat Flow / Component Mole and Heat Flow /
Mass and Heat Flow) map to ``params['mode']`` = mole / mass / heat /
mole_heat / mass_heat. §6.3 is descriptive (it tabulates what each type
conserves but works no numeric example), so each test below asserts exactly
the conservation semantics the book states, including its §6.3.2 use case —
"the modeling of reactors for which ... stoichiometry is unknown. Alkylation
units ... frequently use this application" — reproduced as an alkylation
mass balance (iC4 + C4= -> iC8: mass conserved, moles and species not).
"""
import pytest

from caldyr.core import Component, Flowsheet
from caldyr.core.components_db import molar_mass
from caldyr.unitops import Balance, Mixer


def _fs(component_ids, package="thermo:PR"):
    return Flowsheet(components=[Component(c) for c in component_ids],
                     property_package=package)


# -- ports are param-driven ----------------------------------------------------
def test_param_driven_ports():
    b3 = Balance("B", {"mode": "mole", "n_inlets": 3})
    assert [p.name for p in b3.ports] == ["in1", "in2", "in3", "out1"]
    assert all(p.kind == "material" for p in b3.ports)
    bh = Balance("B", {"mode": "heat", "n_inlets": 1})
    out = next(p for p in bh.ports if p.name == "out1")
    assert out.kind == "energy"


def test_param_errors_are_typed():
    with pytest.raises(ValueError, match="n_inlets"):
        Balance("B", {"mode": "mole", "n_inlets": 0})
    with pytest.raises(ValueError, match="unknown mode"):
        Balance("B", {"mode": "enthalpy"})
    fs = _fs(["methane"])
    fs.add(Balance("B", {"mode": "mass", "n_inlets": 1}))   # no z_out
    fs.feed("F", "B:in1", T=300.0, P=1e5, molar_flow=1.0, z={"methane": 1.0})
    fs.connect("OUT", "B:out1", None)
    with pytest.raises(ValueError, match="z_out"):
        fs.solve()


def test_missing_inlet_is_a_typed_error():
    fs = _fs(["methane"])
    fs.add(Balance("B", {"mode": "mole", "n_inlets": 2}))
    fs.feed("F", "B:in1", T=300.0, P=1e5, molar_flow=1.0, z={"methane": 1.0})
    fs.connect("OUT", "B:out1", None)
    with pytest.raises(ValueError, match="in2"):
        fs.solve()


# -- mode "mole" (book §6.3.1: Component Mole Flow) -----------------------------
def test_mole_mode_conserves_component_moles_and_passes_no_state():
    fs = _fs(["methane", "ethane"])
    fs.add(Balance("B", {"mode": "mole", "n_inlets": 2}))
    fs.feed("F1", "B:in1", T=300.0, P=10e5, molar_flow=2.0,
            z={"methane": 1.0, "ethane": 0.0})
    fs.feed("F2", "B:in2", T=450.0, P=2e5, molar_flow=3.0,
            z={"methane": 0.2, "ethane": 0.8})
    fs.connect("OUT", "B:out1", None)
    rep = fs.solve()
    assert rep.converged
    out = fs.streams["OUT"]
    assert out.molar_flow == pytest.approx(5.0)
    # per-component mole balance: CH4 2.0 + 0.6, C2H6 2.4
    assert out.z["methane"] == pytest.approx(2.6 / 5.0)
    assert out.z["ethane"] == pytest.approx(2.4 / 5.0)
    # "This operation does not pass pressure or temperature" (and no energy
    # balance is made) — the outlet state is deliberately unset.
    assert out.T is None and out.P is None and out.H is None


# -- mode "mass" (book §6.3.2: Mass Flow; alkylation use case) ------------------
def test_mass_mode_alkylation_reactor_of_unknown_stoichiometry():
    """iC4H10 (58.12) + 1-C4H8 (56.11) -> iC8H18 (114.23): the Balance computes
    the product flow from the overall mass balance alone. Mass closes exactly;
    moles (2 -> ~1) and chemical species are NOT conserved — precisely the
    §6.3.2 semantics."""
    comps = ["isobutane", "1-butene", "isooctane"]
    fs = _fs(comps)
    fs.add(Balance("B", {"mode": "mass", "n_inlets": 2,
                         "z_out": {"isooctane": 1.0}}))
    fs.feed("IC4", "B:in1", T=300.0, P=5e5, molar_flow=1.0, z={"isobutane": 1.0})
    fs.feed("C4ENE", "B:in2", T=320.0, P=5e5, molar_flow=1.0, z={"1-butene": 1.0})
    fs.connect("IC8", "B:out1", None)
    rep = fs.solve()
    assert rep.converged
    out = fs.streams["IC8"]

    mass_in = molar_mass("isobutane") + molar_mass("1-butene")     # kg/s
    n_expected = mass_in / molar_mass("isooctane")
    assert out.molar_flow == pytest.approx(n_expected, rel=1e-12)
    # mass conserved...
    assert out.molar_flow * molar_mass("isooctane") == pytest.approx(mass_in, rel=1e-12)
    # ...moles and species are not (2 mol/s in -> ~1 mol/s out, new species).
    assert out.molar_flow == pytest.approx(1.0, abs=0.01)
    assert out.z == {"isooctane": 1.0}
    assert out.T is None and out.P is None      # T/P not passed (book)


# -- mode "heat" (book §6.3.3: Heat Flow) ---------------------------------------
def test_heat_mode_transfers_enthalpy_flow_to_an_energy_stream():
    """"Transfer the enthalpy of a process stream into a second energy
    stream": out1 is an energy port carrying sum(n_i * H_i)."""
    fs = _fs(["water"], package="coolprop:Water")
    fs.add(Balance("B", {"mode": "heat", "n_inlets": 1}))
    fs.feed("F", "B:in1", T=400.0, P=1e5, molar_flow=2.0, z={"water": 1.0})
    fs.connect("QOUT", "B:out1", None)
    rep = fs.solve()
    assert rep.converged
    feed = fs.streams["F"]
    assert rep.duties["QOUT"] == pytest.approx(feed.molar_flow * feed.H, rel=1e-12)
    assert "QOUT" not in fs.streams            # nothing material is passed


# -- mode "mole_heat" (book §6.3.4: Component Mole and Heat Flow) ---------------
def test_mole_heat_mode_balances_material_and_energy_independently():
    """With every inlet known, the §6.3.4 constraints (component moles AND
    total enthalpy conserved) pin the outlet to exactly what an adiabatic
    Mixer produces — the cross-check."""
    def feeds(fs, unit):
        fs.feed("F1", f"{unit}:in1", T=300.0, P=10e5, molar_flow=2.0,
                z={"methane": 1.0, "ethane": 0.0})
        fs.feed("F2", f"{unit}:in2", T=450.0, P=12e5, molar_flow=3.0,
                z={"methane": 0.2, "ethane": 0.8})

    fs_b = _fs(["methane", "ethane"])
    fs_b.add(Balance("B", {"mode": "mole_heat", "n_inlets": 2}))
    feeds(fs_b, "B")
    fs_b.connect("OUT", "B:out1", None)
    assert fs_b.solve().converged

    fs_m = _fs(["methane", "ethane"])
    fs_m.add(Mixer("M"))
    feeds(fs_m, "M")
    fs_m.connect("OUT", "M:out", None)
    assert fs_m.solve().converged

    bal, mix = fs_b.streams["OUT"], fs_m.streams["OUT"]
    assert bal.molar_flow == pytest.approx(mix.molar_flow, rel=1e-12)
    assert bal.T == pytest.approx(mix.T, rel=1e-9)
    assert bal.P == pytest.approx(mix.P)                   # lowest inlet P
    assert bal.H == pytest.approx(mix.H, rel=1e-9)
    for c in ("methane", "ethane"):
        assert bal.z[c] == pytest.approx(mix.z[c], rel=1e-12)


# -- mode "mass_heat" (book §6.3.5: Mass and Heat Flow) --------------------------
def test_mass_heat_mode_conserves_mass_and_energy_but_not_species():
    """n-butane relabeled as isobutane (equal molar mass): total mass AND total
    enthalpy flow are conserved while the species changes. On the engine's
    formation-inclusive basis the conserved enthalpy honestly carries the
    isomerization exotherm (Hf(iC4) < Hf(nC4)), so the outlet leaves hotter."""
    fs = _fs(["n-butane", "isobutane"])
    fs.add(Balance("B", {"mode": "mass_heat", "n_inlets": 1,
                         "z_out": {"isobutane": 1.0}}))
    fs.feed("F", "B:in1", T=350.0, P=1e5, molar_flow=1.0, z={"n-butane": 1.0})
    fs.connect("OUT", "B:out1", None)
    rep = fs.solve()
    assert rep.converged
    feed, out = fs.streams["F"], fs.streams["OUT"]

    # mass balance (isomers: identical molar mass -> identical molar flow)
    assert out.molar_flow == pytest.approx(1.0, rel=1e-12)
    # heat balance: total enthalpy flow conserved on the absolute basis
    assert out.molar_flow * out.H == pytest.approx(feed.molar_flow * feed.H, rel=1e-9)
    # species not conserved; the formation-enthalpy difference shows up as heat
    assert out.z == {"isobutane": 1.0}
    assert out.T > feed.T + 30.0


# -- economics: Balance is a logical op, no equipment ----------------------------
def test_balance_contributes_no_equipment():
    from caldyr.economics.sizing import size_flowsheet
    from caldyr.thermo import make_package

    fs = _fs(["methane"])
    fs.add(Balance("B", {"mode": "mole", "n_inlets": 1}))
    fs.feed("F", "B:in1", T=300.0, P=1e5, molar_flow=1.0, z={"methane": 1.0})
    fs.connect("OUT", "B:out1", None)
    rep = fs.solve()
    pp = make_package(fs.property_package, fs.component_ids)
    assert size_flowsheet(fs, rep, pp) == []
