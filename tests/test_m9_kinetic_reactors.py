"""M9 tests: kinetic reactors (CSTR + PFR) with power-law kinetics.

Validation references:
  * First-order A -> B at constant density and temperature has the classic
    closed forms (H. S. Fogler, *Elements of Chemical Reaction Engineering*,
    4th ed., Ch. 5 — design equations for a first-order reaction; identical in
    O. Levenspiel, *Chemical Reaction Engineering*, 3rd ed., Ch. 5,
    Eqs. 5.11 / 5.17 with eps_A = 0):

        CSTR:  X = k·tau / (1 + k·tau)
        PFR:   X = 1 - exp(-k·tau),         tau = V / Q0

    The test reaction is the gas-phase cis-2-butene -> trans-2-butene
    isomerization (no mole change) at 1 bar / 700 K, where the mixture molar
    volume is composition-independent to ~1e-5 — so the constant-density
    closed form must be reproduced to ~1e-4 (the engine's concentrations come
    from the property package's real molar volume, documented in
    caldyr.unitops.reaction.concentrations).
  * For equal volume and the same kinetics, PFR conversion exceeds CSTR
    conversion for a positive-order reaction (Levenspiel 3e, Ch. 6, Fig. 6.1).
  * Adiabatic operation: exothermic ethylene hydrogenation
    C2H4 + H2 -> C2H6 (dH ~ -137 kJ/mol, CRC) must raise the outlet
    temperature with zero duty and exact total-enthalpy closure (the
    formation-inclusive enthalpy basis carries the heat of reaction).
"""
import math

import pytest

from caldyr.core import Component, Flowsheet
from caldyr.economics.costing import cost_equipment
from caldyr.economics.sizing import size_flowsheet
from caldyr.io import from_dict, to_dict
from caldyr.thermo import make_package
from caldyr.unitops import CSTR, PFR

R_GAS = 8.314462618
T_ISO, P_ISO, V_ISO = 700.0, 1e5, 1.0
K0, EA = 1e7, 110e3                       # k(700 K) = 0.0619 1/s
ISOM = {"stoich": {"cis-2-butene": -1, "trans-2-butene": 1},
        "key": "cis-2-butene", "k0": K0, "Ea": EA}
HYDROG = {"stoich": {"ethylene": -1, "hydrogen": -1, "ethane": 1},
          "key": "ethylene", "k0": 5e4, "Ea": 70e3,
          "orders": {"ethylene": 1, "hydrogen": 1}}


def isom_flowsheet(reactor) -> Flowsheet:
    fs = Flowsheet(components=[Component("cis-2-butene"), Component("trans-2-butene")],
                   property_package="thermo:PR")
    fs.add(reactor)
    fs.feed("F", "R:in1", T=T_ISO, P=P_ISO, molar_flow=1.0,
            z={"cis-2-butene": 1.0, "trans-2-butene": 0.0})
    fs.connect("O", "R:out", None)
    fs.connect("Q", "R:duty", None)
    return fs


def hydrog_flowsheet(reactor) -> Flowsheet:
    fs = Flowsheet(components=[Component("ethylene"), Component("hydrogen"),
                               Component("ethane")],
                   property_package="thermo:PR")
    fs.add(reactor)
    fs.feed("F", "R:in1", T=600.0, P=1e6, molar_flow=2.0,
            z={"ethylene": 0.5, "hydrogen": 0.5, "ethane": 0.0})
    fs.connect("O", "R:out", None)
    fs.connect("Q", "R:duty", None)
    return fs


def conversion(fs, component: str, n_feed: float) -> float:
    o = fs.streams["O"]
    return 1.0 - o.molar_flow * o.z[component] / n_feed


def closed_form_ktau() -> float:
    """k·tau with tau = V/Q0 from the same property package (a property lookup,
    not a reimplementation of the reactor)."""
    pp = make_package("thermo:PR", ["cis-2-butene", "trans-2-butene"])
    k = K0 * math.exp(-EA / (R_GAS * T_ISO))
    q0 = 1.0 * pp.volume(T_ISO, P_ISO, {"cis-2-butene": 1.0})
    return k * V_ISO / q0


# -- closed-form validation (Fogler 4e Ch. 5 / Levenspiel 3e Ch. 5) -----------
def test_cstr_matches_first_order_closed_form():
    ktau = closed_form_ktau()
    fs = isom_flowsheet(CSTR("R", {"V": V_ISO, "T": T_ISO, "reactions": [ISOM]}))
    assert fs.solve().converged
    x_ref = ktau / (1.0 + ktau)
    assert conversion(fs, "cis-2-butene", 1.0) == pytest.approx(x_ref, abs=1e-4)


def test_pfr_matches_first_order_closed_form():
    ktau = closed_form_ktau()
    fs = isom_flowsheet(PFR("R", {"V": V_ISO, "T": T_ISO, "reactions": [ISOM]}))
    assert fs.solve().converged
    x_ref = 1.0 - math.exp(-ktau)
    assert conversion(fs, "cis-2-butene", 1.0) == pytest.approx(x_ref, abs=1e-4)


def test_pfr_beats_cstr_at_equal_volume():
    """Levenspiel 3e Ch. 6: for n > 0 kinetics a PFR always out-converts a CSTR
    of the same volume."""
    fs_c = isom_flowsheet(CSTR("R", {"V": V_ISO, "T": T_ISO, "reactions": [ISOM]}))
    fs_p = isom_flowsheet(PFR("R", {"V": V_ISO, "T": T_ISO, "reactions": [ISOM]}))
    fs_c.solve()
    fs_p.solve()
    assert conversion(fs_p, "cis-2-butene", 1.0) > conversion(fs_c, "cis-2-butene", 1.0)


def test_default_order_is_first_in_key():
    """Omitting 'orders' must mean first order in the key reactant."""
    explicit = dict(ISOM, orders={"cis-2-butene": 1.0})
    implicit = {k: v for k, v in ISOM.items() if k != "orders"}
    fs1 = isom_flowsheet(CSTR("R", {"V": V_ISO, "T": T_ISO, "reactions": [explicit]}))
    fs2 = isom_flowsheet(CSTR("R", {"V": V_ISO, "T": T_ISO, "reactions": [implicit]}))
    fs1.solve()
    fs2.solve()
    assert fs1.streams["O"].z["cis-2-butene"] == pytest.approx(
        fs2.streams["O"].z["cis-2-butene"], rel=1e-12)


def test_cstr_conversion_increases_with_volume():
    def x_at(v):
        fs = isom_flowsheet(CSTR("R", {"V": v, "T": T_ISO, "reactions": [ISOM]}))
        fs.solve()
        return conversion(fs, "cis-2-butene", 1.0)
    assert x_at(2.0) > x_at(1.0) > x_at(0.5)


# -- isothermal duty -----------------------------------------------------------
def test_isothermal_energy_balance_closes_with_duty():
    fs = hydrog_flowsheet(CSTR("R", {"V": 0.005, "T": 600.0, "reactions": [HYDROG]}))
    rep = fs.solve()
    f, o = fs.streams["F"], fs.streams["O"]
    assert rep.duties["Q"] < 0                        # exothermic, held at 600 K
    assert o.molar_flow * o.H == pytest.approx(
        f.molar_flow * f.H + rep.duties["Q"], rel=1e-9)


# -- adiabatic operation ---------------------------------------------------------
@pytest.mark.parametrize("cls", [CSTR, PFR])
def test_adiabatic_exothermic_reactor_heats_up(cls):
    fs = hydrog_flowsheet(cls("R", {"V": 0.05, "reactions": [HYDROG]}))
    rep = fs.solve()
    f, o = fs.streams["F"], fs.streams["O"]
    assert rep.duties["Q"] == pytest.approx(0.0, abs=1e-9)
    assert o.T > f.T + 50.0                           # exothermic temperature rise
    # energy balance: total enthalpy flow conserved
    assert o.molar_flow * o.H == pytest.approx(f.molar_flow * f.H, rel=1e-9)


@pytest.mark.parametrize("cls", [CSTR, PFR])
def test_atom_balance_with_mole_change(cls):
    """C2H4 + H2 -> C2H6 destroys moles; C and H atoms must still balance."""
    fs = hydrog_flowsheet(cls("R", {"V": 0.05, "reactions": [HYDROG]}))
    fs.solve()
    o = fs.streams["O"]
    n = o.molar_flow
    c_atoms = 2 * n * o.z["ethylene"] + 2 * n * o.z["ethane"]
    h_atoms = 4 * n * o.z["ethylene"] + 2 * n * o.z["hydrogen"] + 6 * n * o.z["ethane"]
    assert c_atoms == pytest.approx(2.0, rel=1e-7)    # feed: 1 mol/s C2H4
    assert h_atoms == pytest.approx(6.0, rel=1e-7)    # feed: 4 + 2 mol/s H
    assert n < 2.0                                    # net mole destruction


# -- typed, documented error paths -----------------------------------------------
def test_missing_reactions_raises():
    fs = isom_flowsheet(CSTR("R", {"V": 1.0, "T": T_ISO}))
    with pytest.raises(ValueError, match=r"params\['reactions'\]"):
        fs.solve()


def test_nonpositive_volume_raises():
    fs = isom_flowsheet(PFR("R", {"V": 0.0, "T": T_ISO, "reactions": [ISOM]}))
    with pytest.raises(ValueError, match=r"params\['V'\]"):
        fs.solve()


def test_key_must_be_a_reactant():
    bad = dict(ISOM, key="trans-2-butene")
    fs = isom_flowsheet(CSTR("R", {"V": 1.0, "T": T_ISO, "reactions": [bad]}))
    with pytest.raises(ValueError, match="reactant"):
        fs.solve()


# -- `.flow` round-trip and economics ----------------------------------------------
def test_flow_round_trip_is_exact():
    fs = isom_flowsheet(CSTR("R", {"V": V_ISO, "T": T_ISO, "reactions": [ISOM]}))
    fs.solve()
    doc = to_dict(fs)
    assert doc == to_dict(from_dict(doc))


@pytest.mark.parametrize("cls", [CSTR, PFR])
def test_economics_sizes_on_the_design_volume(cls):
    """Kinetic reactors carry their design volume, so sizing must use it
    directly (vertical vessel), not a residence-time heuristic."""
    fs = isom_flowsheet(cls("R", {"V": V_ISO, "T": T_ISO, "reactions": [ISOM]}))
    rep = fs.solve()
    pp = make_package(fs.property_package, fs.component_ids)
    sizes = size_flowsheet(fs, rep, pp)
    assert len(sizes) == 1
    assert sizes[0].equipment_type == "vessel_vertical"
    assert sizes[0].attribute == pytest.approx(V_ISO)
    assert cost_equipment(sizes[0]).bare_module > 0
