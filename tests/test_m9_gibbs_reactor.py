"""M9 tests: the GibbsReactor (multi-reaction equilibrium by Gibbs minimization
via Cantera's ``equilibrate("TP")`` over gri30.yaml species).

Validation references:
  * Haber-Bosch N2 + 3H2 -> 2NH3 at 673.15 K / 200 bar, stoichiometric feed:
    literature equilibrium gives ~0.34-0.40 mole-fraction NH3 (classic
    Haber-Bosch equilibrium data; e.g. Smith, Van Ness & Abbott, *Intro. to
    Chemical Engineering Thermodynamics* — the same reference test_m2 uses for
    the single-reaction EquilibriumReactor). The Gibbs reactor must agree with
    the EquilibriumReactor within a few mol-% (both use ideal-gas standard
    states; they differ only in the underlying thermochemical data — caldyr's
    formation properties vs Cantera's gri30 NASA polynomials).
  * Steam-methane reforming CH4 + H2O -> CO + 3H2 with the water-gas shift
    CO + H2O -> CO2 + H2 at 1100 K / 1 bar, steam:carbon = 3: a single-reaction
    model cannot produce the CO/CO2 split; at these conditions reforming is
    nearly complete and all of CO/CO2/H2 coexist (Rostrup-Nielsen, *Catalytic
    Steam Reforming*; standard SMR equilibrium behavior). Reforming is
    endothermic, so H2 yield must rise with temperature (Le Chatelier).
"""
import math

import pytest

from caldyr.core import Component, Flowsheet
from caldyr.economics.costing import cost_equipment
from caldyr.economics.sizing import size_flowsheet
from caldyr.io import from_dict, to_dict
from caldyr.thermo import make_package
from caldyr.unitops import CanteraSpeciesError, EquilibriumReactor, GibbsReactor

AMMONIA = {"stoich": {"nitrogen": -1, "hydrogen": -3, "ammonia": 2}, "key": "nitrogen"}


def haber_flowsheet(reactor, *, components=None, z=None, P=2e7) -> Flowsheet:
    comps = components or ["nitrogen", "hydrogen", "ammonia"]
    fs = Flowsheet(components=[Component(c) for c in comps],
                   property_package="thermo:PR")
    fs.add(reactor)
    fs.feed("F", "R:in1", T=673.15, P=P, molar_flow=4.0,
            z=z or {"nitrogen": 0.25, "hydrogen": 0.75, "ammonia": 0.0})
    fs.connect("O", "R:out", None)
    fs.connect("Q", "R:duty", None)
    return fs


def smr_flowsheet(T: float) -> Flowsheet:
    comps = ["methane", "water", "carbon monoxide", "carbon dioxide", "hydrogen"]
    fs = Flowsheet(components=[Component(c) for c in comps],
                   property_package="thermo:PR")
    fs.add(GibbsReactor("R", {"T": T}))
    fs.feed("F", "R:in1", T=800.0, P=1e5, molar_flow=4.0,
            z={"methane": 0.25, "water": 0.75})
    fs.connect("O", "R:out", None)
    fs.connect("Q", "R:duty", None)
    return fs


def atom_flows(stream, atoms_per_component: dict[str, dict[str, float]]) -> dict[str, float]:
    out: dict[str, float] = {}
    for comp, atoms in atoms_per_component.items():
        n = stream.molar_flow * stream.z.get(comp, 0.0)
        for atom, count in atoms.items():
            out[atom] = out.get(atom, 0.0) + n * count
    return out


# -- Haber-Bosch validation --------------------------------------------------
def test_haber_matches_equilibrium_reactor_and_literature():
    fs_g = haber_flowsheet(GibbsReactor("R", {"T": 673.15}))
    fs_e = haber_flowsheet(EquilibriumReactor("R", {"reaction": AMMONIA, "T": 673.15}))
    assert fs_g.solve().converged and fs_e.solve().converged
    y_g = fs_g.streams["O"].z["ammonia"]
    y_e = fs_e.streams["O"].z["ammonia"]
    # within a few mol-% of the single-reaction model, and in the literature band
    assert y_g == pytest.approx(y_e, abs=0.03)
    assert y_g == pytest.approx(0.37, abs=0.06)      # ~0.34-0.40 at 400 C, 200 bar


def test_haber_atom_balance_exact():
    fs = haber_flowsheet(GibbsReactor("R", {"T": 673.15}))
    fs.solve()
    atoms = {"nitrogen": {"N": 2}, "hydrogen": {"H": 2}, "ammonia": {"N": 1, "H": 3}}
    out = atom_flows(fs.streams["O"], atoms)
    assert math.isclose(out["N"], 2.0, rel_tol=1e-9)   # feed: 1 mol/s N2
    assert math.isclose(out["H"], 6.0, rel_tol=1e-9)   # feed: 3 mol/s H2


def test_haber_pressure_trend():
    def y_nh3(P):
        fs = haber_flowsheet(GibbsReactor("R", {"T": 673.15}), P=P)
        fs.solve()
        return fs.streams["O"].z["ammonia"]
    # Le Chatelier: mole-number-reducing reaction favored at higher pressure.
    assert y_nh3(3e7) > y_nh3(2e7) > y_nh3(1e7)


def test_duty_closes_on_caldyr_enthalpy_basis():
    """Composition comes from Cantera but the energy balance is caldyr's: the
    inlet enthalpy plus the reported duty must equal the outlet enthalpy
    computed by the flowsheet's own property package."""
    fs = haber_flowsheet(GibbsReactor("R", {"T": 673.15}))
    rep = fs.solve()
    f, o = fs.streams["F"], fs.streams["O"]
    assert o.molar_flow * o.H == pytest.approx(
        f.molar_flow * f.H + rep.duties["Q"], rel=1e-9)
    assert rep.duties["Q"] < 0          # exothermic held isothermal: remove heat


def test_inert_argon_passes_through():
    """Cantera carries non-reacting species through the minimization — no
    explicit `inerts` parameter needed. Argon flow is conserved exactly."""
    fs = haber_flowsheet(
        GibbsReactor("R", {"T": 673.15}),
        components=["nitrogen", "hydrogen", "ammonia", "argon"],
        z={"nitrogen": 0.225, "hydrogen": 0.675, "ammonia": 0.0, "argon": 0.10})
    fs.solve()
    o = fs.streams["O"]
    assert o.molar_flow * o.z["argon"] == pytest.approx(0.4, rel=1e-9)
    assert o.z["ammonia"] > 0.2          # reaction still proceeds


# -- steam-methane reforming + water-gas shift (multi-reaction) ---------------
def test_smr_produces_full_syngas_slate():
    """A single-reaction equilibrium cannot give the CO/CO2 split; the Gibbs
    minimization must produce CH4, CO, CO2, H2 and H2O simultaneously."""
    fs = smr_flowsheet(1100.0)
    assert fs.solve().converged
    o = fs.streams["O"]
    assert o.z["hydrogen"] > 0.4              # reforming nearly complete at 1100 K
    assert o.z["carbon monoxide"] > 0.01
    assert o.z["carbon dioxide"] > 0.01
    assert o.z["methane"] > 0.0               # trace unconverted CH4
    assert o.z["water"] > 0.1                 # excess steam remains


def test_smr_atom_balances_close():
    fs = smr_flowsheet(1100.0)
    fs.solve()
    atoms = {
        "methane": {"C": 1, "H": 4},
        "water": {"H": 2, "O": 1},
        "carbon monoxide": {"C": 1, "O": 1},
        "carbon dioxide": {"C": 1, "O": 2},
        "hydrogen": {"H": 2},
    }
    out = atom_flows(fs.streams["O"], atoms)
    # feed: 1 mol/s CH4 + 3 mol/s H2O
    assert math.isclose(out["C"], 1.0, rel_tol=1e-8)
    assert math.isclose(out["H"], 10.0, rel_tol=1e-8)
    assert math.isclose(out["O"], 3.0, rel_tol=1e-8)


def test_smr_h2_yield_rises_with_temperature():
    """Reforming is endothermic — hotter shifts equilibrium toward H2."""
    def h2_flow(T):
        fs = smr_flowsheet(T)
        fs.solve()
        o = fs.streams["O"]
        return o.molar_flow * o.z["hydrogen"]
    assert h2_flow(1100.0) > h2_flow(900.0)


# -- typed, documented error paths --------------------------------------------
def test_missing_t_raises_a_clear_error():
    fs = haber_flowsheet(GibbsReactor("R", {}))
    with pytest.raises(ValueError, match=r"params\['T'\] is required"):
        fs.solve()


def test_unmapped_flowing_component_raises_typed_error():
    """n-hexane is not in the gri30 species set; with nonzero flow this must be
    a typed error naming the component, not a KeyError or a silent drop."""
    fs = Flowsheet(components=[Component("methane"), Component("n-hexane")],
                   property_package="thermo:PR")
    fs.add(GibbsReactor("R", {"T": 1000.0}))
    fs.feed("F", "R:in1", T=800.0, P=1e5, molar_flow=1.0,
            z={"methane": 0.5, "n-hexane": 0.5})
    fs.connect("O", "R:out", None)
    fs.connect("Q", "R:duty", None)
    with pytest.raises(CanteraSpeciesError, match="n-hexane"):
        fs.solve()


def test_unmapped_zero_flow_component_passes_through():
    """A flowsheet may carry components the Gibbs reactor cannot map as long as
    they do not flow through it (they cannot participate anyway)."""
    fs = Flowsheet(components=[Component("methane"), Component("water"),
                               Component("hydrogen"), Component("carbon monoxide"),
                               Component("n-hexane")],
                   property_package="thermo:PR")
    fs.add(GibbsReactor("R", {"T": 1100.0}))
    fs.feed("F", "R:in1", T=800.0, P=1e5, molar_flow=4.0,
            z={"methane": 0.25, "water": 0.75, "hydrogen": 0.0,
               "carbon monoxide": 0.0, "n-hexane": 0.0})
    fs.connect("O", "R:out", None)
    fs.connect("Q", "R:duty", None)
    assert fs.solve().converged
    o = fs.streams["O"]
    assert o.z["n-hexane"] == 0.0
    assert o.z["hydrogen"] > 0.4


# -- `.flow` round-trip and economics ------------------------------------------
def test_flow_round_trip_is_exact():
    fs = smr_flowsheet(1100.0)
    fs.solve()
    doc = to_dict(fs)
    assert doc == to_dict(from_dict(doc))


def test_economics_sizes_a_vertical_vessel():
    fs = smr_flowsheet(1100.0)
    rep = fs.solve()
    pp = make_package(fs.property_package, fs.component_ids)
    sizes = size_flowsheet(fs, rep, pp)
    assert len(sizes) == 1
    assert sizes[0].equipment_type == "vessel_vertical"
    assert cost_equipment(sizes[0]).bare_module > 0
