"""M19 tests: the Decanter unit op and the NRTL liquid-liquid (isoactivity)
VLLE flash that gives the activity package three-phase capability.

3-phase (heterogeneous-azeotrope) distillation needs two ingredients: a
property package that can flash a liquid into two liquid phases, and a unit
that settles + (optionally) refluxes one of them. The decant itself is the
property package's ``flash_pt_3p``:

* the cubic-EOS packages (``thermo:PR``/``thermo:SRK``) do it via thermo's
  FlashVLN (two trial EOS liquids) — accurate for the *structure* of a
  water/organic split though a cubic EOS overpredicts the cross-solubility;
* the activity package (``thermo:NRTL``) does it here with a new isoactivity
  liquid-liquid flash (``gamma_i^I x_i^I = gamma_i^II x_i^II``, Rachford-Rice on
  the distribution ratios), which the cubic EOS cannot represent for strongly
  non-ideal polar liquids.

The NRTL-LLE assertions use an illustrative parameter set with a real
miscibility gap (the ChemSep NRTL parameters bundled with thermo are fit to
VLE only and carry essentially no gap for these pairs); the test checks the
*flash is correct* (isoactivity, balances) rather than the exact mutual
solubilities, which a single NRTL set cannot fit jointly with VLE.
"""
import pytest

from caldyr.core import Component, Flowsheet
from caldyr.economics.costing import cost_equipment
from caldyr.economics.sizing import size_flowsheet
from caldyr.io import from_dict, to_dict
from caldyr.thermo import make_package
from caldyr.unitops import Decanter

P_ATM = 101325.0


def _nrtl_gap_pkg(comps):
    """An NRTL package whose binary parameters carry a genuine water/butanol
    miscibility gap (illustrative LLE-focused values; tau_ij = b_ij/T, alpha
    0.4 — they reproduce a physical split: organic ~0.46 water / aqueous ~0.99
    water, near the experimental 0.51/0.98)."""
    from thermo import NRTL
    pp = make_package("thermo:NRTL", comps)
    pp._flasher.liquid.GibbsExcessModel = NRTL(
        T=298.15, xs=[0.5, 0.5], tau_bs=[[0.0, 1300.0], [1100.0, 0.0]],
        alpha_cs=[[0.0, 0.4], [0.4, 0.0]])
    return pp


# -- the NRTL liquid-liquid isoactivity flash ---------------------------------
def test_nrtl_vlle_flash_splits_with_equal_activities():
    pp = _nrtl_gap_pkg(["water", "n-butanol"])
    r = pp.flash_pt_3p(320.0, P_ATM, {"water": 0.5, "n-butanol": 0.5})
    assert r.beta_vapor < 1e-6                       # fully condensed at 320 K
    assert r.beta_light > 0.0 and r.beta_heavy > 0.0  # two liquid phases
    # organic (light) is butanol-rich, aqueous (heavy) is water-rich
    assert r.x_light["n-butanol"] > 0.4 and r.x_heavy["water"] > 0.95
    assert (r.rho_light or 0.0) < (r.rho_heavy or 0.0)
    # isoactivity gamma_i^I x_i^I == gamma_i^II x_i^II to tight tolerance
    gI = pp._gammas(320.0, [r.x_light["water"], r.x_light["n-butanol"]])
    gII = pp._gammas(320.0, [r.x_heavy["water"], r.x_heavy["n-butanol"]])
    for k, i in (("water", 0), ("n-butanol", 1)):
        aI = r.x_light[k] * gI[i]
        aII = r.x_heavy[k] * gII[i]
        assert abs(aI - aII) < 1e-6
    # bulk mole balance: the two phases recombine to the feed
    for c in ("water", "n-butanol"):
        rec = r.beta_light * r.x_light[c] + r.beta_heavy * r.x_heavy[c]
        assert abs(rec - 0.5) < 1e-9
    assert abs(r.beta_light + r.beta_heavy - 1.0) < 1e-9


def test_nrtl_vlle_flash_degrades_to_single_liquid_without_a_gap():
    # The default ChemSep NRTL parameters are VLE-fit (no miscibility gap);
    # the three-phase flash must return ONE liquid, not invent a split.
    pp = make_package("thermo:NRTL", ["ethanol", "water"])
    r = pp.flash_pt_3p(330.0, P_ATM, {"ethanol": 0.5, "water": 0.5})
    assert r.beta_heavy == 0.0
    assert r.beta_light > 0.0 and r.x_light is not None
    # single liquid -> its composition is the feed
    assert abs(r.x_light["ethanol"] - 0.5) < 1e-6


def test_nrtl_vlle_flash_ph_finds_the_same_split():
    pp = _nrtl_gap_pkg(["water", "n-butanol"])
    z = {"water": 0.5, "n-butanol": 0.5}
    pt = pp.flash_pt_3p(320.0, P_ATM, z)
    H = pp.enthalpy(320.0, P_ATM, z)
    ph = pp.flash_ph_3p(P_ATM, H, z)
    assert abs(ph.T - 320.0) < 1.0
    assert ph.beta_heavy > 0.0
    assert abs(ph.x_heavy["water"] - pt.x_heavy["water"]) < 1e-3


# -- the Decanter unit op (on PR, which robustly splits water/organic) ---------
def decanter_fs(reflux_fraction=None, reflux_layer="light", T=320.0):
    """A water/n-butanol heteroazeotrope condensate into a decanter at T."""
    params: dict = {"T": T, "P": P_ATM, "reflux_layer": reflux_layer}
    if reflux_fraction is not None:
        params["reflux_fraction"] = reflux_fraction
    fs = Flowsheet(components=[Component("water"), Component("n-butanol")],
                   property_package="thermo:PR")
    fs.add(Decanter("DEC", params))
    fs.feed("OVHD", "DEC:in1", T=360.0, P=P_ATM, molar_flow=100.0,
            z={"water": 0.6, "n-butanol": 0.4})
    for port in ("liquid_light", "liquid_heavy", "reflux", "product", "vapor"):
        fs.connect(port.upper(), f"DEC:{port}", None)
    fs.connect("Q", "DEC:duty", None)
    return fs


def test_decanter_splits_into_organic_and_aqueous_layers():
    fs = decanter_fs()
    assert fs.solve().converged
    light, heavy = fs.streams["LIQUID_LIGHT"], fs.streams["LIQUID_HEAVY"]
    assert light.molar_flow > 0 and heavy.molar_flow > 0
    # light = organic (butanol-rich, less dense), heavy = aqueous (water-rich)
    assert light.z["n-butanol"] > 0.7
    assert heavy.z["water"] > 0.95
    dec = fs.units["DEC"]
    assert dec.result["rho_light"] < dec.result["rho_heavy"]
    assert dec.result["beta_vapor"] < 1e-6


def test_decanter_mass_balance_is_exact():
    fs = decanter_fs()
    fs.solve()
    feed = fs.streams["OVHD"]
    for c in ("water", "n-butanol"):
        out = sum(fs.streams[s].molar_flow * fs.streams[s].z[c]
                  for s in ("LIQUID_LIGHT", "LIQUID_HEAVY", "REFLUX",
                            "PRODUCT", "VAPOR"))
        assert abs(feed.molar_flow * feed.z[c] - out) < 1e-9


def test_decanter_reflux_split_routes_one_layer():
    fs = decanter_fs(reflux_fraction=0.7, reflux_layer="light")
    fs.solve()
    refl, prod = fs.streams["REFLUX"], fs.streams["PRODUCT"]
    light = fs.streams["LIQUID_LIGHT"]
    # the light layer is split into reflux (70%) + product (30%); its own port empties
    assert light.molar_flow < 1e-9
    assert refl.molar_flow > 0 and prod.molar_flow > 0
    assert abs(refl.molar_flow / (refl.molar_flow + prod.molar_flow) - 0.7) < 1e-9
    # reflux and product carry the SAME (organic) composition
    assert abs(refl.z["n-butanol"] - prod.z["n-butanol"]) < 1e-12
    assert refl.z["n-butanol"] > 0.7
    # the aqueous layer still leaves on its own port
    assert fs.streams["LIQUID_HEAVY"].molar_flow > 0


def test_decanter_reflux_heavy_layer():
    fs = decanter_fs(reflux_fraction=1.0, reflux_layer="heavy")
    fs.solve()
    # all of the heavy (aqueous) layer is refluxed; product empty
    assert fs.streams["REFLUX"].z["water"] > 0.95
    assert fs.streams["PRODUCT"].molar_flow < 1e-9
    assert fs.streams["LIQUID_HEAVY"].molar_flow < 1e-9
    assert fs.streams["LIQUID_LIGHT"].molar_flow > 0   # organic on its own port


def test_decanter_single_liquid_degrades_gracefully():
    # Fully miscible benzene/toluene: one liquid, no second phase, not an error.
    fs = Flowsheet(components=[Component("benzene"), Component("toluene")],
                   property_package="thermo:PR")
    fs.add(Decanter("DEC", {"T": 320.0, "P": P_ATM}))
    fs.feed("F", "DEC:in1", T=350.0, P=P_ATM, molar_flow=10.0,
            z={"benzene": 0.5, "toluene": 0.5})
    for port in ("liquid_light", "liquid_heavy", "reflux", "product", "vapor"):
        fs.connect(port.upper(), f"DEC:{port}", None)
    fs.connect("Q", "DEC:duty", None)
    fs.solve()
    assert fs.streams["LIQUID_LIGHT"].molar_flow > 0
    assert fs.streams["LIQUID_HEAVY"].molar_flow < 1e-9


def test_decanter_requires_temperature():
    fs = Flowsheet(components=[Component("water"), Component("n-butanol")],
                   property_package="thermo:PR")
    fs.add(Decanter("DEC", {"P": P_ATM}))     # no T
    fs.feed("F", "DEC:in1", T=360.0, P=P_ATM, molar_flow=10.0,
            z={"water": 0.6, "n-butanol": 0.4})
    for port in ("liquid_light", "liquid_heavy", "reflux", "product", "vapor"):
        fs.connect(port.upper(), f"DEC:{port}", None)
    fs.connect("Q", "DEC:duty", None)
    with pytest.raises(ValueError, match="T.*required"):
        fs.solve()


def test_decanter_roundtrips_and_costs():
    fs = decanter_fs(reflux_fraction=0.6)
    fs2 = from_dict(to_dict(fs))
    rep = fs2.solve()
    assert rep.converged
    pp = make_package(fs2.property_package, fs2.component_ids)
    sizes = size_flowsheet(fs2, rep, pp)
    # the decanter is sized + costed as a horizontal vessel
    dec_sizes = [s for s in sizes if s.unit_id == "DEC"]
    assert dec_sizes and dec_sizes[0].equipment_type == "vessel_horizontal"
    assert cost_equipment(dec_sizes[0]).bare_module > 0
