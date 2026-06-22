"""M17 tests: Claus sulfur recovery (book §10.2.3).

The Claus process recovers elemental sulfur from the H2S-rich acid gas an amine
unit produces. This module validates the enabling NASA ideal-gas property
package (which carries the sulfur allotropes the cubic EOS cannot), the Claus
reaction furnace (adiabatic) and catalytic converter (isothermal), the sulfur
condenser, and the assembled train.

Validation references (the book reports no stream numbers for §10.2.3, so the
checks are first-principles / against well-established Claus literature):

* **NASA enthalpy basis** — Cantera's NASA polynomials are formation-inclusive;
  H2S and H2O must come out at their formation enthalpies (−20.6, −241.8 kJ/mol;
  e.g. NIST-JANAF), matching caldyr's ``thermo:*`` basis.
* **Thermal stage** — substoichiometric combustion of H2S gives an adiabatic
  flame temperature of ~1000-1400 °C and high-temperature sulfur as the dimer S2
  (Kohl & Nielsen, *Gas Purification* 5e, ch. 8; Gamson & Elkins 1953).
* **Air-demand optimum** — Claus recovery is maximised when the air oxidises
  exactly one-third of the H2S, i.e. the catalytic-section gas sits at the
  stoichiometric H2S:SO2 = 2:1 (the basis of every Claus air-demand controller).
* **Overall recovery** — a thermal stage plus two catalytic converters recovers
  ~96-98 % of feed sulfur at equilibrium (an equilibrium model is a slight
  over-prediction of real, kinetically/sub-dewpoint-limited plants).
* **Atom and energy balances** — exact (sulfur, hydrogen, oxygen, nitrogen) and
  machine-precision energy closure on the whole train.
"""
import pytest

from caldyr.core import Component, Flowsheet
from caldyr.io import from_dict, to_dict
from caldyr.thermo import NasaGasPackage, make_package
from caldyr.thermo.nasa_pkg import NasaSpeciesError
from caldyr.thermo.sulfur import liquid_sulfur_psat, sulfur_hvap_per_atom
from caldyr.unitops import (
    ClausReactor,
    ClausReactorError,
    Heater,
    SulfurCondenser,
    SulfurCondenserError,
)

COMPS = ["hydrogen sulfide", "sulfur dioxide", "S2", "S8", "water", "nitrogen",
         "oxygen", "carbon dioxide", "carbon monoxide", "hydrogen"]

# atoms per component, for balance checks
_ATOMS = {
    "hydrogen sulfide": {"H": 2, "S": 1},
    "sulfur dioxide": {"S": 1, "O": 2},
    "S2": {"S": 2},
    "S8": {"S": 8},
    "water": {"H": 2, "O": 1},
    "nitrogen": {"N": 2},
    "oxygen": {"O": 2},
    "carbon dioxide": {"C": 1, "O": 2},
    "carbon monoxide": {"C": 1, "O": 1},
    "hydrogen": {"H": 2},
}


def _atom_flows(*streams) -> dict[str, float]:
    out: dict[str, float] = {}
    for s in streams:
        for comp, atoms in _ATOMS.items():
            n = s.molar_flow * s.z.get(comp, 0.0)
            for a, k in atoms.items():
                out[a] = out.get(a, 0.0) + n * k
    return out


# ============================ NASA property package ==========================
def test_nasa_enthalpy_is_formation_inclusive():
    pp = make_package("nasa:gas", ["hydrogen sulfide", "water", "nitrogen"])
    # formation enthalpies at 298.15 K (J/mol)
    assert pp.enthalpy(298.15, 1e5, {"hydrogen sulfide": 1.0}) == pytest.approx(
        -20600.0, abs=1500.0)
    assert pp.enthalpy(298.15, 1e5, {"water": 1.0}) == pytest.approx(
        -241800.0, abs=1500.0)


def test_nasa_ideal_gas_volume_and_flash_roundtrip():
    pp = make_package("nasa:claus", COMPS)
    z = {"hydrogen sulfide": 0.3, "water": 0.2, "nitrogen": 0.5}
    R = 8.314462618
    assert pp.volume(600.0, 1.4e5, z) == pytest.approx(R * 600.0 / 1.4e5, rel=1e-12)
    res = pp.flash_pt(600.0, 1.4e5, z)
    assert res.phase == "vapor" and res.vapor_fraction == 1.0
    # PH and PS flashes recover the temperature
    H = pp.enthalpy(600.0, 1.4e5, z)
    assert pp.flash_ph(1.4e5, H, z).T == pytest.approx(600.0, abs=1e-3)
    S = pp.entropy(600.0, 1.4e5, z)
    assert pp.flash_ps(1.4e5, S, z).T == pytest.approx(600.0, abs=1e-3)


def test_nasa_rejects_unmappable_and_pseudo_and_vle():
    with pytest.raises(NasaSpeciesError, match="n-hexane"):
        make_package("nasa:gas", ["n-hexane"])
    pp = make_package("nasa:gas", ["hydrogen sulfide", "S8"])
    # no VLE / liquid in an ideal-gas package
    with pytest.raises(NotImplementedError):
        pp.bubble_point(1e5, {"S8": 1.0})
    with pytest.raises(NotImplementedError):
        pp.k_values(450.0, 1e5, {"S8": 1.0}, {"S8": 1.0})
    assert isinstance(pp, NasaGasPackage)


def test_liquid_sulfur_correlations_anchor_to_boiling_point():
    # thermo's sulfur Antoine reproduces the 717.8 K normal boiling point
    assert liquid_sulfur_psat(717.8) == pytest.approx(1.013e5, rel=0.05)
    # condenser-range vapour pressure is small (tens of Pa) and rises with T
    assert 5.0 < liquid_sulfur_psat(405.0) < 60.0
    assert liquid_sulfur_psat(450.0) > liquid_sulfur_psat(405.0)
    # latent heat ~10 kJ per mol of S atoms near the condensers
    assert sulfur_hvap_per_atom(440.0) == pytest.approx(10400.0, rel=0.2)


# ============================ Claus reaction furnace =========================
def _furnace(o2: float = 45.0, h2s: float = 0.90):
    fs = Flowsheet(components=[Component(c) for c in COMPS],
                   property_package="nasa:claus")
    fs.add(ClausReactor("FURN", {}))                    # adiabatic
    fs.feed("ACID", "FURN:in1", T=320.0, P=1.6e5, molar_flow=100.0,
            z={"hydrogen sulfide": h2s, "carbon dioxide": 1.0 - h2s})
    n2 = o2 * 79.0 / 21.0
    fs.feed("AIR", "FURN:in2", T=480.0, P=1.6e5, molar_flow=o2 + n2,
            z={"oxygen": o2 / (o2 + n2), "nitrogen": n2 / (o2 + n2)})
    fs.connect("g", "FURN:out", None)
    fs.connect("q", "FURN:duty", None)
    return fs


def test_thermal_furnace_flame_temperature_and_s2():
    fs = _furnace()
    rep = fs.solve()
    assert rep.converged
    g = fs.streams["g"]
    assert 1273.0 < g.T < 1673.0          # ~1000-1400 C adiabatic flame
    assert g.z["S2"] > 0.05               # high-T sulfur is the dimer S2
    assert g.z.get("S8", 0.0) < 1e-3      # S8 negligible at flame temperature
    assert rep.duties["q"] == pytest.approx(0.0, abs=1e-3)   # adiabatic


def test_thermal_furnace_atom_balance_exact():
    fs = _furnace()
    fs.solve()
    a_in = _atom_flows(fs.streams["ACID"], fs.streams["AIR"])
    a_out = _atom_flows(fs.streams["g"])
    for atom in ("S", "H", "O", "N"):
        assert a_out[atom] == pytest.approx(a_in[atom], rel=1e-9), atom


# ============================ Catalytic converter ============================
def test_catalytic_converter_makes_s8_and_duty_closes():
    """A 2:1 H2S:SO2 feed over a 530 K isothermal converter forms S8 (the
    low-temperature allotrope), is exothermic (duty removed), and closes its
    own energy balance against the property package."""
    fs = Flowsheet(components=[Component(c) for c in COMPS],
                   property_package="nasa:claus")
    fs.add(ClausReactor("CV", {"T": 530.0}))
    fs.feed("F", "CV:in1", T=530.0, P=1.4e5, molar_flow=100.0,
            z={"hydrogen sulfide": 0.04, "sulfur dioxide": 0.02, "water": 0.30,
               "nitrogen": 0.64})
    fs.connect("o", "CV:out", None)
    fs.connect("q", "CV:duty", None)
    rep = fs.solve()
    assert rep.converged
    o = fs.streams["o"]
    assert o.z["S8"] > 1e-3                       # sulfur formed as S8
    assert o.z["hydrogen sulfide"] < 0.04         # H2S consumed
    assert rep.duties["q"] < 0                    # Claus reaction is exothermic
    # duty closes on the package basis: n_in*h_in + Q == n_out*h_out
    f, oo = fs.streams["F"], fs.streams["o"]
    assert oo.molar_flow * oo.H == pytest.approx(
        f.molar_flow * f.H + rep.duties["q"], rel=1e-9)


# ============================ Sulfur condenser ===============================
def test_sulfur_condenser_balances_and_saturates():
    fs = _furnace()
    fs.add(SulfurCondenser("C1", {"T": 450.0}))
    # rewire furnace -> condenser (rebuild cleanly)
    fs = Flowsheet(components=[Component(c) for c in COMPS],
                   property_package="nasa:claus")
    fs.add(ClausReactor("FURN", {}))
    fs.add(SulfurCondenser("C1", {"T": 450.0}))
    fs.feed("ACID", "FURN:in1", T=320.0, P=1.6e5, molar_flow=100.0,
            z={"hydrogen sulfide": 0.90, "carbon dioxide": 0.10})
    n2 = 45.0 * 79.0 / 21.0
    fs.feed("AIR", "FURN:in2", T=480.0, P=1.6e5, molar_flow=45.0 + n2,
            z={"oxygen": 45.0 / (45.0 + n2), "nitrogen": n2 / (45.0 + n2)})
    fs.connect("g", "FURN:out", "C1:in1")
    fs.connect("q", "FURN:duty", None)
    fs.connect("gas", "C1:gas", None)
    fs.connect("liq", "C1:liquid", None)
    fs.connect("qc", "C1:duty", None)
    rep = fs.solve()
    assert rep.converged
    liq, gas, gin = fs.streams["liq"], fs.streams["gas"], fs.streams["g"]
    # liquid product is elemental sulfur as S8, leaving as a liquid
    assert liq.phase == "liquid" and liq.z["S8"] == pytest.approx(1.0)
    assert liq.molar_flow > 0
    # sulfur-atom balance across the condenser is exact
    s_in = _atom_flows(gin)["S"]
    s_out = _atom_flows(gas)["S"] + liq.molar_flow * 8
    assert s_out == pytest.approx(s_in, rel=1e-9)
    # exit gas is saturated in sulfur (residual S8 set by Psat(T))
    y_s8 = gas.z["S8"]
    assert y_s8 == pytest.approx(liquid_sulfur_psat(450.0) / gas.P, rel=1e-6)
    assert rep.duties["qc"] < 0                   # heat removed (cooling)


def test_condenser_requires_T_and_s8_component():
    fs = Flowsheet(components=[Component(c) for c in COMPS],
                   property_package="nasa:claus")
    fs.add(SulfurCondenser("C", {}))              # no T
    fs.feed("F", "C:in1", T=600.0, P=1.5e5, molar_flow=10.0,
            z={"S2": 0.1, "nitrogen": 0.9})
    fs.connect("g", "C:gas", None)
    fs.connect("l", "C:liquid", None)
    fs.connect("q", "C:duty", None)
    with pytest.raises(SulfurCondenserError, match="params..T.."):
        fs.solve()


def test_condenser_rejects_components_without_s8():
    fs = Flowsheet(components=[Component("S2"), Component("nitrogen")],
                   property_package="nasa:claus")
    fs.add(SulfurCondenser("C", {"T": 450.0}))
    fs.feed("F", "C:in1", T=600.0, P=1.5e5, molar_flow=10.0,
            z={"S2": 0.1, "nitrogen": 0.9})
    fs.connect("g", "C:gas", None)
    fs.connect("l", "C:liquid", None)
    fs.connect("q", "C:duty", None)
    with pytest.raises(SulfurCondenserError, match="requires S8"):
        fs.solve()


# ============================ Full Claus train ===============================
def _claus_train(o2: float = 45.0):
    fs = Flowsheet(components=[Component(c) for c in COMPS],
                   property_package="nasa:claus")
    fs.add(ClausReactor("FURN", {}))
    fs.add(SulfurCondenser("C1", {"T": 450.0}))
    fs.add(Heater("RH1", {"T_out": 530.0}))
    fs.add(ClausReactor("CV1", {"T": 530.0}))
    fs.add(SulfurCondenser("C2", {"T": 440.0}))
    fs.add(Heater("RH2", {"T_out": 480.0}))
    fs.add(ClausReactor("CV2", {"T": 480.0}))
    fs.add(SulfurCondenser("C3", {"T": 405.0}))
    fs.feed("ACID", "FURN:in1", T=320.0, P=1.6e5, molar_flow=100.0,
            z={"hydrogen sulfide": 0.90, "carbon dioxide": 0.10})
    n2 = o2 * 79.0 / 21.0
    fs.feed("AIR", "FURN:in2", T=480.0, P=1.6e5, molar_flow=o2 + n2,
            z={"oxygen": o2 / (o2 + n2), "nitrogen": n2 / (o2 + n2)})
    fs.connect("g0", "FURN:out", "C1:in1")
    fs.connect("q0", "FURN:duty", None)
    fs.connect("L1", "C1:liquid", None)
    fs.connect("qc1", "C1:duty", None)
    fs.connect("g1", "C1:gas", "RH1:in1")
    fs.connect("qr1", "RH1:duty", None)
    fs.connect("g1b", "RH1:out", "CV1:in1")
    fs.connect("qcv1", "CV1:duty", None)
    fs.connect("g2", "CV1:out", "C2:in1")
    fs.connect("L2", "C2:liquid", None)
    fs.connect("qc2", "C2:duty", None)
    fs.connect("g2b", "C2:gas", "RH2:in1")
    fs.connect("qr2", "RH2:duty", None)
    fs.connect("g2c", "RH2:out", "CV2:in1")
    fs.connect("qcv2", "CV2:duty", None)
    fs.connect("g3", "CV2:out", "C3:in1")
    fs.connect("L3", "C3:liquid", None)
    fs.connect("qc3", "C3:duty", None)
    fs.connect("tail", "C3:gas", None)
    return fs


def _recovery(fs) -> float:
    rec = sum(fs.streams[L].molar_flow * 8 for L in ("L1", "L2", "L3"))
    feed = fs.streams["ACID"].molar_flow * fs.streams["ACID"].z["hydrogen sulfide"]
    return rec / feed


def test_full_train_recovery_and_air_demand_optimum():
    """Stoichiometric air (oxidise 1/3 of the H2S: 90·(1/3)·1.5 = 45 mol O2)
    drives the catalytic gas to H2S:SO2 = 2:1 and maximises recovery — the
    Claus air-demand principle. Recovery lands in the literature band for a
    thermal + two-catalytic-bed plant."""
    fs = _claus_train(o2=45.0)
    rep = fs.solve()
    assert rep.converged
    tail = fs.streams["tail"]
    assert tail.z["hydrogen sulfide"] / tail.z["sulfur dioxide"] == pytest.approx(
        2.0, abs=0.1)
    assert 0.95 < _recovery(fs) < 0.99


def test_full_train_air_demand_is_an_optimum():
    # off-stoichiometric air (too little or too much) recovers less sulfur
    r_low = _recovery(_solve(_claus_train(o2=39.0)))
    r_opt = _recovery(_solve(_claus_train(o2=45.0)))
    r_high = _recovery(_solve(_claus_train(o2=51.0)))
    assert r_opt > r_low
    assert r_opt > r_high


def _solve(fs):
    fs.solve()
    return fs


def test_full_train_atom_and_energy_balances_close():
    fs = _claus_train(o2=45.0)
    rep = fs.solve()
    feeds = [fs.streams["ACID"], fs.streams["AIR"]]
    prods = [fs.streams[s] for s in ("L1", "L2", "L3", "tail")]
    # sulfur atom balance
    s_in = _atom_flows(*feeds)["S"]
    s_out = _atom_flows(*prods)["S"]
    assert s_out == pytest.approx(s_in, rel=1e-9)
    # overall energy balance: sum(Hin) + sum(Q) == sum(Hout)
    h_in = sum(s.molar_flow * s.H for s in feeds)
    h_out = sum(s.molar_flow * s.H for s in prods)
    q = sum(rep.duties.values())
    assert h_in + q == pytest.approx(h_out, rel=1e-9)


# ============================ `.flow` round-trip =============================
def test_flow_round_trip_is_exact():
    fs = _claus_train(o2=45.0)
    fs.solve()
    doc = to_dict(fs)
    assert doc == to_dict(from_dict(doc))


def test_furnace_unmapped_flowing_component_raises():
    fs = Flowsheet(components=[Component("hydrogen sulfide"), Component("benzene")],
                   property_package="thermo:PR")    # benzene maps in PR, not nasa
    fs.add(ClausReactor("R", {}))
    fs.feed("F", "R:in1", T=400.0, P=1.5e5, molar_flow=10.0,
            z={"hydrogen sulfide": 0.5, "benzene": 0.5})
    fs.connect("o", "R:out", None)
    fs.connect("q", "R:duty", None)
    with pytest.raises(ClausReactorError, match="benzene"):
        fs.solve()
