"""M16 test: the amine acid-gas PropertyPackage wrapper + a DEA sweetening absorber.

The reactive solubility model (``caldyr.thermo.amine``) is validated against
solubility data in ``tests/test_m16_amine_kent_eisenberg.py``. Here we test the
package that makes it usable by a column — :class:`caldyr.thermo.AmineAcidGasPackage`
(selector ``"amine:DEA"`` / ``"amine:MDEA"``) — and run it through the sum-rates
:class:`caldyr.unitops.absorber.Absorber` on a Hameed 2025 §15.3-style gas
sweetening case.

The book gives no numeric output slate for §15.3, so (per its own emphasis) the
oracle is the underlying physics:
  * the acid-gas K-value *collapses at low loading* (the reactive signature that
    makes amine absorption so favourable) and rises as the solvent saturates;
  * water/amine follow Raoult (amine essentially non-volatile), light gases
    Henry's law (large K, sparingly soluble);
  * the liquid is below the vapour by the heat of absorption (the absorber
    temperature bulge);
  * the absorber **converges, closes its mass balance to machine precision,
    closes its energy balance, removes the acid gas, and lands a rich loading in
    the band the solubility model itself predicts.**
"""
import pytest

from caldyr.core import Stream
from caldyr.thermo import AmineAcidGasPackage, make_package
from caldyr.thermo.amine import acid_gas_partial_pressures
from caldyr.unitops.absorber import Absorber

_COMPS = ["DEA", "water", "CO2", "H2S", "methane"]


def _pkg(amine: str = "DEA", comps: list[str] | None = None) -> AmineAcidGasPackage:
    return make_package(f"amine:{amine}", comps or _COMPS)  # type: ignore[return-value]


# -- make_package wiring + validation -----------------------------------------
def test_selector_builds_the_amine_package():
    for amine in ("DEA", "MDEA"):
        comps = [amine, "water", "CO2", "H2S", "methane"]
        pp = make_package(f"amine:{amine}", comps)
        assert isinstance(pp, AmineAcidGasPackage)
        assert pp.amine == amine


def test_rejects_unknown_amine_and_missing_roles():
    with pytest.raises(ValueError, match="unsupported amine"):
        make_package("amine:MEA", _COMPS)
    with pytest.raises(ValueError, match="no 'water'"):
        AmineAcidGasPackage(["DEA", "CO2", "methane"], "DEA")
    with pytest.raises(ValueError, match="no 'amine'"):
        AmineAcidGasPackage(["water", "CO2", "methane"], "DEA")
    with pytest.raises(ValueError, match="at least one acid gas"):
        AmineAcidGasPackage(["DEA", "water", "methane"], "DEA")


# -- the reactive K-value signature -------------------------------------------
def test_acid_gas_k_collapses_at_low_loading():
    """The whole point of a chemical solvent: at low loading the acid-gas
    K-value is tiny (avid absorption); it rises steeply as the amine saturates."""
    pp = _pkg()
    P = 6.9e6
    lean = {"DEA": 0.062, "water": 0.936, "CO2": 0.0008, "H2S": 0.0004, "methane": 0.0008}
    rich = {"DEA": 0.062, "water": 0.905, "CO2": 0.020, "H2S": 0.010, "methane": 0.003}
    K_lean = pp.k_values(313.15, P, lean, lean)
    K_rich = pp.k_values(313.15, P, rich, rich)
    for gas in ("CO2", "H2S"):
        assert K_lean[gas] < K_rich[gas]      # K rises with loading
        assert K_lean[gas] < 0.1              # avid at low loading


def test_nonreactive_components_follow_raoult_and_henry():
    pp = _pkg()
    P = 6.9e6
    x = {"DEA": 0.062, "water": 0.92, "CO2": 0.01, "H2S": 0.005, "methane": 0.003}
    K = pp.k_values(313.15, P, x, x)
    assert K["methane"] > 100.0               # sparingly soluble light gas
    assert K["DEA"] < 1e-3                     # essentially non-volatile amine
    # water ~ Raoult: Psat(40 C) ~ 7.4 kPa over 6.9 MPa
    assert K["water"] == pytest.approx(7.4e3 / P, rel=0.2)


# -- enthalpy basis: the heat of absorption -----------------------------------
def test_absorbing_acid_gas_releases_the_heat_of_absorption():
    """Moving CO2 from the vapour to the absorbed liquid state must release ~the
    heat of absorption (the liquid reference sits ~ΔH_abs below the gas)."""
    pp = _pkg()
    P, T = 6.9e6, 313.15
    hV_co2 = pp.enthalpy_vapor(T, P, {"CO2": 1.0})
    hL_co2 = pp.enthalpy_liquid(T, P, {"CO2": 1.0})
    assert hL_co2 - hV_co2 == pytest.approx(-70.0e3, rel=0.05)   # DEA CO2 ΔH_abs


# -- flashes ------------------------------------------------------------------
def test_phase_identification_and_ph_roundtrip():
    pp = _pkg()
    P = 6.9e6
    gas = {"methane": 0.94, "CO2": 0.04, "H2S": 0.02}
    solv = {"DEA": 0.062, "water": 0.938}
    assert pp.flash_pt(313.15, P, gas).phase == "vapor"
    assert pp.flash_pt(313.15, P, solv).phase == "liquid"
    # PH flash recovers the temperature its own PT enthalpy was taken at.
    z = {"DEA": 0.06, "water": 0.92, "CO2": 0.012, "H2S": 0.005, "methane": 0.003}
    H = pp.enthalpy(330.0, P, z)
    assert pp.flash_ph(P, H, z).T == pytest.approx(330.0, abs=0.1)


def test_entropy_and_three_phase_are_unsupported():
    pp = _pkg()
    z = {"DEA": 0.06, "water": 0.93, "CO2": 0.01}
    for call in (lambda: pp.entropy(313.15, 6.9e6, z),
                 lambda: pp.flash_ps(6.9e6, 0.0, z),
                 lambda: pp.flash_pt_3p(313.15, 6.9e6, z)):
        with pytest.raises(NotImplementedError):
            call()


# -- the §15.3-style DEA sweetening absorber ----------------------------------
def _sweetening_absorber():
    """A 20-stage DEA absorber on a sour natural gas (Hameed §15.3 proportions:
    ~4.1% CO2, ~1.7% H2S, balance methane) with 28 wt% (~2.7 M) aqueous DEA."""
    pp = _pkg()
    P = 6.9e6                       # ~1000 psia
    z_g = {"methane": 0.9415, "CO2": 0.0413, "H2S": 0.0172}
    gas = Stream(id="sour", components=_COMPS, T=313.15, P=P, molar_flow=350.0,
                 z={c: z_g.get(c, 0.0) for c in _COMPS})
    z_l = {"DEA": 0.062, "water": 0.9355, "CO2": 0.0017, "H2S": 0.0008}
    liq = Stream(id="lean", components=_COMPS, T=313.15, P=P, molar_flow=850.0,
                 z={c: z_l.get(c, 0.0) for c in _COMPS})
    ab = Absorber("ABS", {"n_stages": 20})
    out = ab.solve({"gas_in": gas, "liquid_in": liq}, pp)
    return pp, ab, out, gas, liq


def test_dea_absorber_converges_and_sweetens_the_gas():
    pp, ab, out, gas, liq = _sweetening_absorber()
    d = ab.design
    assert d["energy_residual_rel"] < 1e-6        # energy balance closed
    # Acid gas removed; the inert methane essentially passes through.
    assert d["absorbed"]["CO2"] > 0.95
    assert d["absorbed"]["H2S"] > 0.95
    assert d["absorbed"]["methane"] < 0.02
    # Heat of absorption warms the solvent down the column (temperature bulge).
    assert d["T_bottom"] > d["T_top"] + 5.0


def test_dea_absorber_closes_mass_balance_to_machine_precision():
    pp, ab, out, gas, liq = _sweetening_absorber()
    v, ll = out["vapor_out"], out["liquid_out"]
    F_g = gas.molar_flow
    F_l = liq.molar_flow
    for c in _COMPS:
        inp = F_g * gas.z.get(c, 0.0) + F_l * liq.z.get(c, 0.0)
        outp = v.molar_flow * v.z[c] + ll.molar_flow * ll.z[c]
        assert abs(inp - outp) < 1e-9 * (F_g + F_l)


def test_dea_rich_loading_matches_the_solubility_model():
    """The converged rich loading lands in the DEA band, and the CO2 partial
    pressure the column implies over the rich liquid (y_CO2 * P at the bottom
    stage) is consistent with the underlying solubility model's prediction at
    that loading, temperature and amine molarity."""
    pp, ab, out, gas, liq = _sweetening_absorber()
    d = ab.design
    xb = d["x_bottom"]
    alpha_co2 = xb["CO2"] / xb["DEA"]
    alpha_h2s = xb["H2S"] / xb["DEA"]
    assert 0.15 < alpha_co2 < 0.6                  # physical DEA rich-loading band
    assert 0.02 < alpha_h2s < 0.6

    # 28 wt% DEA ~ 2.7 M, from the converged bottom-stage composition.
    xN = d["x_profile"][-1]
    M = pp._amine_molarity(xN)
    assert M == pytest.approx(2.7, rel=0.25)

    # Self-consistency: the bottom-stage CO2 partial pressure the column carries
    # (y_CO2 * P) equals the solubility model's *coupled* (CO2+H2S competing)
    # prediction at the converged bottom loading, T and molarity.
    T_bot = d["T_profile"][-1]
    p_col_kpa = d["y_profile"][-1]["CO2"] * d["P_profile"][-1] / 1000.0
    p_model_kpa, _ = acid_gas_partial_pressures(
        T_bot, xN["CO2"] / xN["DEA"], xN["H2S"] / xN["DEA"], M, "DEA")
    assert p_col_kpa == pytest.approx(p_model_kpa, rel=0.05)


def test_murphree_efficiency_unlocks_mdea_h2s_selectivity():
    """MDEA's prized H2S-over-CO2 selectivity is *kinetic* — CO2 reacts slowly
    with the tertiary amine. An equilibrium stage absorbs both fully; a Murphree
    efficiency (book §15.3: E_CO2 ~ 0.15, E_H2S ~ 0.8) suppresses CO2 uptake and
    lets H2S through preferentially, with the energy balance still closed."""
    pp = make_package("amine:MDEA", ["MDEA", "water", "CO2", "H2S", "methane"])
    comps = ["MDEA", "water", "CO2", "H2S", "methane"]
    P = 6.9e6
    z_g = {"methane": 0.9415, "CO2": 0.0413, "H2S": 0.0172}
    z_l = {"MDEA": 0.062, "water": 0.9355, "CO2": 0.0017, "H2S": 0.0008}

    def run(murphree=None):
        gas = Stream(id="g", components=comps, T=313.15, P=P, molar_flow=350.0,
                     z={c: z_g.get(c, 0.0) for c in comps})
        liq = Stream(id="l", components=comps, T=313.15, P=P, molar_flow=600.0,
                     z={c: z_l.get(c, 0.0) for c in comps})
        p = {"n_stages": 16}
        if murphree:
            p["murphree"] = murphree
        ab = Absorber("ABS", p)
        ab.solve({"gas_in": gas, "liquid_in": liq}, pp)
        return ab.design

    eq = run()
    assert eq["absorbed"]["CO2"] > 0.95 and eq["absorbed"]["H2S"] > 0.95  # both
    kin = run({"CO2": 0.15, "H2S": 0.8})
    assert kin["energy_residual_rel"] < 1e-6
    assert kin["absorbed"]["CO2"] < eq["absorbed"]["CO2"] - 0.05   # CO2 slips
    assert kin["absorbed"]["H2S"] > kin["absorbed"]["CO2"]          # H2S-selective


def _stream(flows, T, P):
    F = sum(flows.values())
    return Stream(id="s", components=_COMPS, T=T, P=P, molar_flow=F,
                  z={c: flows.get(c, 0.0) / F for c in _COMPS})


def _flows(s):
    return {c: s.molar_flow * s.z[c] for c in _COMPS}


def test_ns_regenerator_strips_the_rich_amine():
    """The amine REGENERATOR (reactive desorption) cannot be solved by sum-rates
    (it limit-cycles) or by the bubble-point method (degenerate hV<hL); the
    Naphtali-Sandholm simultaneous-correction method converges it. A steam
    stripper at ~1.8 bar regenerates a rich DEA: lean loading driven low, acid
    gas to the overhead, mass + energy balances closed."""
    pp = _pkg()
    Prg = 1.8e5
    rich = {"DEA": 52.0, "water": 800.0, "CO2": 14.5, "H2S": 6.0, "methane": 0.05}
    steam = {"water": 140.0}
    reg = Absorber("REGEN", {"n_stages": 8, "method": "naphtali_sandholm",
                             "max_iter": 80})
    out = reg.solve({"liquid_in": _stream(rich, 388.0, Prg),
                     "gas_in": _stream(steam, 396.0, Prg)}, pp)
    d = reg.design
    assert d["energy_residual_rel"] < 1e-6
    lean = _flows(out["liquid_out"])
    assert lean["CO2"] / lean["DEA"] < 0.02          # CO2 stripped out
    assert lean["H2S"] / lean["DEA"] < 0.02          # H2S stripped out
    assert 380.0 < d["T_bottom"] < 410.0             # physical regen temperature
    # mass balance over the regenerator
    feed = {c: rich.get(c, 0.0) + steam.get(c, 0.0) for c in _COMPS}
    acid = _flows(out["vapor_out"])
    for c in _COMPS:
        assert abs(feed[c] - lean[c] - acid[c]) < 1e-7 * sum(feed.values())


def test_amine_recycle_loop_closes():
    """The full §15.3 loop: sweetening absorber -> rich/lean heat-up ->
    NS regenerator -> cooled lean-amine recycle, with water makeup holding the
    circulating water (open-steam stripping removes more water than the steam
    adds, exactly as a real unit's makeup compensates). The recycle tear (lean
    amine) converges and the acid gas is removed then rejected to the overhead."""
    pp = _pkg()
    Pabs, Prg = 6.9e6, 1.8e5
    z_g = {"methane": 0.9415, "CO2": 0.0413, "H2S": 0.0172}
    Fsour, Fsteam, W = 350.0, 140.0, 800.0
    sour = {c: Fsour * z_g.get(c, 0.0) for c in _COMPS}
    lean = {"DEA": 52.0, "water": W, "CO2": 0.05, "H2S": 0.02, "methane": 0.0}

    err = 1.0
    for _ in range(12):
        ab = Absorber("ABS", {"n_stages": 16, "murphree": {"CO2": 0.6}})
        o = ab.solve({"gas_in": _stream(sour, 313.15, Pabs),
                      "liquid_in": _stream(lean, 313.15, Pabs)}, pp)
        rich = _flows(o["liquid_out"])
        reg = Absorber("REG", {"n_stages": 8, "method": "naphtali_sandholm",
                               "max_iter": 80})
        og = reg.solve({"liquid_in": _stream(rich, 388.0, Prg),
                        "gas_in": _stream({"water": Fsteam}, 396.0, Prg)}, pp)
        new = _flows(og["liquid_out"])
        new["water"] = W                       # makeup water closes the balance
        err = sum(abs(new[c] - lean[c]) for c in _COMPS) / sum(lean.values())
        lean = {c: 0.5 * lean[c] + 0.5 * new[c] for c in _COMPS}
        if err < 1e-4:
            break
    assert err < 1e-4                                  # the recycle tear converged
    assert lean["CO2"] / lean["DEA"] < 0.01           # regenerated lean solvent
    assert rich["CO2"] / rich["DEA"] > 0.1            # loaded rich solvent
    # acid gas is removed from the gas and rejected to the regenerator overhead
    assert (1.0 - o["vapor_out"].z["CO2"] * o["vapor_out"].molar_flow
            / sour["CO2"]) > 0.9


def test_regenerator_reflux_condenser_dries_the_acid_gas():
    """The amine regenerator's open-steam overhead is water-rich (most of the
    stripping steam leaves with the acid gas). A partial **reflux condenser** on
    the top stage (``condenser_T``) condenses that water back as internal reflux,
    so the acid-gas product leaves dried — the natural §15.3 refinement. The NS
    solve stays robust (the condenser alone, not combined with a reboiler); the
    overhead temperature pins to ``condenser_T``, the duty is reported, and mass
    balance stays machine-exact."""
    pp = _pkg()
    Prg = 1.8e5
    rich = {"DEA": 52.0, "water": 800.0, "CO2": 14.5, "H2S": 6.0, "methane": 0.0}
    steam = {"water": 140.0}

    wet = Absorber("REGw", {"n_stages": 9, "method": "naphtali_sandholm",
                            "max_iter": 80})
    ow = wet.solve({"liquid_in": _stream(rich, 388.0, Prg),
                    "gas_in": _stream(steam, 396.0, Prg)}, pp)

    dry = Absorber("REGd", {"n_stages": 9, "method": "naphtali_sandholm",
                            "max_iter": 120, "condenser_T": 320.0})
    od = dry.solve({"liquid_in": _stream(rich, 388.0, Prg),
                    "gas_in": _stream(steam, 396.0, Prg)}, pp)
    d = dry.design
    assert d["energy_residual_rel"] < 1e-6
    assert d["T_top"] == pytest.approx(320.0, abs=1e-3)
    assert d["condenser_duty"] < 0.0                       # heat removed
    # the overhead is dramatically drier than the open-steam (no-condenser) case
    assert ow["vapor_out"].z["water"] > 0.7               # wet
    assert od["vapor_out"].z["water"] < 0.2               # dried
    # and the water carried off with the acid gas collapses
    w_wet = ow["vapor_out"].molar_flow * ow["vapor_out"].z["water"]
    w_dry = od["vapor_out"].molar_flow * od["vapor_out"].z["water"]
    assert w_dry < 0.1 * w_wet
    # still a regenerated lean solvent, mass balance machine-exact
    lean = _flows(od["liquid_out"])
    assert lean["CO2"] / lean["DEA"] < 0.01
    feed = {c: rich.get(c, 0.0) + steam.get(c, 0.0) for c in _COMPS}
    acid = _flows(od["vapor_out"])
    for c in _COMPS:
        assert abs(feed[c] - lean[c] - acid[c]) < 1e-7 * sum(feed.values())


def test_more_solvent_absorbs_more():
    """Structural monotonicity: a higher solvent rate removes at least as much
    acid gas (a sanity check that the package drives the column physically)."""
    pp = _pkg()
    P = 6.9e6
    z_g = {"methane": 0.9415, "CO2": 0.0413, "H2S": 0.0172}
    z_l = {"DEA": 0.062, "water": 0.9355, "CO2": 0.0017, "H2S": 0.0008}

    def run(F_l: float) -> float:
        gas = Stream(id="sour", components=_COMPS, T=313.15, P=P, molar_flow=350.0,
                     z={c: z_g.get(c, 0.0) for c in _COMPS})
        liq = Stream(id="lean", components=_COMPS, T=313.15, P=P, molar_flow=F_l,
                     z={c: z_l.get(c, 0.0) for c in _COMPS})
        ab = Absorber("ABS", {"n_stages": 12})
        ab.solve({"gas_in": gas, "liquid_in": liq}, pp)
        return ab.design["absorbed"]["CO2"]

    assert run(700.0) <= run(1100.0) + 1e-9
