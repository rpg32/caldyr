"""M12 tests: the Absorber (sum-rates / Burningham-Otto MESH, no condenser
or reboiler) — gas absorption and stripping.

Validation strategy (sources cited per test):

1. **Kremser closed form** (Hameed, *Chemical Process Simulations using Aspen
   Hysys*, Wiley 2025, eq. (9.1); identical to Seader, Henley & Roper 3e
   ch. 5): for a dilute, near-isothermal absorber the fraction of solute
   absorbed is (A^(N+1) - A)/(A^(N+1) - 1) with A = L/(K V). The MESH
   absorber must reproduce that within a few percent on a system where PR is
   trustworthy (n-pentane from nitrogen into n-decane). Achieved deltas are
   0.3-1.3% over N = 4..10 and A ~ 0.6..1.2.
2. **The book's worked absorber** (Hameed 2025 sec. 9.1: 206 kmol/h of 3%
   SO2 in air + 1.3e5 kg/h water, 20 C / 1 atm, 20 stages, HYSYS-PR). This
   one is reproduced *structurally, not quantitatively*, and the test
   documents exactly why: HYSYS's PR carries a built-in SO2/H2O binary
   interaction parameter that makes K_SO2 ~ 30 at 20 C (hence the book's
   99.07% removal); stock `thermo` PR has no kij for that pair and gives
   K_SO2 ~ 500 (SO2's chemical solubility in water — hydrolysis to H2SO3 —
   cannot be captured by a plain cubic EOS), and the ChemSep NRTL table has
   no SO2/water entry either (ideal fallback, K ~ 3). The dilute-solute slip
   is exponentially sensitive to K, so the 0.0576 kmol/h book number is out
   of reach of any honest stock-thermo run; what *is* reproducible — total
   flows, thermal behavior, machine-exact balances, the structural response
   the book itself highlights (sec. 9.1.2 Note: more stages / more solvent /
   higher P -> more absorption) — is asserted. The Kremser test above
   carries the quantitative burden for the MESH machinery itself.
3. **Stripping** (book sec. 9.2: desorption is the same unit operation with
   the volatile-rich liquid fed at the top and the stripping gas at the
   bottom): the same Absorber unit runs as a stripper with no separate type,
   validated against the Kremser stripping form (book eq. (9.8)).
4. Conservation (machine precision), exact energy closure, `.flow`
   round-trip, typed errors, economics sizing (tower + trays, no
   condenser/reboiler).
"""
import math

import pytest

from caldyr.core import Component, Flowsheet
from caldyr.io import from_dict, to_dict
from caldyr.thermo import make_package
from caldyr.unitops import Absorber, AbsorberError

P_ATM = 101325.0
KMOLH = 1000.0 / 3600.0          # kmol/h -> mol/s


# -- builders -------------------------------------------------------------------
def pentane_absorber(n_stages=6, solvent=40.0, y_c5=0.01, gas=100.0) -> Flowsheet:
    """Dilute oil absorber: n-pentane from nitrogen into n-decane at 1 atm,
    25 C — a system PR handles well (the Kremser benchmark)."""
    fs = Flowsheet(components=[Component("nitrogen"), Component("n-pentane"),
                               Component("n-decane")],
                   property_package="thermo:PR")
    fs.add(Absorber("ABS", {"n_stages": n_stages, "P": P_ATM}))
    fs.feed("GAS", "ABS:gas_in", T=298.15, P=P_ATM, molar_flow=gas,
            z={"nitrogen": 1.0 - y_c5, "n-pentane": y_c5})
    fs.feed("OIL", "ABS:liquid_in", T=298.15, P=P_ATM, molar_flow=solvent,
            z={"n-decane": 1.0})
    fs.connect("GASOUT", "ABS:vapor_out", None)
    fs.connect("LIQOUT", "ABS:liquid_out", None)
    return fs


def so2_absorber(n_stages=20, water_kgh=1.3e5) -> Flowsheet:
    """The book's SO2 absorber (Hameed 2025 sec. 9.1.1): 206 kmol/h of 3 mol%
    SO2 in air, washed with pure water, both at 20 C / 1 atm. 'Air' is not a
    resolvable species in the chemicals database, so it is carried as its
    N2/O2 79/21 split — thermodynamically equivalent for an inert carrier."""
    fs = Flowsheet(components=[Component("water"), Component("nitrogen"),
                               Component("oxygen"), Component("sulfur dioxide")],
                   property_package="thermo:PR")
    fs.add(Absorber("ABS", {"n_stages": n_stages, "P": P_ATM}))
    fs.feed("GAS", "ABS:gas_in", T=293.15, P=P_ATM, molar_flow=206.0 * KMOLH,
            z={"sulfur dioxide": 0.03, "nitrogen": 0.97 * 0.79,
               "oxygen": 0.97 * 0.21})
    fs.feed("WATER", "ABS:liquid_in", T=293.15, P=P_ATM,
            molar_flow=water_kgh / 18.01528 * KMOLH, z={"water": 1.0})
    fs.connect("GASOUT", "ABS:vapor_out", None)
    fs.connect("LIQOUT", "ABS:liquid_out", None)
    return fs


def absorbed_fraction(fs, comp="n-pentane") -> float:
    return fs.units["ABS"].design["absorbed"][comp]


def kremser_fraction(fs, comp, n) -> float:
    """Kremser prediction (book eq. (9.1)) with the effective absorption
    factor A_e = sqrt(A_top * A_bottom) (Edmister's geometric mean) evaluated
    from the *converged* stage states — so the test compares the MESH result
    against the closed form on identical thermodynamics."""
    d = fs.units["ABS"].design
    pp = make_package("thermo:PR", fs.component_ids)
    a = []
    for j in (0, n - 1):
        K = pp.k_values(d["T_profile"][j], d["P_profile"][j],
                        d["x_profile"][j], d["y_profile"][j])[comp]
        a.append(d["L_profile"][j] / (K * d["V_profile"][j]))
    ae = math.sqrt(a[0] * a[1])
    return (ae ** (n + 1) - ae) / (ae ** (n + 1) - 1.0)


# -- 1. Kremser (book eq. 9.1) ----------------------------------------------------
@pytest.mark.parametrize("n,solvent", [(6, 40.0), (10, 40.0), (6, 80.0)])
def test_matches_kremser_dilute_absorber(n, solvent):
    fs = pentane_absorber(n_stages=n, solvent=solvent)
    assert fs.solve().converged
    frac = absorbed_fraction(fs)
    pred = kremser_fraction(fs, "n-pentane", n)
    # Achieved deltas are 0.3-1.3%; Kremser itself assumes constant A, so a
    # few-percent tolerance is the honest bar for a near-isothermal column
    # (the temperature rise here is ~1 K).
    assert frac == pytest.approx(pred, rel=0.04)
    d = fs.units["ABS"].design
    assert max(d["T_profile"]) - min(d["T_profile"]) < 2.0


# -- 2. structure the book highlights (sec. 9.1.2 Note) ----------------------------
def test_more_stages_absorb_more():
    fracs = []
    for n in (4, 8, 12):
        fs = pentane_absorber(n_stages=n)
        assert fs.solve().converged
        fracs.append(absorbed_fraction(fs))
    assert fracs[0] < fracs[1] < fracs[2]


def test_more_solvent_absorbs_more():
    fracs = []
    for solvent in (20.0, 40.0, 80.0):
        fs = pentane_absorber(solvent=solvent)
        assert fs.solve().converged
        fracs.append(absorbed_fraction(fs))
    assert fracs[0] < fracs[1] < fracs[2]


# -- 3. the book's SO2 absorber, structural ----------------------------------------
def test_book_so2_absorber_structure():
    """Hameed 2025 sec. 9.1 (20 stages): book/HYSYS-PR gets Gas_Out =
    204.6 kmol/h at 20.36 C with 0.0576 kmol/h SO2 (99.07% removal). With
    stock PR (no SO2/H2O kij; K_SO2 ~ 500 vs HYSYS's ~30 — see the module
    docstring) the slip is not reproducible, so this test asserts what is:
    convergence, machine-exact balances, the thermal picture (outlet within
    ~1.5 K of the book's 20.36 C — evaporative cooling vs absorption heat),
    water saturation of the exit gas, and a physically ordered profile."""
    fs = so2_absorber()
    rep = fs.solve()
    assert rep.converged
    gas_out, liq_out = fs.streams["GASOUT"], fs.streams["LIQOUT"]
    feed_gas, feed_liq = fs.streams["GAS"], fs.streams["WATER"]

    # Machine-exact component balances.
    for c in fs.component_ids:
        n_in = (feed_gas.molar_flow * feed_gas.z.get(c, 0.0)
                + feed_liq.molar_flow * feed_liq.z.get(c, 0.0))
        n_out = (gas_out.molar_flow * gas_out.z[c]
                 + liq_out.molar_flow * liq_out.z[c])
        assert n_out == pytest.approx(n_in, rel=1e-12, abs=1e-12)
    # Exact energy closure (adiabatic column).
    e_in = feed_gas.molar_flow * feed_gas.H + feed_liq.molar_flow * feed_liq.H
    e_out = gas_out.molar_flow * gas_out.H + liq_out.molar_flow * liq_out.H
    assert e_out == pytest.approx(e_in, rel=1e-12)
    assert fs.units["ABS"].design["energy_residual_rel"] < 1e-9

    # Thermal picture: both products stay within ~1.5 K of the book's 20.36 C
    # gas outlet (the feeds are at 20 C; little heat moves at this dilution).
    assert gas_out.T == pytest.approx(20.36 + 273.15, abs=1.5)
    assert liq_out.T == pytest.approx(293.15, abs=1.5)
    # The exit gas leaves water-saturated (the book's 204.6 kmol/h vs
    # 206 kmol/h fed already nets SO2 removal against water pickup).
    assert gas_out.z["water"] > 0.01
    # Total gas product within a few % of the book's 204.6 kmol/h (ours is
    # higher because the SO2 stays in the gas: 209.7 kmol/h, +2.5%).
    assert gas_out.molar_flow == pytest.approx(204.6 * KMOLH, rel=0.04)
    # Some SO2 is absorbed and the response direction is right, but with
    # K_SO2 ~ 500 the absorbed fraction is ~7%, not the book's 99% — the
    # honest stock-thermo result (see module docstring).
    frac = fs.units["ABS"].design["absorbed"]["sulfur dioxide"]
    assert 0.0 < frac < 0.5


# -- 4. stripping: the same unit, mirrored (book sec. 9.2) --------------------------
def test_stripper_is_the_same_unit():
    """Gas desorption (book sec. 9.2): feed the volatile-rich liquid at the
    top ('liquid_in') and the stripping gas at the bottom ('gas_in') — the
    same Absorber unit strips. Checked against the Kremser stripping form,
    book eq. (9.8): fraction stripped = (S^(N+1) - S)/(S^(N+1) - 1) with
    S = K V / L."""
    n = 8
    fs = Flowsheet(components=[Component("nitrogen"), Component("n-pentane"),
                               Component("n-decane")],
                   property_package="thermo:PR")
    fs.add(Absorber("STR", {"n_stages": n, "P": P_ATM}))
    # Volatile-rich liquid on top, stripping nitrogen at the bottom.
    fs.feed("RICH", "STR:liquid_in", T=298.15, P=P_ATM, molar_flow=40.0,
            z={"n-decane": 0.98, "n-pentane": 0.02})
    fs.feed("N2", "STR:gas_in", T=298.15, P=P_ATM, molar_flow=60.0,
            z={"nitrogen": 1.0})
    fs.connect("LEAN", "STR:liquid_out", None)
    fs.connect("OFFGAS", "STR:vapor_out", None)
    assert fs.solve().converged

    rich, lean = fs.streams["RICH"], fs.streams["LEAN"]
    c5_in = rich.molar_flow * rich.z["n-pentane"]
    c5_out = lean.molar_flow * lean.z["n-pentane"]
    stripped = 1.0 - c5_out / c5_in
    assert stripped > 0.5                       # the stripping gas works

    # Kremser stripping prediction with the effective S = K V / L
    # (geometric mean of top and bottom stage states).
    d = fs.units["STR"].design
    pp = make_package("thermo:PR", fs.component_ids)
    s = []
    for j in (0, n - 1):
        K = pp.k_values(d["T_profile"][j], d["P_profile"][j],
                        d["x_profile"][j], d["y_profile"][j])["n-pentane"]
        s.append(K * d["V_profile"][j] / d["L_profile"][j])
    se = math.sqrt(s[0] * s[1])
    pred = (se ** (n + 1) - se) / (se ** (n + 1) - 1.0)
    # Evaporative cooling makes this less isothermal than the absorber case;
    # 7% is the honest Kremser tolerance here (achieved ~2-4%).
    assert stripped == pytest.approx(pred, rel=0.07)


# -- 5. result integrity, IO, errors -------------------------------------------------
def test_profiles_present_and_consistent():
    fs = pentane_absorber()
    assert fs.solve().converged
    d = fs.units["ABS"].design
    n = d["n_stages"]
    for key in ("T_profile", "P_profile", "L_profile", "V_profile",
                "x_profile", "y_profile"):
        assert isinstance(d[key], list) and len(d[key]) == n
    for row in d["x_profile"] + d["y_profile"]:
        assert sum(row.values()) == pytest.approx(1.0, abs=1e-9)
    # The traffic honors the feeds: V_N+1-equivalent is the gas feed entering
    # stage N, L_0-equivalent the solvent; products match the profiles.
    assert d["V_top"] == d["V_profile"][0]
    assert d["L_bottom"] == pytest.approx(d["L_profile"][-1], rel=1e-6)
    assert d["N"] == float(n)


def test_flow_round_trip_is_exact():
    fs = pentane_absorber()
    fs.solve()
    doc = to_dict(fs)
    assert doc == to_dict(from_dict(doc))
    fs2 = from_dict(doc)
    fs2.solve()
    assert fs2.streams["GASOUT"].z["n-pentane"] == pytest.approx(
        fs.streams["GASOUT"].z["n-pentane"], rel=1e-9)


def test_missing_inlet_raises():
    fs = Flowsheet(components=[Component("nitrogen"), Component("n-pentane"),
                               Component("n-decane")],
                   property_package="thermo:PR")
    fs.add(Absorber("ABS", {"n_stages": 4}))
    fs.feed("GAS", "ABS:gas_in", T=298.15, P=P_ATM, molar_flow=100.0,
            z={"nitrogen": 0.99, "n-pentane": 0.01})
    fs.connect("GASOUT", "ABS:vapor_out", None)
    fs.connect("LIQOUT", "ABS:liquid_out", None)
    with pytest.raises(AbsorberError, match="liquid_in"):
        fs.solve()


def test_bad_params_raise():
    for bad in (None, 0, -3):
        fs = pentane_absorber()
        fs.units["ABS"].params["n_stages"] = bad
        with pytest.raises(AbsorberError, match="n_stages"):
            fs.solve()
    fs = pentane_absorber()
    fs.units["ABS"].params["dP_stage"] = -10.0
    with pytest.raises(AbsorberError, match="dP_stage"):
        fs.solve()


# -- 6. thermo: the new per-phase calls are consistent with bubble_point --------------
def test_k_values_match_bubble_point_at_saturation():
    """pp.k_values (phi-phi at arbitrary T,P) must agree with the K-values
    implied by the saturated VF=0 flash when evaluated exactly at the bubble
    point — same EOS, two code paths."""
    pp = make_package("thermo:PR", ["benzene", "toluene"])
    x = {"benzene": 0.5, "toluene": 0.5}
    res = pp.bubble_point(P_ATM, x)
    assert res.y is not None and res.x is not None
    k_flash = {c: res.y[c] / res.x[c] for c in x}
    k_phi = pp.k_values(res.T, P_ATM, x, res.y)
    for c in x:
        assert k_phi[c] == pytest.approx(k_flash[c], rel=1e-6)
    # Per-phase enthalpies at saturation match the flash's phase enthalpies.
    assert pp.enthalpy_liquid(res.T, P_ATM, x) == pytest.approx(
        res.H_liquid, rel=1e-9)
    assert pp.enthalpy_vapor(res.T, P_ATM, res.y) == pytest.approx(
        res.H_vapor, rel=1e-9)


# -- 7. economics: tower + trays, no condenser/reboiler -------------------------------
def test_absorber_is_sized_as_tower_plus_trays():
    import math as _math

    from caldyr.economics.sizing import (
        ABSORBER_TRAY_EFFICIENCY, SizingOptions, size_flowsheet,
    )

    fs = pentane_absorber(n_stages=8)
    rep = fs.solve()
    pp = make_package("thermo:PR", fs.component_ids)
    sizes = size_flowsheet(fs, rep, pp, SizingOptions())
    by_id = {s.unit_id: s for s in sizes}
    assert set(by_id) == {"ABS", "ABS.trays"}
    assert by_id["ABS"].equipment_type == "vessel_vertical"
    assert by_id["ABS.trays"].equipment_type == "tray_sieve"
    # 8 ideal stages at the (low) absorption efficiency of 0.4 -> 20 trays.
    assert by_id["ABS.trays"].quantity == _math.ceil(8 / ABSORBER_TRAY_EFFICIENCY)
    assert by_id["ABS"].diameter_m > 0.1
