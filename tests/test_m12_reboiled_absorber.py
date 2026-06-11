"""M12 tests: the ReboiledAbsorber (stripping tower).

Validated against the worked example of Hameed, *Chemical Process
Simulations using Aspen Hysys*, Wiley 2025, sec. 9.3.5: a 400 kmol/h
saturated-liquid feed of 55 mol% n-pentane / 45 mol% n-heptane fed to the
top of a stripping tower (HYSYS: 8 stages + reboiler -> n_stages = 9 here),
top pressure 101.3 kPa, reboiler 110 kPa, Peng-Robinson.

Book results (figs. 9.36-9.38):
  * Ovhd Prod Rate spec = 220 kmol/h  ->  x(nC5) in vapor = 0.899,
    x(nC7) in bottoms = 0.876, boilup ratio ~ 1.093.
  * 99% nC5 recovery overhead       ->  Vap 242.4 kmol/h at 52.62 C with
    x(nC5) = 0.8985; bottoms 157.6 kmol/h at 99.32 C with x(nC7) = 0.9860.

PR is excellent for this alkane pair, and the achieved agreement is
essentially exact: compositions within 0.0005 mole fraction, temperatures
within 0.1 K, rates by construction (well inside the 2-10% HYSYS-vs-thermo
tolerance band). Conservation, spec round-trips, `.flow` IO, typed errors,
and economics sizing (tower + trays + reboiler, no condenser) follow.
"""
import pytest

from caldyr.core import Component, Flowsheet
from caldyr.io import from_dict, to_dict
from caldyr.thermo import make_package
from caldyr.unitops import AbsorberError, ReboiledAbsorber

KMOLH = 1000.0 / 3600.0          # kmol/h -> mol/s
Z_FEED = {"n-pentane": 0.55, "n-heptane": 0.45}


def feed_bubble_T() -> float:
    pp = make_package("thermo:PR", ["n-pentane", "n-heptane"])
    bub, _ = pp.bubble_dew(110000.0, Z_FEED)
    return bub


def book_tower(**spec) -> Flowsheet:
    """The book's stripping tower: 9 stages incl. the reboiler, saturated
    liquid feed at the top, 101.3 kPa top / 110 kPa reboiler."""
    fs = Flowsheet(components=[Component("n-pentane"), Component("n-heptane")],
                   property_package="thermo:PR")
    params = {"n_stages": 9, "P": 101300.0, "P_bottom": 110000.0}
    params.update(spec)
    fs.add(ReboiledAbsorber("RA", params))
    fs.feed("FEED", "RA:feed", T=feed_bubble_T(), P=110000.0,
            molar_flow=400.0 * KMOLH, z=dict(Z_FEED))
    fs.connect("VAP", "RA:vapor_out", None)
    fs.connect("BOT", "RA:bottoms", None)
    fs.connect("QR", "RA:reboiler_duty", None)
    return fs


# -- 1. the book's two operating points ---------------------------------------------
def test_book_overhead_rate_spec():
    """Ovhd Prod Rate = 220 kmol/h (book fig. 9.36): x(nC5)_vap = 0.899,
    x(nC7)_btm = 0.876, boilup ~ 1.093. Achieved: 0.8989 / 0.8765 / 1.092."""
    fs = book_tower(vapor_rate=220.0 * KMOLH)
    assert fs.solve().converged
    vap, bot = fs.streams["VAP"], fs.streams["BOT"]
    d = fs.units["RA"].design
    assert vap.molar_flow == pytest.approx(220.0 * KMOLH, rel=1e-12)
    assert bot.molar_flow == pytest.approx(180.0 * KMOLH, rel=1e-12)
    assert vap.z["n-pentane"] == pytest.approx(0.899, abs=0.005)
    assert bot.z["n-heptane"] == pytest.approx(0.876, abs=0.005)
    assert d["boilup_ratio"] == pytest.approx(1.093, rel=0.03)
    assert d["Q_reboiler"] > 0.0


def test_book_99pct_recovery_point():
    """The book's 99%-recovery answer (figs. 9.37-9.38): Vap 242.4 kmol/h at
    52.62 C, x(nC5) = 0.8985; bottoms 157.6 kmol/h at 99.32 C,
    x(nC7) = 0.9860. Specifying the book's converged overhead rate must
    reproduce its compositions, temperatures and the 99% recovery."""
    fs = book_tower(vapor_rate=242.4 * KMOLH)
    assert fs.solve().converged
    vap, bot = fs.streams["VAP"], fs.streams["BOT"]
    assert vap.z["n-pentane"] == pytest.approx(0.8985, abs=0.005)
    assert bot.z["n-heptane"] == pytest.approx(0.9860, abs=0.005)
    assert vap.T - 273.15 == pytest.approx(52.62, abs=1.0)
    assert bot.T - 273.15 == pytest.approx(99.32, abs=1.0)
    recovery = (vap.molar_flow * vap.z["n-pentane"]
                / (400.0 * KMOLH * Z_FEED["n-pentane"]))
    assert recovery == pytest.approx(0.99, abs=0.005)


# -- 2. the three heat-input specs are consistent -------------------------------------
def test_boilup_and_duty_specs_round_trip():
    fs = book_tower(vapor_rate=220.0 * KMOLH)
    fs.solve()
    d = fs.units["RA"].design

    fs_b = book_tower(boilup_ratio=d["boilup_ratio"])
    fs_b.solve()
    assert fs_b.streams["VAP"].molar_flow == pytest.approx(
        220.0 * KMOLH, rel=1e-4)

    fs_q = book_tower(reboiler_duty=d["Q_reboiler"])
    fs_q.solve()
    assert fs_q.streams["VAP"].molar_flow == pytest.approx(
        220.0 * KMOLH, rel=1e-4)
    assert fs_q.streams["VAP"].z["n-pentane"] == pytest.approx(
        fs.streams["VAP"].z["n-pentane"], rel=1e-4)


def test_more_boilup_strips_more():
    recov = []
    for s in (0.8, 1.1, 1.5):
        fs = book_tower(boilup_ratio=s)
        assert fs.solve().converged
        vap = fs.streams["VAP"]
        recov.append(vap.molar_flow * vap.z["n-pentane"])
    assert recov[0] < recov[1] < recov[2]


# -- 3. conservation and result integrity ---------------------------------------------
def test_balances_close_exactly():
    fs = book_tower(vapor_rate=220.0 * KMOLH)
    rep = fs.solve()
    feed, vap, bot = fs.streams["FEED"], fs.streams["VAP"], fs.streams["BOT"]
    for c in feed.components:
        n_in = feed.molar_flow * feed.z[c]
        n_out = vap.molar_flow * vap.z[c] + bot.molar_flow * bot.z[c]
        assert n_out == pytest.approx(n_in, rel=1e-12, abs=1e-12)
    # F h_F + Q_reb = V h_V + B h_B — exact (Q_reb closes the balance).
    e_in = feed.molar_flow * feed.H + rep.duties["QR"]
    e_out = vap.molar_flow * vap.H + bot.molar_flow * bot.H
    assert e_out == pytest.approx(e_in, rel=1e-12)
    assert rep.duties["QR"] > 0.0
    # The independent reboiler-stage balance agrees to the MESH tolerance.
    assert fs.units["RA"].design["energy_residual_rel"] < 1e-3


def test_profiles_and_design_keys():
    fs = book_tower(vapor_rate=220.0 * KMOLH)
    assert fs.solve().converged
    d = fs.units["RA"].design
    n = d["n_stages"]
    assert n == 9
    assert d["N"] == 8.0                    # trays (the reboiler is no tray)
    for key in ("T_profile", "P_profile", "L_profile", "V_profile",
                "x_profile", "y_profile"):
        assert isinstance(d[key], list) and len(d[key]) == n
    T = d["T_profile"]
    assert all(a < b for a, b in zip(T, T[1:]))      # hotter going down
    assert d["P_profile"][0] == pytest.approx(101300.0)
    assert d["P_profile"][-1] == pytest.approx(110000.0)
    assert d["V_top"] == pytest.approx(220.0 * KMOLH, rel=1e-12)


def test_flow_round_trip_is_exact():
    fs = book_tower(vapor_rate=220.0 * KMOLH)
    fs.solve()
    doc = to_dict(fs)
    assert doc == to_dict(from_dict(doc))
    fs2 = from_dict(doc)
    fs2.solve()
    assert fs2.streams["VAP"].z["n-pentane"] == pytest.approx(
        fs.streams["VAP"].z["n-pentane"], rel=1e-9)


# -- 4. typed errors -------------------------------------------------------------------
def test_spec_validation():
    with pytest.raises(AbsorberError, match="exactly one"):
        book_tower().solve()
    with pytest.raises(AbsorberError, match="exactly one"):
        book_tower(vapor_rate=10.0, boilup_ratio=1.0).solve()
    with pytest.raises(AbsorberError, match="strictly between"):
        book_tower(vapor_rate=500.0 * KMOLH).solve()
    with pytest.raises(AbsorberError, match="n_stages"):
        book_tower(vapor_rate=10.0, n_stages=1).solve()
    with pytest.raises(AbsorberError, match="feed_stage"):
        book_tower(vapor_rate=10.0, feed_stage=9).solve()
    with pytest.raises(AbsorberError, match="not both"):
        book_tower(vapor_rate=10.0, dP_stage=100.0).solve()


# -- 5. economics: tower + trays + reboiler --------------------------------------------
def test_sized_as_tower_trays_reboiler():
    import math

    from caldyr.economics import TEAConfig, analyze
    from caldyr.economics.sizing import STRIPPER_TRAY_EFFICIENCY

    fs = book_tower(vapor_rate=220.0 * KMOLH)
    rep = fs.solve()
    res = analyze(fs, rep, TEAConfig(product_component="n-heptane",
                                     product_min_fraction=0.8,
                                     prices_per_kg={"n-heptane": 0.8,
                                                    "n-pentane": 0.6}))
    by_id = {s.unit_id: s for s in res.sizes}
    assert set(by_id) == {"RA", "RA.trays", "RA.reboiler"}
    assert by_id["RA.trays"].quantity == math.ceil(8 / STRIPPER_TRAY_EFFICIENCY)
    assert by_id["RA.reboiler"].utility is not None
    costs = {c.unit_id: c for c in res.costs}
    assert costs["RA"].bare_module > 0
