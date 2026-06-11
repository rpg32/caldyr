"""M13 tests: the liquid-liquid ExtractionColumn (counter-current LLE cascade).

Validation strategy (sources cited per test):

1. **The book's extractor** (Hameed, *Chemical Process Simulations using
   Aspen Hysys*, Wiley 2025, sec. 9.4): 100 kg/h of 50/50 wt% acetone/water
   fed at the top, 150 kg/h of 3-methylhexane at the bottom, 8 trays, 1 bar,
   25 C. This polar system is reproduced *structurally, not quantitatively*,
   and the tests document exactly why: the book runs HYSYS-PRSV and our LLE
   flash is stock-`thermo` PR — for acetone between water and a C7 alkane a
   cubic EOS expels essentially *all* non-water species from the aqueous
   phase (K_acetone ~ 1e5 organic/aqueous, where NRTL/experiment put the
   distribution coefficient at O(1), e.g. Seader, Henley & Roper,
   *Separation Process Principles* 3e, ch. 8 acetone/water/MIBK type
   systems). The book's own headline result — the bottom water stream leaves
   with x(water) = 1.0 (sec. 9.4.2 step 9) — *is* reproduced exactly,
   because HYSYS-PRSV behaves the same way there; the acetone recovery and
   the extract loading are asserted as structure (solute transfers to the
   solvent phase; phases stay ordered top/bottom), not as NRTL-quantitative
   compositions.
2. **Kremser closed form on a PR-friendly system** (Seader 3e ch. 8,
   eq. (8-13)/(5-48): fraction extracted = (E^(N+1) - E)/(E^(N+1) - 1) with
   extraction factor E = K S/F): methanol distributed between water and
   toluene — a hydrocarbon/water split PR represents robustly (finite
   K_methanol ~ O(10) organic/aqueous). The cascade must reproduce the
   closed form evaluated on its *own* converged stage K-values (Edmister
   geometric-mean E of the end stages, same construction as the absorber's
   Kremser test); achieved deltas are 1.7-2.3% on the E < 1 cases and
   < 0.1% on the E > 1 cases, over N = 3..6 and E ~ 0.14..6 (Kremser itself
   assumes constant E, so a few percent is the honest bar).
3. Structure the theory demands (Seader 3e ch. 8): recovery rises with
   stages and with solvent rate; the ascending stream's solute loading grows
   monotonically bottom -> top.
4. Conservation to machine precision, exact energy closure through the duty
   port, `.flow` round-trip, full stage profiles on ``unit.design``,
   economics sizing (tower + trays, no condenser/reboiler), typed errors
   (missing specs, NRTL package, vaporizing stage, fully miscible feeds).
"""
import math

import pytest

from caldyr.core import Component, Flowsheet
from caldyr.economics.sizing import SizingOptions, size_flowsheet
from caldyr.io import from_dict, to_dict
from caldyr.thermo import make_package
from caldyr.unitops import ExtractionColumn, ExtractionColumnError

P_ATM = 101325.0
KMOLH = 1000.0 / 3600.0          # kmol/h -> mol/s
MW = {"acetone": 58.0791, "water": 18.01528, "3-methylhexane": 100.2019}


# -- builders -------------------------------------------------------------------
def book_extractor(n_stages=8, solvent_kgh=150.0) -> Flowsheet:
    """The book's acetone extractor (Hameed 2025 sec. 9.4.2): 100 kg/h of
    50/50 wt% acetone/water (the heavy liquid, fed at the top) against
    3-methylhexane (the light solvent, fed at the bottom), 8 trays, 1 bar,
    top temperature 25 C."""
    f_acetone = 50.0 / MW["acetone"] * KMOLH          # mol/s
    f_water = 50.0 / MW["water"] * KMOLH
    fs = Flowsheet(components=[Component("acetone"), Component("water"),
                               Component("3-methylhexane")],
                   property_package="thermo:PR")
    fs.add(ExtractionColumn("EXT", {"n_stages": n_stages, "T": 298.15,
                                    "P": 1e5}))
    fs.feed("FEED", "EXT:feed_in", T=298.15, P=1e5,
            molar_flow=f_acetone + f_water,
            z={"acetone": f_acetone / (f_acetone + f_water),
               "water": f_water / (f_acetone + f_water)})
    fs.feed("SOLV", "EXT:solvent_in", T=298.15, P=1e5,
            molar_flow=solvent_kgh / MW["3-methylhexane"] * KMOLH,
            z={"3-methylhexane": 1.0})
    fs.connect("EXTRACT", "EXT:extract_out", None)
    fs.connect("RAFF", "EXT:raffinate_out", None)
    fs.connect("Q", "EXT:duty", None)
    return fs


def methanol_extractor(n_stages, solvent, x_meoh=0.02, feed=100.0,
                       pkg="thermo:PR") -> Flowsheet:
    """The quantitative (Kremser) system: toluene extracting methanol out of
    water at 25 C, 1 atm. The aqueous feed is the heavy liquid (top); the
    toluene solvent is the light liquid (bottom) and leaves overhead as the
    extract. Under PR the methanol distribution coefficient is finite
    (K ~ 5-20 organic/aqueous depending on loading), so the extraction
    factor E = K S/F sweeps the interesting O(0.1)..O(10) range with
    sensible solvent rates."""
    fs = Flowsheet(components=[Component("methanol"), Component("toluene"),
                               Component("water")],
                   property_package=pkg)
    fs.add(ExtractionColumn("EXT", {"n_stages": n_stages, "T": 298.15,
                                    "P": P_ATM}))
    fs.feed("FEED", "EXT:feed_in", T=298.15, P=P_ATM, molar_flow=feed,
            z={"water": 1.0 - x_meoh, "methanol": x_meoh})
    fs.feed("SOLV", "EXT:solvent_in", T=298.15, P=P_ATM, molar_flow=solvent,
            z={"toluene": 1.0})
    fs.connect("EXTRACT", "EXT:extract_out", None)
    fs.connect("RAFF", "EXT:raffinate_out", None)
    fs.connect("Q", "EXT:duty", None)
    return fs


def kremser_extracted(fs, n) -> float:
    """Kremser prediction (Seader 3e eq. (8-13)) with the effective
    extraction factor E_e = sqrt(E_top * E_bottom) (Edmister geometric mean)
    evaluated from the *converged* stage states: K_j is the actual
    extract/raffinate composition ratio of the two streams leaving stage j
    (they are in equilibrium), and E_j = K_j * (ascending flow)/(descending
    flow) there — so the closed form is compared on identical
    thermodynamics, exactly like the absorber's Kremser test."""
    d = fs.units["EXT"].design
    es = []
    for j in (0, n - 1):
        k = (d["x_extract_profile"][j]["methanol"]
             / d["x_raffinate_profile"][j]["methanol"])
        es.append(k * d["E_profile"][j] / d["R_profile"][j])
    ee = math.sqrt(es[0] * es[1])
    return (ee ** (n + 1) - ee) / (ee ** (n + 1) - 1.0)


def total_in(fs, c) -> float:
    return sum(s.molar_flow * s.z.get(c, 0.0)
               for sid, s in fs.streams.items() if sid in ("FEED", "SOLV"))


# -- 1. the book's extractor (sec. 9.4) -------------------------------------------
def test_book_9_4_acetone_extractor():
    """The book's headline check (sec. 9.4.2 step 9): 'the mole fraction for
    water [in the bottom stream] is 1' — reproduced exactly. The acetone
    goes overhead with the solvent (the Rich-Sol stream feeding the book's
    T-101 solvent-recovery column). Quantitative caveat per the module
    docstring: PR makes the acetone transfer total (recovery 1.0), which
    matches the book's PRSV raffinate but is *not* an NRTL-grade LLE."""
    fs = book_extractor()
    assert fs.solve().converged
    e, r = fs.streams["EXTRACT"], fs.streams["RAFF"]

    # Book: the bottoms heavy-liquid water stream is pure water.
    assert r.z["water"] == pytest.approx(1.0, abs=1e-3)
    assert r.z["acetone"] == pytest.approx(0.0, abs=1e-3)
    assert r.z["3-methylhexane"] == pytest.approx(0.0, abs=1e-3)

    # The solute transferred to the solvent phase: the extract carries
    # essentially all the acetone and all the solvent.
    d = fs.units["EXT"].design
    assert d["recovery"]["acetone"] > 0.99
    assert d["recovery"]["3-methylhexane"] > 0.99
    assert e.z["acetone"] > 0.3                      # acetone-loaded extract
    assert e.z["3-methylhexane"] > 0.5

    # Phases ordered: both products liquid at the operating T, P.
    assert e.phase == "liquid" and r.phase == "liquid"
    assert e.T == pytest.approx(298.15) and r.T == pytest.approx(298.15)

    # Total-flow sanity against the feeds (2.41 / 2.72 kmol/h split).
    assert e.molar_flow + r.molar_flow == pytest.approx(
        fs.streams["FEED"].molar_flow + fs.streams["SOLV"].molar_flow,
        rel=1e-12)


def test_book_orientation_densities_reported():
    """The book's rule (sec. 9.4.2 step 4): the heavier liquid enters at the
    top. The model reports both feed densities as diagnostics (not policed —
    see the module docstring on PR liquid densities)."""
    fs = book_extractor()
    fs.solve()
    d = fs.units["EXT"].design
    assert d["rho_feed_in"] is not None and d["rho_feed_in"] > 0.0
    assert d["rho_solvent_in"] is not None and d["rho_solvent_in"] > 0.0
    # PR's aqueous-feed density does beat the alkane's, even if its absolute
    # water density is poor.
    assert d["rho_feed_in"] > d["rho_solvent_in"]


# -- 2. Kremser closed form (Seader 3e ch. 8) ---------------------------------------
@pytest.mark.parametrize("n,solvent", [(3, 5.0), (6, 3.0), (3, 30.0),
                                       (5, 30.0)])
def test_matches_kremser_dilute_extraction(n, solvent):
    fs = methanol_extractor(n, solvent)
    assert fs.solve().converged
    actual = fs.units["EXT"].design["recovery"]["methanol"]
    pred = kremser_extracted(fs, n)
    # Achieved deltas are 1.7-2.3% on the E < 1 cases and < 0.1% on the
    # E > 1 cases; Kremser assumes a constant extraction factor, so a few
    # percent is the honest tolerance.
    assert actual == pytest.approx(pred, rel=0.05)


# -- 3. structure: stages and solvent rate (Seader 3e ch. 8) -------------------------
def test_recovery_rises_with_stages():
    rec = []
    for n in (1, 3, 5):
        fs = methanol_extractor(n, 30.0)
        assert fs.solve().converged
        rec.append(fs.units["EXT"].design["recovery"]["methanol"])
    assert rec[0] < rec[1] < rec[2]
    assert rec[2] > 0.999                  # E ~ 6: near-total recovery


def test_recovery_rises_with_solvent_rate():
    rec = []
    for s in (5.0, 10.0, 30.0):
        fs = methanol_extractor(3, s)
        assert fs.solve().converged
        rec.append(fs.units["EXT"].design["recovery"]["methanol"])
    assert rec[0] < rec[1] < rec[2]


def test_profiles_published_and_solute_loading_monotone():
    n = 4
    fs = methanol_extractor(n, 10.0)
    assert fs.solve().converged
    d = fs.units["EXT"].design
    for key in ("E_profile", "R_profile", "T_profile", "P_profile",
                "x_extract_profile", "x_raffinate_profile"):
        assert len(d[key]) == n
    assert all(v > 0.0 for v in d["E_profile"])
    assert all(v > 0.0 for v in d["R_profile"])
    # The ascending (solvent) stream loads up monotonically bottom -> top.
    solute_up = [d["E_profile"][j] * d["x_extract_profile"][j]["methanol"]
                 for j in range(n)]
    assert all(solute_up[j] > solute_up[j + 1] for j in range(n - 1))
    # The descending (aqueous) stream is stripped monotonically top -> bottom.
    solute_dn = [d["R_profile"][j] * d["x_raffinate_profile"][j]["methanol"]
                 for j in range(n)]
    assert all(solute_dn[j] > solute_dn[j + 1] for j in range(n - 1))


# -- 4. conservation, energy, round-trip, sizing --------------------------------------
def test_mass_balance_machine_exact():
    fs = book_extractor()
    assert fs.solve().converged
    e, r = fs.streams["EXTRACT"], fs.streams["RAFF"]
    for c in fs.component_ids:
        n_out = e.molar_flow * e.z[c] + r.molar_flow * r.z[c]
        assert n_out == pytest.approx(total_in(fs, c), rel=1e-12, abs=1e-12)


def test_energy_balance_closes_through_duty():
    fs = methanol_extractor(4, 10.0)
    rep = fs.solve()
    assert rep.converged
    e, r = fs.streams["EXTRACT"], fs.streams["RAFF"]
    f, s = fs.streams["FEED"], fs.streams["SOLV"]
    lhs = e.molar_flow * e.H + r.molar_flow * r.H
    rhs = f.molar_flow * f.H + s.molar_flow * s.H + rep.duties["Q"]
    assert lhs == pytest.approx(rhs, rel=1e-12)


def test_flow_round_trip():
    fs = book_extractor()
    fs.solve()
    doc = to_dict(fs)
    assert doc == to_dict(from_dict(doc))
    fs2 = from_dict(doc)
    fs2.solve()
    for sid in ("EXTRACT", "RAFF"):
        assert fs2.streams[sid].molar_flow == pytest.approx(
            fs.streams[sid].molar_flow, rel=1e-9)
        for c in fs.component_ids:
            assert fs2.streams[sid].z[c] == pytest.approx(
                fs.streams[sid].z[c], rel=1e-9, abs=1e-12)


def test_economics_sizer_tower_and_trays():
    fs = book_extractor()
    rep = fs.solve()
    pp = make_package("thermo:PR", fs.component_ids)
    sizes = size_flowsheet(fs, rep, pp, SizingOptions())
    ext = [s for s in sizes if s.unit_id.startswith("EXT")]
    assert {s.equipment_type for s in ext} == {"vessel_vertical", "tray_sieve"}
    tower = next(s for s in ext if s.equipment_type == "vessel_vertical")
    trays = next(s for s in ext if s.equipment_type == "tray_sieve")
    assert math.isfinite(tower.attribute) and tower.attribute > 0.0
    assert tower.diameter_m is not None and tower.diameter_m > 0.0
    # 8 ideal stages at the 0.25 liquid-liquid stage efficiency -> 32 trays.
    assert trays.quantity == 32
    assert math.isfinite(trays.attribute) and trays.attribute > 0.0
    # No condenser/reboiler/utility items for an extractor.
    assert all(s.utility is None for s in ext)


# -- 5. typed errors -------------------------------------------------------------------
def test_missing_temperature_raises():
    fs = book_extractor()
    del fs.units["EXT"].params["T"]
    with pytest.raises(ExtractionColumnError, match=r"params\['T'\]"):
        fs.solve()


def test_missing_or_bad_stage_count_raises():
    fs = book_extractor()
    del fs.units["EXT"].params["n_stages"]
    with pytest.raises(ExtractionColumnError, match="n_stages"):
        fs.solve()
    fs = book_extractor()
    fs.units["EXT"].params["n_stages"] = 0
    with pytest.raises(ExtractionColumnError, match=">= 1"):
        fs.solve()


def test_nrtl_package_raises_clear_error():
    """The NRTL activity package has no three-phase flash; the column must
    say so instead of failing inside the property package."""
    fs = methanol_extractor(3, 10.0, pkg="thermo:NRTL")
    with pytest.raises(ExtractionColumnError, match="three-phase"):
        fs.solve()


def test_vaporizing_stage_raises():
    """At 390 K and 1 bar the acetone/water/alkane stage boils — a
    liquid-liquid extractor must refuse, not silently route a vapor."""
    fs = book_extractor()
    fs.units["EXT"].params["T"] = 390.0
    with pytest.raises(ExtractionColumnError, match="vaporizes"):
        fs.solve()


def test_fully_miscible_feeds_raise():
    """Benzene 'extracted' with toluene: one liquid phase everywhere. The
    column must diagnose miscibility, not report a convergence failure."""
    fs = Flowsheet(components=[Component("benzene"), Component("toluene")],
                   property_package="thermo:PR")
    fs.add(ExtractionColumn("EXT", {"n_stages": 3, "T": 298.15, "P": P_ATM}))
    fs.feed("FEED", "EXT:feed_in", T=298.15, P=P_ATM, molar_flow=10.0,
            z={"benzene": 1.0})
    fs.feed("SOLV", "EXT:solvent_in", T=298.15, P=P_ATM, molar_flow=10.0,
            z={"toluene": 1.0})
    fs.connect("E", "EXT:extract_out", None)
    fs.connect("R", "EXT:raffinate_out", None)
    fs.connect("Q", "EXT:duty", None)
    with pytest.raises(ExtractionColumnError, match="miscible"):
        fs.solve()


def test_missing_inlet_raises():
    pp = make_package("thermo:PR", ["acetone", "water", "3-methylhexane"])
    unit = ExtractionColumn("EXT", {"n_stages": 3, "T": 298.15, "P": 1e5})
    with pytest.raises(ExtractionColumnError, match="feed_in"):
        unit.solve({}, pp)
