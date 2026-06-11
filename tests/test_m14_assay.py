"""M14 tests: petroleum assay characterization (caldyr.assay) and
pseudo-components in the cubic-EOS thermo layer.

References:
* Hameed, K.W., "Chemical Process Simulations using Aspen HYSYS", Wiley 2025,
  ch. 10 ("Refinery Process"): the sec. 10.2.1 worked crude assay (figs.
  10.22-10.31) — TBP / molecular-weight / API-gravity curves, light ends,
  bulk MW 300 & 48.75 API, and the 28+8+2 = 38-hypocomponent cut scheme.
  (The same assay as the classic Aspen HYSYS oil-characterization tutorial.)
* Riazi, M.R., "Characterization and Properties of Petroleum Fractions",
  ASTM MNL50 (2005): D86->TBP interconversion (Table 3.5, the Riazi-Daubert /
  API method), MW / Tc / Pc / omega / Cp correlations as cited in
  caldyr/assay.py.
* Pure-component anchor data (n-decane, n-octane): Poling, Prausnitz &
  O'Connell, "The Properties of Gases and Liquids" 5e, App. A (also the
  `chemicals` databank): used to check each correlation's accuracy class.
"""
import warnings

import pytest

from caldyr.assay import (
    acentric_factor,
    api_to_sg,
    astm_d86_to_tbp,
    characterize_assay,
    kesler_lee_cp_ig,
    kesler_lee_pc,
    kesler_lee_tc,
    riazi_daubert_mw,
    riazi_daubert_mw_1980,
    riazi_daubert_pc,
    riazi_daubert_tc,
    sg_to_api,
    watson_k,
)
from caldyr.core import (
    Component,
    Flowsheet,
    is_pseudo_component,
    pseudo_constants,
    resolve_component,
)
from caldyr.core.components_db import molar_mass
from caldyr.io import from_dict, to_dict
from caldyr.thermo import make_package
from caldyr.unitops import FlashDrum, Heater, ShortcutColumn


def degF(f: float) -> float:
    return (f - 32.0) / 1.8 + 273.15


# -- the book's worked assay (Hameed 2025 sec. 10.2.1) -------------------------
# TBP curve, LV% vs F (fig. 10.25); MW curve (fig. 10.26); API gravity curve
# (fig. 10.27); light ends in LV% (fig. 10.22, total 1.13%); bulk MW 300 and
# standard density 48.75 API (fig. 10.23).
BOOK_TBP = [(0.0, 80.0), (10.0, 255.0), (20.0, 349.0), (30.0, 430.0),
            (40.0, 527.0), (50.0, 635.0), (60.0, 751.0), (70.0, 915.0),
            (80.0, 1095.0), (90.0, 1277.0), (98.0, 1410.0)]      # (LV%, F)
BOOK_MW = [(0.0, 68.0), (10.0, 119.0), (20.0, 150.0), (30.0, 182.0),
           (40.0, 225.0), (50.0, 282.0), (60.0, 350.0), (70.0, 456.0),
           (80.0, 585.0), (90.0, 713.0), (98.0, 838.0)]          # (LV%, g/mol)
BOOK_API = [(13.0, 63.28), (33.0, 54.86), (57.0, 45.91),
            (74.0, 38.21), (91.0, 26.01)]                        # (LV%, API)
BOOK_LIGHT_ENDS = {"isobutane": 0.19, "n-butane": 0.11,
                   "isopentane": 0.37, "n-pentane": 0.46}        # LV%
BOOK_BULK_MW = 300.0          # g/mol (user bulk input the book enters)
BOOK_BULK_API = 48.75
# Cut scheme (fig. 10.31 / step 20): start 100 F; to 800 F in 28 cuts of 25 F,
# 800-1200 F in 8 cuts of 50 F, 1200-1400 F in 2 cuts of 100 F -> 38 hypos.
BOOK_CUT_POINTS_F = ([100.0 + 25.0 * i for i in range(29)]
                     + [800.0 + 50.0 * i for i in range(1, 9)]
                     + [1200.0 + 100.0 * i for i in range(1, 3)])


def book_assay(**overrides):
    kwargs = dict(
        kind="TBP",
        api_gravity=BOOK_BULK_API,
        sg_curve=[(v, api_to_sg(a)) for v, a in BOOK_API],
        mw_curve=BOOK_MW,
        cut_points=[degF(f) for f in BOOK_CUT_POINTS_F],
        light_ends=dict(BOOK_LIGHT_ENDS),
    )
    kwargs.update(overrides)
    return characterize_assay([(v, degF(f)) for v, f in BOOK_TBP], **kwargs)


# =============================================================================
# 1. correlation unit checks against literature anchors
# =============================================================================
# n-decane: Tb = 447.30 K, SG = 0.7342, MW = 142.28 g/mol, Tc = 617.7 K,
# Pc = 21.1 bar, omega = 0.4923 (Poling 5e App. A; Riazi 2005 Table 2.1 uses
# the same values when demonstrating these very correlations).
N_DECANE = dict(Tb=447.30, SG=0.7342, MW=142.28, Tc=617.7, Pc=21.1e5, omega=0.4923)
# n-octane: Tb = 398.8 K, SG = 0.7070, MW = 114.23, Tc = 568.7 K, Pc = 24.9 bar
N_OCTANE = dict(Tb=398.8, SG=0.7070, MW=114.23, Tc=568.7, Pc=24.9e5, omega=0.3996)


@pytest.mark.parametrize("ref", [N_DECANE, N_OCTANE], ids=["n-decane", "n-octane"])
def test_riazi_daubert_mw_vs_literature(ref):
    # The extended (1987 / API TDB 2B2.1) form is the working correlation
    # (~2-4% AAD in range); the original 1980 form is looser (~5-7% here).
    assert riazi_daubert_mw(ref["Tb"], ref["SG"]) * 1000 == pytest.approx(ref["MW"], rel=0.05)
    assert riazi_daubert_mw_1980(ref["Tb"], ref["SG"]) * 1000 == pytest.approx(ref["MW"], rel=0.08)


@pytest.mark.parametrize("ref", [N_DECANE, N_OCTANE], ids=["n-decane", "n-octane"])
def test_critical_constants_vs_literature(ref):
    # Tc: both Kesler-Lee and Riazi-Daubert land within ~1%.
    assert kesler_lee_tc(ref["Tb"], ref["SG"]) == pytest.approx(ref["Tc"], rel=0.01)
    assert riazi_daubert_tc(ref["Tb"], ref["SG"]) == pytest.approx(ref["Tc"], rel=0.01)
    # Pc is the soft constant of characterization; ~5% is the accuracy class.
    assert kesler_lee_pc(ref["Tb"], ref["SG"]) == pytest.approx(ref["Pc"], rel=0.05)
    assert riazi_daubert_pc(ref["Tb"], ref["SG"]) == pytest.approx(ref["Pc"], rel=0.06)


@pytest.mark.parametrize("ref", [N_DECANE, N_OCTANE], ids=["n-decane", "n-octane"])
def test_acentric_factor_vs_literature(ref):
    # Lee-Kesler omega evaluated at the *experimental* Tc/Pc reproduces the
    # databank omega to <1% — the correlation itself is that good; the larger
    # errors in practice come from feeding it estimated Tc/Pc.
    kw = watson_k(ref["Tb"], ref["SG"])
    omega = acentric_factor(ref["Tb"], ref["Tc"], ref["Pc"], kw)
    assert omega == pytest.approx(ref["omega"], rel=0.01)


def test_kesler_lee_cp_vs_databank_n_decane():
    # Ideal-gas Cp of n-decane from the Kesler-Lee petroleum correlation vs
    # the thermo databank correlation: within ~3% across 300-500 K.
    from thermo import ChemicalConstantsPackage

    _, props = ChemicalConstantsPackage.from_IDs(["decane"])
    kw = watson_k(N_DECANE["Tb"], N_DECANE["SG"])
    a0, a1, a2 = kesler_lee_cp_ig(N_DECANE["MW"] / 1000.0, kw)
    for T in (300.0, 400.0, 500.0):
        cp_corr = a0 + a1 * T + a2 * T * T
        cp_ref = float(props.HeatCapacityGases[0](T))
        assert cp_corr == pytest.approx(cp_ref, rel=0.03), f"T={T}"


def test_watson_k_and_api_round_trip():
    assert api_to_sg(sg_to_api(0.85)) == pytest.approx(0.85, rel=1e-12)
    # n-decane is a textbook Kw ~ 12.7 paraffin
    assert watson_k(N_DECANE["Tb"], N_DECANE["SG"]) == pytest.approx(12.67, abs=0.05)


# =============================================================================
# 2. ASTM D86 -> TBP conversion vs a published worked pair
# =============================================================================
def test_d86_to_tbp_published_naphtha():
    """Riazi-Daubert (API) D86->TBP on a naphtha fraction. Coefficients:
    Riazi 2005 Table 3.5 (Riazi & Daubert 1986 / API TDB 3A1.1). Reference
    output: the worked conversion of this exact fraction in Humadi,
    'Properties of Petroleum & Natural Gas' lecture 11 (Tikrit Univ. College
    of Petroleum Processes Eng.), computed with the same published table —
    agreement to <=0.8 K is their spreadsheet rounding."""
    d86_C = [(0.0, 138.8), (10.0, 149.6), (30.0, 158.8), (50.0, 165.8),
             (70.0, 169.9), (90.0, 178.1), (95.0, 180.4)]
    tbp_pub_K = [382.40, 405.36, 425.86, 440.38, 448.59, 460.25, 460.80]
    out = astm_d86_to_tbp([(v, t + 273.15) for v, t in d86_C])
    assert [v for v, _ in out] == [0.0, 10.0, 30.0, 50.0, 70.0, 90.0, 95.0]
    for (_, t_calc), t_ref in zip(out, tbp_pub_K):
        assert t_calc == pytest.approx(t_ref, abs=0.8)


def test_d86_to_tbp_shape():
    # TBP is a sharper separation than D86: front of the curve drops, back
    # rises (Riazi 2005 sec. 3.1.2; the lecture's fig. on p. 8 shows exactly
    # this crossover near the 50% point).
    d86 = [(0.0, 412.0), (10.0, 423.0), (30.0, 432.0), (50.0, 439.0),
           (70.0, 443.0), (90.0, 451.0), (95.0, 454.0)]
    tbp = dict(astm_d86_to_tbp(d86))
    d86d = dict(d86)
    assert tbp[0.0] < d86d[0.0]
    assert tbp[10.0] < d86d[10.0]
    assert tbp[90.0] > d86d[90.0]


def test_characterize_from_d86_curve():
    d86 = [(0.0, 412.0), (10.0, 423.0), (30.0, 432.0), (50.0, 439.0),
           (70.0, 443.0), (90.0, 451.0), (95.0, 454.0)]
    res = characterize_assay(d86, kind="ASTM_D86", api_gravity=55.0, n_cuts=10)
    tbp = dict(astm_d86_to_tbp(d86))
    assert res.kind == "ASTM_D86"
    assert len(res.cuts) == 10
    # cuts live inside the converted-TBP span, not the raw D86 span
    assert res.cuts[0].Tb >= tbp[0.0] - 1e-6
    assert res.cuts[-1].Tb <= tbp[95.0] + 1e-6
    # input curve preserved for plotting alongside the working TBP curve
    assert res.curves["input_T_K"] == [t for _, t in sorted(d86)]


# =============================================================================
# 3. the book's worked assay (Hameed 2025 sec. 10.2.1)
# =============================================================================
def test_book_assay_cut_count_matches_hysys():
    res = book_assay()
    # (800-IBP)/25 -> 28, +8, +2 = 38 hypocomponents (book, step 20-2b)
    assert len(res.cuts) == 38


def test_book_assay_cut_nbps():
    res = book_assay()
    # NBPs are mid-volume boiling points: monotone, inside their cut range,
    # and the cut ranges tile the book's boundaries.
    for cut in res.cuts:
        assert cut.T_lo - 1e-9 <= cut.Tb <= cut.T_hi + 1e-9
    nbps = [c.Tb for c in res.cuts]
    assert nbps == sorted(nbps)
    # first range ends at 125 F, last reaches the curve's 1410 F end point
    assert res.cuts[0].T_hi == pytest.approx(degF(125.0), abs=1e-6)
    assert res.cuts[-1].T_hi == pytest.approx(degF(1410.0), abs=0.5)
    # interior boundaries are exactly the book's 25/50/100 F grid
    assert res.cuts[10].T_lo == pytest.approx(degF(100.0 + 25.0 * 10), abs=1e-9)
    assert res.cuts[28].T_lo == pytest.approx(degF(800.0), abs=1e-9)
    assert res.cuts[36].T_lo == pytest.approx(degF(1200.0), abs=1e-9)


def test_book_assay_cut_mws_vs_book_curve():
    """With the book's MW curve supplied (the HYSYS 'Molecular Wt. Curve:
    Dependent' input), cut MWs interpolate it directly. Check the table
    values are hit at their own mid-volume locations."""
    res = book_assay()
    from scipy.interpolate import PchipInterpolator

    mw_curve = PchipInterpolator([v for v, _ in BOOK_MW], [m for _, m in BOOK_MW])
    for cut in res.cuts:
        assert cut.MW * 1000.0 == pytest.approx(
            float(mw_curve(min(cut.vol_mid_pct, 98.0))), rel=0.01), cut.id


def test_book_assay_correlation_mw_vs_book_curve():
    """Without an MW curve, MW comes from Riazi-Daubert(Tb, SG). Against the
    book's measured MW curve the correlation achieves, over the 33 cuts inside
    its stated validity (Tb <= 850 K — through the vacuum-gas-oil range):
    mean |d| = 3.7%, max = 16% (at the very edge of validity; mid-barrel cuts
    sit at 2-5%) — inside the 5-15% class agreement expected between
    different simulators' characterization sets. The heavy-resid
    extrapolation beyond 850 K drifts to +47% and is excluded as out of
    correlation range (documented in caldyr.assay; supply mw_curve for
    resid-accurate work)."""
    res = book_assay(mw_curve=None)
    from scipy.interpolate import PchipInterpolator

    mw_curve = PchipInterpolator([v for v, _ in BOOK_MW], [m for _, m in BOOK_MW])
    deltas = []
    for cut in res.cuts:
        if cut.Tb > 850.0:
            continue
        book = float(mw_curve(min(cut.vol_mid_pct, 98.0)))
        deltas.append(abs(cut.MW * 1000.0 - book) / book)
    assert len(deltas) >= 25
    assert max(deltas) < 0.17
    assert sum(deltas) / len(deltas) < 0.05


def test_book_assay_cut_sgs_follow_api_curve():
    res = book_assay()
    # SG interpolates the book's API curve (converted): increasing with Tb,
    # spanning the curve's API 63.28 -> 26.01 (SG 0.726 -> 0.898).
    sgs = [c.SG for c in res.cuts]
    assert sgs == sorted(sgs)
    assert min(sgs) == pytest.approx(api_to_sg(63.28), abs=0.002)
    assert max(sgs) == pytest.approx(api_to_sg(26.01), abs=0.002)


def test_book_assay_bulk_properties():
    res = book_assay()
    # Volume-blended gravity of the characterized curves vs the book's bulk
    # input 48.75 API: Delta ~ 1.6 API (the curves and the bulk value are not
    # perfectly consistent in the book's data; HYSYS resolves this by forcing
    # the bulk to win — "matching Whole Crude properties takes precedence",
    # sec. 10.1.1.3. We report what the curves imply.)
    assert res.bulk["API"] == pytest.approx(BOOK_BULK_API, abs=2.0)
    # Blend MW implied by the book's own MW + density curves is ~234 g/mol,
    # 22% below the book's bulk input of 300 — same data inconsistency, kept
    # honest here rather than silently rescaled.
    assert res.bulk["MW"] * 1000.0 == pytest.approx(234.0, rel=0.05)
    assert res.bulk["MW"] * 1000.0 == pytest.approx(BOOK_BULK_MW, rel=0.25)
    # Whole-crude Watson K ~ 11.5-13.5 for a paraffinic-intermediate crude
    assert 11.0 < res.bulk["watson_k"] < 13.5
    # the constant-Watson-K path (no density curve) reproduces the bulk API
    res_kw = book_assay(sg_curve=None)
    assert res_kw.bulk["API"] == pytest.approx(BOOK_BULK_API, abs=0.6)


def test_book_assay_fractions_close():
    res = book_assay()
    vol = sum(c.vol_frac for c in res.cuts) + sum(
        f["vol_frac"] for f in res.light_end_fractions.values())
    mass = sum(c.mass_frac for c in res.cuts) + sum(
        f["mass_frac"] for f in res.light_end_fractions.values())
    mol = sum(c.mole_frac for c in res.cuts) + sum(
        f["mole_frac"] for f in res.light_end_fractions.values())
    assert vol == pytest.approx(1.0, abs=1e-12)
    assert mass == pytest.approx(1.0, abs=1e-12)
    assert mol == pytest.approx(1.0, abs=1e-12)
    # light ends were entered as 1.13 LV% with the cut region covering the
    # remaining 1.13-98 LV% of the curve -> 1.13/98 of the accounted volume
    le_vol = sum(f["vol_frac"] for f in res.light_end_fractions.values())
    assert le_vol == pytest.approx(1.13 / 98.0, rel=1e-3)


def test_book_assay_curves_are_plot_ready():
    res = book_assay()
    assert len(res.curves["vol_pct"]) == len(res.curves["tbp_K"]) == 101
    assert res.curves["tbp_K"][0] == pytest.approx(degF(80.0), abs=15.0)
    assert res.curves["tbp_K"][-1] == pytest.approx(degF(1410.0), abs=0.5)
    assert res.curves["input_vol_pct"] == [v for v, _ in BOOK_TBP]


# =============================================================================
# 4. pseudo-components in the thermo layer
# =============================================================================
def test_pseudo_component_registers_and_resolves():
    comp = Component(id="NBP_TEST_A", pseudo={
        "MW": 0.120, "Tb": 423.15, "SG": 0.78, "Tc": 620.0, "Pc": 2.5e6,
        "omega": 0.35, "Cp_ig": [10.0, 0.45, -1.5e-4]})
    assert is_pseudo_component("NBP_TEST_A")
    assert pseudo_constants("NBP_TEST_A")["Tb"] == 423.15
    assert molar_mass("NBP_TEST_A") == pytest.approx(0.120)
    resolved = resolve_component("NBP_TEST_A")
    assert resolved.pseudo is not None and resolved.pseudo["MW"] == 0.120
    assert comp.pseudo == resolved.pseudo


def test_pseudo_component_missing_constant_raises():
    with pytest.raises(ValueError, match="missing required constant"):
        Component(id="NBP_BAD", pseudo={"MW": 0.1, "Tb": 400.0})


def test_pseudo_package_flash_no_formation_warnings():
    """A PR package mixing a databank species with pseudo cuts: PT/PH flashes
    work, are consistent, and raise NO 'no formation enthalpy' warnings —
    pseudo Hf/Gf are 0 by design (reactions on pseudos unsupported), not
    missing data."""
    Component(id="NBP_TEST_L", pseudo={
        "MW": 0.100, "Tb": 380.0, "SG": 0.73, "Tc": 585.0, "Pc": 2.8e6,
        "omega": 0.30, "Cp_ig": [5.0, 0.40, -1.2e-4]})
    Component(id="NBP_TEST_H", pseudo={
        "MW": 0.250, "Tb": 600.0, "SG": 0.85, "Tc": 810.0, "Pc": 1.4e6,
        "omega": 0.70, "Cp_ig": [12.0, 0.95, -3.0e-4]})
    with warnings.catch_warnings():
        warnings.filterwarnings("error", message=".*formation.*")
        pp = make_package("thermo:PR", ["propane", "NBP_TEST_L", "NBP_TEST_H"])
    z = {"propane": 0.2, "NBP_TEST_L": 0.5, "NBP_TEST_H": 0.3}
    res = pp.flash_pt(450.0, 3e5, z)
    assert res.phase == "VLE"
    assert 0.0 < res.vapor_fraction < 1.0
    # lights concentrate in the vapor, heavies in the liquid
    assert res.y["propane"] > res.x["propane"]
    assert res.y["NBP_TEST_H"] < res.x["NBP_TEST_H"]
    # PH flash round-trips the PT enthalpy
    back = pp.flash_ph(3e5, res.H, z)
    assert back.T == pytest.approx(450.0, abs=1e-5)
    # enthalpy balance basis: pure-pseudo enthalpy is finite and Cp-driven
    h1 = pp.enthalpy(400.0, 3e5, {"NBP_TEST_H": 1.0})
    h2 = pp.enthalpy(410.0, 3e5, {"NBP_TEST_H": 1.0})
    assert h2 > h1


def test_pseudo_rejected_by_nrtl_and_coolprop():
    Component(id="NBP_TEST_N", pseudo={
        "MW": 0.150, "Tb": 450.0, "SG": 0.80, "Tc": 650.0, "Pc": 2.2e6,
        "omega": 0.45, "Cp_ig": [8.0, 0.6, -2e-4]})
    with pytest.raises(ValueError, match="pseudo"):
        make_package("thermo:NRTL", ["water", "NBP_TEST_N"])
    with pytest.raises(ValueError, match="pseudo"):
        make_package("coolprop:Water", ["NBP_TEST_N"])


def test_recharacterized_pseudo_does_not_reuse_stale_flasher():
    """Same id, new constants -> the flasher cache must rebuild (keyed on the
    constants, not just the id)."""
    base = {"MW": 0.140, "SG": 0.78, "Tc": 640.0, "Pc": 2.0e6,
            "omega": 0.42, "Cp_ig": [8.0, 0.55, -1.8e-4]}
    Component(id="NBP_TEST_R", pseudo={**base, "Tb": 440.0})
    pp1 = make_package("thermo:PR", ["propane", "NBP_TEST_R"])
    assert pp1._flasher.constants.Tcs[1] == pytest.approx(640.0)
    _, dew1 = pp1.bubble_dew(1e5, {"propane": 0.5, "NBP_TEST_R": 0.5})
    Component(id="NBP_TEST_R", pseudo={**base, "Tb": 460.0, "Tc": 680.0})
    pp2 = make_package("thermo:PR", ["propane", "NBP_TEST_R"])
    assert pp2._flasher.constants.Tcs[1] == pytest.approx(680.0)
    # the dew point is set by the heavy component: a heavier re-characterized
    # cut must shift it up (it cannot if a stale flasher were reused)
    _, dew2 = pp2.bubble_dew(1e5, {"propane": 0.5, "NBP_TEST_R": 0.5})
    assert dew2 > dew1 + 5.0


# =============================================================================
# 5. .flow round-trip
# =============================================================================
def test_flow_load_with_pseudo_component_entries():
    """A `.flow` document whose component entries carry `pseudo` loads through
    Component(**c) (from_dict's existing path), registers the constants, and
    solves under PR."""
    doc = {
        "schema": "caldyr.flow/1",
        "components": [
            {"id": "propane", "name": "propane"},
            {"id": "NBP_TEST_F", "name": "NBP_TEST_F", "pseudo": {
                "MW": 0.180, "Tb": 500.0, "SG": 0.82, "Tc": 700.0,
                "Pc": 1.8e6, "omega": 0.5, "Cp_ig": [9.0, 0.7, -2.4e-4]}},
        ],
        "property_package": "thermo:PR",
        "units": [{"id": "H", "type": "Heater", "params": {"T_out": 420.0}}],
        "streams": [
            {"id": "f", "from": None, "to": "H:in1",
             "spec": {"T": 350.0, "P": 4e5, "molar_flow": 10.0,
                      "z": {"propane": 0.4, "NBP_TEST_F": 0.6}}},
            {"id": "p", "from": "H:out", "to": None},
        ],
    }
    fs = from_dict(doc)
    assert fs.components[1].pseudo is not None
    rep = fs.solve()
    assert rep.converged
    assert fs.streams["p"].T == pytest.approx(420.0)


def test_flow_roundtrip_in_session():
    """to_dict -> from_dict -> solve works within a session: the registry
    still holds the constants even though to_dict (in the off-limits io/
    module) does not yet serialize the `pseudo` field — the documented
    1-line io gap (see the M14 report)."""
    Component(id="NBP_TEST_S", pseudo={
        "MW": 0.160, "Tb": 470.0, "SG": 0.80, "Tc": 670.0, "Pc": 1.9e6,
        "omega": 0.46, "Cp_ig": [8.5, 0.62, -2.1e-4]})
    fs = Flowsheet(components=[resolve_component("propane"),
                               resolve_component("NBP_TEST_S")],
                   property_package="thermo:PR")
    fs.add(Heater("H", {"T_out": 410.0}))
    fs.connect("f", None, "H:in1")
    fs.connect("p", "H:out", None)
    fs.feed("f", "H:in1", T=350.0, P=4e5, molar_flow=5.0,
            z={"propane": 0.5, "NBP_TEST_S": 0.5})
    assert fs.solve().converged
    fs2 = from_dict(to_dict(fs))
    assert fs2.solve().converged
    assert fs2.streams["p"].T == pytest.approx(410.0)


# =============================================================================
# 6. end-to-end: crude characterized -> heated -> flashed -> column
# =============================================================================
@pytest.fixture(scope="module")
def crude_result():
    # 12 cuts keeps the 16-component PR flashes fast while spanning the barrel
    return book_assay(cut_points=None, n_cuts=12)


def test_crude_heat_and_flash_sensible_split(crude_result):
    """The book's flowsheet front end (sec. 10.2.1 steps 27-29): crude at
    450 F / 75 psia, fired to 650 F (10 psi drop), flashed near tower-feed
    pressure. The split must be physically sensible: substantial vaporization,
    every cut's K-value monotone decreasing with NBP, light ends overhead,
    resid in the liquid, balances closed."""
    res = crude_result
    fs = Flowsheet(components=res.components(), property_package="thermo:PR")
    fs.add(Heater("FURNACE", {"T_out": degF(650.0), "dP": 68947.6}))   # 10 psi
    fs.add(FlashDrum("FLASH", {"P": 448159.0}))                        # 65 psia
    fs.connect("crude", None, "FURNACE:in1")
    fs.connect("hot", "FURNACE:out", "FLASH:in1")
    fs.connect("vap", "FLASH:vapor", None)
    fs.connect("liq", "FLASH:liquid", None)
    fs.feed("crude", "FURNACE:in1", T=degF(450.0), P=517107.0,         # 75 psia
            molar_flow=100.0, z=res.mole_fractions())
    rep = fs.solve()
    assert rep.converged

    vap, liq = fs.streams["vap"], fs.streams["liq"]
    assert vap.molar_flow + liq.molar_flow == pytest.approx(100.0, rel=1e-8)
    # component balances
    for cid in fs.component_ids:
        in_i = 100.0 * fs.streams["crude"].z[cid]
        out_i = vap.molar_flow * vap.z.get(cid, 0.0) + liq.molar_flow * liq.z.get(cid, 0.0)
        assert out_i == pytest.approx(in_i, rel=1e-6, abs=1e-10)
    # substantial but partial vaporization at 650 F / 65 psia
    assert 0.3 < vap.molar_flow / 100.0 < 0.95
    # K-values monotone in boiling point (lighter pseudos to vapor)
    ks = [(vap.z.get(c.id, 0.0) + 1e-300) / (liq.z.get(c.id, 0.0) + 1e-300)
          for c in res.cuts]
    assert all(k1 >= k2 for k1, k2 in zip(ks, ks[1:]))
    assert ks[0] > 1.0 > ks[-1]
    # light ends go overhead almost completely
    for le in res.light_end_fractions:
        z_f = fs.streams["crude"].z[le]
        rec = vap.molar_flow * vap.z.get(le, 0.0) / (100.0 * z_f)
        assert rec > 0.90, le      # iC5 lands at ~0.94 at 65 psia / 650 F
    # the heaviest cut stays >99% in the liquid
    heavy = res.cuts[-1].id
    rec_heavy = liq.molar_flow * liq.z[heavy] / (100.0 * fs.streams["crude"].z[heavy])
    assert rec_heavy > 0.99


def test_crude_shortcut_column_naphtha_split(crude_result):
    """A FUG shortcut column making a naphtha-type cut on the characterized
    crude (LK/HK = adjacent ~190/250 C cuts). Converges and meets its
    recovery specs — pseudo-components behave as ordinary EOS species."""
    res = crude_result
    ids = [c.id for c in res.cuts]
    lk, hk = ids[2], ids[3]
    fs = Flowsheet(components=res.components(), property_package="thermo:PR")
    fs.add(ShortcutColumn("COL", {
        "light_key": lk, "heavy_key": hk,
        "recovery_light": 0.95, "recovery_heavy": 0.95,
        "P": 2e5, "partial_condenser": True}))
    fs.connect("crude", None, "COL:in1")
    fs.connect("dist", "COL:distillate", None)
    fs.connect("btms", "COL:bottoms", None)
    fs.feed("crude", "COL:in1", T=degF(450.0), P=3e5, molar_flow=100.0,
            z=res.mole_fractions())
    rep = fs.solve()
    assert rep.converged
    dist, btms = fs.streams["dist"], fs.streams["btms"]
    z_feed = fs.streams["crude"].z
    rec_lk = dist.molar_flow * dist.z.get(lk, 0.0) / (100.0 * z_feed[lk])
    rec_hk = btms.molar_flow * btms.z.get(hk, 0.0) / (100.0 * z_feed[hk])
    assert rec_lk == pytest.approx(0.95, abs=0.01)
    assert rec_hk == pytest.approx(0.95, abs=0.01)
    # non-keys distribute monotonically: everything lighter than the LK goes
    # overhead, everything heavier than the HK goes down
    for cid in ids[:2]:
        rec = dist.molar_flow * dist.z.get(cid, 0.0) / (100.0 * z_feed[cid])
        assert rec > 0.99, cid
    for cid in ids[5:]:
        rec = btms.molar_flow * btms.z.get(cid, 0.0) / (100.0 * z_feed[cid])
        assert rec > 0.99, cid
    design = fs.units["COL"].design
    assert design["N"] > design["N_min"] > 1.0
    assert design["R"] > design["R_min"] > 0.0


# =============================================================================
# 7. input validation
# =============================================================================
def test_curve_validation_errors():
    with pytest.raises(ValueError, match="at least 5 points"):
        characterize_assay([(0, 300.0), (50, 400.0), (100, 500.0)], api_gravity=40)
    pts = [(0, 300.0), (10, 350.0), (30, 340.0), (50, 420.0), (90, 500.0)]
    with pytest.raises(ValueError, match="strictly increasing"):
        characterize_assay(pts, api_gravity=40)
    with pytest.raises(ValueError, match="kind"):
        characterize_assay([(0, 300.0), (10, 350.0), (30, 380.0), (50, 420.0),
                            (90, 500.0)], kind="D1160", api_gravity=40)


def test_light_ends_validation():
    pts = [(0, 300.0), (10, 350.0), (30, 380.0), (50, 420.0), (90, 500.0)]
    with pytest.raises(ValueError, match="light ends"):
        characterize_assay(pts, api_gravity=40, light_ends={"propane": 100.0})
    with pytest.raises(ValueError, match="light end"):
        characterize_assay(pts, api_gravity=40, light_ends={"propane": -1.0})
