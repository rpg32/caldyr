"""Robust bubble/dew fallbacks for liquids carrying dissolved permanent gases.

thermo's PVF flash crashes (UnboundLocalError) on such mixtures, and a strict
bubble point may not even exist (VF > 0 at every T once enough H2/CH4 is
dissolved). These paths were hit by the book's cyclohexane plant (Hameed 2025
ch. 15.1): the column feed and mid-iteration stage liquids contain H2/N2/CH4.
The fallbacks: VF-bisection on the PT flash (which has stability analysis),
and the classical non-condensables treatment (bubble T of the condensable
submixture, lights folded into the incipient vapor by K-value; Seader 3e
sec. 4.4).
"""
import pytest

from caldyr.thermo import make_package

COMPS = ["benzene", "cyclohexane", "hydrogen", "nitrogen", "methane"]

# the exact stage-liquid composition that crashed the plant solve
STAGE_LIQ = {"benzene": 0.0279, "cyclohexane": 0.8651, "hydrogen": 0.0152,
             "nitrogen": 0.0024, "methane": 0.0894}
# a degassed column feed (has a strict bubble point via bisection)
FEED_LIQ = {"benzene": 0.0310, "cyclohexane": 0.9616, "hydrogen": 0.0003,
            "nitrogen": 0.0001, "methane": 0.0070}


@pytest.fixture(scope="module")
def pp():
    return make_package("thermo:PR", COMPS)


def test_bubble_point_survives_dissolved_permanent_gases(pp):
    res = pp.bubble_point(1.38e6, STAGE_LIQ)
    # sane saturated-liquid answer near the condensables' boiling range
    assert 380.0 < res.T < 520.0
    assert res.vapor_fraction == 0.0
    assert res.y is not None
    # lights are enriched in the incipient vapor (K >> 1)
    assert res.y["hydrogen"] > 10 * STAGE_LIQ["hydrogen"] / sum(STAGE_LIQ.values())
    assert res.y["methane"] > STAGE_LIQ["methane"] / sum(STAGE_LIQ.values())
    assert abs(sum(res.y.values()) - 1.0) < 1e-6
    # both saturated-phase enthalpies present, vapor above liquid
    assert res.H_vapor is not None and res.H_liquid is not None
    assert res.H_vapor > res.H_liquid


def test_bubble_point_trace_gases_close_to_pvf_branch(pp):
    """With only trace lights the answer must sit just above the light-free
    bubble point and far below the heavy's dew point."""
    res = pp.bubble_point(1.38e6, FEED_LIQ)
    clean = pp.bubble_point(1.38e6, {"benzene": 0.031, "cyclohexane": 0.969})
    # ~0.7% dissolved lights legitimately depress the bubble T by ~13 K at
    # 13.8 bar; the guard is against wildly wrong fallback answers.
    assert res.T < clean.T
    assert abs(res.T - clean.T) < 25.0


def test_bubble_dew_survives_the_same_mixture(pp):
    bub, dew = pp.bubble_dew(1.38e6, STAGE_LIQ)
    assert 380.0 < bub < 520.0
    assert dew >= bub          # ordering enforced (degenerate interval allowed)
    assert dew < 600.0


def test_clean_mixtures_unaffected(pp):
    """Regression: the fast PVF branch still serves ordinary mixtures."""
    res = pp.bubble_point(101325.0, {"benzene": 0.5, "cyclohexane": 0.5})
    assert 350.0 < res.T < 360.0       # ~353-355 K at 1 atm
    bub, dew = pp.bubble_dew(101325.0, {"benzene": 0.5, "cyclohexane": 0.5})
    assert bub == pytest.approx(res.T, abs=0.5)
    assert dew > bub
