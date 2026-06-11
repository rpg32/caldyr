"""M11 tests: pinch analysis / heat-integration targeting.

Validation reference: the classic four-stream problem of **Kemp, *Pinch
Analysis and Process Integration*, 2e (2007), Table 2.1 / Sec. 2.3** (the
Linnhoff "User Guide" example; also Smith, *Chemical Process Design and
Integration*, 2005, Ch. 16):

    stream  type   Ts (C)  Tt (C)  CP (kW/K)   Q (kW)
    1       cold     20      135      2.0        230
    2       hot     170       60      3.0        330
    3       cold     80      140      4.0        240
    4       hot     150       30      1.5        180

At dT_min = 10 K the published problem-table targets are QH_min = 20 kW,
QC_min = 60 kW, with the pinch at a shifted temperature of 85 C (hot-side
90 C / cold-side 80 C).

Flowsheet-level check: the ammonia synthesis loop (examples/04) heats the
reactor feed (PREHEAT, cold demand) and cools the reactor effluent over an
overlapping temperature range (COOL + the exothermic isothermal reactor, hot
availability), so the heat-recovery potential must be positive and the
minimum hot utility no larger than the heating bought today.
"""
import pytest

from caldyr.analysis import ThermalStream, pinch_analysis, pinch_from_streams
from caldyr.core import Component, Flowsheet
from caldyr.unitops import EquilibriumReactor, FlashDrum, Heater, Mixer, Splitter

C = 273.15   # K offset

# Kemp 2e Table 2.1, in SI (K, W). Hot/cold is inferred from the T direction.
KEMP_STREAMS = [
    {"Tin": 20.0 + C, "Tout": 135.0 + C, "Q": 230e3},    # cold, CP 2 kW/K
    {"Tin": 170.0 + C, "Tout": 60.0 + C, "Q": 330e3},    # hot,  CP 3 kW/K
    {"Tin": 80.0 + C, "Tout": 140.0 + C, "Q": 240e3},    # cold, CP 4 kW/K
    {"Tin": 150.0 + C, "Tout": 30.0 + C, "Q": 180e3},    # hot,  CP 1.5 kW/K
]


def test_kemp_four_stream_targets_match_published_values():
    res = pinch_from_streams(KEMP_STREAMS, dt_min=10.0)
    assert res.qh_min == pytest.approx(20e3, rel=1e-9)        # Kemp 2e: 20 kW
    assert res.qc_min == pytest.approx(60e3, rel=1e-9)        # Kemp 2e: 60 kW
    assert res.pinch_T_shifted == pytest.approx(85.0 + C)     # Kemp 2e: 85 C
    assert res.pinch_T_hot == pytest.approx(90.0 + C)
    assert res.pinch_T_cold == pytest.approx(80.0 + C)
    # Current (un-integrated) utility use and the recovery headline.
    assert res.current_hot_utility == pytest.approx(470e3)    # 230 + 240 kW
    assert res.current_cold_utility == pytest.approx(510e3)   # 330 + 180 kW
    assert res.heat_recovery_potential == pytest.approx(450e3, rel=1e-9)


def test_overall_energy_balance_identity():
    """QC_min - QH_min == (total hot) - (total cold), always."""
    res = pinch_from_streams(KEMP_STREAMS, dt_min=10.0)
    assert res.qc_min - res.qh_min == pytest.approx(510e3 - 470e3, rel=1e-9)


def test_targets_tighten_as_dt_min_shrinks():
    loose = pinch_from_streams(KEMP_STREAMS, dt_min=20.0)
    base = pinch_from_streams(KEMP_STREAMS, dt_min=10.0)
    tight = pinch_from_streams(KEMP_STREAMS, dt_min=0.0)
    assert loose.qh_min >= base.qh_min >= tight.qh_min
    assert loose.qc_min >= base.qc_min >= tight.qc_min


def test_composite_curves_are_plot_ready():
    res = pinch_from_streams(KEMP_STREAMS, dt_min=10.0)
    # Hot composite spans the total hot duty; cold curve is offset by QC_min
    # so the curves align at dT_min.
    assert res.hot_composite[0][0] == pytest.approx(0.0)
    assert res.hot_composite[-1][0] == pytest.approx(510e3)
    assert res.cold_composite[0][0] == pytest.approx(res.qc_min)
    assert res.cold_composite[-1][0] == pytest.approx(res.qc_min + 470e3)
    for curve in (res.hot_composite, res.cold_composite):
        temps = [t for _, t in curve]
        enth = [h for h, _ in curve]
        assert temps == sorted(temps)                  # monotone in T
        assert enth == sorted(enth)                    # cumulative H


def test_thermal_stream_objects_and_isothermal_widening():
    streams = [
        ThermalStream(T_in=400.0, T_out=400.0, Q=1e5, kind="cold"),   # reboiler-like
        ThermalStream(T_in=450.0, T_out=450.0, Q=1e5, kind="hot"),    # condenser-like
    ]
    res = pinch_from_streams(streams, dt_min=10.0)
    # The hot duty sits 50 K above the cold one: full recovery is possible,
    # so both utility targets are zero (the duties match exactly).
    assert res.qh_min == pytest.approx(0.0, abs=1e-6)
    assert res.qc_min == pytest.approx(0.0, abs=1e-6)
    assert res.heat_recovery_potential == pytest.approx(1e5, rel=1e-9)
    # Isothermal duties were widened to 1 K segments.
    assert all(abs(s.T_in - s.T_out) == pytest.approx(1.0) for s in res.streams)


def test_isothermal_dict_uses_duty_sign_for_kind():
    # Engine convention: Q > 0 enters the process (cold demand), Q < 0 leaves.
    res = pinch_from_streams([{"Tin": 400.0, "Tout": 400.0, "Q": 5e4}], dt_min=10.0)
    assert res.streams[0].kind == "cold"
    assert res.qh_min == pytest.approx(5e4)
    assert res.qc_min == pytest.approx(0.0, abs=1e-9)
    assert res.pinch_T_shifted is None                  # threshold problem


def test_negligible_and_empty_inputs():
    empty = pinch_from_streams([], dt_min=10.0)
    assert empty.qh_min == empty.qc_min == 0.0
    assert empty.pinch_T_shifted is None
    tiny = pinch_from_streams([{"Tin": 300.0, "Tout": 400.0, "Q": 1e-9}])
    assert tiny.qh_min == 0.0 and not tiny.streams


def test_bad_inputs_raise():
    with pytest.raises(ValueError, match="dt_min"):
        pinch_from_streams(KEMP_STREAMS, dt_min=-1.0)
    with pytest.raises(ValueError, match="Tin"):
        pinch_from_streams([{"T": 300.0, "Q": 1e5}])
    with pytest.raises(ValueError, match="kind"):
        pinch_from_streams([ThermalStream(300.0, 400.0, 1e5, kind="warm")])


# -- flowsheet level ---------------------------------------------------------
AMMONIA = {"stoich": {"nitrogen": -1, "hydrogen": -3, "ammonia": 2}, "key": "nitrogen"}


def ammonia_loop() -> Flowsheet:
    """The examples/04 Haber-Bosch loop: PREHEAT (cold) vs RXN+COOL (hot)."""
    fs = Flowsheet(
        components=[Component(c) for c in ("nitrogen", "hydrogen", "ammonia", "argon")],
        property_package="thermo:PR")
    fs.add(Mixer("MIX", {"dP": 0.0}))
    fs.add(Heater("PREHEAT", {"T_out": 673.15}))
    fs.add(EquilibriumReactor("RXN", {"reaction": AMMONIA, "T": 673.15}))
    fs.add(Heater("COOL", {"T_out": 250.0}))
    fs.add(FlashDrum("SEP", {"T": 250.0, "P": 2e7}))
    fs.add(Splitter("SPLIT", {"split": 0.90}))
    fs.feed("MAKEUP", "MIX:in1", T=300.0, P=2e7, molar_flow=100.0,
            z={"nitrogen": 0.2475, "hydrogen": 0.7425, "ammonia": 0.0, "argon": 0.01})
    fs.connect("S1", "MIX:out", "PREHEAT:in1")
    fs.connect("S2", "PREHEAT:out", "RXN:in1")
    fs.connect("S3", "RXN:out", "COOL:in1")
    fs.connect("S4", "COOL:out", "SEP:in1")
    fs.connect("PRODUCT", "SEP:liquid", None)
    fs.connect("VAP", "SEP:vapor", "SPLIT:in1")
    fs.connect("RECYCLE", "SPLIT:out1", "MIX:in2")
    fs.connect("PURGE", "SPLIT:out2", None)
    for u in ("PREHEAT", "RXN", "COOL", "SEP"):
        fs.connect(f"Q_{u}", f"{u}:duty", None)
    return fs


def test_ammonia_loop_has_nonzero_recovery_potential():
    fs = ammonia_loop()
    report = fs.solve(tol=1e-7, max_iter=400)
    assert report.converged
    res = pinch_analysis(fs, report, dt_min=10.0)

    by_unit = {s.unit_id: s for s in res.streams}
    # PREHEAT is the only heating demand; COOL and the exothermic isothermal
    # reactor reject heat over an overlapping temperature range.
    assert by_unit["PREHEAT"].kind == "cold"
    assert by_unit["COOL"].kind == "hot"
    assert by_unit["RXN"].kind == "hot"
    assert res.current_hot_utility == pytest.approx(report.duties["Q_PREHEAT"], rel=1e-9)

    # Sanity: the target can never exceed what we buy today, and the overlap
    # means real recovery is on the table.
    assert res.qh_min <= res.current_hot_utility + 1e-6
    assert res.heat_recovery_potential > 0.0
    assert res.heat_recovery_potential == \
        pytest.approx(res.current_hot_utility - res.qh_min, rel=1e-12)
    # Same overall identity as the stream-level form.
    assert res.qc_min - res.qh_min == \
        pytest.approx(res.current_cold_utility - res.current_hot_utility, rel=1e-6)
