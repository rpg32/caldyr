"""M5 acceptance tests: the FastAPI bridge (headless, via TestClient).

The DoD — build + solve + cost a flowsheet — is exercised over HTTP exactly as
the browser will: POST a ``.flow`` document, get back resolved streams, costs,
and an optimization result. No physics lives in the transport layer; these tests
confirm the engine is faithfully exposed.
"""
import math

import pytest
from fastapi.testclient import TestClient

from api.main import app
from caldyr.io import to_dict
from caldyr.unitops import EquilibriumReactor, FlashDrum, Heater, Mixer, Splitter

from caldyr.core import Component, Flowsheet  # noqa: E402  (after api path bootstrap)

client = TestClient(app)
AMMONIA = {"stoich": {"nitrogen": -1, "hydrogen": -3, "ammonia": 2}, "key": "nitrogen"}


def mixer_heater_flow() -> dict:
    return {
        "schema": "caldyr.flow/1",
        "components": [{"id": "water"}, {"id": "ethanol"}],
        "property_package": "thermo:PR",
        "units": [
            {"id": "MIX", "type": "Mixer", "params": {"dP": 0.0}},
            {"id": "H", "type": "Heater", "params": {"T_out": 350.0, "dP": 0.0}},
        ],
        "streams": [
            {"id": "S1", "from": None, "to": "MIX:in1",
             "spec": {"T": 298.15, "P": 101325, "molar_flow": 10,
                      "z": {"water": 0.6, "ethanol": 0.4}}},
            {"id": "S2", "from": None, "to": "MIX:in2",
             "spec": {"T": 320.0, "P": 101325, "molar_flow": 5, "z": {"water": 1.0}}},
            {"id": "S3", "from": "MIX:out", "to": "H:in1"},
            {"id": "S4", "from": "H:out", "to": None},
            {"id": "Q", "from": "H:duty", "to": None},
        ],
    }


def flash_recycle_flow() -> dict:
    fs = Flowsheet(components=[Component("n-pentane"), Component("n-octane")],
                   property_package="thermo:PR")
    fs.add(Mixer("MIX", {"dP": 0.0}))
    fs.add(FlashDrum("FL", {"T": 360.0, "P": 101325.0}))
    fs.add(Splitter("SP", {"split": 0.6}))
    fs.feed("FEED", "MIX:in1", T=330.0, P=101325.0, molar_flow=10.0,
            z={"n-pentane": 0.5, "n-octane": 0.5})
    fs.connect("MIXOUT", "MIX:out", "FL:in1")
    fs.connect("VAP", "FL:vapor", None)
    fs.connect("LIQ", "FL:liquid", "SP:in1")
    fs.connect("RECY", "SP:out1", "MIX:in2")
    fs.connect("BOT", "SP:out2", None)
    fs.connect("Q", "FL:duty", None)
    return to_dict(fs)


def ammonia_flow() -> dict:
    fs = Flowsheet(components=[Component(c) for c in ("nitrogen", "hydrogen", "ammonia", "argon")],
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
    return to_dict(fs)


# -- metadata --------------------------------------------------------------
def test_health_and_metadata():
    assert client.get("/health").json()["status"] == "ok"
    types = {u["type"] for u in client.get("/unit-types").json()}
    assert {"Mixer", "Heater", "Flash", "EquilibriumReactor"} <= types
    # every unit-type advertises its ports for the palette
    heater = next(u for u in client.get("/unit-types").json() if u["type"] == "Heater")
    assert {"in1", "out", "duty"} == {p["name"] for p in heater["ports"]}
    assert "thermo:NRTL" in {p["id"] for p in client.get("/property-packages").json()}


# -- solve -----------------------------------------------------------------
@pytest.mark.parametrize("backend", ["sequential", "equation_oriented"])
def test_solve_mixer_heater(backend):
    r = client.post("/solve", json={"flow": mixer_heater_flow(), "backend": backend})
    assert r.status_code == 200
    body = r.json()
    assert body["report"]["converged"]
    assert body["streams"]["S4"]["T"] == pytest.approx(350.0, abs=1e-3)
    assert body["report"]["duties"]["Q"] == pytest.approx(298188.0, rel=1e-4)


def test_solve_recycle_both_backends_agree():
    seq = client.post("/solve", json={"flow": flash_recycle_flow(), "backend": "sequential"}).json()
    eo = client.post("/solve", json={"flow": flash_recycle_flow(),
                                     "backend": "equation_oriented"}).json()
    for sid in ("VAP", "BOT", "RECY"):
        assert math.isclose(seq["streams"][sid]["molar_flow"],
                            eo["streams"][sid]["molar_flow"], rel_tol=1e-5)
    assert eo["report"]["tear_streams"] == []


def test_solve_rejects_malformed_flow():
    assert client.post("/solve", json={"flow": {"schema": "nope"}}).status_code == 422


# -- cost ------------------------------------------------------------------
def test_cost_ammonia_loop():
    r = client.post("/cost", json={"flow": ammonia_flow(),
                                   "config": {"product_component": "ammonia"}})
    assert r.status_code == 200
    d = r.json()
    assert d["capital"]["tci"] > d["capital"]["isbl"] > 0
    assert d["opex"]["total"] > 0
    assert d["profitability"]["lcop"] > 0
    assert {e["unit_id"] for e in d["equipment"]} == {"PREHEAT", "RXN", "COOL", "SEP"}
    assert d["tornado"][0]["variable"].startswith("feed price")    # dominant lever


# -- optimize --------------------------------------------------------------
def test_optimize_min_duty_with_recovery_constraint():
    req = {
        "flow": flash_recycle_flow(),
        "objective": {"sense": "min", "metric": {"type": "duty", "stream": "Q"}},
        "design_vars": [{"unit_id": "FL", "param": "T", "lower": 340.0,
                         "upper": 370.0, "initial": 360.0}],
        "constraints": [{"metric": {"type": "component_rate", "stream": "VAP",
                                    "component": "n-pentane"}, "op": ">=", "value": 4.2}],
    }
    d = client.post("/optimize", json=req).json()
    assert d["success"]
    assert 340.0 <= d["design"]["FL.T"] <= 370.0
    vap = d["streams"]["VAP"]
    assert vap["molar_flow"] * vap["z"]["n-pentane"] == pytest.approx(4.2, abs=3e-3)


# -- round-trip ------------------------------------------------------------
def test_flow_roundtrip_is_stable():
    flow = ammonia_flow()
    once = client.post("/flow/roundtrip", json=flow).json()
    twice = client.post("/flow/roundtrip", json=once).json()
    assert once == twice


# -- phase envelope (UI plots) ----------------------------------------------
def test_envelope_for_feed_stream():
    d = client.post("/envelope", json={
        "flow": mixer_heater_flow(), "stream": "S1", "n": 12,
    }).json()
    assert d["stream"] == "S1"
    assert len(d["points"]) >= 8                     # most of the grid feasible
    for p in d["points"]:
        # dew point is never below the bubble point
        assert p["T_dew"] >= p["T_bubble"] - 1e-6
    # higher pressure -> higher bubble point (monotone for this mixture)
    assert d["points"][-1]["T_bubble"] > d["points"][0]["T_bubble"]


def test_envelope_unknown_stream_is_422():
    r = client.post("/envelope", json={"flow": mixer_heater_flow(), "stream": "NOPE"})
    assert r.status_code == 422


# -- monte-carlo over /cost ---------------------------------------------------
def test_cost_with_monte_carlo_samples():
    d = client.post("/cost", json={
        "flow": ammonia_flow(),
        "config": {"product_component": "ammonia"},
        "tornado": False,
        "monte_carlo": 60,
    }).json()
    mc = d["monte_carlo"]
    assert mc["n"] == 60 and len(mc["lcop_samples"]) == 60
    assert mc["lcop"]["p10"] <= mc["lcop"]["p50"] <= mc["lcop"]["p90"]
    # distribution brackets the deterministic value
    assert mc["lcop"]["p10"] < d["profitability"]["lcop"] < mc["lcop"]["p90"]


def test_solve_report_includes_residual_history():
    d = client.post("/solve", json={"flow": flash_recycle_flow()}).json()
    h = d["report"]["history"]
    assert len(h) == d["report"]["iterations"] >= 2
    assert h[-1] < h[0]                              # residual fell


# -- analysis tools: property table / relief / pinch -------------------------
def npentane_flow() -> dict:
    """A single n-pentane stream, the book's §2.1.4 Property Table case."""
    return {
        "schema": "caldyr.flow/1",
        "components": [{"id": "n-pentane"}],
        "property_package": "thermo:PR",
        "units": [{"id": "H", "type": "Heater", "params": {"T_out": 350.0, "dP": 0.0}}],
        "streams": [
            {"id": "F", "from": None, "to": "H:in1",
             "spec": {"T": 300.0, "P": 1.2e6, "molar_flow": 10, "z": {"n-pentane": 1.0}}},
            {"id": "O", "from": "H:out", "to": None},
            {"id": "Q", "from": "H:duty", "to": None},
        ],
    }


def test_property_table_grid_and_trends():
    d = client.post("/property-table", json={
        "flow": npentane_flow(), "stream": "F",
        "T": [500.0, 550.0, 600.0], "P": [1.2e6, 1.6e6],
        "props": ["mass_density", "vapor_fraction"],
    }).json()
    assert d["T"] == [500.0, 550.0, 600.0] and d["P"] == [1.2e6, 1.6e6]
    rho = d["values"]["mass_density"]
    assert len(rho) == 3 and len(rho[0]) == 2          # (n_T, n_P) grid
    # density falls with T (down a column) and rises with P (across a row) — the
    # qualitative shape of the book's Fig. 2.20.
    assert rho[0][0] > rho[2][0]
    assert rho[0][1] > rho[0][0]
    assert d["failures"] == []


def test_property_table_explicit_z_overrides_stream():
    d = client.post("/property-table", json={
        "flow": npentane_flow(), "z": {"n-pentane": 1.0},
        "T": [550.0], "P": [1.5e6], "props": ["mass_density"],
    }).json()
    assert d["values"]["mass_density"][0][0] > 0.0


def test_property_table_needs_a_composition_source():
    r = client.post("/property-table", json={
        "flow": npentane_flow(), "T": [500.0], "P": [1.2e6]})
    assert r.status_code == 422


def test_relief_vapor_sizes_and_selects_orifice():
    d = client.post("/relief", json={
        "phase": "vapor", "W": 10.0, "T": 500.0, "M": 0.072,
        "Z": 0.9, "k": 1.05, "P1": 1.2e6,
    }).json()
    assert d["phase"] == "vapor" and d["critical"] is True
    assert d["area_m2"] > 0.0 and d["area_cm2"] == pytest.approx(d["area_m2"] * 1e4)
    assert d["orifice"] in list("DEFGHJKLMNPQRT")     # an API 526 letter
    assert 0.0 < d["capacity_used"] <= 1.0            # fits the selected orifice


def test_relief_vapor_derives_M_from_stream():
    """M omitted -> derived from the flow/stream composition."""
    d = client.post("/relief", json={
        "phase": "vapor", "W": 10.0, "T": 500.0, "k": 1.05, "P1": 1.2e6,
        "flow": npentane_flow(), "stream": "F",
    }).json()
    # n-pentane M ~ 0.0721 kg/mol -> a sane area.
    assert d["area_m2"] > 0.0 and d["orifice"] is not None


def test_relief_liquid_path():
    d = client.post("/relief", json={
        "phase": "liquid", "W": 5.0, "rho": 600.0, "P1": 1.2e6, "P2": 1.0e5,
    }).json()
    assert d["phase"] == "liquid" and d["critical"] is None
    assert d["area_m2"] > 0.0


def test_relief_bad_inputs_are_handled():
    # vapor missing T/k -> 422
    assert client.post("/relief", json={
        "phase": "vapor", "W": 10.0, "M": 0.072, "P1": 1.2e6}).status_code == 422
    # subcritical backpressure -> the physics raises -> 400
    r = client.post("/relief", json={
        "phase": "vapor", "W": 10.0, "T": 500.0, "M": 0.072, "Z": 0.9,
        "k": 1.4, "P1": 1.2e6, "backpressure": 1.1e6})
    assert r.status_code == 400


def test_pinch_targets_and_composites():
    d = client.post("/pinch", json={"flow": ammonia_flow(), "dt_min": 10.0}).json()
    assert d["dt_min"] == 10.0
    assert d["qh_min"] >= 0.0 and d["qc_min"] >= 0.0
    assert d["heat_recovery_potential"] >= -1e-6
    # composite curves are plot-ready (H, T) point lists when both kinds exist
    assert isinstance(d["hot_composite"], list)
    assert isinstance(d["cold_composite"], list)
    assert {s["kind"] for s in d["streams"]} <= {"hot", "cold"}
