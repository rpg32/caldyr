"""M6 acceptance tests: the typed AI tool layer.

The DoD — "build me an ammonia loop and cost it" produces a solved, costed
flowsheet — is a *capability* of the tool layer: an agent only sequences these
tool calls. So the headline test drives the same tools deterministically (no LLM)
and asserts the flowsheet solves and costs. The live agent (agent.run) needs an
API key and is not exercised here.
"""
import math

import pytest

from caldyr.ai import (
    AgentSession,
    LLMResponse,
    OllamaBackend,
    OpenAIBackend,
    ToolCall,
    anthropic_tools,
    dispatch,
    make_backend,
    run,
)
from caldyr.ai.tools import TOOLS

AMMONIA_RECIPE = [
    ("new_flowsheet", {"components": ["nitrogen", "hydrogen", "ammonia", "argon"],
                       "property_package": "thermo:PR"}),
    ("add_unit", {"id": "MIX", "type": "Mixer", "params": {"dP": 0.0}}),
    ("add_unit", {"id": "PREHEAT", "type": "Heater", "params": {"T_out": 673.15}}),
    ("add_unit", {"id": "RXN", "type": "EquilibriumReactor",
                  "params": {"reaction": {"stoich": {"nitrogen": -1, "hydrogen": -3, "ammonia": 2},
                                          "key": "nitrogen"}, "T": 673.15}}),
    ("add_unit", {"id": "COOL", "type": "Heater", "params": {"T_out": 250.0}}),
    ("add_unit", {"id": "SEP", "type": "Flash", "params": {"T": 250.0, "P": 2e7}}),
    ("add_unit", {"id": "SPLIT", "type": "Splitter", "params": {"split": 0.9}}),
    ("add_feed", {"id": "MAKEUP", "to": "MIX:in1", "T": 300.0, "P": 2e7, "molar_flow": 100.0,
                  "z": {"nitrogen": 0.2475, "hydrogen": 0.7425, "ammonia": 0.0, "argon": 0.01}}),
    ("connect", {"id": "S1", "from": "MIX:out", "to": "PREHEAT:in1"}),
    ("connect", {"id": "S2", "from": "PREHEAT:out", "to": "RXN:in1"}),
    ("connect", {"id": "S3", "from": "RXN:out", "to": "COOL:in1"}),
    ("connect", {"id": "S4", "from": "COOL:out", "to": "SEP:in1"}),
    ("connect", {"id": "PRODUCT", "from": "SEP:liquid", "to": None}),
    ("connect", {"id": "VAP", "from": "SEP:vapor", "to": "SPLIT:in1"}),
    ("connect", {"id": "RECYCLE", "from": "SPLIT:out1", "to": "MIX:in2"}),
    ("connect", {"id": "PURGE", "from": "SPLIT:out2", "to": None}),
]

PENTANE_FLASH = [
    ("new_flowsheet", {"components": ["n-pentane", "n-octane"]}),
    ("add_unit", {"id": "MIX", "type": "Mixer", "params": {"dP": 0.0}}),
    ("add_unit", {"id": "FL", "type": "Flash", "params": {"T": 360.0, "P": 101325.0}}),
    ("add_unit", {"id": "SP", "type": "Splitter", "params": {"split": 0.6}}),
    ("add_feed", {"id": "FEED", "to": "MIX:in1", "T": 330.0, "P": 101325.0,
                  "molar_flow": 10.0, "z": {"n-pentane": 0.5, "n-octane": 0.5}}),
    ("connect", {"id": "MIXOUT", "from": "MIX:out", "to": "FL:in1"}),
    ("connect", {"id": "VAP", "from": "FL:vapor", "to": None}),
    ("connect", {"id": "LIQ", "from": "FL:liquid", "to": "SP:in1"}),
    ("connect", {"id": "RECY", "from": "SP:out1", "to": "MIX:in2"}),
    ("connect", {"id": "BOT", "from": "SP:out2", "to": None}),
]


def _run(session, calls):
    out = {}
    for name, args in calls:
        out = dispatch(session, name, args)
        assert out["ok"], (name, out)
    return out


# -- tool schemas ----------------------------------------------------------
def test_tool_schemas_are_well_formed():
    names = {t.name for t in TOOLS}
    assert {"new_flowsheet", "add_unit", "add_feed", "connect", "solve", "cost"} <= names
    for t in anthropic_tools():
        assert t["name"] and t["description"]
        assert t["input_schema"]["type"] == "object"


def test_list_tools_describe_the_engine():
    s = AgentSession()
    ut = dispatch(s, "list_unit_types", {})
    assert any(u["type"] == "EquilibriumReactor" for u in ut["unit_types"])
    pp = dispatch(s, "list_property_packages", {})
    assert any(p["id"] == "thermo:NRTL" for p in pp["property_packages"])


def test_list_property_packages_serves_the_full_registry():
    """The tool serves thermo.AVAILABLE_PACKAGES — not a hand-copied subset."""
    s = AgentSession()
    ids = {p["id"] for p in dispatch(s, "list_property_packages", {})["property_packages"]}
    assert {"thermo:PR", "thermo:SRK", "thermo:NRTL", "thermo:UNIFAC",
            "thermo:UNIFAC-LLE", "coolprop:Water", "amine:DEA", "amine:MDEA",
            "amine:MEA", "nasa:gas"} <= ids


def test_list_unit_types_single_type_full_docs():
    s = AgentSession()
    out = dispatch(s, "list_unit_types", {"type": "RigorousColumn"})
    assert out["ok"] and out["type"] == "RigorousColumn"
    assert "n_stages" in out["doc"]                   # full params documentation
    assert any(p["name"] == "distillate" for p in out["ports"])
    assert dispatch(s, "list_unit_types", {"type": "Nope"})["ok"] is False


def test_system_prompt_covers_every_unit_type():
    from caldyr.ai.agent import SYSTEM_PROMPT
    from caldyr.unitops import REGISTRY
    missing = [name for name in REGISTRY if name not in SYSTEM_PROMPT]
    assert not missing, f"SYSTEM_PROMPT lacks port docs for {missing}"


# -- the DoD: build + solve + cost via tools -------------------------------
def test_agent_tools_build_solve_cost_ammonia_loop():
    s = AgentSession()
    _run(s, AMMONIA_RECIPE)

    solved = dispatch(s, "solve", {"tol": 1e-7})
    assert solved["converged"]
    assert solved["streams"]                       # a populated stream table
    assert "stream_table" in solved

    costed = dispatch(s, "cost", {"product_component": "ammonia"})
    assert costed["lcop"] > 0
    assert costed["capital"]["tci"] > 0
    assert costed["annual_production_kg"] > 0
    # LCOP matches the hand-built example flowsheet (duty sizing is robust to the
    # tools not wiring duty ports).
    assert costed["lcop"] == pytest.approx(0.690, abs=0.01)


def test_export_flow_roundtrips_through_io():
    from caldyr.io import from_dict
    s = AgentSession()
    _run(s, AMMONIA_RECIPE)
    flow = dispatch(s, "export_flow", {})["flow"]
    fs = from_dict(flow)                            # the engine can re-load it
    rep = fs.solve(tol=1e-7)
    assert rep.converged


def test_cost_auto_solves_if_needed():
    s = AgentSession()
    _run(s, AMMONIA_RECIPE)
    # no explicit solve — cost should solve first
    costed = dispatch(s, "cost", {"product_component": "ammonia"})
    assert costed["ok"] and costed["lcop"] > 0


def test_optimize_tool():
    s = AgentSession()
    _run(s, PENTANE_FLASH)
    res = dispatch(s, "optimize", {
        "objective": {"sense": "min", "metric": {"type": "duty", "stream": "FL.duty"}},
        "design_vars": [{"unit_id": "FL", "param": "T", "lower": 340.0, "upper": 370.0,
                         "initial": 360.0}],
        "constraints": [{"metric": {"type": "component_rate", "stream": "VAP",
                                    "component": "n-pentane"}, "op": ">=", "value": 4.2}],
    })
    assert res["ok"] and res["success"]
    assert 340.0 <= res["design"]["FL.T"] <= 370.0


# -- capability-parity tools (Phase 1 of the AI-centric plan) ----------------
def test_remove_unit_takes_attached_streams_and_invalidates():
    s = AgentSession()
    _run(s, PENTANE_FLASH)
    dispatch(s, "solve", {})
    out = dispatch(s, "remove_unit", {"id": "SP"})
    assert out["ok"] and set(out["removed_streams"]) == {"LIQ", "RECY", "BOT"}
    fs = s.flowsheet
    assert "SP" not in fs.units
    assert {c.stream_id for c in fs.connections} == {"FEED", "MIXOUT", "VAP"}
    assert s.report is None                        # structural edit -> stale solve
    assert dispatch(s, "remove_unit", {"id": "GHOST"})["ok"] is False


def test_remove_stream_then_flowsheet_still_solves():
    s = AgentSession()
    _run(s, PENTANE_FLASH)
    out = dispatch(s, "remove_stream", {"id": "RECY"})   # break the recycle
    assert out["ok"]
    assert "RECY" not in {c.stream_id for c in s.flowsheet.connections}
    assert dispatch(s, "remove_stream", {"id": "RECY"})["ok"] is False
    solved = dispatch(s, "solve", {})                    # now acyclic
    assert solved["ok"] and solved["converged"]


def test_set_param_edits_one_value():
    s = AgentSession()
    _run(s, PENTANE_FLASH)
    dispatch(s, "solve", {})
    out = dispatch(s, "set_param", {"unit_id": "SP", "param": "split", "value": 0.8})
    assert out["ok"] and s.flowsheet.units["SP"].params["split"] == 0.8
    assert s.report is None                        # param edit -> stale solve
    assert dispatch(s, "solve", {})["converged"]
    # null clears the param back to the unit's default
    ok = dispatch(s, "set_param", {"unit_id": "MIX", "param": "dP", "value": None})
    assert ok["ok"] and "dP" not in s.flowsheet.units["MIX"].params
    assert dispatch(s, "set_param",
                    {"unit_id": "GHOST", "param": "T", "value": 1.0})["ok"] is False


def test_pinch_analysis_tool():
    s = AgentSession()
    _run(s, AMMONIA_RECIPE)
    out = dispatch(s, "pinch_analysis", {"dt_min": 10.0})   # auto-solves
    assert out["ok"] and out["dt_min"] == 10.0
    assert out["qh_min"] >= 0 and out["qc_min"] >= 0
    assert out["heat_recovery_potential"] >= 0
    kinds = {t["kind"] for t in out["streams"]}
    assert kinds == {"hot", "cold"}                # PRE heats (cold), COOL cools (hot)
    assert "minimum utilities" in out["summary"]


def test_sweep_parameter_tool_leaves_session_untouched():
    s = AgentSession()
    _run(s, PENTANE_FLASH)
    dispatch(s, "solve", {})
    report_before = s.report
    out = dispatch(s, "sweep_parameter", {
        "unit_id": "FL", "param": "T", "start": 340.0, "stop": 370.0, "steps": 4,
        "metric": {"type": "component_rate", "stream": "VAP", "component": "n-pentane"}})
    assert out["ok"] and len(out["points"]) == 4
    assert out["points"][0]["value"] == 340.0 and out["points"][-1]["value"] == 370.0
    ys = [p["metric"] for p in out["points"] if p["metric"] is not None]
    assert len(ys) >= 3 and ys == sorted(ys)       # hotter flash boils off more pentane
    # a sensitivity study is read-only: param and solve state are untouched
    assert s.flowsheet.units["FL"].params["T"] == 360.0
    assert s.report is report_before
    bad = dispatch(s, "sweep_parameter", {
        "unit_id": "FL", "param": "T", "start": 340.0, "stop": 370.0, "steps": 99,
        "metric": {"type": "flow", "stream": "VAP"}})
    assert bad["ok"] is False and "steps" in bad["error"]
    typo = dispatch(s, "sweep_parameter", {
        "unit_id": "FL", "param": "T", "start": 340.0, "stop": 370.0,
        "metric": {"type": "flow", "stream": "VAPOUR"}})
    assert typo["ok"] is False and "VAPOUR" in typo["error"]


def test_property_table_tool():
    s = AgentSession()
    _run(s, PENTANE_FLASH)
    out = dispatch(s, "property_table", {
        "z": {"n-pentane": 0.5, "n-octane": 0.5},
        "T": [300.0, 400.0, 500.0], "P": 101325.0,
        "props": ["mass_density", "vapor_fraction"]})
    assert out["ok"] and len(out["T"]) == 3 and len(out["P"]) == 1
    vf = [row[0] for row in out["values"]["vapor_fraction"]]
    assert vf[0] < vf[-1]                          # heating vaporizes the mixture
    # composition can come from a stream instead of an explicit z
    from_stream = dispatch(s, "property_table",
                           {"stream": "FEED", "T": 350.0, "P": 101325.0})
    assert from_stream["ok"] and "vapor_fraction" in from_stream["values"]
    # helpful errors: no composition; foreign component
    assert dispatch(s, "property_table",
                    {"T": 300.0, "P": 1e5})["ok"] is False
    foreign = dispatch(s, "property_table",
                       {"z": {"benzene": 1.0}, "T": 300.0, "P": 1e5})
    assert foreign["ok"] is False and "benzene" in foreign["error"]


# -- robustness the LLM relies on ------------------------------------------
def test_ports_are_validated_with_helpful_errors():
    s = AgentSession()
    _run(s, [("new_flowsheet", {"components": ["water"]}),
             ("add_unit", {"id": "MIX", "type": "Mixer", "params": {}})])
    bad = dispatch(s, "connect", {"id": "x", "from": "MIX:out", "to": "MIX:nope"})
    assert bad["ok"] is False and "no port" in bad["error"]
    missing = dispatch(s, "connect", {"id": "y", "from": "GHOST:out", "to": None})
    assert missing["ok"] is False and "no unit" in missing["error"]


def test_feed_composition_is_normalized():
    s = AgentSession()
    _run(s, [("new_flowsheet", {"components": ["nitrogen", "hydrogen", "argon"]}),
             ("add_unit", {"id": "H", "type": "Heater", "params": {"T_out": 400.0}})])
    out = dispatch(s, "add_feed", {"id": "F", "to": "H:in1", "T": 300.0, "P": 1e5,
                                   "molar_flow": 10.0,
                                   "z": {"nitrogen": 0.5, "hydrogen": 0.5, "argon": 0.01}})
    assert out["ok"]                                  # 1.01 tolerated (normalized)
    z = s.flowsheet.streams["F"].z
    assert math.isclose(sum(z.values()), 1.0, rel_tol=1e-9)


def test_structural_edit_invalidates_solve_so_cost_resolves():
    """Connecting a product after solving must re-solve before costing — the bug
    that blocked the live agent."""
    s = AgentSession()
    no_product = [c for c in AMMONIA_RECIPE if c[1].get("id") != "PRODUCT"]
    _run(s, no_product)
    dispatch(s, "solve", {"tol": 1e-7})
    # product not connected yet -> cost should report the actionable error
    err = dispatch(s, "cost", {"product_component": "ammonia"})
    assert err["ok"] is False and "no product stream" in err["error"]
    # connect the product, then cost: it re-solves and succeeds
    dispatch(s, "connect", {"id": "PRODUCT", "from": "SEP:liquid", "to": None})
    assert s.report is None                           # edit invalidated the solve
    ok = dispatch(s, "cost", {"product_component": "ammonia"})
    assert ok["ok"] and ok["lcop"] > 0


# -- backend selection + agent loop (no network) ---------------------------
def test_make_backend_selection():
    assert isinstance(make_backend("ollama"), OllamaBackend)
    assert isinstance(make_backend("openai"), OpenAIBackend)
    assert make_backend("ollama").name == "ollama"          # local by default
    with pytest.raises(ValueError):
        make_backend("not-a-provider")


def test_mcp_server_constructs_with_the_tools():
    from caldyr.ai.mcp_server import build_server
    server, session = build_server()
    assert server.name == "caldyr"
    assert isinstance(session, AgentSession)


class _FakeBackend:
    """Replays a fixed tool-call script, so the agent loop is testable offline."""
    name = "fake"
    model = "scripted"

    def __init__(self, recipe):
        self._calls = [ToolCall(id=f"c{i}", name=n, arguments=a)
                       for i, (n, a) in enumerate(recipe)]
        self._i = 0

    def complete(self, system, turns, tools):
        if self._i < len(self._calls):
            call = self._calls[self._i]
            self._i += 1
            return LLMResponse(text="", tool_calls=[call])
        return LLMResponse(text="Done: built, solved, and costed.", tool_calls=[])


def test_agent_loop_drives_tools_with_a_fake_backend():
    recipe = AMMONIA_RECIPE + [
        ("solve", {"tol": 1e-7}),
        ("cost", {"product_component": "ammonia"}),
    ]
    res = run("build and cost an ammonia loop", backend=_FakeBackend(recipe))
    assert res.session.report.converged
    assert res.session.tea.profitability.lcop > 0
    assert "Done" in res.final_text
    assert res.tool_calls[-1] == "cost"


# -- error handling --------------------------------------------------------
def test_errors_are_returned_not_raised():
    s = AgentSession()
    assert dispatch(s, "no_such_tool", {})["ok"] is False
    assert dispatch(s, "add_unit", {"id": "X", "type": "Mixer"})["ok"] is False  # no flowsheet
    dispatch(s, "new_flowsheet", {"components": ["water"]})
    bad = dispatch(s, "add_unit", {"id": "X", "type": "Nonexistent"})
    assert bad["ok"] is False and "unknown unit type" in bad["error"]
    # a heater with no inlet stream is an engine error, surfaced (not raised)
    dispatch(s, "add_unit", {"id": "H", "type": "Heater", "params": {"T_out": 350.0}})
    dispatch(s, "connect", {"id": "O", "from": "H:out", "to": None})
    assert dispatch(s, "solve", {})["ok"] is False
