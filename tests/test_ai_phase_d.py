"""Phase D acceptance: AI diagnostics tools, the incremental ChatAgent, and the
solver's iteration callback. The LLM is faked (deterministic transcript) so
these run offline; the live path is exercised manually against Ollama.
"""
from __future__ import annotations

import json

from caldyr.ai.chat import ChatAgent
from caldyr.ai.llm import LLMResponse, ToolCall
from caldyr.ai.session import AgentSession
from caldyr.ai.tools import dispatch
from caldyr.core import Component, Flowsheet
from caldyr.io import to_dict
from caldyr.unitops import FlashDrum, Mixer, Splitter


def flash_recycle() -> Flowsheet:
    fs = Flowsheet(components=[Component("n-pentane"), Component("n-octane")],
                   property_package="thermo:PR")
    fs.add(Mixer("MIX", {"dP": 0.0}))
    fs.add(FlashDrum("FL", {"T": 360.0, "P": 101325.0}))
    fs.add(Splitter("SP", {"split": 0.6}))
    fs.feed("FEED", "MIX:in1", T=330.0, P=101325.0, molar_flow=10.0,
            z={"n-pentane": 0.5, "n-octane": 0.5})
    fs.connect("M1", "MIX:out", "FL:in1")
    fs.connect("VAP", "FL:vapor", None)
    fs.connect("LIQ", "FL:liquid", "SP:in1")
    fs.connect("RECY", "SP:out1", "MIX:in2")
    fs.connect("BOT", "SP:out2", None)
    return fs


# -- diagnostics tools --------------------------------------------------------
def test_describe_flowsheet_tool():
    s = AgentSession()
    s.flowsheet = flash_recycle()
    out = dispatch(s, "describe_flowsheet", {})
    assert out["ok"]
    assert {u["id"] for u in out["units"]} == {"MIX", "FL", "SP"}
    assert [f["stream"] for f in out["feeds"]] == ["FEED"]
    assert {p["stream"] for p in out["products"]} == {"VAP", "BOT"}
    assert "Not solved yet" in out["summary"]


def test_explain_convergence_before_and_after_solve():
    s = AgentSession()
    s.flowsheet = flash_recycle()
    out = dispatch(s, "explain_convergence", {})
    assert out["ok"] and out["converged"] is None      # graceful pre-solve answer

    solve_out = dispatch(s, "solve", {})
    assert solve_out["ok"]
    out2 = dispatch(s, "explain_convergence", {})
    assert out2["ok"] and out2["converged"] is True
    assert out2["tear_streams"] == ["RECY"]
    assert len(out2["residual_history"]) == out2["iterations"]
    assert "Converged" in out2["summary"]


# -- iteration callback -------------------------------------------------------
def test_on_iteration_callback_streams_progress():
    fs = flash_recycle()
    seen: list[tuple[int, float]] = []
    report = fs.solve(on_iteration=lambda i, r: seen.append((i, r)))
    assert report.converged
    assert [i for i, _ in seen] == list(range(1, len(seen) + 1))
    assert len(seen) == report.iterations
    assert seen[-1][1] < seen[0][1]


# -- ChatAgent with a scripted backend ---------------------------------------
class FakeBackend:
    """Replays a fixed transcript: edit SP.split, solve, then summarize."""
    name = "fake"
    model = "scripted"

    def __init__(self) -> None:
        self.step = 0

    def complete(self, system: str, turns: list[dict], tools: list[dict]) -> LLMResponse:
        self.step += 1
        if self.step == 1:
            return LLMResponse(text="Updating the splitter.", tool_calls=[
                ToolCall("c1", "add_unit",
                         {"id": "SP", "type": "Splitter", "params": {"split": 0.8}}),
            ])
        if self.step == 2:
            # add_unit upserts: re-adding 'SP' replaces it (how chat edits params)
            last = json.loads(turns[-1]["content"])
            assert last["ok"] is True
            return LLMResponse(tool_calls=[ToolCall("c2", "solve", {})])
        return LLMResponse(text="Done — solved.")


class FakeParityBackend:
    """Drives the Phase-1 parity tools: set_param, then remove_unit, then solve."""
    name = "fake"
    model = "scripted"

    def __init__(self) -> None:
        self.step = 0

    def complete(self, system: str, turns: list[dict], tools: list[dict]) -> LLMResponse:
        self.step += 1
        if self.step == 1:
            return LLMResponse(text="Retuning the flash, dropping the splitter.", tool_calls=[
                ToolCall("c1", "set_param", {"unit_id": "FL", "param": "T", "value": 355.0}),
                ToolCall("c2", "remove_unit", {"id": "SP"}),
            ])
        if self.step == 2:
            last = json.loads(turns[-1]["content"])
            assert last["ok"] is True and set(last["removed_streams"]) == {"LIQ", "RECY", "BOT"}
            return LLMResponse(tool_calls=[
                ToolCall("c3", "connect", {"id": "LIQ", "from": "FL:liquid", "to": None}),
                ToolCall("c4", "solve", {}),
            ])
        return LLMResponse(text="Done — once-through now.")


def test_chat_agent_edits_via_parity_tools():
    """set_param / remove_unit flow through the chat loop and land in the doc."""
    agent = ChatAgent(backend=FakeParityBackend())
    agent.load_flow(to_dict(flash_recycle()))

    events: list[dict] = []
    out = agent.send("run the flash at 355 K once-through, no recycle",
                     on_event=events.append)

    assert all(e["ok"] for e in events if e["type"] == "tool_result")
    assert out["text"] == "Done — once-through now."
    assert not any(u["id"] == "SP" for u in out["flow"]["units"])
    sids = {s["id"] for s in out["flow"]["streams"]}
    assert "RECY" not in sids and "BOT" not in sids
    fl = next(u for u in out["flow"]["units"] if u["id"] == "FL")
    assert fl["params"]["T"] == 355.0
    assert "solved" in out["flow"]                  # the re-solve converged and rode back


def test_remove_prunes_logical_ops_and_tear_guesses():
    """Removal keeps the flowsheet solvable: logical ops and solver hints that
    reference the removed unit/streams go with it."""
    fs = flash_recycle()
    fs.logical.append({"type": "adjust", "vary": ["SP", "split"],
                       "spec": {"type": "flow", "stream": "VAP"}, "target": 5.0})
    fs.logical.append({"type": "set", "target": ["FL", "T"],
                       "source": ["FL", "P"], "multiplier": 1.0})
    fs.solver_hints["tear_guesses"] = {"RECY": {"molar_flow": 2.0}}

    removed = fs.remove_unit("SP")
    assert set(removed) == {"LIQ", "RECY", "BOT"}
    assert fs.solver_hints["tear_guesses"] == {}    # RECY guess went with the stream
    assert fs.logical == [{"type": "set", "target": ["FL", "T"],
                           "source": ["FL", "P"], "multiplier": 1.0}]

    fs2 = flash_recycle()
    fs2.logical.append({"type": "adjust", "vary": ["FL", "T"],
                        "spec": {"type": "flow", "stream": "VAP"}, "target": 5.0})
    fs2.remove_stream("VAP")                        # the adjust's observed stream
    assert fs2.logical == []


def test_chat_agent_syncs_canvas_flow_and_returns_doc():
    agent = ChatAgent(backend=FakeBackend())
    agent.load_flow(to_dict(flash_recycle()))

    events: list[dict] = []
    out = agent.send("change the split to 0.8 and solve",
                     on_event=events.append)

    kinds = [e["type"] for e in events]
    assert kinds.count("tool_call") == 2
    assert kinds.count("tool_result") == 2
    assert all(e["ok"] for e in events if e["type"] == "tool_result")

    assert out["text"] == "Done — solved."
    assert out["flow"] is not None
    # the edit landed in the returned doc
    sp = next(u for u in out["flow"]["units"] if u["id"] == "SP")
    assert sp["params"]["split"] == 0.8
    # conversation persists across messages
    assert any(t["role"] == "user" for t in agent.turns)
    # solved state rides back in the exported doc
    assert "solved" in out["flow"]
