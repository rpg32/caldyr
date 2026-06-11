"""Working state for an AI agent building a flowsheet through tool calls.

A session is the equivalent of a REPL: tools mutate the flowsheet it holds, then
solve and cost it. Keeping the state here (rather than re-sending the whole
flowsheet with every tool call) lets the agent build incrementally.
"""
from __future__ import annotations

from dataclasses import dataclass

from ..core import Flowsheet


@dataclass
class AgentSession:
    flowsheet: Flowsheet | None = None
    report: object = None          # last SolveReport
    tea: object = None             # last TEAResult

    def require_flowsheet(self) -> Flowsheet:
        if self.flowsheet is None:
            raise ValueError("no flowsheet yet — call new_flowsheet first")
        return self.flowsheet

    def require_solved(self):
        if self.report is None:
            raise ValueError("flowsheet not solved yet — call solve first")
        return self.report
