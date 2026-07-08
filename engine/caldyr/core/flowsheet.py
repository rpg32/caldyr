from dataclasses import dataclass, field
from typing import Any

from .component import Component
from .stream import Stream
from .unitop import UnitOp


@dataclass(frozen=True)
class Connection:
    """A stream edge: (from_unit, from_port) -> (to_unit, to_port).
    A None endpoint denotes a flowsheet boundary (feed or product)."""
    stream_id: str
    from_unit: str | None
    from_port: str | None
    to_unit: str | None
    to_port: str | None


@dataclass
class Flowsheet:
    components: list[Component] = field(default_factory=list)
    property_package: str = "thermo:PR"
    units: dict[str, UnitOp] = field(default_factory=dict)
    connections: list[Connection] = field(default_factory=list)
    # Stream states keyed by id. Feeds carry a full spec up front; internal and
    # product streams are filled in by the solver. Also serves as the `solved`
    # cache for `.flow` IO.
    streams: dict[str, Stream] = field(default_factory=dict)
    # Flowsheet-level logical operations (Set/Adjust dicts; see
    # caldyr.solver.logical) — persisted in `.flow` under "logical".
    logical: list[dict] = field(default_factory=list)
    # Solver hints (tear-stream guesses / tolerance override; see
    # SequentialModularSolver) — persisted in `.flow` under "solver_hints".
    solver_hints: dict = field(default_factory=dict)
    # The most recent SolveReport (set by solve(); not serialized).
    last_report: Any = None

    @property
    def component_ids(self) -> list[str]:
        return [c.id for c in self.components]

    def add(self, unit: UnitOp) -> UnitOp:
        self.units[unit.id] = unit
        return unit

    def connect(self, stream_id, src, dst) -> None:
        """src/dst are "UNIT:port" strings or None for a boundary."""
        fu, fp = (src.split(":") + [None])[:2] if src else (None, None)
        tu, tp = (dst.split(":") + [None])[:2] if dst else (None, None)
        self.connections.append(Connection(stream_id, fu, fp, tu, tp))

    def remove_stream(self, stream_id: str) -> None:
        """Remove one stream edge (feed, internal, or product) and its state.
        Logical ops and tear guesses referencing it are pruned so the flowsheet
        stays solvable."""
        if (stream_id not in self.streams
                and all(c.stream_id != stream_id for c in self.connections)):
            raise ValueError(f"no stream {stream_id!r}; streams so far: "
                             f"{sorted({c.stream_id for c in self.connections})}")
        self.connections = [c for c in self.connections if c.stream_id != stream_id]
        self.streams.pop(stream_id, None)
        self.logical = [op for op in self.logical
                        if (op.get("spec") or {}).get("stream") != stream_id]
        (self.solver_hints.get("tear_guesses") or {}).pop(stream_id, None)

    def remove_unit(self, unit_id: str) -> list[str]:
        """Remove a unit plus every stream attached to it (their other endpoint
        becomes an unconnected port, as in a graph editor). Logical ops varying
        or reading the unit's params are pruned. Returns the removed stream ids."""
        if unit_id not in self.units:
            raise ValueError(f"no unit {unit_id!r}; units so far: {list(self.units)}")
        del self.units[unit_id]
        dead = list(dict.fromkeys(c.stream_id for c in self.connections
                                  if unit_id in (c.from_unit, c.to_unit)))
        for sid in dead:
            self.remove_stream(sid)

        def _refs_unit(op: dict) -> bool:
            return any(isinstance(ref := op.get(role), (list, tuple))
                       and ref and ref[0] == unit_id
                       for role in ("target", "source", "vary"))

        self.logical = [op for op in self.logical if not _refs_unit(op)]
        return dead

    def feed(self, stream_id, dst, *, T=None, P=None, molar_flow=None, z=None) -> Stream:
        """Declare a boundary feed stream into ``dst`` ("UNIT:port") and register
        its spec. Returns the created Stream."""
        s = Stream(
            id=stream_id,
            components=list(self.component_ids),
            T=T, P=P, molar_flow=molar_flow, z=dict(z or {}),
        )
        self.streams[stream_id] = s
        self.connect(stream_id, None, dst)
        return s

    def solve(self, backend: str = "sequential", *, tol: float = 1e-6,
              max_iter: int = 200, method: str = "wegstein", on_iteration=None):
        """Build the selected property package and solve via the chosen backend.

        ``"sequential"`` — sequential-modular (acyclic sweep; recycles torn and
        converged with Wegstein). ``"equation_oriented"`` — assemble all unit
        equations and solve them simultaneously (no tearing). Both share physics
        and agree within tolerance. Mutates ``self.streams`` with the result.

        Logical ops in ``self.logical`` are honored in a documented order:
        **Sets** are applied (target_param = k*source + b) before each inner
        solve; **Adjusts** then wrap the whole solve in an outer root find (so
        they work under either backend). See :mod:`caldyr.solver.logical`. The
        report is also stashed on ``self.last_report``.
        """
        from ..solver.logical import solve_flowsheet

        return solve_flowsheet(self, backend, tol=tol, max_iter=max_iter, method=method,
                               on_iteration=on_iteration)
