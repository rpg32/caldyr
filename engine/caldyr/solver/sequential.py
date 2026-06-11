from dataclasses import dataclass, field

from ..core import EnergyStream, Stream
from . import recycle


@dataclass
class SolveReport:
    converged: bool
    iterations: int = 0
    residual: float | None = None
    tol: float = 1e-6
    method: str = "direct"
    order: list[str] = field(default_factory=list)
    tear_streams: list[str] = field(default_factory=list)
    history: list[float] = field(default_factory=list)   # residual per iteration
    duties: dict[str, float] = field(default_factory=dict)   # energy-stream id -> W
    messages: list[str] = field(default_factory=list)


class CycleError(Exception):
    """Raised when topological ordering is impossible after tearing (should not
    happen for a valid flowsheet — indicates a tear-selection bug)."""


class SequentialModularSolver:
    """Sequential-modular solver.

    Acyclic flowsheets solve in a single topological sweep. Flowsheets with
    recycles are torn (DFS back edges) and converged by direct substitution,
    accelerated with Wegstein once two iterates are available. Never fails
    silently: a non-converged run returns ``converged=False`` with the residual
    history rather than raising.
    """

    def __init__(self, tol: float = 1e-6, max_iter: int = 200, method: str = "wegstein",
                 on_iteration=None) -> None:
        self.tol = tol
        self.max_iter = max_iter
        self.method = method  # "wegstein" or "direct"
        # Optional progress hook called as on_iteration(iteration, residual)
        # after every recycle sweep (used by the API's WebSocket channel).
        self.on_iteration = on_iteration

    def solve(self, flowsheet, pp) -> SolveReport:
        # Solver hints (the Recycle-block equivalent, persisted in `.flow` under
        # "solver_hints"): per-solve tolerance override and tear-stream guesses.
        hints = getattr(flowsheet, "solver_hints", None) or {}
        if hints.get("tear_tolerance") is not None:
            self.tol = float(hints["tear_tolerance"])

        tears = recycle.find_tears(flowsheet)
        order = self._topological_order(flowsheet, skip=set(tears))

        # Boundary feeds are constant; resolve them once (fill H/phase/VF).
        for conn in flowsheet.connections:
            if conn.from_unit is None and conn.to_unit is not None:
                self._resolve_feed(flowsheet.streams[conn.stream_id], pp)

        if not tears:
            duties = self._sweep(flowsheet, order, pp)
            return SolveReport(
                converged=True, iterations=1, residual=0.0, tol=self.tol,
                method="direct", order=order, duties=duties,
                messages=[f"acyclic; solved {len(order)} units in one sweep"],
            )

        return self._solve_recycle(flowsheet, order, tears, pp,
                                   tear_guesses=hints.get("tear_guesses") or {})

    # -- recycle loop ------------------------------------------------------
    def _solve_recycle(self, flowsheet, order, tears, pp,
                       tear_guesses: dict | None = None) -> SolveReport:
        comps = flowsheet.component_ids
        guess_msgs: list[str] = []
        for tid in tears:
            hint = (tear_guesses or {}).get(tid)
            if hint is not None:                # user-provided starting guess
                flowsheet.streams[tid] = self._seed_from_hint(tid, hint, comps)
                guess_msgs.append(f"tear {tid!r} seeded from solver_hints guess")
            elif tid not in flowsheet.streams:  # default: empty recycle
                flowsheet.streams[tid] = self._seed_tear(tid, comps, pp)

        history: list[float] = []
        duties: dict[str, float] = {}
        x_prev: dict[str, list[float]] | None = None
        g_prev: dict[str, list[float]] | None = None
        converged = False
        used_method = "direct"
        residual = None

        for it in range(1, self.max_iter + 1):
            guess = {tid: recycle.tear_vector(flowsheet.streams[tid], comps) for tid in tears}
            duties = self._sweep(flowsheet, order, pp)
            computed = {tid: recycle.tear_vector(flowsheet.streams[tid], comps) for tid in tears}
            comp_streams = {tid: flowsheet.streams[tid] for tid in tears}

            residual = recycle.residual(guess, computed)
            history.append(residual)
            if self.on_iteration is not None:
                self.on_iteration(it, residual)
            if residual < self.tol:
                converged = True
                break

            # Build the next guess (direct substitution, or Wegstein if enabled
            # and a previous iterate exists), then write it back for the next sweep.
            use_wegstein = (
                self.method == "wegstein" and x_prev is not None and g_prev is not None
            )
            used_method = "wegstein" if use_wegstein else "direct"
            for tid in tears:
                if use_wegstein and x_prev is not None and g_prev is not None:
                    nxt = recycle.wegstein_update(
                        x_prev[tid], g_prev[tid], guess[tid], computed[tid]
                    )
                else:
                    nxt = computed[tid]
                flowsheet.streams[tid] = recycle.stream_from_vector(
                    tid, nxt, comps, P=comp_streams[tid].P, fallback_z=comp_streams[tid].z
                )
            x_prev, g_prev = guess, computed

        msg = (
            f"recycle converged in {len(history)} iterations (residual {residual:.2e})"
            if converged
            else f"NOT converged after {self.max_iter} iterations (residual {residual:.2e})"
        )
        return SolveReport(
            converged=converged, iterations=len(history), residual=residual,
            tol=self.tol, method=used_method, order=order, tear_streams=list(tears),
            history=history, duties=duties, messages=guess_msgs + [msg],
        )

    # -- one sweep ---------------------------------------------------------
    def _sweep(self, flowsheet, order, pp) -> dict[str, float]:
        duties: dict[str, float] = {}
        for uid in order:
            unit = flowsheet.units[uid]
            inlets = self._gather_inlets(flowsheet, uid)
            outlets = unit.solve(inlets, pp)
            self._store_outlets(flowsheet, uid, outlets, duties)
        return duties

    # -- internals ---------------------------------------------------------
    @staticmethod
    def _seed_from_hint(tid: str, hint: dict, comps: list[str]) -> Stream:
        """Build a tear-stream starting guess from a `.flow` solver-hints entry
        ``{T, P, molar_flow, z}`` (missing fields get neutral defaults)."""
        z = {c: float(v) for c, v in (hint.get("z") or {}).items()}
        if not z or sum(z.values()) <= 0.0:
            z = {c: 1.0 / len(comps) for c in comps}
        return Stream(
            id=tid, components=list(comps),
            T=float(hint.get("T", 298.15)), P=float(hint.get("P", 101325.0)),
            molar_flow=float(hint.get("molar_flow", 0.0)), z=z,
        )

    @staticmethod
    def _seed_tear(tid: str, comps: list[str], pp) -> Stream:
        """Start a recycle empty: zero flow, equimolar placeholder composition.
        Zero flow means mixers ignore it on the first sweep, so the loop builds
        the recycle up from the genuine first-pass result."""
        z = {c: 1.0 / len(comps) for c in comps}
        return Stream(id=tid, components=list(comps), T=298.15, P=101325.0,
                      molar_flow=0.0, z=z)

    @staticmethod
    def _resolve_feed(stream: Stream, pp) -> None:
        stream.validate()
        z = stream.normalized_z()
        if stream.T is None or stream.P is None:
            raise ValueError(f"feed {stream.id!r} needs both T and P specified")
        res = pp.flash_pt(stream.T, stream.P, z)
        stream.H = res.H
        stream.phase = res.phase
        stream.vapor_fraction = res.vapor_fraction

    @staticmethod
    def _gather_inlets(flowsheet, uid: str) -> dict[str, Stream]:
        inlets: dict[str, Stream] = {}
        for conn in flowsheet.connections:
            if conn.to_unit == uid:
                stream = flowsheet.streams.get(conn.stream_id)
                if stream is None:
                    raise ValueError(
                        f"unit {uid!r} inlet stream {conn.stream_id!r} is undefined "
                        f"(upstream not solved or feed spec missing)"
                    )
                inlets[conn.to_port] = stream
        return inlets

    @staticmethod
    def _store_outlets(flowsheet, uid: str, outlets, duties: dict[str, float]) -> None:
        port_to_stream = {
            conn.from_port: conn.stream_id
            for conn in flowsheet.connections
            if conn.from_unit == uid
        }
        for port, obj in outlets.items():
            if isinstance(obj, EnergyStream):
                sid = port_to_stream.get(port, obj.id)
                obj.id = sid
                if obj.duty is None:
                    raise ValueError(f"energy stream {sid!r} produced without a duty")
                duties[sid] = obj.duty
                continue
            sid = port_to_stream.get(port)
            if sid is None:
                continue  # produced but unconnected (e.g. an unused outlet)
            obj.id = sid
            flowsheet.streams[sid] = obj

    @staticmethod
    def _topological_order(flowsheet, skip: set[str]) -> list[str]:
        """Kahn topological sort over internal edges, ignoring torn streams."""
        units = list(flowsheet.units)
        succ: dict[str, list[str]] = {u: [] for u in units}
        indeg: dict[str, int] = {u: 0 for u in units}
        for conn in flowsheet.connections:
            if conn.stream_id in skip:
                continue
            if conn.from_unit in indeg and conn.to_unit in indeg:
                succ[conn.from_unit].append(conn.to_unit)
                indeg[conn.to_unit] += 1

        ready = [u for u in units if indeg[u] == 0]
        order: list[str] = []
        while ready:
            u = ready.pop(0)
            order.append(u)
            for v in succ[u]:
                indeg[v] -= 1
                if indeg[v] == 0:
                    ready.append(v)

        if len(order) != len(units):
            stuck = sorted(set(units) - set(order))
            raise CycleError(
                f"could not order units {stuck} after tearing; tear selection failed"
            )
        return order
