"""Equation-oriented (simultaneous) solve backend.

Where the sequential-modular solver propagates stream-by-stream and tears
recycles, this backend assembles *all* the flowsheet's equations into one system
and solves them at once with a Newton method (scipy). Recycles need no tearing —
they are just more equations.

The trick that keeps the physics identical to the sequential solver: the residual
for each unit is ``x_outlet − unit.solve(x_inlets)``. We reuse the very same
``unit.solve`` and ``PropertyPackage`` (real PR/NRTL flashes), so the two
backends converge to the same answer — which is exactly what the M4 cross-check
asserts. A production EO backend (IDAES/Pyomo + IPOPT, with algebraic property
models and analytic derivatives) can replace the scipy core behind this same
interface; it is intentionally pluggable.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import root

from ..core import EnergyStream, Stream
from .sequential import SolveReport
from .recycle import find_tears

_T_SCALE = 100.0       # K
_P_SCALE = 1e5         # Pa (-> pressure residuals in ~bar)


@dataclass
class _Layout:
    """Maps the global variable vector to/from computed-stream states."""
    components: list[str]
    stream_ids: list[str]                 # computed material streams (the unknowns)
    flow_scale: float
    nc: int = field(init=False)

    def __post_init__(self):
        self.nc = len(self.components)

    @property
    def block(self) -> int:
        return self.nc + 2                # [n_1..n_C, T, P] per stream

    @property
    def size(self) -> int:
        return self.block * len(self.stream_ids)

    def scales(self) -> np.ndarray:
        per = np.array([self.flow_scale] * self.nc + [_T_SCALE, _P_SCALE])
        return np.tile(per, len(self.stream_ids))

    def pack(self, streams: dict[str, Stream]) -> np.ndarray:
        out = np.empty(self.size)
        for i, sid in enumerate(self.stream_ids):
            s = streams[sid]
            z = s.normalized_z() if sum(s.z.values()) > 0 else {}
            n = s.molar_flow or 0.0
            base = i * self.block
            for j, c in enumerate(self.components):
                out[base + j] = n * z.get(c, 0.0)
            out[base + self.nc] = s.T if s.T is not None else 298.15
            out[base + self.nc + 1] = s.P if s.P is not None else 101325.0
        return out

    def unpack_one(self, x: np.ndarray, i: int) -> Stream:
        base = i * self.block
        flows = x[base: base + self.nc]
        total = float(flows.sum())
        if total > 1e-12:
            z = {c: float(flows[j]) / total for j, c in enumerate(self.components)}
        else:
            total, z = 0.0, {c: 1.0 / self.nc for c in self.components}
        return Stream(id=self.stream_ids[i], components=list(self.components),
                      T=float(x[base + self.nc]), P=float(x[base + self.nc + 1]),
                      molar_flow=total, z=z)


class EquationOrientedSolver:
    """Solve a whole flowsheet's equations simultaneously (scipy Newton)."""

    def __init__(self, tol: float = 1e-8, max_nfev: int | None = None) -> None:
        self.tol = tol
        self.max_nfev = max_nfev

    def solve(self, flowsheet, pp) -> SolveReport:
        layout = self._layout(flowsheet, pp)

        # Resolve boundary feeds once (fill H/phase) — they are fixed inputs.
        for conn in flowsheet.connections:
            if conn.from_unit is None and conn.to_unit is not None:
                self._resolve_feed(flowsheet.streams[conn.stream_id], pp)

        self._warm_start(flowsheet, layout, pp)
        x0 = layout.pack(flowsheet.streams) / layout.scales()

        scales = layout.scales()

        def residual(x_scaled: np.ndarray) -> np.ndarray:
            x = x_scaled * scales
            # Write the trial states into the flowsheet (H unset -> recomputed).
            for i, sid in enumerate(layout.stream_ids):
                flowsheet.streams[sid] = layout.unpack_one(x, i)
            target = x.copy()
            computed = np.empty_like(x)
            for i, sid in enumerate(layout.stream_ids):
                computed[i * layout.block:(i + 1) * layout.block] = \
                    self._computed_block(flowsheet, sid, layout, pp)
            return (computed - target) / scales

        sol = root(residual, x0, method="hybr",
                   options={"xtol": self.tol, **({"maxfev": self.max_nfev}
                                                 if self.max_nfev else {})})

        # Commit the solution and a final pass to populate duties / resolved H.
        x = sol.x * scales
        for i, sid in enumerate(layout.stream_ids):
            flowsheet.streams[sid] = layout.unpack_one(x, i)
        duties = self._finalize(flowsheet, pp)

        res = float(np.max(np.abs(sol.fun))) if sol.fun is not None else float("nan")
        converged = bool(sol.success) and res < max(self.tol * 100, 1e-6)
        return SolveReport(
            converged=converged, iterations=int(sol.nfev), residual=res,
            tol=self.tol, method="equation_oriented",
            order=list(flowsheet.units), tear_streams=[],
            duties=duties,
            messages=[f"simultaneous solve: {sol.nfev} residual evals; {sol.message}"],
        )

    # -- helpers -----------------------------------------------------------
    def _layout(self, flowsheet, pp) -> _Layout:
        sids = []
        for conn in flowsheet.connections:
            if conn.from_unit is None:
                continue
            unit = flowsheet.units[conn.from_unit]
            port = next((p for p in unit.ports if p.name == conn.from_port), None)
            if port is not None and port.kind == "material":
                sids.append(conn.stream_id)
        feeds = [flowsheet.streams[c.stream_id] for c in flowsheet.connections
                 if c.from_unit is None and c.stream_id in flowsheet.streams]
        flow_scale = max((s.molar_flow or 0.0) for s in feeds) if feeds else 1.0
        return _Layout(list(flowsheet.component_ids), sids, max(flow_scale, 1.0))

    def _computed_block(self, flowsheet, sid: str, layout: _Layout, pp) -> np.ndarray:
        """The (n_i, T, P) a unit *computes* for the stream it produces."""
        conn = next(c for c in flowsheet.connections if c.stream_id == sid)
        unit = flowsheet.units[conn.from_unit]
        inlets = self._gather_inlets(flowsheet, conn.from_unit)
        outlets = unit.solve(inlets, pp)
        s = outlets[conn.from_port]
        z = s.normalized_z() if sum(s.z.values()) > 0 else {}
        n = s.molar_flow or 0.0
        block = np.empty(layout.block)
        for j, c in enumerate(layout.components):
            block[j] = n * z.get(c, 0.0)
        block[layout.nc] = s.T
        block[layout.nc + 1] = s.P
        return block

    @staticmethod
    def _gather_inlets(flowsheet, uid: str) -> dict[str, Stream]:
        inlets = {}
        for conn in flowsheet.connections:
            if conn.to_unit == uid and conn.stream_id in flowsheet.streams:
                inlets[conn.to_port] = flowsheet.streams[conn.stream_id]
        return inlets

    def _warm_start(self, flowsheet, layout: _Layout, pp) -> None:
        """One forward sweep (tear streams seeded empty) for a sane initial guess."""
        tears = set(find_tears(flowsheet))
        feed = next((flowsheet.streams[c.stream_id] for c in flowsheet.connections
                     if c.from_unit is None and c.stream_id in flowsheet.streams), None)
        for sid in layout.stream_ids:
            if sid not in flowsheet.streams:
                z = (feed.normalized_z() if feed else
                     {c: 1.0 / layout.nc for c in layout.components})
                flowsheet.streams[sid] = Stream(
                    id=sid, components=list(layout.components),
                    T=feed.T if feed else 300.0, P=feed.P if feed else 101325.0,
                    molar_flow=0.0 if sid in tears else (feed.molar_flow if feed else 1.0),
                    z=z)
        from .sequential import SequentialModularSolver
        order = SequentialModularSolver._topological_order(flowsheet, skip=tears)
        for uid in order:
            inlets = self._gather_inlets(flowsheet, uid)
            try:
                outlets = flowsheet.units[uid].solve(inlets, pp)
            except Exception:
                continue
            self._store(flowsheet, uid, outlets)

    @staticmethod
    def _store(flowsheet, uid, outlets) -> None:
        port_to_stream = {c.from_port: c.stream_id for c in flowsheet.connections
                          if c.from_unit == uid}
        for port, obj in outlets.items():
            sid = port_to_stream.get(port)
            if sid is not None and isinstance(obj, Stream):
                obj.id = sid
                flowsheet.streams[sid] = obj

    def _finalize(self, flowsheet, pp) -> dict[str, float]:
        duties: dict[str, float] = {}
        for uid, unit in flowsheet.units.items():
            inlets = self._gather_inlets(flowsheet, uid)
            outlets = unit.solve(inlets, pp)
            port_to_stream = {c.from_port: c.stream_id for c in flowsheet.connections
                              if c.from_unit == uid}
            for port, obj in outlets.items():
                sid = port_to_stream.get(port)
                if isinstance(obj, EnergyStream):
                    if obj.duty is not None:
                        duties[sid or obj.id] = obj.duty   # default id if unwired
                elif sid is not None:
                    obj.id = sid
                    flowsheet.streams[sid] = obj
        return duties

    @staticmethod
    def _resolve_feed(stream: Stream, pp) -> None:
        stream.validate()
        z = stream.normalized_z()
        res = pp.flash_pt(stream.T, stream.P, z)
        stream.H, stream.phase, stream.vapor_fraction = res.H, res.phase, res.vapor_fraction
