"""Recycle-convergence support for the sequential-modular solver.

Tear-stream selection (DFS back edges) + a per-variable Wegstein accelerator over
direct substitution. The tear "state vector" for a stream is its per-component
molar flows followed by temperature; pressure is carried by direct substitution
(it is fixed by the producing unit and converges in one sweep).
"""
from __future__ import annotations

from ..core import Stream

# Wegstein acceleration factor is clipped to this range. q = 0 is pure direct
# substitution; negative q accelerates a monotonic approach; the lower bound
# caps the step so a noisy slope estimate cannot destabilize the iteration.
Q_MIN, Q_MAX = -5.0, 0.0
_FLOOR = 1e-9


def find_tears(flowsheet) -> list[str]:
    """Return a set of stream ids whose removal makes the unit graph acyclic.

    DFS over units; an edge to a unit currently on the recursion stack is a back
    edge, and the stream it carries is torn. Not guaranteed minimal, but always
    valid (the torn graph is a DAG). Deterministic in unit insertion order.
    """
    succ: dict[str, list[tuple[str, str]]] = {u: [] for u in flowsheet.units}
    for conn in flowsheet.connections:
        if conn.from_unit in succ and conn.to_unit in succ:
            succ[conn.from_unit].append((conn.to_unit, conn.stream_id))

    WHITE, GRAY, BLACK = 0, 1, 2
    color = {u: WHITE for u in flowsheet.units}
    tears: list[str] = []
    tear_set: set[str] = set()

    def visit(u: str) -> None:
        color[u] = GRAY
        for v, sid in succ[u]:
            if sid in tear_set:
                continue
            if color[v] == GRAY:                 # back edge -> tear this stream
                if sid not in tear_set:
                    tear_set.add(sid)
                    tears.append(sid)
            elif color[v] == WHITE:
                visit(v)
        color[u] = BLACK

    for u in flowsheet.units:
        if color[u] == WHITE:
            visit(u)
    return tears


def tear_vector(stream: Stream, components: list[str]) -> list[float]:
    """[n_c for each component] + [T] for a tear stream."""
    n = stream.molar_flow or 0.0
    total = sum(stream.z.values())
    z = {k: v / total for k, v in stream.z.items()} if total > 0 else {}
    T = stream.T if stream.T is not None else 298.15
    return [n * z.get(c, 0.0) for c in components] + [T]


def stream_from_vector(
    sid: str, vec: list[float], components: list[str], P: float, fallback_z: dict[str, float]
) -> Stream:
    """Rebuild a tear-stream guess from an updated state vector. Enthalpy is left
    unset so the consuming unit re-resolves it from (T, P, z)."""
    *flows, T = vec
    total = sum(flows)
    if total > 1e-12:
        z = {c: f / total for c, f in zip(components, flows)}
    else:
        total, z = 0.0, dict(fallback_z)
    return Stream(id=sid, components=list(components), T=T, P=P, molar_flow=total, z=z)


def residual(guess: dict[str, list[float]], computed: dict[str, list[float]]) -> float:
    """Max relative change across every tear variable, normalized by the larger
    of the two magnitudes so a 0 -> finite first step reads as ~1, not ~1e9."""
    worst = 0.0
    for tid in guess:
        for a, b in zip(computed[tid], guess[tid]):
            worst = max(worst, abs(a - b) / (max(abs(a), abs(b)) + _FLOOR))
    return worst


def wegstein_update(
    x_prev: list[float], g_prev: list[float], x: list[float], g: list[float]
) -> list[float]:
    """One per-variable Wegstein step. x is the current guess, g = f(x) the swept
    result; *_prev are the previous iteration's pair."""
    out = []
    for xp, gp, xk, gk in zip(x_prev, g_prev, x, g):
        dx = xk - xp
        if abs(dx) < _FLOOR:
            q = 0.0                              # no slope info -> direct substitution
        else:
            s = (gk - gp) / dx
            denom = s - 1.0
            q = 0.0 if abs(denom) < _FLOOR else max(Q_MIN, min(Q_MAX, s / denom))
        out.append(q * xk + (1.0 - q) * gk)
    return out
