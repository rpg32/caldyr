// Stable stream numbering shared by the PFD stream flags and the stream table,
// so a given wire shows the same "#" on the canvas and in the table. Numbering
// follows flow order: feed streams first, then downstream by topological rank of
// the source node, tie-broken by the stream id's natural number (creation order).
import type { Edge } from "@xyflow/react";
import type { CaldyrNode } from "../flow";

/** Topological rank of each node (source-to-target). Boundary feeds and any node
 *  with no incoming wire rank 0; nodes on a recycle cycle keep a stable rank from
 *  their remaining in-degree order rather than blowing up. */
function topoRank(nodes: CaldyrNode[], edges: Edge[]): Map<string, number> {
  const indeg = new Map<string, number>(nodes.map((n) => [n.id, 0]));
  const succ = new Map<string, string[]>(nodes.map((n) => [n.id, []]));
  for (const e of edges) {
    if (!indeg.has(e.source) || !indeg.has(e.target)) continue;
    indeg.set(e.target, (indeg.get(e.target) ?? 0) + 1);
    succ.get(e.source)!.push(e.target);
  }
  const rank = new Map<string, number>();
  let queue = nodes.filter((n) => (indeg.get(n.id) ?? 0) === 0).map((n) => n.id);
  let r = 0;
  const seen = new Set<string>();
  while (queue.length) {
    const next: string[] = [];
    for (const id of queue) {
      if (seen.has(id)) continue;
      seen.add(id);
      rank.set(id, r);
      for (const t of succ.get(id) ?? []) {
        indeg.set(t, (indeg.get(t) ?? 1) - 1);
        if ((indeg.get(t) ?? 0) <= 0) next.push(t);
      }
    }
    queue = next;
    r++;
  }
  // Anything left unranked (inside a cycle) sorts after the acyclic part.
  for (const n of nodes) if (!rank.has(n.id)) rank.set(n.id, r);
  return rank;
}

/** Natural (numeric) key of a stream id, e.g. "S12" -> 12; falls back to 0. */
function idNumber(id: string): number {
  const m = id.match(/(\d+)/);
  return m ? Number(m[1]) : 0;
}

/** Map every edge id to a 1-based stream number in stable flow order. */
export function streamNumbers(nodes: CaldyrNode[], edges: Edge[]): Map<string, number> {
  const rank = topoRank(nodes, edges);
  const ordered = [...edges].sort((a, b) => {
    const ra = rank.get(a.source) ?? 0;
    const rb = rank.get(b.source) ?? 0;
    if (ra !== rb) return ra - rb;
    const na = idNumber(a.id);
    const nb = idNumber(b.id);
    if (na !== nb) return na - nb;
    return a.id < b.id ? -1 : a.id > b.id ? 1 : 0;
  });
  const out = new Map<string, number>();
  ordered.forEach((e, i) => out.set(e.id, i + 1));
  return out;
}
