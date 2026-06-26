// Cheap client-side flowsheet checks that catch wiring mistakes before a solve,
// so the user gets a clear, node-targeted message instead of an engine exception.
import type { Edge } from "@xyflow/react";
import type { CaldyrNode } from "../flow";

/** Unit nodes that have at least one material inlet port but no incoming stream
 *  at all — i.e. a unit dropped on the canvas and never fed. */
export function unfedUnits(
  nodes: CaldyrNode[], edges: Edge[],
): { id: string; label: string }[] {
  const fed = new Set(edges.map((e) => e.target));
  return nodes
    .filter((n) => n.data.kind === "unit"
      && n.data.ports.some((p) => p.direction === "inlet" && p.kind === "material")
      && !fed.has(n.id))
    .map((n) => ({ id: n.id, label: n.data.label }));
}
