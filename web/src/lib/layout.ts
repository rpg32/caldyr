// Auto-layout via ELK's layered algorithm (left-to-right, like a PFD reads).
// ELK is ~1.4 MB, so it loads lazily on first use.
import type { Edge } from "@xyflow/react";
import type { CaldyrNode } from "../flow";

let elkPromise: Promise<InstanceType<typeof import("elkjs").default>> | null = null;
const getElk = () => {
  elkPromise ??= import("elkjs/lib/elk.bundled.js").then((m) => new m.default());
  return elkPromise;
};

const OPTIONS: Record<string, string> = {
  "elk.algorithm": "layered",
  "elk.direction": "RIGHT",
  "elk.layered.spacing.nodeNodeBetweenLayers": "110",
  "elk.spacing.nodeNode": "70",
  "elk.layered.nodePlacement.strategy": "NETWORK_SIMPLEX",
  "elk.layered.considerModelOrder.strategy": "NODES_AND_EDGES",
  // Break recycle cycles at the edge that runs against flowsheet order (the
  // tear), not somewhere on the main forward path.
  "elk.layered.cycleBreaking.strategy": "MODEL_ORDER",
};

export async function autoLayout(
  nodes: CaldyrNode[],
  edges: Edge[],
): Promise<Map<string, { x: number; y: number }>> {
  const graph = {
    id: "root",
    layoutOptions: OPTIONS,
    children: nodes.map((n) => ({
      id: n.id,
      width: n.measured?.width ?? 140,
      height: n.measured?.height ?? 54,
    })),
    edges: edges.map((e) => ({ id: e.id, sources: [e.source], targets: [e.target] })),
  };
  const elk = await getElk();
  const res = await elk.layout(graph);
  const out = new Map<string, { x: number; y: number }>();
  for (const child of res.children ?? []) {
    out.set(child.id, { x: Math.round(child.x ?? 0), y: Math.round(child.y ?? 0) });
  }
  return out;
}
