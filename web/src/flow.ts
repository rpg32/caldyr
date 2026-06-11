// Canvas <-> .flow document conversion. This is what makes flowsheet-as-canvas
// and flowsheet-as-code the same object (the engine's io schema). UI-only state
// (boundary node positions, product, backend) rides in `meta.ui`, which the
// engine ignores on load.
import type { Edge, Node } from "@xyflow/react";
import type { UiMeta } from "./lib/persist";
import type { FlowDoc, NodeData, Port, UnitType } from "./types";

interface FlowUnit {
  id: string;
  type: string;
  params: Record<string, unknown>;
  xy?: [number, number];
}
interface FlowStream {
  id: string;
  from: string | null;
  to: string | null;
  spec?: Record<string, unknown>;
}

export type CaldyrNode = Node<NodeData>;

const firstPort = (ports: Port[], dir: "inlet" | "outlet"): string =>
  ports.find((p) => p.direction === dir)?.name ?? (dir === "inlet" ? "in1" : "out");

const feedSpec = (n: CaldyrNode, components: string[]) => {
  const p = n.data.params as Record<string, unknown>;
  // Drop mole fractions for components no longer in the flowsheet.
  const z = Object.fromEntries(
    Object.entries((p.z as Record<string, number>) ?? {})
      .filter(([c]) => components.includes(c)),
  );
  return { T: p.T, P: p.P, molar_flow: p.molar_flow, z };
};

/** Serialize the canvas into a .flow document the engine can solve. */
export interface FlowExtras {
  logical?: Record<string, unknown>[];
  solver_hints?: Record<string, unknown>;
}

export function canvasToFlow(
  nodes: CaldyrNode[],
  edges: Edge[],
  components: string[],
  propertyPackage: string,
  ui?: UiMeta,
  extras?: FlowExtras,
): FlowDoc {
  const byId = new Map(nodes.map((n) => [n.id, n]));
  const streams: FlowStream[] = [];
  const usedOutlets = new Set<string>();
  const boundaryXy: Record<string, [number, number]> = {};

  for (const e of edges) {
    const src = byId.get(e.source);
    const tgt = byId.get(e.target);
    if (!src || !tgt) continue;
    const fromPort = e.sourceHandle ?? firstPort(src.data.ports, "outlet");
    const toPort = e.targetHandle ?? firstPort(tgt.data.ports, "inlet");

    const stream: FlowStream = { id: e.id, from: null, to: null };
    if (src.data.kind === "feed") {
      stream.spec = feedSpec(src, components);
      boundaryXy[e.id] = [Math.round(src.position.x), Math.round(src.position.y)];
    } else {
      stream.from = `${src.id}:${fromPort}`;
      usedOutlets.add(stream.from);
    }
    if (tgt.data.kind === "product") {
      stream.to = null;
      boundaryXy[e.id] = [Math.round(tgt.position.x), Math.round(tgt.position.y)];
    } else {
      stream.to = `${tgt.id}:${toPort}`;
    }
    streams.push(stream);
  }

  const units: FlowUnit[] = [];
  for (const n of nodes) {
    if (n.data.kind !== "unit") continue;
    units.push({
      id: n.id,
      type: n.data.unitType as string,
      params: n.data.params,
      xy: [Math.round(n.position.x), Math.round(n.position.y)],
    });
    // Any unconnected outlet (product or duty) becomes a boundary stream so the
    // engine reports it.
    for (const p of n.data.ports) {
      if (p.direction !== "outlet") continue;
      const key = `${n.id}:${p.name}`;
      if (!usedOutlets.has(key)) streams.push({ id: `${n.id}_${p.name}`, from: key, to: null });
    }
  }

  const doc: FlowDoc = {
    schema: "caldyr.flow/1",
    meta: { ui: { ...ui, boundary_xy: boundaryXy } },
    components: components.map((id) => ({ id })),
    property_package: propertyPackage,
    units,
    streams,
  };
  if (extras?.logical?.length) doc.logical = extras.logical;
  if (extras?.solver_hints && Object.keys(extras.solver_hints).length) {
    doc.solver_hints = extras.solver_hints;
  }
  return doc;
}

export interface CanvasState {
  nodes: CaldyrNode[];
  edges: Edge[];
  components: string[];
  propertyPackage: string;
  ui: UiMeta;
  extras: FlowExtras;
}

/** Build canvas nodes/edges from a .flow document (loading an example/file). */
export function flowToCanvas(flow: FlowDoc, unitTypes: UnitType[]): CanvasState {
  const portsOf = new Map(unitTypes.map((u) => [u.type, u.ports]));
  const units = (flow.units ?? []) as FlowUnit[];
  const streams = (flow.streams ?? []) as FlowStream[];
  const components = ((flow.components ?? []) as { id: string }[]).map((c) => c.id);
  const meta = (flow.meta ?? {}) as { ui?: UiMeta };
  const ui: UiMeta = meta.ui ?? {};
  const boundaryXy = ui.boundary_xy ?? {};

  const nodes: CaldyrNode[] = [];
  const edges: Edge[] = [];

  units.forEach((u, i) => {
    nodes.push({
      id: u.id,
      type: "caldyr",
      position: { x: u.xy?.[0] ?? 80 + i * 200, y: u.xy?.[1] ?? 120 },
      data: {
        kind: "unit", label: u.id, unitType: u.type,
        ports: portsOf.get(u.type) ?? [], params: u.params ?? {},
      },
    });
  });

  const portKind = (endpoint: string | null): string | undefined => {
    if (!endpoint) return undefined;
    const [uid, port] = endpoint.split(":");
    return portsOf.get(units.find((u) => u.id === uid)?.type ?? "")
      ?.find((p) => p.name === port)?.kind;
  };

  let feedN = 0;
  let prodN = 0;
  for (const s of streams) {
    let source = s.from?.split(":")[0];
    let sourceHandle = s.from?.split(":")[1];
    let target = s.to?.split(":")[0];
    let targetHandle = s.to?.split(":")[1];

    if (s.from === null && s.spec) {
      const id = `FEED_${s.id}`;
      const anchor = nodes.find((n) => n.id === target);
      const xy = boundaryXy[s.id];
      nodes.push({
        id, type: "caldyr",
        position: xy
          ? { x: xy[0], y: xy[1] }
          : { x: (anchor?.position.x ?? 200) - 180, y: (anchor?.position.y ?? 120) + feedN++ * 90 },
        data: {
          kind: "feed", label: s.id,
          ports: [{ name: "out", direction: "outlet", kind: "material" }],
          params: s.spec,
        },
      });
      source = id;
      sourceHandle = "out";
    }
    if (s.to === null) {
      if (portKind(s.from) === "energy") continue; // duty sinks are implied
      const id = `PROD_${s.id}`;
      const anchor = nodes.find((n) => n.id === source);
      const xy = boundaryXy[s.id];
      nodes.push({
        id, type: "caldyr",
        position: xy
          ? { x: xy[0], y: xy[1] }
          : { x: (anchor?.position.x ?? 200) + 200, y: (anchor?.position.y ?? 120) + prodN++ * 90 },
        data: {
          kind: "product", label: s.id,
          ports: [{ name: "in", direction: "inlet", kind: "material" }], params: {},
        },
      });
      target = id;
      targetHandle = "in";
    }
    if (source && target) {
      edges.push({ id: s.id, source, target, sourceHandle, targetHandle, label: s.id });
    }
  }

  return {
    nodes, edges, components,
    propertyPackage: (flow.property_package as string) ?? "thermo:PR",
    ui,
    extras: {
      logical: (flow.logical as Record<string, unknown>[] | undefined) ?? [],
      solver_hints: (flow.solver_hints as Record<string, unknown> | undefined) ?? {},
    },
  };
}
