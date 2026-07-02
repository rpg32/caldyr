import {
  Background, Controls, MiniMap, Panel, ReactFlow, SelectionMode, ViewportPortal,
  useReactFlow,
  type Edge, type EdgeTypes, type Node, type NodeTypes, type OnSelectionChangeParams,
} from "@xyflow/react";
import { useCallback, useEffect, useMemo } from "react";
import type { CaldyrNode } from "../flow";
import { streamNumbers } from "../lib/streamNumbers";
import { useStore } from "../store";
import { ChatPanel } from "./ChatPanel";
import { GroupNode, type GroupNodeType } from "./GroupNode";
import { SignalEdge, type SignalEdgeType } from "./SignalEdge";
import { StreamEdge, type StreamEdgeType } from "./StreamEdge";
import { UnitNode } from "./UnitNode";

const nodeTypes: NodeTypes = { caldyr: UnitNode, caldyrGroup: GroupNode };
const edgeTypes: EdgeTypes = { stream: StreamEdge, signal: SignalEdge };

const MINIMAP_COLORS: Record<string, string> = {
  feed: "#22c55e",
  product: "#94a3b8",
  unit: "#475569",
};

/** Fit the view whenever the store bumps fitNonce (e.g. after auto-layout). */
function FitWatcher() {
  const fitNonce = useStore((s) => s.fitNonce);
  const { fitView } = useReactFlow();
  useEffect(() => {
    if (fitNonce > 0) void fitView({ duration: 350, padding: 0.15 });
  }, [fitNonce, fitView]);
  return null;
}

/** Alignment guide lines, drawn in flow coordinates while dragging. */
function HelperLinesOverlay() {
  const lines = useStore((s) => s.helperLines);
  if (lines.horizontal === undefined && lines.vertical === undefined) return null;
  return (
    <ViewportPortal>
      {lines.vertical !== undefined && (
        <div className="helper-line" style={{
          left: lines.vertical, top: -10000, width: 1, height: 20000,
        }} />
      )}
      {lines.horizontal !== undefined && (
        <div className="helper-line" style={{
          top: lines.horizontal, left: -10000, height: 1, width: 20000,
        }} />
      )}
    </ViewportPortal>
  );
}

function Legend() {
  const colorMode = useStore((s) => s.colorMode);
  const solveRes = useStore((s) => s.solveRes);
  if (colorMode === "none") return null;
  return (
    <Panel position="bottom-left" className="legend">
      {colorMode === "phase" ? (
        <>
          <span><i style={{ background: "#f97316" }} /> vapor</span>
          <span><i style={{ background: "#3b82f6" }} /> liquid</span>
          <span><i style={{ background: "#a855f7" }} /> two-phase</span>
          <span><i style={{ background: "var(--energy)" }} /> duty</span>
        </>
      ) : (
        <>
          <span><i style={{ background: "#3b82f6" }} /> cold</span>
          <span className="legend-grad" />
          <span><i style={{ background: "#ef4444" }} /> hot</span>
        </>
      )}
      {!solveRes && <span className="legend-note">solve to color streams</span>}
    </Panel>
  );
}

export function CanvasView() {
  const nodes = useStore((s) => s.nodes);
  const edges = useStore((s) => s.edges);
  const solveRes = useStore((s) => s.solveRes);
  const colorMode = useStore((s) => s.colorMode);
  const pinnedStreams = useStore((s) => s.pinnedStreams);
  const viewMode = useStore((s) => s.viewMode);
  const groups = useStore((s) => s.groups);
  const logical = useStore((s) => s.logical);
  const onNodesChange = useStore((s) => s.onNodesChange);
  const onEdgesChange = useStore((s) => s.onEdgesChange);
  const onConnect = useStore((s) => s.onConnect);
  const setSelection = useStore((s) => s.setSelection);
  const setTab = useStore((s) => s.setTab);
  const togglePin = useStore((s) => s.togglePin);
  const theme = useStore((s) => s.theme);

  // Which edges are duty/work connections (either endpoint is an energy port)?
  const energyEdges = useMemo(() => {
    const byId = new Map(nodes.map((n) => [n.id, n]));
    const kindOf = (node: CaldyrNode | undefined, handle: string | null | undefined) =>
      node?.data.ports.find((p) => p.name === handle)?.kind;
    return new Set(edges
      .filter((e) =>
        kindOf(byId.get(e.source), e.sourceHandle) === "energy"
        || kindOf(byId.get(e.target), e.targetHandle) === "energy")
      .map((e) => e.id));
  }, [nodes, edges]);

  // Solved temperature range for the heat-map mode.
  const [tMin, tMax] = useMemo(() => {
    const ts = Object.values(solveRes?.streams ?? {})
      .map((s) => s.T).filter((t): t is number => t != null);
    return ts.length ? [Math.min(...ts), Math.max(...ts)] : [0, 1];
  }, [solveRes]);

  // collapsed-group bookkeeping: member id -> group id
  const collapsedOwner = useMemo(() => {
    const m = new Map<string, string>();
    for (const g of groups) {
      if (g.collapsed) for (const id of g.members) m.set(id, g.id);
    }
    return m;
  }, [groups]);

  const renderNodes: Node[] = useMemo(() => {
    const out: Node[] = nodes.map((n) =>
      collapsedOwner.has(n.id) ? { ...n, hidden: true } : n);
    for (const g of groups) {
      if (!g.collapsed) continue;
      out.push({
        id: g.id,
        type: "caldyrGroup",
        position: { x: g.xy[0], y: g.xy[1] },
        data: { label: g.label, count: g.members.length },
      } satisfies GroupNodeType);
    }
    return out;
  }, [nodes, groups, collapsedOwner]);

  // Stable stream numbers for the PFD flags (shared with the stream table).
  const streamNos = useMemo(() => streamNumbers(nodes, edges), [nodes, edges]);

  const styledEdges: Edge[] = useMemo(() => {
    const bfd = viewMode === "bfd";
    const out: Edge[] = [];
    for (const e of edges) {
      const energy = energyEdges.has(e.id);
      if (bfd && energy) continue;                      // BFD hides duty lines
      const srcGrp = collapsedOwner.get(e.source);
      const tgtGrp = collapsedOwner.get(e.target);
      if (srcGrp && tgtGrp) continue;                   // internal to collapsed group(s)
      const state = solveRes?.streams[e.id];
      out.push({
        ...e,
        // reroute boundary edges of collapsed groups onto the group block
        source: srcGrp ?? e.source,
        sourceHandle: srcGrp ? "out" : e.sourceHandle,
        target: tgtGrp ?? e.target,
        targetHandle: tgtGrp ? "in" : e.targetHandle,
        type: "stream" as const,
        data: {
          state,
          energy,
          colorMode,
          tNorm: state?.T != null && tMax > tMin
            ? (state.T - tMin) / (tMax - tMin) : undefined,
          pinned: pinnedStreams.includes(e.id),
          plain: bfd,
          pfd: !bfd,
          streamNo: streamNos.get(e.id),
        },
      } satisfies StreamEdgeType);
    }
    // P&ID view: instrumentation overlay from logical ops
    if (viewMode === "pid") {
      const byId = new Map(nodes.map((n) => [n.id, n]));
      const edgeById = new Map(edges.map((e) => [e.id, e]));
      let ac = 0;
      let st = 0;
      for (const op of logical) {
        let from: string | undefined;
        let to: string | undefined;
        let tag = "";
        let detail = "";
        if (op.type === "adjust") {
          const vary = op.vary as [string, string] | undefined;
          const spec = op.spec as { stream?: string; type?: string } | undefined;
          // anchor at the measured stream's source unit (its target may be a
          // product node, which has no source handle to hang an edge on)
          const measured = edgeById.get(spec?.stream ?? "");
          from = measured ? measured.source : undefined;
          to = vary?.[0];
          tag = `AC-${++ac}`;
          detail = `Adjust ${vary?.join(".")} until ${spec?.type} of ${spec?.stream} = ${op.value}`;
        } else if (op.type === "set") {
          const src = op.source as [string, string] | undefined;
          const tgt = op.target as [string, string] | undefined;
          from = src?.[0];
          to = tgt?.[0];
          tag = `SET-${++st}`;
          detail = `${tgt?.join(".")} = ${op.multiplier ?? 1} × ${src?.join(".")} + ${op.offset ?? 0}`;
        }
        if (!from || !to || !byId.has(from) || !byId.has(to)) continue;
        out.push({
          id: `__sig_${tag}`,
          source: collapsedOwner.get(from) ?? from,
          target: collapsedOwner.get(to) ?? to,
          type: "signal" as const,
          selectable: false,
          deletable: false,
          data: { tag, detail },
        } satisfies SignalEdgeType);
      }
    }
    return out;
  }, [edges, nodes, solveRes, energyEdges, colorMode, tMin, tMax, pinnedStreams,
      viewMode, logical, collapsedOwner, streamNos]);

  const onSelectionChange = useCallback(
    ({ nodes: ns, edges: es }: OnSelectionChangeParams) => {
      if (ns.length) setSelection({ kind: "node", id: ns[0].id });
      else if (es.length) setSelection({ kind: "edge", id: es[0].id });
      else setSelection(null);
    },
    [setSelection],
  );

  return (
    <main className="relative min-w-0">
      <ChatPanel />
      <ReactFlow
        nodes={renderNodes}
        edges={styledEdges}
        onNodesChange={onNodesChange}
        onEdgesChange={onEdgesChange}
        onConnect={onConnect}
        onSelectionChange={onSelectionChange}
        onNodeClick={() => setTab("params")}
        onEdgeClick={() => setTab("params")}
        onEdgeDoubleClick={(_, e) => togglePin(e.id)}
        nodeTypes={nodeTypes}
        edgeTypes={edgeTypes}
        deleteKeyCode={["Backspace", "Delete"]}
        colorMode={theme}
        snapToGrid
        snapGrid={[8, 8]}
        selectionMode={SelectionMode.Partial}
        fitView
      >
        <Background gap={16} />
        <Controls />
        <MiniMap pannable zoomable
          nodeColor={(n) => MINIMAP_COLORS[(n as CaldyrNode).data.kind] ?? "#475569"} />
        <FitWatcher />
        <HelperLinesOverlay />
        <Legend />
      </ReactFlow>
    </main>
  );
}
