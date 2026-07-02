// Custom edge: smoothstep path, stream-name label, phase/temperature coloring,
// hover tooltip with the solved state, and a pinnable data callout
// (double-click an edge to pin/unpin).
import {
  BaseEdge, EdgeLabelRenderer, getSmoothStepPath, type Edge, type EdgeProps,
} from "@xyflow/react";
import { X } from "lucide-react";
import { useState } from "react";
import { compositionRows, fmtFrac, streamMassFlow } from "../lib/composition";
import { defaultUnit, fmtDim } from "../lib/units";
import { useStore, type ColorMode } from "../store";
import type { StreamState } from "../types";

export interface StreamEdgeData extends Record<string, unknown> {
  state?: StreamState;
  energy?: boolean;
  colorMode: ColorMode;
  tNorm?: number; // 0..1 within the solved temperature range
  pinned?: boolean;
  plain?: boolean; // BFD view: no labels/callouts
  pfd?: boolean;   // PFD/P&ID view: show a numbered stream flag instead of a name
  streamNo?: number; // stable stream number (see lib/streamNumbers)
}
export type StreamEdgeType = Edge<StreamEdgeData, "stream">;

const PHASE_COLORS: Record<string, string> = {
  vapor: "#f97316",
  liquid: "#3b82f6",
  VLE: "#a855f7",
  supercritical: "#14b8a6",
};

/** Blue -> red lerp for the temperature heat-map mode. */
function tempColor(t: number): string {
  const a = [59, 130, 246];   // blue-500
  const b = [239, 68, 68];    // red-500
  const c = a.map((v, i) => Math.round(v + (b[i] - v) * Math.max(0, Math.min(1, t))));
  return `rgb(${c[0]},${c[1]},${c[2]})`;
}

function edgeColor(data: StreamEdgeData | undefined, selected: boolean | undefined): string {
  if (data?.energy) return "var(--energy)";
  if (data?.colorMode === "phase" && data.state?.phase) {
    return PHASE_COLORS[data.state.phase] ?? "var(--line)";
  }
  if (data?.colorMode === "temperature" && data.tNorm !== undefined) {
    return tempColor(data.tNorm);
  }
  return selected ? "var(--accent)" : "var(--edge)";
}

function Callout({ id, state, pinned }: { id: string; state: StreamState; pinned: boolean }) {
  const togglePin = useStore((s) => s.togglePin);
  const unitSet = useStore((s) => s.unitSet);
  const mw = useStore((s) => s.solveRes?.molar_mass);
  const massFlow = streamMassFlow(state.z, state.molar_flow, mw);
  return (
    <div className="edge-callout nodrag nopan">
      <div className="edge-callout-head">
        <b>{id}</b>
        {pinned && (
          <button onClick={() => togglePin(id)} aria-label={`Unpin ${id}`} title="Unpin">
            <X size={11} />
          </button>
        )}
      </div>
      <div>{fmtDim("temperature", state.T, unitSet, 1)} {defaultUnit("temperature", unitSet)}
        {" · "}{fmtDim("pressure", state.P, unitSet, 1)} {defaultUnit("pressure", unitSet)}</div>
      <div>{fmtDim("molar_flow", state.molar_flow, unitSet, 2)} {defaultUnit("molar_flow", unitSet)}
        {" · "}{state.phase ?? "?"}
        {state.vapor_fraction != null && state.phase === "VLE"
          ? ` (VF ${state.vapor_fraction.toFixed(2)})` : ""}
      </div>
      {massFlow != null && (
        <div>{fmtDim("mass_flow", massFlow, unitSet, 2)} {defaultUnit("mass_flow", unitSet)}</div>
      )}
      <CompositionLines state={state} />
    </div>
  );
}

/** Top components by mole fraction in the pinned/hover callout. */
function CompositionLines({ state }: { state: StreamState }) {
  const rows = compositionRows(state.z, state.molar_flow);
  if (!rows.length) return null;
  const shown = rows.slice(0, 4);
  const rest = rows.length - shown.length;
  return (
    <div className="edge-callout-comp">
      {shown.map((r) => (
        <div key={r.comp}>{r.comp} {fmtFrac(r.frac)}</div>
      ))}
      {rest > 0 && <div>+{rest} more</div>}
    </div>
  );
}

export function StreamEdge({
  id, sourceX, sourceY, targetX, targetY, sourcePosition, targetPosition,
  data, selected, label,
}: EdgeProps<StreamEdgeType>) {
  const [hovered, setHovered] = useState(false);
  const [path, labelX, labelY] = getSmoothStepPath({
    sourceX, sourceY, targetX, targetY, sourcePosition, targetPosition,
    borderRadius: 2, // near-orthogonal corners for a classic PFD look
  });
  const color = edgeColor(data, selected);
  const showCallout = data?.state && (data.pinned || hovered);
  // A directional arrowhead that follows the edge colour. Colours vary per edge
  // (phase/temperature map to concrete rgb; others to a theme CSS var), so each
  // edge carries its own marker rather than sharing a fixed palette.
  const markerId = `stream-arrow-${id}`;
  // Material streams in PFD/P&ID get a numbered diamond flag; the stream name
  // moves to the hover/pin callout. Energy streams keep the plain text label.
  const flagged = data?.pfd && !data?.energy && data?.streamNo != null;

  return (
    <g onMouseEnter={() => setHovered(true)} onMouseLeave={() => setHovered(false)}>
      <marker id={markerId} viewBox="0 0 10 10" refX="8.5" refY="5"
        markerWidth="6" markerHeight="6" orient="auto">
        <path d="M0 0 L10 5 L0 10 z" style={{ fill: color }} />
      </marker>
      <BaseEdge
        id={id}
        path={path}
        markerEnd={`url(#${markerId})`}
        style={{
          stroke: color,
          strokeWidth: selected || hovered ? 2.2 : 1.5,
          strokeDasharray: data?.energy ? "5 4" : undefined,
        }}
      />
      <EdgeLabelRenderer>
        {!data?.plain && (flagged ? (
          <div
            className={`stream-flag${selected ? " selected" : ""}`}
            style={{ transform: `translate(-50%, -50%) translate(${labelX}px, ${labelY}px)` }}
            title={String(label ?? id)}
          >
            <span>{data!.streamNo}</span>
          </div>
        ) : (
          <div
            className={`edge-label${selected ? " selected" : ""}`}
            style={{ transform: `translate(-50%, -50%) translate(${labelX}px, ${labelY}px)` }}
          >
            {String(label ?? id)}
          </div>
        ))}
        {!data?.plain && showCallout && (
          <div
            className="edge-callout-anchor"
            style={{ transform: `translate(-50%, 0) translate(${labelX}px, ${labelY + 10}px)` }}
          >
            <Callout id={id} state={data!.state!} pinned={!!data!.pinned} />
          </div>
        )}
      </EdgeLabelRenderer>
    </g>
  );
}
