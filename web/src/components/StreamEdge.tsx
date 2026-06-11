// Custom edge: smoothstep path, stream-name label, phase/temperature coloring,
// hover tooltip with the solved state, and a pinnable data callout
// (double-click an edge to pin/unpin).
import {
  BaseEdge, EdgeLabelRenderer, getSmoothStepPath, type Edge, type EdgeProps,
} from "@xyflow/react";
import { X } from "lucide-react";
import { useState } from "react";
import { useStore, type ColorMode } from "../store";
import type { StreamState } from "../types";

export interface StreamEdgeData extends Record<string, unknown> {
  state?: StreamState;
  energy?: boolean;
  colorMode: ColorMode;
  tNorm?: number; // 0..1 within the solved temperature range
  pinned?: boolean;
  plain?: boolean; // BFD view: no labels/callouts
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

const fmt = (x: number | null | undefined, d = 1) =>
  x == null ? "—" : x.toLocaleString(undefined, { maximumFractionDigits: d });

function Callout({ id, state, pinned }: { id: string; state: StreamState; pinned: boolean }) {
  const togglePin = useStore((s) => s.togglePin);
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
      <div>{fmt(state.T, 1)} K · {fmt(state.P != null ? state.P / 1000 : null, 0)} kPa</div>
      <div>{fmt(state.molar_flow, 2)} mol/s · {state.phase ?? "?"}
        {state.vapor_fraction != null && state.phase === "VLE"
          ? ` (VF ${state.vapor_fraction.toFixed(2)})` : ""}
      </div>
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
    borderRadius: 6,
  });
  const color = edgeColor(data, selected);
  const showCallout = data?.state && (data.pinned || hovered);

  return (
    <g onMouseEnter={() => setHovered(true)} onMouseLeave={() => setHovered(false)}>
      <BaseEdge
        id={id}
        path={path}
        style={{
          stroke: color,
          strokeWidth: selected || hovered ? 2.2 : 1.5,
          strokeDasharray: data?.energy ? "5 4" : undefined,
        }}
      />
      <EdgeLabelRenderer>
        {!data?.plain && (
          <div
            className={`edge-label${selected ? " selected" : ""}`}
            style={{ transform: `translate(-50%, -50%) translate(${labelX}px, ${labelY}px)` }}
          >
            {String(label ?? id)}
          </div>
        )}
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
