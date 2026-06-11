// P&ID-view instrumentation overlay edge: a dashed "signal" line with an
// instrument bubble, synthesized from the flowsheet's logical ops (Adjust →
// controller, Set → setpoint link). Not selectable/deletable — edit the ops in
// the flowsheet panel.
import {
  BaseEdge, EdgeLabelRenderer, getBezierPath, type Edge, type EdgeProps,
} from "@xyflow/react";

export interface SignalEdgeData extends Record<string, unknown> {
  tag: string;        // e.g. "AC-1" (adjust) or "SET-1"
  detail: string;     // tooltip text
}
export type SignalEdgeType = Edge<SignalEdgeData, "signal">;

export function SignalEdge({
  id, sourceX, sourceY, targetX, targetY, sourcePosition, targetPosition, data,
}: EdgeProps<SignalEdgeType>) {
  const [path, labelX, labelY] = getBezierPath({
    sourceX, sourceY, targetX, targetY, sourcePosition, targetPosition,
    curvature: 0.4,
  });
  return (
    <>
      <BaseEdge
        id={id}
        path={path}
        style={{ stroke: "var(--warn)", strokeWidth: 1.2, strokeDasharray: "3 4" }}
      />
      <EdgeLabelRenderer>
        <div
          className="instrument-bubble"
          style={{ transform: `translate(-50%, -50%) translate(${labelX}px, ${labelY}px)` }}
          title={data?.detail}
        >
          {data?.tag}
        </div>
      </EdgeLabelRenderer>
    </>
  );
}
