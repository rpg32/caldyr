// A collapsed sub-flowsheet block: stands in for its member nodes on the
// canvas. Double-click to expand; the flowsheet panel lists groups too.
import { Handle, Position, type Node, type NodeProps } from "@xyflow/react";
import { Boxes } from "lucide-react";
import { useStore } from "../store";

export interface GroupNodeData extends Record<string, unknown> {
  label: string;
  count: number;
}
export type GroupNodeType = Node<GroupNodeData, "caldyrGroup">;

export function GroupNode({ id, data, selected }: NodeProps<GroupNodeType>) {
  const toggle = useStore((s) => s.toggleGroupCollapse);
  return (
    <div
      className={`node node-group${selected ? " selected" : ""}`}
      onDoubleClick={(e) => {
        e.stopPropagation();
        toggle(id);
      }}
      title="Double-click to expand"
    >
      <Handle id="in" type="target" position={Position.Left}
        style={{ background: "var(--accent)", width: 9, height: 9 }} />
      <span className="node-glyph"><Boxes size={24} aria-hidden /></span>
      <span className="node-text">
        <div className="node-title">{data.label}</div>
        <div className="node-sub">{data.count} units</div>
      </span>
      <Handle id="out" type="source" position={Position.Right}
        style={{ background: "var(--accent)", width: 9, height: 9 }} />
    </div>
  );
}
