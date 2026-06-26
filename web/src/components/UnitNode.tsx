import { Handle, Position, type NodeProps } from "@xyflow/react";
import { useEffect, useRef, useState } from "react";
import type { CaldyrNode } from "../flow";
import { useStore } from "../store";
import type { Port } from "../types";
import { glyphFor } from "./glyphs";

// One handle per port: inlets on the left, outlets on the right, spaced evenly.
// Energy ports (duties/work) render in amber so they read as utility connections.
function PortHandles({ ports, dir }: { ports: Port[]; dir: "inlet" | "outlet" }) {
  const side = ports.filter((p) => p.direction === dir);
  const position = dir === "inlet" ? Position.Left : Position.Right;
  const type = dir === "inlet" ? "target" : "source";
  return (
    <>
      {side.map((p, i) => (
        <Handle
          key={p.name}
          id={p.name}
          type={type}
          position={position}
          style={{
            top: `${((i + 1) / (side.length + 1)) * 100}%`,
            background: p.kind === "energy" ? "var(--energy)" : "var(--accent)",
            width: 9, height: 9,
          }}
          title={`${p.name} (${p.kind})`}
        />
      ))}
    </>
  );
}

/** Node title that turns into an inline rename input on double-click. */
function NodeTitle({ id, label }: { id: string; label: string }) {
  const renameNode = useStore((s) => s.renameNode);
  const [editing, setEditing] = useState(false);
  const [draft, setDraft] = useState(label);
  const inputRef = useRef<HTMLInputElement>(null);

  useEffect(() => {
    if (editing) inputRef.current?.select();
  }, [editing]);

  if (!editing) {
    return (
      <div
        className="node-title"
        title="Double-click to rename"
        onDoubleClick={(e) => {
          e.stopPropagation();
          setDraft(label);
          setEditing(true);
        }}
      >
        {label}
      </div>
    );
  }
  const finish = (apply: boolean) => {
    setEditing(false);
    if (apply) renameNode(id, draft);
  };
  return (
    <input
      ref={inputRef}
      className="node-rename nodrag"
      value={draft}
      onChange={(e) => setDraft(e.target.value)}
      onBlur={() => finish(true)}
      onKeyDown={(e) => {
        e.stopPropagation();
        if (e.key === "Enter") finish(true);
        if (e.key === "Escape") finish(false);
      }}
      onClick={(e) => e.stopPropagation()}
      aria-label="Rename node"
    />
  );
}

export function UnitNode({ id, data, selected }: NodeProps<CaldyrNode>) {
  const bfd = useStore((s) => s.viewMode) === "bfd";
  const invalid = useStore((s) => s.invalidNodes.includes(id));
  const cls = `node node-${data.kind}${selected ? " selected" : ""}`
    + `${bfd ? " node-bfd" : ""}${invalid ? " node-invalid" : ""}`;
  const subtitle =
    data.kind === "unit" ? data.unitType
    : data.kind === "feed" ? "feed"
    : "product";
  return (
    <div className={cls}>
      <PortHandles ports={data.ports} dir="inlet" />
      {!bfd && <span className="node-glyph">{glyphFor(data.kind, data.unitType)}</span>}
      <span className="node-text">
        <NodeTitle id={id} label={data.label} />
        {!bfd && <div className="node-sub">{subtitle}</div>}
      </span>
      <PortHandles ports={data.ports} dir="outlet" />
    </div>
  );
}
