import { Handle, Position, type NodeProps } from "@xyflow/react";
import { useEffect, useRef, useState } from "react";
import type { CaldyrNode } from "../flow";
import { useStore } from "../store";
import type { Port } from "../types";
import { humanizeType } from "../lib/format";
import { portAnchors, symbolFor, type Anchor } from "./pfdSymbols";

// BFD view: one handle per port, inlets left / outlets right, spaced evenly.
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

// PFD view: handles pinned to the symbol's physical nozzles (see pfdSymbols). The
// side↔offset is turned into a Position + left/top so edges route orthogonally.
function AnchoredHandles({ ports, anchors }: { ports: Port[]; anchors: Map<string, Anchor> }) {
  return (
    <>
      {ports.map((p) => {
        const a = anchors.get(p.name);
        if (!a) return null;
        const type = p.direction === "inlet" ? "target" : "source";
        const along = `${a.offset * 100}%`;
        const style: React.CSSProperties = {
          background: p.kind === "energy" ? "var(--energy)" : "var(--accent)",
          width: 8, height: 8,
        };
        if (a.position === Position.Left || a.position === Position.Right) style.top = along;
        else style.left = along;
        return (
          <Handle key={p.name} id={p.name} type={type} position={a.position}
            style={style} title={`${p.name} (${p.kind})`} />
        );
      })}
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

function subtitleOf(data: CaldyrNode["data"]): string {
  // Spaced words for CamelCase engine types ("Equilibrium Reactor"); boundary
  // nodes just read "feed" / "product".
  return data.kind === "unit" ? humanizeType(data.unitType ?? "") : data.kind;
}

export function UnitNode({ id, data, selected }: NodeProps<CaldyrNode>) {
  const bfd = useStore((s) => s.viewMode) === "bfd";
  const invalid = useStore((s) => s.invalidNodes.includes(id));

  // BFD keeps the labelled card (correct for block diagrams); PFD/P&ID draw the
  // equipment symbol as the node body with the tag + type below it.
  if (bfd) {
    const cls = `node node-${data.kind} node-bfd${selected ? " selected" : ""}`
      + `${invalid ? " node-invalid" : ""}`;
    return (
      <div className={cls}>
        <PortHandles ports={data.ports} dir="inlet" />
        <span className="node-text">
          <NodeTitle id={id} label={data.label} />
        </span>
        <PortHandles ports={data.ports} dir="outlet" />
      </div>
    );
  }

  const { w, h, body } = symbolFor(data.kind, data.unitType);
  const anchors = portAnchors(data.kind, data.unitType, data.ports);
  const cls = `pfd-node kind-${data.kind}${selected ? " selected" : ""}`
    + `${invalid ? " invalid" : ""}`;
  return (
    <div className={cls}>
      <div className="pfd-symbol" style={{ width: w, height: h }}>
        {body}
        <AnchoredHandles ports={data.ports} anchors={anchors} />
      </div>
      <div className="pfd-label" style={{ maxWidth: Math.max(w, 96) }}>
        <NodeTitle id={id} label={data.label} />
        <div className="node-sub">{subtitleOf(data)}</div>
      </div>
    </div>
  );
}
