// Alignment helper lines: while dragging a node, snap it to other nodes'
// edges/centers when within a small threshold, and report the guide lines to
// draw. Adapted from the xyflow helper-lines example, trimmed to our needs.
import type { NodePositionChange } from "@xyflow/react";
import type { CaldyrNode } from "../flow";

export interface HelperLines {
  horizontal?: number; // y in flow coordinates
  vertical?: number;   // x in flow coordinates
}

const THRESHOLD = 6;

interface Bounds { left: number; right: number; top: number; bottom: number; cx: number; cy: number }

const boundsOf = (x: number, y: number, w: number, h: number): Bounds => ({
  left: x, right: x + w, top: y, bottom: y + h, cx: x + w / 2, cy: y + h / 2,
});

const sizeOf = (n: CaldyrNode) => ({
  w: n.measured?.width ?? 140,
  h: n.measured?.height ?? 54,
});

/** Mutates `change.position` to snap, returns the lines to render. */
export function applyHelperLines(
  change: NodePositionChange,
  nodes: CaldyrNode[],
): HelperLines {
  const lines: HelperLines = {};
  if (!change.position) return lines;
  const moving = nodes.find((n) => n.id === change.id);
  if (!moving) return lines;

  const { w, h } = sizeOf(moving);
  const a = boundsOf(change.position.x, change.position.y, w, h);
  let bestV = THRESHOLD;
  let bestH = THRESHOLD;

  for (const other of nodes) {
    if (other.id === change.id) continue;
    const { w: ow, h: oh } = sizeOf(other);
    const b = boundsOf(other.position.x, other.position.y, ow, oh);

    // candidate vertical alignments: (my edge/center) vs (their edge/center)
    const vPairs: [number, number, number][] = [
      [a.left, b.left, b.left],
      [a.left, b.right, b.right],
      [a.right, b.left, b.left - w],
      [a.right, b.right, b.right - w],
      [a.cx, b.cx, b.cx - w / 2],
    ];
    for (const [mine, theirs, snapX] of vPairs) {
      const d = Math.abs(mine - theirs);
      if (d < bestV) {
        bestV = d;
        change.position.x = snapX;
        lines.vertical = theirs;
      }
    }

    const hPairs: [number, number, number][] = [
      [a.top, b.top, b.top],
      [a.top, b.bottom, b.bottom],
      [a.bottom, b.top, b.top - h],
      [a.bottom, b.bottom, b.bottom - h],
      [a.cy, b.cy, b.cy - h / 2],
    ];
    for (const [mine, theirs, snapY] of hPairs) {
      const d = Math.abs(mine - theirs);
      if (d < bestH) {
        bestH = d;
        change.position.y = snapY;
        lines.horizontal = theirs;
      }
    }
  }
  return lines;
}
