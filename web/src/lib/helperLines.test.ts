import { describe, expect, it } from "vitest";
import type { NodePositionChange } from "@xyflow/react";
import type { CaldyrNode } from "../flow";
import { applyHelperLines } from "./helperLines";

const node = (id: string, x: number, y: number): CaldyrNode => ({
  id,
  type: "caldyr",
  position: { x, y },
  measured: { width: 140, height: 54 },
  data: { kind: "unit", label: id, unitType: "Mixer", ports: [], params: {} },
});

const change = (id: string, x: number, y: number): NodePositionChange => ({
  id, type: "position", position: { x, y }, dragging: true,
});

describe("applyHelperLines", () => {
  it("snaps to another node's left edge when within threshold", () => {
    const nodes = [node("A", 100, 100), node("B", 400, 300)];
    const c = change("A", 403, 100); // 3px off B's left edge
    const lines = applyHelperLines(c, nodes);
    expect(c.position!.x).toBe(400);
    expect(lines.vertical).toBe(400);
  });

  it("snaps to another node's top edge (same-size nodes)", () => {
    const nodes = [node("A", 100, 100), node("B", 400, 300)];
    const c = change("A", 100, 305); // 5px below B's top
    const lines = applyHelperLines(c, nodes);
    expect(c.position!.y).toBe(300);
    expect(lines.horizontal).toBe(300);
  });

  it("snaps centers when sizes differ", () => {
    const small: CaldyrNode = { ...node("A", 100, 100), measured: { width: 140, height: 40 } };
    const nodes = [small, node("B", 400, 300)];
    // B center y = 327; A at y=305 → center 325 (2px off, closer than any edge)
    const c = change("A", 100, 305);
    const lines = applyHelperLines(c, nodes);
    expect(c.position!.y).toBe(307); // 327 - 40/2
    expect(lines.horizontal).toBe(327);
  });

  it("does nothing when out of range", () => {
    const nodes = [node("A", 100, 100), node("B", 400, 300)];
    const c = change("A", 200, 150);
    const lines = applyHelperLines(c, nodes);
    expect(c.position).toEqual({ x: 200, y: 150 });
    expect(lines.vertical).toBeUndefined();
    expect(lines.horizontal).toBeUndefined();
  });
});
