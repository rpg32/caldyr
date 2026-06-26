import type { Edge } from "@xyflow/react";
import { describe, expect, it } from "vitest";
import type { CaldyrNode } from "../flow";
import { unfedUnits } from "./validate";

const unit = (id: string): CaldyrNode => ({
  id, type: "unit", position: { x: 0, y: 0 },
  data: {
    kind: "unit", unitType: "Flash", label: id, params: {},
    ports: [
      { name: "in1", direction: "inlet", kind: "material" },
      { name: "vapor", direction: "outlet", kind: "material" },
    ],
  },
} as unknown as CaldyrNode);

const edge = (source: string, target: string): Edge =>
  ({ id: `${source}-${target}`, source, target } as Edge);

describe("unfedUnits", () => {
  it("flags a unit with a material inlet and no incoming stream", () => {
    const r = unfedUnits([unit("A"), unit("B")], [edge("A", "B")]);
    expect(r.map((u) => u.id)).toEqual(["A"]); // A is unfed; B is fed by A
  });

  it("returns nothing when every unit is fed", () => {
    expect(unfedUnits([unit("B")], [edge("FEED", "B")])).toEqual([]);
  });

  it("ignores feed/product boundary nodes (no material inlet ports)", () => {
    const feed = {
      id: "F", type: "feed", position: { x: 0, y: 0 },
      data: { kind: "feed", label: "F", params: {}, ports: [
        { name: "out", direction: "outlet", kind: "material" }] },
    } as unknown as CaldyrNode;
    expect(unfedUnits([feed], [])).toEqual([]);
  });
});
