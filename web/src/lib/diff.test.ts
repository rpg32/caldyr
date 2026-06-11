import { describe, expect, it } from "vitest";
import type { FlowDoc } from "../types";
import { diffFlows, mergePositions } from "./diff";

const base = (): FlowDoc => ({
  schema: "caldyr.flow/1",
  components: [{ id: "a" }, { id: "b" }],
  property_package: "thermo:PR",
  units: [
    { id: "H", type: "Heater", params: { T_out: 350 }, xy: [100, 50] },
    { id: "F", type: "Flash", params: { T: 300, P: 1e5 }, xy: [300, 50] },
  ],
  streams: [
    { id: "S1", from: null, to: "H:in1", spec: { T: 298, P: 1e5, molar_flow: 1, z: { a: 1 } } },
    { id: "S2", from: "H:out", to: "F:in1" },
  ],
});

describe("diffFlows", () => {
  it("reports no changes for identical documents", () => {
    expect(diffFlows(base(), base()).count).toBe(0);
  });

  it("detects added units, param changes, and new streams", () => {
    const next = base();
    (next.units as Record<string, unknown>[]).push({ id: "P", type: "Pump", params: {} });
    (next.units as { params: Record<string, unknown> }[])[0].params = { T_out: 400 };
    (next.streams as Record<string, unknown>[]).push({ id: "S3", from: "F:vapor", to: null });
    const d = diffFlows(base(), next);
    expect(d.addedUnits).toEqual([{ id: "P", type: "Pump" }]);
    expect(d.changedParams).toEqual([{ unit: "H", key: "T_out", from: 350, to: 400 }]);
    expect(d.addedStreams).toEqual(["S3"]);
    expect(d.count).toBe(3);
  });

  it("detects removals, rewires, and component changes", () => {
    const next = base();
    next.units = (next.units as { id: string }[]).filter((u) => u.id !== "F");
    next.streams = [
      { id: "S1", from: null, to: "H:in1" },
      { id: "S2", from: "H:out", to: null },  // rewired (was F:in1)
    ];
    next.components = [{ id: "a" }];
    const d = diffFlows(base(), next);
    expect(d.removedUnits).toEqual(["F"]);
    expect(d.rewiredStreams).toEqual(["S2"]);
    expect(d.removedComponents).toEqual(["b"]);
  });
});

describe("mergePositions", () => {
  it("keeps existing positions and staggers new units below", () => {
    const next = base();
    (next.units as Record<string, unknown>[]).push({ id: "P", type: "Pump", params: {} });
    delete (next.units as { xy?: number[] }[])[0].xy; // AI drops coordinates
    const merged = mergePositions(base(), next);
    const units = merged.units as { id: string; xy?: number[] }[];
    expect(units.find((u) => u.id === "H")!.xy).toEqual([100, 50]); // restored
    expect(units.find((u) => u.id === "F")!.xy).toEqual([300, 50]); // kept
    const p = units.find((u) => u.id === "P")!;
    expect(p.xy).toBeDefined();
    expect(p.xy![1]).toBeGreaterThan(50); // staggered below existing rows
  });

  it("carries the canvas UI meta into the proposed doc", () => {
    const cur = { ...base(), meta: { ui: { product: "b", color_mode: "phase" } } };
    const merged = mergePositions(cur, base());
    expect((merged.meta as { ui: { product: string } }).ui.product).toBe("b");
  });
});
