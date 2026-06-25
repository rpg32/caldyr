// Round-trip tests for canvas <-> .flow conversion — the invariant that keeps
// flowsheet-as-canvas and flowsheet-as-code the same object.
import { describe, expect, it } from "vitest";
import { ammoniaLoop } from "./examples";
import { canvasToFlow, flowToCanvas } from "./flow";
import type { FlowDoc, UnitType } from "./types";

// Port catalog mirroring the engine's /unit-types for the types used in tests.
const UNIT_TYPES: UnitType[] = [
  { type: "Mixer", doc: "", ports: [
    { name: "in1", direction: "inlet", kind: "material" },
    { name: "in2", direction: "inlet", kind: "material" },
    { name: "out", direction: "outlet", kind: "material" },
  ]},
  { type: "Heater", doc: "", ports: [
    { name: "in1", direction: "inlet", kind: "material" },
    { name: "out", direction: "outlet", kind: "material" },
    { name: "duty", direction: "outlet", kind: "energy" },
  ]},
  { type: "EquilibriumReactor", doc: "", ports: [
    { name: "in1", direction: "inlet", kind: "material" },
    { name: "out", direction: "outlet", kind: "material" },
    { name: "duty", direction: "outlet", kind: "energy" },
  ]},
  { type: "Flash", doc: "", ports: [
    { name: "in1", direction: "inlet", kind: "material" },
    { name: "vapor", direction: "outlet", kind: "material" },
    { name: "liquid", direction: "outlet", kind: "material" },
    { name: "duty", direction: "outlet", kind: "energy" },
  ]},
  { type: "Splitter", doc: "", ports: [
    { name: "in1", direction: "inlet", kind: "material" },
    { name: "out1", direction: "outlet", kind: "material" },
    { name: "out2", direction: "outlet", kind: "material" },
  ]},
];

interface FlowStreamish {
  id: string; from: string | null; to: string | null; spec?: Record<string, unknown>;
}
const streamsOf = (doc: FlowDoc) => doc.streams as FlowStreamish[];
const unitsOf = (doc: FlowDoc) =>
  doc.units as { id: string; type: string; params: Record<string, unknown>; xy?: number[] }[];

describe("flowToCanvas", () => {
  it("creates a node per unit plus boundary feed/product nodes", () => {
    const c = flowToCanvas(ammoniaLoop, UNIT_TYPES);
    const kinds = c.nodes.map((n) => n.data.kind);
    expect(kinds.filter((k) => k === "unit")).toHaveLength(6);
    expect(kinds.filter((k) => k === "feed")).toHaveLength(1);   // MAKEUP
    expect(kinds.filter((k) => k === "product")).toHaveLength(2); // PRODUCT, PURGE
  });

  it("connects edges to the right handles and labels them with the stream id", () => {
    const c = flowToCanvas(ammoniaLoop, UNIT_TYPES);
    const recycle = c.edges.find((e) => e.id === "RECYCLE")!;
    expect(recycle.source).toBe("SPLIT");
    expect(recycle.sourceHandle).toBe("out1");
    expect(recycle.target).toBe("MIX");
    expect(recycle.targetHandle).toBe("in2");
    expect(recycle.label).toBe("RECYCLE");
  });

  it("preserves component list and property package", () => {
    const c = flowToCanvas(ammoniaLoop, UNIT_TYPES);
    expect(c.components).toEqual(["nitrogen", "hydrogen", "ammonia", "argon"]);
    expect(c.propertyPackage).toBe("thermo:PR");
  });
});

describe("canvasToFlow ∘ flowToCanvas round-trip", () => {
  const c = flowToCanvas(ammoniaLoop, UNIT_TYPES);
  const doc = canvasToFlow(c.nodes, c.edges, c.components, c.propertyPackage);

  it("preserves every unit with type, params and position", () => {
    const orig = unitsOf(ammoniaLoop);
    const round = unitsOf(doc);
    expect(round.map((u) => u.id).sort()).toEqual(orig.map((u) => u.id).sort());
    for (const u of orig) {
      const r = round.find((x) => x.id === u.id)!;
      expect(r.type).toBe(u.type);
      expect(r.params).toEqual(u.params);
      expect(r.xy).toEqual(u.xy);
    }
  });

  it("preserves connectivity of every original stream", () => {
    const orig = streamsOf(ammoniaLoop);
    const round = streamsOf(doc);
    for (const s of orig) {
      const r = round.find((x) => x.id === s.id);
      // Unwired duty streams are regenerated under engine-default ids; skip them.
      if (s.from?.includes(":duty")) continue;
      expect(r, `stream ${s.id} survives the round-trip`).toBeDefined();
      expect(r!.from).toBe(s.from);
      expect(r!.to).toBe(s.to);
    }
  });

  it("carries the feed spec through unchanged", () => {
    const orig = streamsOf(ammoniaLoop).find((s) => s.id === "MAKEUP")!;
    const round = streamsOf(doc).find((s) => s.id === "MAKEUP")!;
    expect(round.spec).toEqual(orig.spec);
  });

  it("emits boundary streams for unconnected outlets (duties)", () => {
    const round = streamsOf(doc);
    // Each Heater/Reactor/Flash duty port is unwired on the canvas → boundary stream.
    for (const id of ["PREHEAT_duty", "RXN_duty", "COOL_duty", "SEP_duty"]) {
      const s = round.find((x) => x.id === id);
      expect(s, `expected boundary duty stream ${id}`).toBeDefined();
      expect(s!.to).toBeNull();
    }
  });

  it("filters stale components out of feed compositions", () => {
    const trimmed = canvasToFlow(c.nodes, c.edges, ["nitrogen", "hydrogen", "ammonia"],
      c.propertyPackage);
    const feed = streamsOf(trimmed).find((s) => s.id === "MAKEUP")!;
    expect(Object.keys(feed.spec!.z as Record<string, number>)).not.toContain("argon");
  });

  it("stores boundary node positions in meta.ui and restores them", () => {
    const meta = (doc.meta as { ui: { boundary_xy: Record<string, [number, number]> } }).ui;
    expect(meta.boundary_xy.MAKEUP).toBeDefined();
    const again = flowToCanvas(doc, UNIT_TYPES);
    const feedNode = again.nodes.find((n) => n.data.kind === "feed")!;
    expect([Math.round(feedNode.position.x), Math.round(feedNode.position.y)])
      .toEqual(meta.boundary_xy.MAKEUP);
  });

  it("round-trips per-field unit overrides through meta.ui", () => {
    const withUnits = canvasToFlow(c.nodes, c.edges, c.components, c.propertyPackage,
      { unit_overrides: { "PREHEAT:T_out": "°F", "MAKEUP:P": "bar" } });
    const meta = (withUnits.meta as { ui: { unit_overrides?: Record<string, string> } }).ui;
    expect(meta.unit_overrides).toEqual({ "PREHEAT:T_out": "°F", "MAKEUP:P": "bar" });
    const back = flowToCanvas(withUnits, UNIT_TYPES);
    expect(back.ui.unit_overrides).toEqual({ "PREHEAT:T_out": "°F", "MAKEUP:P": "bar" });
  });

  it("is a fixed point: a second round-trip yields the identical document", () => {
    const c2 = flowToCanvas(doc, UNIT_TYPES);
    const doc2 = canvasToFlow(c2.nodes, c2.edges, c2.components, c2.propertyPackage);
    expect(doc2).toEqual(doc);
  });
});
