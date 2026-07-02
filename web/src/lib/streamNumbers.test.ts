import type { Edge } from "@xyflow/react";
import { describe, expect, it } from "vitest";
import type { CaldyrNode } from "../flow";
import { streamNumbers } from "./streamNumbers";

const node = (id: string, kind: "unit" | "feed" | "product" = "unit"): CaldyrNode => ({
  id, type: "caldyr", position: { x: 0, y: 0 },
  data: { kind, label: id, unitType: kind === "unit" ? "Heater" : undefined, ports: [], params: {} },
} as unknown as CaldyrNode);

const edge = (id: string, source: string, target: string): Edge =>
  ({ id, source, target } as Edge);

describe("streamNumbers", () => {
  it("numbers streams in flow order: feed first, then downstream", () => {
    // FEED -> A -> B, with the wires deliberately created out of order.
    const nodes = [node("A"), node("B"), node("FEED", "feed"), node("P", "product")];
    const edges = [
      edge("S2", "A", "B"),
      edge("S1", "FEED", "A"),
      edge("S3", "B", "P"),
    ];
    const nos = streamNumbers(nodes, edges);
    expect(nos.get("S1")).toBe(1); // FEED source ranks first
    expect(nos.get("S2")).toBe(2);
    expect(nos.get("S3")).toBe(3);
  });

  it("breaks ties on source rank by the id's natural number", () => {
    // Two feeds into a mixer: both source nodes rank 0, so id order decides.
    const nodes = [node("F1", "feed"), node("F2", "feed"), node("M")];
    const edges = [edge("S10", "F2", "M"), edge("S2", "F1", "M")];
    const nos = streamNumbers(nodes, edges);
    expect(nos.get("S2")).toBe(1);
    expect(nos.get("S10")).toBe(2);
  });

  it("stays stable and total on a recycle loop (no blow-up)", () => {
    // A -> B -> A cycle plus a feed; every edge still gets a distinct number.
    const nodes = [node("FEED", "feed"), node("A"), node("B")];
    const edges = [
      edge("S1", "FEED", "A"),
      edge("S2", "A", "B"),
      edge("S3", "B", "A"), // recycle
    ];
    const nos = streamNumbers(nodes, edges);
    expect(new Set(nos.values()).size).toBe(3);
    expect([...nos.values()].sort()).toEqual([1, 2, 3]);
    expect(nos.get("S1")).toBe(1);
  });
});
