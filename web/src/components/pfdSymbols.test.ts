import { Position } from "@xyflow/react";
import { describe, expect, it } from "vitest";
import type { Port } from "../types";
import { portAnchors } from "./pfdSymbols";

const P = (name: string, direction: "inlet" | "outlet", kind: "material" | "energy" = "material"): Port =>
  ({ name, direction, kind });

describe("portAnchors", () => {
  it("places a flash drum's nozzles physically (feed left, vapor up, liquid down, duty bottom)", () => {
    const ports = [
      P("in1", "inlet"), P("vapor", "outlet"), P("liquid", "outlet"),
      P("duty", "outlet", "energy"),
    ];
    const a = portAnchors("unit", "Flash", ports);
    expect(a.get("in1")!.position).toBe(Position.Left);
    expect(a.get("vapor")!.position).toBe(Position.Top);
    expect(a.get("liquid")!.position).toBe(Position.Bottom);
    expect(a.get("duty")!.position).toBe(Position.Bottom);
  });

  it("puts column products on the right, near top (overhead) and bottom (bottoms)", () => {
    const ports = [
      P("in1", "inlet"), P("distillate", "outlet"), P("bottoms", "outlet"),
      P("condenser_duty", "outlet", "energy"), P("reboiler_duty", "outlet", "energy"),
    ];
    const a = portAnchors("unit", "ShortcutColumn", ports);
    expect(a.get("distillate")).toEqual({ position: Position.Right, offset: 0.1 });
    expect(a.get("bottoms")).toEqual({ position: Position.Right, offset: 0.9 });
    expect(a.get("condenser_duty")!.position).toBe(Position.Top);
    expect(a.get("reboiler_duty")!.position).toBe(Position.Bottom);
  });

  it("computes RigorousColumn feed anchors from the actual (variable) port list", () => {
    const ports = [
      P("in1", "inlet"), P("in2", "inlet"),
      P("distillate", "outlet"), P("bottoms", "outlet"), P("side1", "outlet"),
      P("condenser_duty", "outlet", "energy"), P("reboiler_duty", "outlet", "energy"),
    ];
    const a = portAnchors("unit", "RigorousColumn", ports);
    expect(a.get("in1")!.position).toBe(Position.Left);
    expect(a.get("in2")!.position).toBe(Position.Left);
    expect(a.get("in1")!.offset).toBeLessThan(a.get("in2")!.offset); // top-to-bottom
    expect(a.get("side1")!.position).toBe(Position.Right);
  });

  it("falls back to evenly-spaced left/right/bottom for an unknown type", () => {
    const ports = [
      P("in1", "inlet"), P("in2", "inlet"),
      P("out", "outlet"), P("duty", "outlet", "energy"),
    ];
    const a = portAnchors("unit", "SomeFutureUnit", ports);
    expect(a.get("in1")).toEqual({ position: Position.Left, offset: 1 / 3 });
    expect(a.get("in2")).toEqual({ position: Position.Left, offset: 2 / 3 });
    expect(a.get("out")).toEqual({ position: Position.Right, offset: 0.5 });
    expect(a.get("duty")).toEqual({ position: Position.Bottom, offset: 0.5 });
  });

  it("assigns every port an anchor even when a table omits one", () => {
    // Heater's table lacks an imagined extra outlet; it must still be wired.
    const ports = [
      P("in1", "inlet"), P("out", "outlet"), P("duty", "outlet", "energy"),
      P("vent", "outlet"),
    ];
    const a = portAnchors("unit", "Heater", ports);
    expect([...a.keys()].sort()).toEqual(["duty", "in1", "out", "vent"]);
    expect(a.get("vent")).toBeDefined();
  });

  it("anchors boundary feed/product on the flow axis", () => {
    expect(portAnchors("feed", undefined, [P("out", "outlet")]).get("out"))
      .toEqual({ position: Position.Right, offset: 0.5 });
    expect(portAnchors("product", undefined, [P("in", "inlet")]).get("in"))
      .toEqual({ position: Position.Left, offset: 0.5 });
  });
});
