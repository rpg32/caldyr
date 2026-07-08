import { describe, expect, it, vi } from "vitest";
import type { UnitType } from "../types";
import { buildCommands, filterCommands, fuzzyScore, type CommandContext } from "./commands";

const unitTypes: UnitType[] = [
  { type: "Flash", doc: "", ports: [] },
  { type: "HeatExchanger", doc: "", ports: [] },
];

function ctx(): CommandContext {
  return {
    unitTypes,
    solve: vi.fn(), cost: vi.fn(), runAutoLayout: vi.fn(), groupSelection: vi.fn(),
    newFlowsheet: vi.fn(), openFile: vi.fn(), saveFile: vi.fn(),
    toggleProjects: vi.fn(), toggleScenarios: vi.fn(), toggleChat: vi.fn(),
    toggleGlossary: vi.fn(), toggleTutorials: vi.fn(), toggleSettings: vi.fn(),
    setViewMode: vi.fn(), setColorMode: vi.fn(), setUnitSet: vi.fn(), setTab: vi.fn(),
    addUnitAtCenter: vi.fn(),
  };
}

describe("buildCommands", () => {
  it("includes core actions and one Add command per unit type", () => {
    const cmds = buildCommands(ctx());
    const ids = new Set(cmds.map((c) => c.id));
    expect(ids.has("solve")).toBe(true);
    expect(ids.has("cost")).toBe(true);
    expect(ids.has("view-pfd")).toBe(true);
    expect(ids.has("tab-streams")).toBe(true);
    expect(ids.has("add-Flash")).toBe(true);
    expect(ids.has("add-HeatExchanger")).toBe(true);
    expect(cmds.filter((c) => c.section === "Add unit")).toHaveLength(2);
  });

  it("dispatches the right store action when a command runs", () => {
    const c = ctx();
    const cmds = buildCommands(c);
    cmds.find((x) => x.id === "solve")!.run();
    expect(c.solve).toHaveBeenCalledOnce();

    cmds.find((x) => x.id === "view-bfd")!.run();
    expect(c.setViewMode).toHaveBeenCalledWith("bfd");

    cmds.find((x) => x.id === "add-Flash")!.run();
    expect(c.addUnitAtCenter).toHaveBeenCalledWith(unitTypes[0]);
  });
});

describe("filterCommands", () => {
  const cmds = buildCommands(ctx());

  it("returns the full list in registry order for an empty query", () => {
    expect(filterCommands(cmds, "")).toEqual(cmds);
    expect(filterCommands(cmds, "   ")).toEqual(cmds);
  });

  it("ranks a direct title hit first", () => {
    expect(filterCommands(cmds, "solve")[0].id).toBe("solve");
    expect(filterCommands(cmds, "add flash")[0].id).toBe("add-Flash");
  });

  it("matches keywords, not just titles", () => {
    // 'lcop' only appears in the Cost command's keywords.
    expect(filterCommands(cmds, "lcop").map((c) => c.id)).toContain("cost");
  });

  it("excludes non-subsequence queries", () => {
    expect(filterCommands(cmds, "zzqx")).toHaveLength(0);
  });
});

describe("fuzzyScore", () => {
  it("returns null for a non-subsequence", () => {
    expect(fuzzyScore("Solve", "xyz")).toBeNull();
  });
  it("scores a contiguous substring above a scattered subsequence", () => {
    const contiguous = fuzzyScore("Save file", "save")!;
    const scattered = fuzzyScore("Solve and validate everything", "save")!;
    expect(contiguous).toBeGreaterThan(scattered);
  });
  it("treats an empty query as a neutral match", () => {
    expect(fuzzyScore("anything", "")).toBe(0);
  });
});
