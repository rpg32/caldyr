// Command registry for the Ctrl+K palette. Kept UI-free and pure so the list and
// the fuzzy filter are unit-testable: buildCommands() takes the store actions it
// needs plus the unit-op catalog and returns a flat, dispatchable command list.
import type { UnitType } from "../types";
import { UNIT_SETS, type UnitSet } from "./units";
import type { ColorMode, Tab, ViewMode } from "../store";

export type CommandSection = "Run" | "File" | "View" | "Go to" | "Add unit" | "Help";

export interface Command {
  id: string;
  title: string;
  section: CommandSection;
  run: () => void;
  hint?: string;        // keybinding hint shown on the right (e.g. "Ctrl+S")
  keywords?: string;    // extra search terms not in the title
}

// The slice of store state + actions the palette drives. A structural subset so
// `useStore.getState()` satisfies it without extra plumbing.
export interface CommandContext {
  unitTypes: UnitType[];
  solve: () => void;
  cost: () => void;
  runAutoLayout: () => void;
  groupSelection: () => void;
  newFlowsheet: () => void;
  openFile: () => void;
  saveFile: () => void;
  toggleProjects: () => void;
  toggleScenarios: () => void;
  toggleChat: () => void;
  toggleGlossary: () => void;
  toggleTutorials: () => void;
  toggleSettings: () => void;
  setViewMode: (m: ViewMode) => void;
  setColorMode: (m: ColorMode) => void;
  setUnitSet: (u: UnitSet) => void;
  setTab: (t: Tab) => void;
  addUnitAtCenter: (t: UnitType) => void;
}

const VIEWS: { id: ViewMode; label: string }[] = [
  { id: "bfd", label: "BFD" }, { id: "pfd", label: "PFD" }, { id: "pid", label: "P&ID" },
];
const COLORS: ColorMode[] = ["none", "phase", "temperature"];
const TABS: { id: Tab; label: string }[] = [
  { id: "params", label: "Params" }, { id: "streams", label: "Streams" },
  { id: "economics", label: "Econ" }, { id: "optimize", label: "Opt" },
  { id: "study", label: "Study" }, { id: "calc", label: "Calc" },
  { id: "tools", label: "Tools" },
];

/** Flat, ordered command list. Order here is the default (empty-query) order. */
export function buildCommands(c: CommandContext): Command[] {
  const cmds: Command[] = [
    { id: "solve", title: "Solve", section: "Run", run: c.solve, keywords: "run converge" },
    { id: "cost", title: "Cost", section: "Run", run: c.cost, keywords: "economics tea lcop" },
    { id: "arrange", title: "Arrange (auto-layout)", section: "Run", run: c.runAutoLayout, keywords: "layout elk tidy" },
    { id: "group", title: "Group selection", section: "Run", run: c.groupSelection, keywords: "collapse block" },

    { id: "new", title: "New flowsheet", section: "File", run: c.newFlowsheet, keywords: "reset clear" },
    { id: "open", title: "Open .flow file", section: "File", run: c.openFile, hint: "Ctrl+O" },
    { id: "save", title: "Save .flow file", section: "File", run: c.saveFile, hint: "Ctrl+S" },
    { id: "projects", title: "Projects & templates", section: "File", run: c.toggleProjects, keywords: "template gallery load" },
    { id: "cases", title: "Cases / scenarios", section: "File", run: c.toggleScenarios, keywords: "scenario compare" },

    ...VIEWS.map((v): Command => ({
      id: `view-${v.id}`, title: `View: ${v.label}`, section: "View",
      run: () => c.setViewMode(v.id), keywords: "diagram bfd pfd pid",
    })),
    ...COLORS.map((m): Command => ({
      id: `color-${m}`, title: `Color streams: ${m}`, section: "View",
      run: () => c.setColorMode(m), keywords: "phase temperature heatmap",
    })),
    ...UNIT_SETS.map((u): Command => ({
      id: `units-${u}`, title: `Units: ${u}`, section: "View",
      run: () => c.setUnitSet(u as UnitSet), keywords: "unit system si field",
    })),

    ...TABS.map((t): Command => ({
      id: `tab-${t.id}`, title: `Go to ${t.label}`, section: "Go to",
      run: () => c.setTab(t.id), keywords: "panel tab inspector",
    })),

    { id: "tutorials", title: "Guided tutorials", section: "Help", run: c.toggleTutorials, keywords: "tutorial learn walkthrough guide getting started" },
    { id: "glossary", title: "Open Glossary", section: "Help", run: c.toggleGlossary, keywords: "terms lcop reflux pinch" },
    { id: "settings", title: "Open Settings", section: "Help", run: c.toggleSettings, keywords: "cost assumptions prices" },
    { id: "copilot", title: "Toggle Copilot", section: "Help", run: c.toggleChat, keywords: "ai chat assistant" },
  ];

  // One "Add <UnitType>" command per palette unit op (drops at canvas center).
  for (const t of c.unitTypes) {
    cmds.push({
      id: `add-${t.type}`, title: `Add ${t.type}`, section: "Add unit",
      run: () => c.addUnitAtCenter(t), keywords: `insert new unit ${t.description ?? ""}`,
    });
  }
  return cmds;
}

/** Subsequence fuzzy score of `query` against `text` (both lowercased). Returns
 *  null when not a subsequence; higher is a better match. Rewards contiguous
 *  runs, word-boundary starts, and exact substring hits. */
export function fuzzyScore(text: string, query: string): number | null {
  if (!query) return 0;
  const t = text.toLowerCase();
  const q = query.toLowerCase();
  if (t.includes(q)) {
    // Strong bonus for a contiguous hit, extra if it's at a word boundary.
    const at = t.indexOf(q);
    const boundary = at === 0 || /\W|_/.test(t[at - 1]) ? 6 : 0;
    return 100 + boundary + q.length - at * 0.1;
  }
  let score = 0;
  let ti = 0;
  let prev = -2;
  for (const ch of q) {
    const found = t.indexOf(ch, ti);
    if (found === -1) return null;
    score += 1;
    if (found === prev + 1) score += 3;                     // contiguous run
    if (found === 0 || /\W|_|\s/.test(t[found - 1])) score += 2; // word start
    prev = found;
    ti = found + 1;
  }
  return score;
}

/** Filter + rank commands for a query. Empty query keeps registry order. */
export function filterCommands(commands: Command[], query: string): Command[] {
  const q = query.trim();
  if (!q) return commands;
  const scored: { cmd: Command; score: number; idx: number }[] = [];
  commands.forEach((cmd, idx) => {
    const hay = `${cmd.title} ${cmd.keywords ?? ""}`;
    const score = fuzzyScore(hay, q);
    if (score !== null) scored.push({ cmd, score, idx });
  });
  scored.sort((a, b) => (b.score - a.score) || (a.idx - b.idx)); // stable on ties
  return scored.map((x) => x.cmd);
}
