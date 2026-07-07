// File save/load and localStorage autosave. Saved files are engine-compatible
// .flow JSON; UI-only state (boundary node positions, product, backend) rides
// in `meta.ui`, which the engine's from_dict ignores.

import type { FlowDoc } from "../types";

const AUTOSAVE_KEY = "caldyr.autosave.v1";
const THEME_KEY = "caldyr.theme";

export interface UiMeta {
  product?: string;
  backend?: string;
  boundary_xy?: Record<string, [number, number]>;
  color_mode?: string;
  pinned_streams?: string[];
  view_mode?: string;
  unit_overrides?: Record<string, string>;
  cost_config?: Record<string, unknown>;
  scenarios?: Record<string, unknown>[];
  calcs?: { name: string; expr: string }[];
  groups?: {
    id: string; label: string; members: string[];
    collapsed: boolean; xy: [number, number];
  }[];
}

// -- structural gate for untrusted .flow JSON ---------------------------------
// Shared files, old autosaves, and hand-edited localStorage all reach the store
// through JSON.parse; the canvas then indexes into units/streams/components.
// Reject anything whose containers aren't the expected shapes so a crafted or
// corrupt doc fails here (null) instead of throwing mid-render.

const isRecord = (v: unknown): v is Record<string, unknown> =>
  typeof v === "object" && v !== null && !Array.isArray(v);

function recordList(v: unknown): Record<string, unknown>[] | null {
  if (v === undefined) return [];
  if (!Array.isArray(v) || !v.every(isRecord)) return null;
  return v as Record<string, unknown>[];
}

export function sanitizeFlowDoc(v: unknown): FlowDoc | null {
  if (!isRecord(v)) return null;
  if (typeof v.schema !== "string" || !v.schema.startsWith("caldyr.flow/")) return null;
  const units = recordList(v.units);
  const streams = recordList(v.streams);
  const components = recordList(v.components);
  if (!units || !streams || !components) return null;
  for (const u of units) {
    if (typeof u.id !== "string" || typeof u.type !== "string") return null;
    if (u.params !== undefined && !isRecord(u.params)) return null;
    if (u.xy !== undefined &&
        !(Array.isArray(u.xy) && u.xy.every((n) => typeof n === "number"))) {
      delete u.xy;
    }
  }
  for (const s of streams) {
    if (typeof s.id !== "string") return null;
    if (s.from !== undefined && s.from !== null && typeof s.from !== "string") return null;
    if (s.to !== undefined && s.to !== null && typeof s.to !== "string") return null;
    if (s.spec !== undefined && s.spec !== null && !isRecord(s.spec)) return null;
  }
  for (const c of components) {
    if (typeof c.id !== "string") return null;
  }
  if (v.meta !== undefined && !isRecord(v.meta)) delete v.meta;
  const meta = v.meta as Record<string, unknown> | undefined;
  if (meta && meta.ui !== undefined && !isRecord(meta.ui)) delete meta.ui;
  const ui = meta?.ui as Record<string, unknown> | undefined;
  if (ui && ui.calcs !== undefined) {
    const ok = Array.isArray(ui.calcs) && ui.calcs.every(
      (r) => isRecord(r) && typeof r.name === "string" && typeof r.expr === "string");
    if (!ok) delete ui.calcs;
  }
  return v as FlowDoc;
}

export function downloadFlow(doc: FlowDoc, filename = "flowsheet.flow"): void {
  const blob = new Blob([JSON.stringify(doc, null, 2)], { type: "application/json" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

export function pickFlowFile(): Promise<FlowDoc | null> {
  return new Promise((resolve, reject) => {
    const input = document.createElement("input");
    input.type = "file";
    input.accept = ".flow,.json,application/json";
    input.onchange = () => {
      const file = input.files?.[0];
      if (!file) return resolve(null);
      file.text()
        .then((text) => {
          const doc = sanitizeFlowDoc(JSON.parse(text));
          if (!doc) {
            throw new Error("not a valid caldyr .flow document");
          }
          resolve(doc);
        })
        .catch(reject);
    };
    // Cancelling the picker fires no event; resolve(null) only on empty change.
    input.click();
  });
}

export function autosave(doc: FlowDoc): void {
  try {
    localStorage.setItem(AUTOSAVE_KEY, JSON.stringify(doc));
  } catch {
    /* quota/private mode: autosave is best-effort */
  }
}

export function loadAutosave(): FlowDoc | null {
  try {
    const raw = localStorage.getItem(AUTOSAVE_KEY);
    return raw ? sanitizeFlowDoc(JSON.parse(raw)) : null;
  } catch {
    return null;
  }
}

export function clearAutosave(): void {
  try {
    localStorage.removeItem(AUTOSAVE_KEY);
  } catch { /* ignore */ }
}

export function loadPref(key: string): string | null {
  try {
    return localStorage.getItem(`caldyr.${key}`);
  } catch {
    return null;
  }
}

export function savePref(key: string, value: string): void {
  try {
    localStorage.setItem(`caldyr.${key}`, value);
  } catch { /* ignore */ }
}

// -- named projects (localStorage) -------------------------------------------
export interface SavedProject {
  name: string;
  savedAt: string; // ISO
  doc: FlowDoc;
}

const PROJECTS_KEY = "caldyr.projects.v1";

export function listProjects(): SavedProject[] {
  try {
    const raw: unknown = JSON.parse(localStorage.getItem(PROJECTS_KEY) ?? "[]");
    if (!Array.isArray(raw)) return [];
    const out: SavedProject[] = [];
    for (const p of raw) {
      if (typeof p !== "object" || p === null) continue;
      const { name, savedAt, doc } = p as Record<string, unknown>;
      const clean = sanitizeFlowDoc(doc);
      if (typeof name !== "string" || !clean) continue;
      out.push({ name, savedAt: typeof savedAt === "string" ? savedAt : "", doc: clean });
    }
    return out;
  } catch {
    return [];
  }
}

export function saveProject(name: string, doc: FlowDoc): SavedProject[] {
  const projects = listProjects().filter((p) => p.name !== name);
  projects.unshift({ name, savedAt: new Date().toISOString(), doc });
  try {
    localStorage.setItem(PROJECTS_KEY, JSON.stringify(projects.slice(0, 50)));
  } catch { /* quota */ }
  return listProjects();
}

export function deleteProject(name: string): SavedProject[] {
  const projects = listProjects().filter((p) => p.name !== name);
  try {
    localStorage.setItem(PROJECTS_KEY, JSON.stringify(projects));
  } catch { /* ignore */ }
  return projects;
}

export type Theme = "dark" | "light";

export function loadTheme(): Theme {
  return localStorage.getItem(THEME_KEY) === "light" ? "light" : "dark";
}

export function saveTheme(t: Theme): void {
  try {
    localStorage.setItem(THEME_KEY, t);
  } catch { /* ignore */ }
}
