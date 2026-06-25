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
  calcs?: { name: string; expr: string }[];
  groups?: {
    id: string; label: string; members: string[];
    collapsed: boolean; xy: [number, number];
  }[];
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
          const doc = JSON.parse(text) as FlowDoc;
          if (typeof doc.schema !== "string" || !doc.schema.startsWith("caldyr.flow/")) {
            throw new Error("not a caldyr .flow document (missing schema)");
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
    return raw ? (JSON.parse(raw) as FlowDoc) : null;
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
    return JSON.parse(localStorage.getItem(PROJECTS_KEY) ?? "[]") as SavedProject[];
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
