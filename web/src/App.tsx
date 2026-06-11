import { useEffect, useRef } from "react";
import { CanvasView } from "./components/CanvasView";
import { Inspector } from "./components/Inspector";
import { Palette } from "./components/Palette";
import { ProjectsDialog } from "./components/ProjectsDialog";
import { Toasts } from "./components/Toasts";
import { Toolbar } from "./components/Toolbar";
import { Tour } from "./components/Tour";
import { useStore } from "./store";

/** True when the event target is a text-entry element (don't steal its keys). */
function inTextInput(e: KeyboardEvent): boolean {
  const t = e.target as HTMLElement | null;
  if (!t) return false;
  return t.tagName === "INPUT" || t.tagName === "TEXTAREA" || t.tagName === "SELECT"
    || t.isContentEditable;
}

function useKeyboardShortcuts() {
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if (!(e.ctrlKey || e.metaKey) || inTextInput(e)) return;
      const s = useStore.getState();
      const k = e.key.toLowerCase();
      if (k === "z" && e.shiftKey) { e.preventDefault(); s.redo(); }
      else if (k === "z") { e.preventDefault(); s.undo(); }
      else if (k === "y") { e.preventDefault(); s.redo(); }
      else if (k === "c") { s.copySelection(); }
      else if (k === "v") { e.preventDefault(); s.paste(); }
      else if (k === "d") { e.preventDefault(); s.duplicateSelection(); }
      else if (k === "s") { e.preventDefault(); s.saveFile(); }
      else if (k === "o") { e.preventDefault(); void s.openFile(); }
      else if (k === "a") { e.preventDefault(); s.selectAll(); }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, []);
}

/** Debounced localStorage autosave whenever the flowsheet changes. */
function useAutosave() {
  const timer = useRef<number | undefined>(undefined);
  useEffect(() => {
    const unsub = useStore.subscribe((state, prev) => {
      if (
        state.nodes === prev.nodes && state.edges === prev.edges
        && state.components === prev.components
        && state.propertyPackage === prev.propertyPackage
        && state.product === prev.product
        && state.logical === prev.logical
        && state.solverHints === prev.solverHints
        && state.groups === prev.groups
        && state.viewMode === prev.viewMode
        && state.calcs === prev.calcs
      ) return;
      window.clearTimeout(timer.current);
      timer.current = window.setTimeout(() => useStore.getState().autosaveNow(), 800);
    });
    return () => {
      unsub();
      window.clearTimeout(timer.current);
    };
  }, []);
}

/** Drag handle that resizes the inspector panel. */
function PanelResizer() {
  const setInspectorWidth = useStore((s) => s.setInspectorWidth);
  const onPointerDown = (e: React.PointerEvent) => {
    e.preventDefault();
    const onMove = (ev: PointerEvent) =>
      setInspectorWidth(window.innerWidth - ev.clientX);
    const onUp = () => {
      window.removeEventListener("pointermove", onMove);
      window.removeEventListener("pointerup", onUp);
    };
    window.addEventListener("pointermove", onMove);
    window.addEventListener("pointerup", onUp);
  };
  return (
    <div
      role="separator"
      aria-orientation="vertical"
      aria-label="Resize panel"
      className="w-1.5 cursor-col-resize bg-panel transition-colors hover:bg-accent/40"
      onPointerDown={onPointerDown}
    />
  );
}

export function App() {
  const init = useStore((s) => s.init);
  const theme = useStore((s) => s.theme);
  const inspectorWidth = useStore((s) => s.inspectorWidth);

  useEffect(() => {
    document.documentElement.dataset.theme = theme;
  }, [theme]);

  useEffect(() => { void init(); }, [init]);
  useKeyboardShortcuts();
  useAutosave();

  return (
    <div className="flex h-screen flex-col bg-bg text-text">
      <Toolbar />
      <div
        className="grid min-h-0 flex-1"
        style={{ gridTemplateColumns: `150px 1fr 6px ${inspectorWidth}px` }}
      >
        <Palette />
        <CanvasView />
        <PanelResizer />
        <Inspector />
      </div>
      <Toasts />
      <ProjectsDialog />
      <Tour />
    </div>
  );
}
