import {
  BookOpen, Bot, Boxes, Calculator, CircleHelp, Command as CommandIcon, FilePlus2,
  FolderKanban, FolderOpen, GraduationCap, Layers, MoreHorizontal, Moon, Network,
  Play, Redo2, Save, Settings, Sun, Undo2,
} from "lucide-react";
import { useEffect, useRef, useState } from "react";
import { UNIT_SETS, type UnitSet } from "../lib/units";
import { useStore, type ColorMode, type ViewMode } from "../store";
import { relaunchTour } from "./Tour";
import { Button } from "./ui";

const VIEWS: { id: ViewMode; label: string; title: string }[] = [
  { id: "bfd", label: "BFD", title: "Block flow diagram — plain blocks, no duty lines" },
  { id: "pfd", label: "PFD", title: "Process flow diagram — full detail" },
  { id: "pid", label: "P&ID", title: "PFD + control/logical-op instrumentation overlay" },
];

/** Vertical rule separating toolbar clusters. */
function Sep() {
  return <span className="mx-0.5 h-5 w-px shrink-0 bg-line" aria-hidden />;
}

/** Track a container's width so clusters can collapse into the overflow menu
 *  instead of wrapping onto a second row. */
function useWidth<T extends HTMLElement>() {
  const ref = useRef<T>(null);
  const [w, setW] = useState(0);
  useEffect(() => {
    const el = ref.current;
    if (!el) return;
    const ro = new ResizeObserver((entries) => setW(entries[0].contentRect.width));
    ro.observe(el);
    setW(el.getBoundingClientRect().width);
    return () => ro.disconnect();
  }, []);
  return [ref, w] as const;
}

const selectCls = "rounded-md border border-line bg-panel2 px-2 py-1.5 text-text";

function ColorBySelect() {
  const colorMode = useStore((s) => s.colorMode);
  const setColorMode = useStore((s) => s.setColorMode);
  return (
    <label className="flex items-center gap-1.5 text-muted">
      <span className="text-[11px]">color by</span>
      <select className={selectCls} value={colorMode}
        onChange={(e) => setColorMode(e.target.value as ColorMode)}
        title="Color streams by solved phase or temperature">
        <option value="none">none</option>
        <option value="phase">phase</option>
        <option value="temperature">temperature</option>
      </select>
    </label>
  );
}

function UnitsSelect() {
  const unitSet = useStore((s) => s.unitSet);
  const setUnitSet = useStore((s) => s.setUnitSet);
  return (
    <label className="flex items-center gap-1.5 text-muted">
      <span className="text-[11px]">units</span>
      <select className={selectCls} value={unitSet}
        onChange={(e) => setUnitSet(e.target.value as UnitSet)}
        title="Unit system for all inputs and outputs (engine stays SI internally)">
        {UNIT_SETS.map((u) => <option key={u} value={u}>{u}</option>)}
      </select>
    </label>
  );
}

function ViewToggle() {
  const viewMode = useStore((s) => s.viewMode);
  const setViewMode = useStore((s) => s.setViewMode);
  return (
    <div className="flex shrink-0 overflow-hidden rounded-md border border-line" role="radiogroup"
      aria-label="Diagram view">
      {VIEWS.map((v) => (
        <button key={v.id} role="radio" aria-checked={viewMode === v.id} title={v.title}
          className={`cursor-pointer border-0 px-2 py-1.5 text-[12px] transition-colors ${
            viewMode === v.id ? "bg-accent/20 text-accent" : "bg-panel2 text-muted hover:text-text"
          }`}
          onClick={() => setViewMode(v.id)}>
          {v.label}
        </button>
      ))}
    </div>
  );
}

/** A single action row inside the ⋯ overflow menu. */
function MenuItem({ icon, label, onClick }: { icon: React.ReactNode; label: string; onClick: () => void }) {
  return (
    <button className="flex w-full items-center gap-2 rounded-md px-2 py-1.5 text-left text-[13px] text-text hover:bg-panel2"
      onClick={onClick}>
      <span className="text-muted">{icon}</span>{label}
    </button>
  );
}

/** ⋯ overflow menu: always holds Arrange + Group; also absorbs whichever clusters
 *  don't fit at the current width (view controls, then Projects/Cases). */
function Overflow({ level }: { level: number }) {
  const s = useStore();
  const [open, setOpen] = useState(false);
  const ref = useRef<HTMLDivElement>(null);
  useEffect(() => {
    if (!open) return;
    const onDown = (e: MouseEvent) => {
      if (ref.current && !ref.current.contains(e.target as Node)) setOpen(false);
    };
    const onKey = (e: KeyboardEvent) => e.key === "Escape" && setOpen(false);
    window.addEventListener("mousedown", onDown);
    window.addEventListener("keydown", onKey);
    return () => {
      window.removeEventListener("mousedown", onDown);
      window.removeEventListener("keydown", onKey);
    };
  }, [open]);
  const act = (fn: () => void) => () => { setOpen(false); fn(); };
  return (
    <div className="relative" ref={ref}>
      <Button variant="ghost" icon={<MoreHorizontal size={16} />} aria-label="More actions"
        title="More actions" aria-expanded={open} onClick={() => setOpen((o) => !o)} />
      {open && (
        <div role="menu"
          className="absolute right-0 top-full z-40 mt-1 w-56 rounded-lg border border-line bg-panel p-1 shadow-2xl">
          <MenuItem icon={<Network size={14} />} label="Arrange (auto-layout)"
            onClick={act(() => void s.runAutoLayout())} />
          <MenuItem icon={<Boxes size={14} />} label="Group selection"
            onClick={act(s.groupSelection)} />
          {level >= 1 && (
            <div className="mt-1 border-t border-line px-2 pb-1 pt-2">
              <ColorBySelect /><div className="h-1.5" /><UnitsSelect />
            </div>
          )}
          {level >= 2 && (
            <div className="mt-1 border-t border-line pt-1">
              <MenuItem icon={<FolderKanban size={14} />} label="Projects & templates"
                onClick={act(s.toggleProjects)} />
              <MenuItem icon={<Layers size={14} />} label="Cases / scenarios"
                onClick={act(s.toggleScenarios)} />
            </div>
          )}
        </div>
      )}
    </div>
  );
}

export function Toolbar() {
  const s = useStore();
  const [ref, w] = useWidth<HTMLElement>();
  // Progressive collapse: narrower widths fold more clusters into the ⋯ menu
  // rather than wrapping the toolbar onto extra rows. 0 = everything inline.
  const level = w === 0 ? 0 : w < 900 ? 2 : w < 1120 ? 1 : 0;

  return (
    <header ref={ref}
      className="flex min-w-0 flex-nowrap items-center gap-x-1.5 overflow-hidden border-b border-line bg-panel px-3 py-2">
      <strong className="shrink-0 text-[15px] text-accent">Caldyr</strong>

      <Sep />
      {/* File */}
      <Button variant="ghost" icon={<FilePlus2 size={14} />} title="New flowsheet"
        onClick={s.newFlowsheet}>New</Button>
      <Button variant="ghost" icon={<FolderOpen size={14} />} title="Open .flow file (Ctrl+O)"
        onClick={() => void s.openFile()}>Open</Button>
      <Button variant="ghost" icon={<Save size={14} />} title="Save .flow file (Ctrl+S)"
        onClick={s.saveFile}>Save</Button>
      <Button variant="ghost" icon={<Undo2 size={14} />} title="Undo (Ctrl+Z)"
        disabled={!s.past.length} onClick={s.undo} aria-label="Undo" />
      <Button variant="ghost" icon={<Redo2 size={14} />} title="Redo (Ctrl+Y)"
        disabled={!s.future.length} onClick={s.redo} aria-label="Redo" />

      <Sep />
      {/* Run */}
      <label className="flex items-center gap-1.5 text-muted">
        <span className="sr-only">Solver backend</span>
        <select className={selectCls} value={s.backend}
          onChange={(e) => s.setBackend(e.target.value)} title="Solver backend">
          <option value="sequential">sequential-modular</option>
          <option value="equation_oriented">equation-oriented</option>
        </select>
      </label>
      <Button variant="primary" icon={<Play size={14} />} busy={s.busy === "solve"}
        disabled={!s.engineReady || s.busy !== null} onClick={() => void s.solve()}>
        Solve
      </Button>
      <Button variant="primary" icon={<Calculator size={14} />} busy={s.busy === "cost"}
        disabled={!s.engineReady || s.busy !== null} onClick={() => void s.cost()}>
        Cost
      </Button>

      <Sep />
      {/* View */}
      {level < 1 && <ColorBySelect />}
      {level < 1 && <UnitsSelect />}
      <ViewToggle />

      {/* Overflow (Arrange, Group, + collapsed clusters) */}
      <Overflow level={level} />

      {/* Flexible middle: absorbs slack and truncates first when space is tight. */}
      <span className="mx-1 ml-auto min-w-0 truncate text-muted" role="status">{s.status}</span>

      {/* Right cluster */}
      <Button variant="ghost" icon={<CommandIcon size={13} />} title="Command palette (Ctrl+K)"
        aria-label="Command palette" onClick={s.toggleCommandPalette} />
      {level < 2 && (
        <>
          <Button variant="ghost" icon={<FolderKanban size={14} />}
            title="Projects & templates" onClick={s.toggleProjects}>Projects</Button>
          <Button variant="ghost" icon={<Layers size={14} />}
            title="Cases / scenarios — save, switch, compare" onClick={s.toggleScenarios}>Cases</Button>
        </>
      )}
      <Button variant={s.chatOpen ? "primary" : "ghost"} icon={<Bot size={14} />}
        title="AI copilot chat" onClick={s.toggleChat}>Copilot</Button>
      <Button variant="ghost" aria-label="Guided tutorials" title="Guided tutorials — load a flowsheet and follow the step-by-step guide"
        icon={<GraduationCap size={14} />} onClick={s.toggleTutorials} />
      <Button variant="ghost" aria-label="Glossary" title="Glossary of terms (LCOP, reflux, pinch…)"
        icon={<BookOpen size={14} />} onClick={s.toggleGlossary} />
      <Button variant="ghost" aria-label="Cost assumptions settings" title="Cost assumptions (prices, factors, sizing)"
        icon={<Settings size={14} />} onClick={s.toggleSettings} />
      <Button variant="ghost" aria-label="Show the quick tour" title="Quick tour"
        icon={<CircleHelp size={14} />} onClick={relaunchTour} />
      <Button variant="ghost" aria-label="Toggle color theme" title="Toggle light/dark"
        icon={s.theme === "dark" ? <Sun size={14} /> : <Moon size={14} />}
        onClick={s.toggleTheme} />
    </header>
  );
}
