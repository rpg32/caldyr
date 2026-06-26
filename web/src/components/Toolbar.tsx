import {
  Bot, Boxes, Calculator, CircleHelp, FilePlus2, FolderKanban, FolderOpen, Moon,
  Network, Play, Redo2, Save, Settings, Sun, Undo2,
} from "lucide-react";
import { UNIT_SETS, type UnitSet } from "../lib/units";
import { useStore, type ColorMode, type ViewMode } from "../store";
import { relaunchTour } from "./Tour";
import { Button } from "./ui";

const VIEWS: { id: ViewMode; label: string; title: string }[] = [
  { id: "bfd", label: "BFD", title: "Block flow diagram — plain blocks, no duty lines" },
  { id: "pfd", label: "PFD", title: "Process flow diagram — full detail" },
  { id: "pid", label: "P&ID", title: "PFD + control/logical-op instrumentation overlay" },
];

function Sep() {
  return <span className="mx-1 h-5 w-px shrink-0 bg-line" aria-hidden />;
}

export function Toolbar() {
  const s = useStore();
  return (
    <header className="flex flex-wrap items-center gap-x-2 gap-y-1.5 border-b border-line bg-panel px-3 py-2">
      <strong className="text-[15px] text-accent">Caldyr</strong>

      <Sep />
      <Button variant="ghost" icon={<FilePlus2 size={14} />} title="New flowsheet"
        onClick={s.newFlowsheet}>New</Button>
      <Button variant="ghost" icon={<FolderOpen size={14} />} title="Open .flow file (Ctrl+O)"
        onClick={() => void s.openFile()}>Open</Button>
      <Button variant="ghost" icon={<Save size={14} />} title="Save .flow file (Ctrl+S)"
        onClick={s.saveFile}>Save</Button>

      <Sep />
      <Button variant="ghost" icon={<Undo2 size={14} />} title="Undo (Ctrl+Z)"
        disabled={!s.past.length} onClick={s.undo} aria-label="Undo" />
      <Button variant="ghost" icon={<Redo2 size={14} />} title="Redo (Ctrl+Y)"
        disabled={!s.future.length} onClick={s.redo} aria-label="Redo" />

      <Sep />
      <label className="flex items-center gap-1.5 text-muted">
        <span className="sr-only">Solver backend</span>
        <select
          className="rounded-md border border-line bg-panel2 px-2 py-1.5 text-text"
          value={s.backend}
          onChange={(e) => s.setBackend(e.target.value)}
          title="Solver backend"
        >
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
      <Button variant="ghost" icon={<Network size={14} />} title="Auto-arrange the flowsheet"
        onClick={() => void s.runAutoLayout()}>Arrange</Button>
      <Button variant="ghost" icon={<Boxes size={14} />}
        title="Group the selected nodes into a collapsible block"
        onClick={s.groupSelection}>Group</Button>
      <label className="flex items-center gap-1.5 text-muted">
        <span className="text-[11px]">color by</span>
        <select
          className="rounded-md border border-line bg-panel2 px-2 py-1.5 text-text"
          value={s.colorMode}
          onChange={(e) => s.setColorMode(e.target.value as ColorMode)}
          title="Color streams by solved phase or temperature"
        >
          <option value="none">none</option>
          <option value="phase">phase</option>
          <option value="temperature">temperature</option>
        </select>
      </label>
      <label className="flex items-center gap-1.5 text-muted">
        <span className="text-[11px]">units</span>
        <select
          className="rounded-md border border-line bg-panel2 px-2 py-1.5 text-text"
          value={s.unitSet}
          onChange={(e) => s.setUnitSet(e.target.value as UnitSet)}
          title="Unit system for all inputs and outputs (engine stays SI internally)"
        >
          {UNIT_SETS.map((u) => <option key={u} value={u}>{u}</option>)}
        </select>
      </label>
      <div className="flex overflow-hidden rounded-md border border-line" role="radiogroup"
        aria-label="Diagram view">
        {VIEWS.map((v) => (
          <button key={v.id} role="radio" aria-checked={s.viewMode === v.id} title={v.title}
            className={`cursor-pointer border-0 px-2 py-1.5 text-[12px] transition-colors ${
              s.viewMode === v.id ? "bg-accent/20 text-accent" : "bg-panel2 text-muted hover:text-text"
            }`}
            onClick={() => s.setViewMode(v.id)}>
            {v.label}
          </button>
        ))}
      </div>

      <Sep />
      <Button variant="ghost" icon={<FolderKanban size={14} />}
        title="Projects & templates" onClick={s.toggleProjects}>
        Projects
      </Button>

      <Sep />
      <Button variant={s.chatOpen ? "primary" : "ghost"} icon={<Bot size={14} />}
        title="AI copilot chat" onClick={s.toggleChat}>
        Copilot
      </Button>

      <span className="ml-auto truncate text-muted" role="status">{s.status}</span>
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
