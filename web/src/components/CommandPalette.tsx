// Ctrl+K command palette: fuzzy-searchable list of every toolbar/canvas action
// plus "Add <unit>" for each palette unit op. Opened from the global shortcut
// (see App), closed with Esc; arrows + Enter navigate. The command registry and
// fuzzy filter live in lib/commands (unit-tested); this is just the UI shell.
import { useReactFlow } from "@xyflow/react";
import { useEffect, useMemo, useRef, useState } from "react";
import { buildCommands, filterCommands, type Command } from "../lib/commands";
import { useStore } from "../store";
import type { UnitType } from "../types";

export function CommandPalette() {
  const open = useStore((s) => s.commandPaletteOpen);
  const setOpen = useStore((s) => s.setCommandPalette);
  const rf = useReactFlow();
  const [query, setQuery] = useState("");
  const [active, setActive] = useState(0);
  const inputRef = useRef<HTMLInputElement>(null);
  const listRef = useRef<HTMLDivElement>(null);

  // Rebuild per open: store actions are stable and unitTypes settle after init.
  // "Add <unit>" drops at the centre of the visible canvas via the live viewport.
  const commands = useMemo(() => {
    const s = useStore.getState();
    const addUnitAtCenter = (t: UnitType) => {
      const pane = document.querySelector<HTMLElement>(".react-flow");
      const r = pane?.getBoundingClientRect();
      const at = r
        ? rf.screenToFlowPosition({ x: r.left + r.width / 2, y: r.top + r.height / 2 })
        : undefined;
      s.addUnit("unit", t, at);
    };
    return buildCommands({ ...s, addUnitAtCenter });
  }, [open, rf]);

  const results = useMemo(() => filterCommands(commands, query), [commands, query]);

  useEffect(() => {
    if (open) {
      setQuery("");
      setActive(0);
      // focus after paint so the input is mounted
      const id = window.setTimeout(() => inputRef.current?.focus(), 0);
      return () => window.clearTimeout(id);
    }
  }, [open]);
  useEffect(() => setActive(0), [query]);
  useEffect(() => {
    listRef.current?.querySelector('[data-active="true"]')
      ?.scrollIntoView({ block: "nearest" });
  }, [active, results]);

  if (!open) return null;

  const run = (cmd?: Command) => {
    if (!cmd) return;
    setOpen(false);
    cmd.run();
  };

  const onKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === "Escape") { e.preventDefault(); setOpen(false); }
    else if (e.key === "ArrowDown") { e.preventDefault(); setActive((a) => Math.min(a + 1, results.length - 1)); }
    else if (e.key === "ArrowUp") { e.preventDefault(); setActive((a) => Math.max(a - 1, 0)); }
    else if (e.key === "Enter") { e.preventDefault(); run(results[active]); }
  };

  return (
    <div className="fixed inset-0 z-50 flex items-start justify-center bg-black/40 pt-[12vh]"
      onClick={() => setOpen(false)}>
      <div role="dialog" aria-modal="true" aria-label="Command palette"
        className="flex max-h-[70vh] w-[560px] max-w-[92vw] flex-col overflow-hidden rounded-xl border border-line bg-panel shadow-2xl"
        onClick={(e) => e.stopPropagation()} onKeyDown={onKeyDown}>
        <input
          ref={inputRef}
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          placeholder="Type a command…  (Esc to close)"
          aria-label="Command"
          className="w-full border-b border-line bg-transparent px-4 py-3 text-[14px] text-text outline-none"
        />
        <div ref={listRef} className="min-h-0 flex-1 overflow-auto py-1">
          {results.length === 0 && (
            <div className="px-4 py-6 text-center text-[12.5px] text-muted">No matching commands</div>
          )}
          {results.map((cmd, i) => (
            <button
              key={cmd.id}
              data-active={i === active}
              className={`flex w-full items-center gap-2.5 px-4 py-1.5 text-left text-[13px] ${
                i === active ? "bg-accent/15 text-accent" : "text-text hover:bg-panel2"
              }`}
              onMouseMove={() => setActive(i)}
              onClick={() => run(cmd)}
            >
              <span className="w-14 shrink-0 text-[10px] uppercase tracking-wide text-muted">{cmd.section}</span>
              <span className="min-w-0 flex-1 truncate">{cmd.title}</span>
              {cmd.hint && (
                <kbd className="shrink-0 rounded border border-line bg-panel2 px-1.5 py-0.5 text-[10px] text-muted">{cmd.hint}</kbd>
              )}
            </button>
          ))}
        </div>
      </div>
    </div>
  );
}
