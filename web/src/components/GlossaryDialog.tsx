// Searchable glossary / cheat-sheet of the process + techno-economic jargon the
// app uses. Opened from the toolbar; terms also appear inline as <Term> tooltips.
import { Search, X } from "lucide-react";
import { useEffect, useMemo, useRef, useState } from "react";
import { GLOSSARY } from "../lib/glossary";
import { useStore } from "../store";

export function GlossaryDialog() {
  const open = useStore((s) => s.glossaryOpen);
  const toggle = useStore((s) => s.toggleGlossary);
  const [q, setQ] = useState("");
  const inputRef = useRef<HTMLInputElement>(null);

  useEffect(() => {
    if (!open) return;
    const onKey = (e: KeyboardEvent) => e.key === "Escape" && toggle();
    window.addEventListener("keydown", onKey);
    inputRef.current?.focus();
    return () => window.removeEventListener("keydown", onKey);
  }, [open, toggle]);

  const rows = useMemo(() => {
    const needle = q.trim().toLowerCase();
    const list = [...GLOSSARY].sort((a, b) => a.term.localeCompare(b.term));
    if (!needle) return list;
    return list.filter((e) =>
      e.term.toLowerCase().includes(needle)
      || e.def.toLowerCase().includes(needle)
      || (e.aka ?? []).some((a) => a.toLowerCase().includes(needle)));
  }, [q]);

  if (!open) return null;

  return (
    <div className="fixed inset-0 z-40 flex items-center justify-center bg-black/40" onClick={toggle}>
      <div role="dialog" aria-modal="true" aria-label="Glossary"
        className="flex max-h-[82vh] w-[560px] flex-col rounded-xl border border-line bg-panel p-4 shadow-2xl"
        onClick={(e) => e.stopPropagation()}>
        <div className="mb-2 flex items-center">
          <b className="text-[15px]">Glossary</b>
          <button className="ml-auto cursor-pointer text-muted hover:text-text"
            onClick={toggle} aria-label="Close glossary"><X size={16} /></button>
        </div>
        <div className="mb-2 flex items-center gap-1.5 rounded-md border border-line bg-panel2 px-2">
          <Search size={13} className="text-muted" aria-hidden />
          <input ref={inputRef} value={q} onChange={(e) => setQ(e.target.value)}
            placeholder="Search terms (LCOP, reflux, pinch…)" aria-label="Search glossary"
            className="w-full bg-transparent py-1.5 text-text outline-none" />
        </div>
        <div className="overflow-auto">
          {rows.length === 0 && <div className="p-3 text-muted">No match for “{q}”.</div>}
          <dl className="m-0">
            {rows.map((e) => (
              <div key={e.term} className="border-b border-line/60 py-1.5">
                <dt className="text-[13px] font-semibold text-text">
                  {e.term}
                  {e.aka?.length ? <span className="ml-2 text-[11px] font-normal text-muted">{e.aka.join(", ")}</span> : null}
                </dt>
                <dd className="m-0 text-[12px] leading-relaxed text-muted">{e.def}</dd>
              </div>
            ))}
          </dl>
        </div>
      </div>
    </div>
  );
}
