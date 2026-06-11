// Projects & templates: save/load named flowsheets (localStorage) and start
// from a curated template gallery.
import { FilePlus2, FolderOpen, Save, Trash2, X } from "lucide-react";
import { useEffect, useRef, useState } from "react";
import { useStore } from "../store";
import { TEMPLATES } from "../templates";
import { Button, PanelTitle } from "./ui";

export function ProjectsDialog() {
  const open = useStore((s) => s.projectsOpen);
  const toggle = useStore((s) => s.toggleProjects);
  const projects = useStore((s) => s.projects);
  const saveCurrent = useStore((s) => s.saveCurrentProject);
  const loadProject = useStore((s) => s.loadProject);
  const removeProject = useStore((s) => s.removeProject);
  const loadFlowDoc = useStore((s) => s.loadFlowDoc);
  const setProduct = useStore((s) => s.setProduct);
  const [name, setName] = useState("");
  const dialogRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    if (!open) return;
    const onKey = (e: KeyboardEvent) => e.key === "Escape" && toggle();
    window.addEventListener("keydown", onKey);
    dialogRef.current?.querySelector("input")?.focus();
    return () => window.removeEventListener("keydown", onKey);
  }, [open, toggle]);

  if (!open) return null;

  return (
    <div className="fixed inset-0 z-40 flex items-center justify-center bg-black/40"
      onClick={toggle}>
      <div ref={dialogRef} role="dialog" aria-modal="true" aria-label="Projects and templates"
        className="max-h-[80vh] w-[560px] overflow-auto rounded-xl border border-line bg-panel p-4 shadow-2xl"
        onClick={(e) => e.stopPropagation()}>
        <div className="mb-2 flex items-center">
          <b className="text-[15px]">Projects & templates</b>
          <button className="ml-auto cursor-pointer text-muted hover:text-text"
            onClick={toggle} aria-label="Close dialog"><X size={16} /></button>
        </div>

        <PanelTitle>Save current flowsheet</PanelTitle>
        <div className="flex items-center gap-1.5">
          <input
            className="min-w-0 flex-1 rounded-md border border-line bg-panel2 px-2 py-1.5 text-text"
            placeholder="project name"
            value={name}
            aria-label="Project name"
            onChange={(e) => setName(e.target.value)}
            onKeyDown={(e) => {
              if (e.key === "Enter" && name.trim()) {
                saveCurrent(name);
                setName("");
              }
            }}
          />
          <Button variant="primary" icon={<Save size={13} />} disabled={!name.trim()}
            onClick={() => { saveCurrent(name); setName(""); }}>
            Save
          </Button>
        </div>

        <PanelTitle>Saved projects</PanelTitle>
        {projects.length === 0 && (
          <div className="text-[12px] text-muted">Nothing saved yet.</div>
        )}
        {projects.map((p) => (
          <div key={p.name}
            className="my-1 flex items-center gap-2 rounded-md border border-line bg-panel2/60 px-2.5 py-1.5">
            <span className="min-w-0 flex-1 truncate">{p.name}</span>
            <span className="text-[11px] text-muted">
              {new Date(p.savedAt).toLocaleString()}
            </span>
            <Button variant="ghost" icon={<FolderOpen size={12} />}
              onClick={() => loadProject(p.name)} aria-label={`Open ${p.name}`}>open</Button>
            <button className="cursor-pointer p-1 text-muted hover:text-bad"
              onClick={() => removeProject(p.name)} aria-label={`Delete ${p.name}`}>
              <Trash2 size={13} />
            </button>
          </div>
        ))}

        <PanelTitle>Templates</PanelTitle>
        <div className="grid grid-cols-2 gap-2">
          {TEMPLATES.map((t) => (
            <button key={t.name}
              className="cursor-pointer rounded-lg border border-line bg-panel2 p-2.5 text-left transition-colors hover:border-accent"
              onClick={() => {
                loadFlowDoc(t.flow, `template "${t.name}"`);
                setProduct(t.product);
                toggle();
              }}>
              <div className="mb-0.5 flex items-center gap-1.5 font-semibold">
                <FilePlus2 size={13} className="text-accent" aria-hidden />{t.name}
              </div>
              <div className="text-[11.5px] leading-snug text-muted">{t.blurb}</div>
            </button>
          ))}
        </div>
      </div>
    </div>
  );
}
