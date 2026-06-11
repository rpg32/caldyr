// First-run onboarding: a short sequence of coachmarks anchored to the UI.
// Shows once (flag in localStorage); re-launchable from the help button.
import { X } from "lucide-react";
import { useEffect, useState } from "react";
import { loadPref, savePref } from "../lib/persist";
import { Button } from "./ui";

interface Step {
  title: string;
  body: string;
  anchor: string; // CSS selector to highlight (best-effort)
}

const STEPS: Step[] = [
  {
    title: "Build",
    body: "Add a Feed, unit operations, and Products from the palette, then drag "
      + "between ports to connect streams. Double-click anything to rename it.",
    anchor: "aside.border-r",
  },
  {
    title: "Configure",
    body: "Click a unit or stream to edit it in the inspector. With nothing "
      + "selected, the Params tab holds the component list, property package, "
      + "and logical ops (Set/Adjust).",
    anchor: "aside.border-l",
  },
  {
    title: "Solve & analyze",
    body: "Solve runs the engine (watch live convergence in the status bar); "
      + "Cost runs the full techno-economic analysis. The Streams/Econ/Opt/Study "
      + "tabs hold tables, plots, optimization, and parameter sweeps.",
    anchor: "header",
  },
  {
    title: "Views",
    body: "Switch BFD / PFD / P&ID for block, process, and instrumentation "
      + "views. Color streams by phase or temperature; Arrange auto-lays-out "
      + "the flowsheet; Group collapses sections into blocks.",
    anchor: "[role=radiogroup]",
  },
  {
    title: "Copilot",
    body: "The Copilot chats with a local LLM that can build, edit, solve, and "
      + "explain flowsheets. Every edit it proposes arrives as a diff you "
      + "accept or reject.",
    anchor: "header",
  },
];

const SEEN_KEY = "tour_seen";

export function Tour() {
  const [step, setStep] = useState<number | null>(
    loadPref(SEEN_KEY) ? null : 0,
  );

  // The toolbar's help button re-launches the tour without a reload.
  useEffect(() => {
    const h = () => setStep(0);
    window.addEventListener("caldyr:tour", h);
    return () => window.removeEventListener("caldyr:tour", h);
  }, []);

  useEffect(() => {
    if (step === null) return;
    const el = document.querySelector(STEPS[step].anchor);
    el?.classList.add("tour-highlight");
    return () => el?.classList.remove("tour-highlight");
  }, [step]);

  if (step === null) return null;
  const s = STEPS[step];
  const done = () => {
    savePref(SEEN_KEY, "1");
    setStep(null);
  };

  return (
    <div className="fixed bottom-6 left-1/2 z-50 w-[420px] -translate-x-1/2 rounded-xl border border-accent/60 bg-panel p-3.5 shadow-2xl"
      role="dialog" aria-label="Onboarding tour">
      <div className="mb-1 flex items-center gap-2">
        <b>{s.title}</b>
        <span className="text-[11px] text-muted">{step + 1} / {STEPS.length}</span>
        <button className="ml-auto cursor-pointer text-muted hover:text-text"
          onClick={done} aria-label="Dismiss tour"><X size={14} /></button>
      </div>
      <p className="mb-2 text-[12.5px] leading-relaxed text-muted">{s.body}</p>
      <div className="flex items-center gap-2">
        {step > 0 && (
          <Button variant="ghost" onClick={() => setStep(step - 1)}>Back</Button>
        )}
        <span className="ml-auto" />
        <Button variant="ghost" onClick={done}>Skip</Button>
        {step < STEPS.length - 1 ? (
          <Button variant="primary" onClick={() => setStep(step + 1)}>Next</Button>
        ) : (
          <Button variant="primary" onClick={done}>Done</Button>
        )}
      </div>
    </div>
  );
}

export function relaunchTour(): void {
  window.dispatchEvent(new Event("caldyr:tour"));
}
