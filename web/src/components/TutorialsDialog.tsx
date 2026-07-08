// Guided tutorials launcher: each entry loads its starting flowsheet onto the
// canvas and opens the matching step-by-step guide, so a new user can follow
// along in the app instead of only reading the docs.
import { ExternalLink, GraduationCap, Play, X } from "lucide-react";
import { useEffect, useRef } from "react";
import { useStore } from "../store";
import { TEMPLATES } from "../templates";
import { Button, PanelTitle } from "./ui";

const DOCS_BASE = "https://rpg32.github.io/caldyr";

interface Tutorial {
  title: string;
  goal: string;
  time: string;
  template: string;   // must match a TEMPLATES entry name
  slug: string;       // docs page slug
}

const TUTORIALS: Tutorial[] = [
  {
    title: "1 · Your first flowsheet: ammonia loop",
    goal: "Build a Haber-Bosch synthesis loop with recycle and purge, solve the "
      + "tear stream, then size and cost it to a levelized cost of ammonia.",
    time: "~15 min",
    template: "Ammonia loop",
    slug: "tutorial-ammonia-loop",
  },
  {
    title: "2 · Design and cost a distillation column",
    goal: "Take a benzene/toluene split from a shortcut (FUG) design through a "
      + "rigorous MESH check, then price the reflux-ratio trade-off.",
    time: "~20 min",
    template: "Benzene/toluene column",
    slug: "tutorial-distillation-tea",
  },
  {
    title: "3 · Optimize against an economic objective",
    goal: "Let the solver pick design variables that minimize cost, then quantify "
      + "the answer with tornado sensitivity and Monte-Carlo uncertainty.",
    time: "~20 min",
    template: "Ammonia loop",
    slug: "tutorial-optimization",
  },
];

export function TutorialsDialog() {
  const open = useStore((s) => s.tutorialsOpen);
  const toggle = useStore((s) => s.toggleTutorials);
  const loadFlowDoc = useStore((s) => s.loadFlowDoc);
  const setProduct = useStore((s) => s.setProduct);
  const dialogRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    if (!open) return;
    const onKey = (e: KeyboardEvent) => e.key === "Escape" && toggle();
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [open, toggle]);

  if (!open) return null;

  const guideUrl = (slug: string) => `${DOCS_BASE}/${slug}/`;
  const openGuide = (slug: string) =>
    window.open(guideUrl(slug), "_blank", "noopener,noreferrer");

  const start = (t: Tutorial) => {
    const tpl = TEMPLATES.find((x) => x.name === t.template);
    if (tpl) {
      loadFlowDoc(tpl.flow, `tutorial "${t.title}"`);
      setProduct(tpl.product);
    }
    openGuide(t.slug);
    toggle();
  };

  return (
    <div className="fixed inset-0 z-40 flex items-center justify-center bg-black/40"
      onClick={toggle}>
      <div ref={dialogRef} role="dialog" aria-modal="true" aria-label="Guided tutorials"
        className="max-h-[80vh] w-[560px] overflow-auto rounded-xl border border-line bg-panel p-4 shadow-2xl"
        onClick={(e) => e.stopPropagation()}>
        <div className="mb-1 flex items-center gap-2">
          <GraduationCap size={16} className="text-accent" aria-hidden />
          <b className="text-[15px]">Guided tutorials</b>
          <button className="ml-auto cursor-pointer text-muted hover:text-text"
            onClick={toggle} aria-label="Close dialog"><X size={16} /></button>
        </div>
        <p className="mb-1 text-[12px] leading-snug text-muted">
          Each tutorial loads its starting flowsheet onto the canvas and opens the
          step-by-step guide in a new tab. Work along in the app as you read.
        </p>

        <PanelTitle>Start a tutorial</PanelTitle>
        <div className="flex flex-col gap-2">
          {TUTORIALS.map((t) => (
            <div key={t.slug}
              className="rounded-lg border border-line bg-panel2 p-3">
              <div className="mb-0.5 flex items-baseline gap-2">
                <span className="font-semibold">{t.title}</span>
                <span className="text-[11px] text-muted">{t.time}</span>
              </div>
              <p className="mb-2 text-[12px] leading-snug text-muted">{t.goal}</p>
              <div className="flex items-center gap-2">
                <Button variant="primary" icon={<Play size={13} />}
                  onClick={() => start(t)}>
                  Load flowsheet &amp; open guide
                </Button>
                <Button variant="ghost" icon={<ExternalLink size={12} />}
                  onClick={() => openGuide(t.slug)}
                  title="Open the guide without changing the canvas">
                  Guide only
                </Button>
              </div>
            </div>
          ))}
        </div>

        <PanelTitle>More</PanelTitle>
        <a className="inline-flex items-center gap-1.5 text-[12.5px] text-accent hover:underline"
          href={`${DOCS_BASE}/tutorials/`} target="_blank" rel="noopener noreferrer">
          <ExternalLink size={12} aria-hidden /> Full example catalog and docs
        </a>
      </div>
    </div>
  );
}
