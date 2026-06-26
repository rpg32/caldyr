import { ArrowRightFromLine, ArrowRightToLine } from "lucide-react";
import { useState, type MouseEvent } from "react";
import { useStore } from "../store";
import { PanelTitle } from "./ui";

interface Tip { text: string; top: number; left: number; maxH: number }

export function Palette() {
  const unitTypes = useStore((s) => s.unitTypes);
  const addUnit = useStore((s) => s.addUnit);
  const [tip, setTip] = useState<Tip | null>(null);

  const item =
    "w-full rounded-md border border-line bg-panel2 px-2 py-1.5 text-left text-[13px] " +
    "cursor-pointer hover:border-accent transition-colors";

  // A full, wrapping tooltip positioned just right of the hovered item. `fixed`
  // (viewport-relative) so it isn't clipped by the palette's own scroll overflow.
  const show = (text: string) => (e: MouseEvent<HTMLButtonElement>) => {
    if (!text) return;
    const r = e.currentTarget.getBoundingClientRect();
    setTip({ text, top: r.top, left: r.right + 8, maxH: window.innerHeight - r.top - 12 });
  };
  const hide = () => setTip(null);

  return (
    <aside className="flex flex-col gap-1 overflow-auto border-r border-line bg-panel p-2">
      <PanelTitle>Boundaries</PanelTitle>
      <button className={item} onClick={() => addUnit("feed")}
        onMouseEnter={show("Material feed with a temperature, pressure, molar-flow and composition spec — the boundary source of the flowsheet.")}
        onMouseLeave={hide}>
        <ArrowRightFromLine size={12} className="mr-1 inline text-ok" aria-hidden />
        Feed
      </button>
      <button className={item} onClick={() => addUnit("product")}
        onMouseEnter={show("Product sink (boundary outlet). The stream feeding it is reported by the engine after a solve; pick one as the costing product.")}
        onMouseLeave={hide}>
        <ArrowRightToLine size={12} className="mr-1 inline text-muted" aria-hidden />
        Product
      </button>

      <PanelTitle>Unit ops</PanelTitle>
      {unitTypes.map((t) => (
        <button key={t.type} className={item} onClick={() => addUnit("unit", t)}
          onMouseEnter={show(t.description || t.doc)} onMouseLeave={hide}>
          {t.type}
        </button>
      ))}

      {tip && (
        <div
          className="pointer-events-none fixed z-50 w-72 overflow-auto rounded-md border border-line bg-panel p-2.5 text-[11px] leading-relaxed text-text shadow-2xl"
          style={{ top: tip.top, left: tip.left, maxHeight: tip.maxH }}
          role="tooltip"
        >
          {tip.text}
        </div>
      )}
    </aside>
  );
}
