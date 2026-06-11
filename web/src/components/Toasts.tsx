import { AlertCircle, CheckCircle2, Info, X } from "lucide-react";
import { useEffect } from "react";
import { useStore, type Toast } from "../store";

const ICONS = {
  error: <AlertCircle size={15} className="shrink-0 text-bad" aria-hidden />,
  success: <CheckCircle2 size={15} className="shrink-0 text-ok" aria-hidden />,
  info: <Info size={15} className="shrink-0 text-accent" aria-hidden />,
};

function ToastItem({ toast }: { toast: Toast }) {
  const dismiss = useStore((s) => s.dismissToast);
  useEffect(() => {
    const t = setTimeout(() => dismiss(toast.id), 6000);
    return () => clearTimeout(t);
  }, [toast.id, dismiss]);
  return (
    <div
      role="alert"
      className="pointer-events-auto flex max-w-[380px] items-start gap-2 rounded-lg border border-line bg-panel px-3 py-2 shadow-lg"
    >
      {ICONS[toast.kind]}
      <span className="min-w-0 break-words text-[13px]">{toast.msg}</span>
      <button
        className="ml-auto cursor-pointer text-muted hover:text-text"
        onClick={() => dismiss(toast.id)}
        aria-label="Dismiss notification"
      >
        <X size={13} />
      </button>
    </div>
  );
}

export function Toasts() {
  const toasts = useStore((s) => s.toasts);
  if (!toasts.length) return null;
  return (
    <div className="pointer-events-none fixed bottom-4 right-4 z-50 flex flex-col gap-2">
      {toasts.map((t) => <ToastItem key={t.id} toast={t} />)}
    </div>
  );
}
