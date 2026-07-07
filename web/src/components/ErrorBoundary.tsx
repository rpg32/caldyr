// Top-level crash containment. The app restores the autosaved session on every
// boot, so a document that throws during render would otherwise white-screen
// the UI on every launch until the user manually clears localStorage. The
// fallback offers exactly that escape hatch.
import { Component, type ErrorInfo, type ReactNode } from "react";
import { clearAutosave } from "../lib/persist";

interface Props {
  children: ReactNode;
}

interface State {
  error: Error | null;
}

export class ErrorBoundary extends Component<Props, State> {
  state: State = { error: null };

  static getDerivedStateFromError(error: Error): State {
    return { error };
  }

  componentDidCatch(error: Error, info: ErrorInfo): void {
    console.error("caldyr crashed:", error, info.componentStack);
  }

  render(): ReactNode {
    if (!this.state.error) return this.props.children;
    return (
      <div style={{
        display: "flex", flexDirection: "column", alignItems: "center",
        justifyContent: "center", height: "100vh", gap: "12px",
        fontFamily: "system-ui, sans-serif", background: "#0b0e14", color: "#e6e6e6",
        padding: "24px", textAlign: "center",
      }}>
        <h1 style={{ fontSize: "20px", margin: 0 }}>Caldyr hit an unexpected error</h1>
        <p style={{ margin: 0, opacity: 0.8, maxWidth: "48em" }}>
          {this.state.error.message}
        </p>
        <div style={{ display: "flex", gap: "8px", marginTop: "8px" }}>
          <button
            onClick={() => location.reload()}
            style={{
              padding: "8px 16px", borderRadius: "6px", border: "1px solid #444",
              background: "#1a1f2b", color: "#e6e6e6", cursor: "pointer",
            }}
          >
            Reload
          </button>
          <button
            onClick={() => { clearAutosave(); location.reload(); }}
            style={{
              padding: "8px 16px", borderRadius: "6px", border: "1px solid #444",
              background: "#2b1a1a", color: "#e6e6e6", cursor: "pointer",
            }}
          >
            Clear saved session &amp; reload
          </button>
        </div>
        <p style={{ margin: 0, opacity: 0.6, fontSize: "13px", maxWidth: "48em" }}>
          "Clear saved session" discards the autosaved flowsheet (named projects
          are kept). If this keeps happening, please open an issue with the
          message above.
        </p>
      </div>
    );
  }
}
