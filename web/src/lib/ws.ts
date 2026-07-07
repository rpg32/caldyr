// Minimal WebSocket helpers for the API's /ws/* channels (proxied by Vite).

export const wsUrl = (path: string): string => {
  const proto = location.protocol === "https:" ? "wss" : "ws";
  return `${proto}://${location.host}/api${path}`;
};

/** Open a socket, send one JSON message, and stream parsed events to the
 * callback until `done` returns true for one of them (then resolve with it). */
export function requestOverWs<T extends { type: string }>(
  path: string,
  payload: unknown,
  onEvent: (e: T) => void,
  isFinal: (e: T) => boolean,
): Promise<T> {
  return new Promise((resolve, reject) => {
    const ws = new WebSocket(wsUrl(path));
    ws.onerror = () => reject(new Error("WebSocket connection failed"));
    ws.onclose = (ev) => {
      if (!ev.wasClean) reject(new Error("WebSocket closed unexpectedly"));
    };
    ws.onopen = () => ws.send(JSON.stringify(payload));
    ws.onmessage = (m) => {
      let e: T;
      try {
        e = JSON.parse(m.data as string) as T;
      } catch {
        ws.close();
        reject(new Error("malformed message from the engine"));
        return;
      }
      onEvent(e);
      if (isFinal(e)) {
        ws.close();
        resolve(e);
      }
    };
  });
}

/** A persistent chat socket: send() per user message, events stream back. */
export class ChatSocket {
  private ws: WebSocket | null = null;

  connect(onEvent: (e: Record<string, unknown> & { type: string }) => void,
          onClose: () => void): Promise<void> {
    return new Promise((resolve, reject) => {
      const ws = new WebSocket(wsUrl("/ws/chat"));
      ws.onopen = () => {
        this.ws = ws;
        resolve();
      };
      ws.onerror = () => reject(new Error("chat socket failed to connect"));
      ws.onmessage = (m) => {
        let e: (Record<string, unknown> & { type: string }) | null = null;
        try {
          e = JSON.parse(m.data as string);
        } catch {
          return; // malformed frame: drop it rather than crash the handler
        }
        if (e && typeof e === "object") onEvent(e);
      };
      ws.onclose = () => {
        this.ws = null;
        onClose();
      };
    });
  }

  get connected(): boolean {
    return this.ws?.readyState === WebSocket.OPEN;
  }

  send(payload: unknown): void {
    if (!this.connected) throw new Error("chat socket is not connected");
    this.ws!.send(JSON.stringify(payload));
  }

  close(): void {
    this.ws?.close();
  }
}
