/// <reference types="vitest/config" />
import tailwindcss from "@tailwindcss/vite";
import react from "@vitejs/plugin-react";
import { defineConfig } from "vite";

// The engine API runs on :8753; proxy /api there so the browser hits one origin.
export default defineConfig({
  plugins: [react(), tailwindcss()],
  server: {
    // Distinct port so it doesn't clash with other local Vite apps.
    port: 5273,
    strictPort: true,
    proxy: {
      "/api": {
        target: "http://localhost:8753",
        changeOrigin: true,
        ws: true, // /ws/solve and /ws/chat ride the same proxy
        rewrite: (p) => p.replace(/^\/api/, ""),
      },
    },
  },
  preview: { port: 5274, strictPort: true },
  test: {
    environment: "node",
    include: ["src/**/*.test.ts"],
  },
});
