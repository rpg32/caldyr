import { defineConfig } from "@playwright/test";

// E2E smoke tests. Boots the engine API and the Vite dev server (reusing them
// if already running), then drives the app in Chromium.
export default defineConfig({
  testDir: "./e2e",
  timeout: 60_000,
  use: {
    baseURL: "http://localhost:5273",
  },
  webServer: [
    {
      command: "python -m uvicorn api.main:app --port 8753",
      cwd: "..",
      url: "http://localhost:8753/health",
      reuseExistingServer: true,
      timeout: 30_000,
    },
    {
      command: "npm run dev",
      url: "http://localhost:5273",
      reuseExistingServer: true,
      timeout: 30_000,
    },
  ],
});
