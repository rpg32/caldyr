// Smoke test: the M5 definition-of-done (build + solve + cost in the browser)
// must keep holding as the UI evolves.
import { expect, test, type Page } from "@playwright/test";

test.beforeEach(async ({ page }) => {
  // skip the first-run tour in the fresh Playwright profile
  await page.addInitScript(() => localStorage.setItem("caldyr.tour_seen", "1"));
});

async function loadAmmoniaTemplate(page: Page) {
  await expect(page.getByRole("button", { name: "Mixer", exact: true })).toBeVisible();
  await page.getByRole("button", { name: "Projects" }).click();
  await page.getByRole("button", { name: /Ammonia loop/ }).click();
  await expect(page.locator(".node")).toHaveCount(9);
}

test("load template, solve, and cost in the browser", async ({ page }) => {
  await page.goto("/");
  await loadAmmoniaTemplate(page);

  // Solve converges.
  await page.getByRole("button", { name: "Solve" }).click();
  await expect(page.getByRole("status")).toContainText("Solved", { timeout: 30_000 });
  await expect(page.getByText("converged")).toBeVisible();

  // Cost reports an LCOP.
  await page.getByRole("button", { name: "Cost" }).click();
  await expect(page.getByRole("status")).toContainText("LCOP", { timeout: 30_000 });
  await expect(page.getByText("EQUIPMENT (INSTALLED)")).toBeVisible();
});

test("undo restores a deleted unit", async ({ page }) => {
  await page.goto("/");
  await loadAmmoniaTemplate(page);

  await page.getByRole("button", { name: "Mixer", exact: true }).click();
  await expect(page.locator(".node")).toHaveCount(10);

  await page.getByRole("button", { name: "Undo" }).click();
  await expect(page.locator(".node")).toHaveCount(9);
});

test("view modes and projects dialog", async ({ page }) => {
  await page.goto("/");
  await loadAmmoniaTemplate(page);

  // BFD strips glyphs and edge labels
  await page.getByRole("radio", { name: "BFD" }).click();
  await expect(page.locator(".node-bfd")).toHaveCount(9);
  await expect(page.locator(".edge-label")).toHaveCount(0);
  await page.getByRole("radio", { name: "PFD" }).click();
  await expect(page.locator(".node-bfd")).toHaveCount(0);

  // save + reopen a named project
  await page.getByRole("button", { name: "Projects" }).click();
  const dialog = page.getByRole("dialog", { name: "Projects and templates" });
  await dialog.getByLabel("Project name").fill("e2e-project");
  await dialog.getByRole("button", { name: "Save", exact: true }).click();
  await expect(dialog.getByText("e2e-project")).toBeVisible();
});
