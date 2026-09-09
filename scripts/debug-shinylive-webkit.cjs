const { createRequire } = require("node:module");

const requirePlaywright = createRequire(
  "/Users/shaheedazaad/Projects/weekli_monorepo/node_modules/playwright-core/package.json",
);
const { chromium, webkit } = requirePlaywright("playwright-core");

(async () => {
  const engine = process.env.BROWSER_ENGINE === "chromium" ? chromium : webkit;
  const browser = await engine.launch({ headless: true });
  const page = await browser.newPage();
  const errors = [];
  const failedResponses = [];

  page.on("console", (message) => {
    if (message.type() === "error") errors.push(message.text());
  });
  page.on("pageerror", (error) => errors.push(error.stack || `${error.name}: ${error.message}`));
  page.on("response", (response) => {
    if (!response.ok()) failedResponses.push(`${response.status()} ${response.url()}`);
  });

  if (process.env.MINIMAL_APP === "1") {
    await page.route("**/app.json", (route) =>
      route.fulfill({
        contentType: "application/json",
        body: JSON.stringify([
          {
            name: "app.R",
            type: "text",
            content: "library(shiny)\nshinyApp(fluidPage('MINIMAL_OK'), function(input, output, session) {})",
          },
        ]),
      }),
    );
  }

  await page.goto("http://127.0.0.1:8765/", { waitUntil: "domcontentloaded" });
  await page.waitForTimeout(15_000);

  const loading = await page.locator(".loading-wrapper").count();
  const appFrames = page.frames().length - 1;
  const exactError = errors.some((error) =>
    error.includes("Maximum call stack size exceeded") ||
    error.includes("Cannot destructure property 'type'"),
  );

  console.log(JSON.stringify({ loading, appFrames, errors, failedResponses }, null, 2));
  await browser.close();

  if (loading > 0 || exactError) process.exitCode = 1;
})().catch((error) => {
  console.error(error);
  process.exitCode = 2;
});
