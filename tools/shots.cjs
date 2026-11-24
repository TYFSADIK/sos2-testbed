// Screenshots of the 3D viewer and the Grafana dashboard (headless Chromium).
// usage: node tools/shots.cjs <playwright-core dir> <out dir>
const path = require("path");
const { chromium } = require(path.resolve(process.argv[2]));
const out = process.argv[3];
(async () => {
  const b = await chromium.launch({ args: ["--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"] });
  const pg = await b.newPage({ viewport: { width: 1600, height: 900 } });
  const errs = []; pg.on("pageerror", (e) => errs.push(e.message));
  await pg.goto("http://127.0.0.1:8000/", { waitUntil: "load", timeout: 120000 });
  await pg.waitForFunction(() => window.__sos && window.__sos.loaded, null, { timeout: 120000 });
  await pg.waitForTimeout(6000);
  await pg.screenshot({ path: `${out}/viewer-globe.png` });
  // the GEO belt in data-age colours, with a Canadian satellite selected
  await pg.click('[data-c="age"]');
  await pg.evaluate(() => document.querySelector("#geo tr.click") && document.querySelector("#geo tr.click").click());
  await pg.waitForTimeout(5000);
  await pg.screenshot({ path: `${out}/viewer-geo.png` });
  // a conjunction
  await pg.evaluate(() => document.querySelector("#conj tr.click") && document.querySelector("#conj tr.click").click());
  await pg.waitForTimeout(7000);
  await pg.screenshot({ path: `${out}/viewer-conjunction.png` });
  const g = await b.newPage({ viewport: { width: 1600, height: 1250 } });
  await g.goto("http://127.0.0.1:3000/d/sos-health?orgId=1&kiosk", { waitUntil: "networkidle", timeout: 120000 });
  await g.waitForTimeout(15000);
  await g.screenshot({ path: `${out}/grafana.png`, fullPage: true, timeout: 120000 });
  console.log(JSON.stringify({ errors: errs.slice(0, 5) }));
  await b.close();
})();
