"use strict";
Cesium.Ion.defaultAccessToken = "";
const C = Cesium;
const viewer = new C.Viewer("globe", {
  baseLayer: C.ImageryLayer.fromProviderAsync(C.TileMapServiceImageryProvider.fromUrl(C.buildModuleUrl("Assets/Textures/NaturalEarthII"))),
  baseLayerPicker: false, geocoder: false, homeButton: false, sceneModePicker: false, navigationHelpButton: false,
  animation: false, timeline: false, fullscreenButton: false, infoBox: false, selectionIndicator: false,
});
const scene = viewer.scene;
scene.skyAtmosphere.show = true; scene.globe.enableLighting = false;
viewer.camera.setView({ destination: C.Cartesian3.fromDegrees(-85, 25, 95000000) });
const S = { mode: "regime", pos: null, sel: null, pts: scene.primitives.add(new C.PointPrimitiveCollection()), byNorad: {} };
window.__sos = S;
const PAL = {
  regime: [["LEO", "#5ea2ff"], ["MEO", "#62d39a"], ["GEO", "#ffcf5a"], ["HEO", "#ff8a65"], ["GTO", "#c792ea"]],
  type: [["Payload", "#5ea2ff"], ["Rocket body", "#ffb454"], ["Debris", "#ff6b6b"], ["Unknown", "#9aa4b2"]],
  age: [["< 1 day", "#62d39a"], ["1-3 days", "#ffcf5a"], ["3-7 days", "#ff8a65"], ["> 7 days", "#ff4d6d"]],
};
const tracks = new C.CustomDataSource("tracks"); viewer.dataSources.add(tracks);

function colour(i) {
  const p = S.pos;
  let k;
  if (S.mode === "regime") k = p.regime[i];
  else if (S.mode === "type") k = p.type[i];
  else { const a = p.age[i]; k = a < 1 ? 0 : a < 3 ? 1 : a < 7 ? 2 : 3; }
  return C.Color.fromCssColorString(PAL[S.mode][k][1]);
}
function legend() {
  document.getElementById("legend").innerHTML = PAL[S.mode].map(([n, c]) => `<span><i class="dot" style="background:${c}"></i>${n}</span>`).join("");
}
async function refresh() {
  const p = await (await fetch("/api/positions")).json();
  S.pos = p;
  document.getElementById("n").textContent = p.n.toLocaleString();
  document.getElementById("clock").textContent = p.t.slice(11, 16);
  const seen = new Set();
  for (let i = 0; i < p.n; i++) {
    const n = p.norad[i]; seen.add(n);
    const pos = C.Cartesian3.fromDegrees(p.lon[i], p.lat[i], p.alt[i] * 1000);
    let pt = S.byNorad[n];
    if (!pt) pt = S.byNorad[n] = S.pts.add({ position: pos, pixelSize: 3, id: { norad: n } });
    pt.position = pos; pt.color = colour(i);
    pt.pixelSize = n === S.sel ? 9 : (p.regime[i] === 2 ? 4 : 3);
    pt.outlineWidth = n === S.sel ? 2 : 0; pt.outlineColor = C.Color.WHITE;
  }
  S.loaded = true;
}
async function select(n) {
  S.sel = n;
  const o = await (await fetch(`/api/object/${n}`)).json();
  tracks.entities.removeAll();
  tracks.entities.add({ polyline: { positions: C.Cartesian3.fromDegreesArrayHeights(o.track.flatMap(([x, y, h]) => [x, y, h * 1000])),
    width: 2, material: C.Color.fromCssColorString("#ffffff").withAlpha(0.8), arcType: C.ArcType.NONE } });
  const sc = o.satcat || {}, u = o.uncertainty_1sigma_km;
  const fresh = o.age_days < 1 ? "#62d39a" : o.age_days < 3 ? "#ffcf5a" : "#ff6b6b";
  document.getElementById("sel").innerHTML = `<div class="card"><b>${o.name}</b> <span class="dim">#${o.norad} · ${o.intdes}</span>
    <div>${sc.type || "?"} · owner ${sc.owner || "?"} · launched ${sc.launch || "?"} · status ${sc.status || "?"}</div>
    <div>${o.regime}: ${o.orbit.perigee_km} x ${o.orbit.apogee_km} km, ${o.orbit.inc_deg}°, ${o.orbit.period_min} min</div>
    <div>Element set age: <b style="color:${fresh}">${o.age_days} days</b> (epoch ${o.epoch.slice(0, 16)}Z)</div>
    <div>1-sigma now: radial ${u.radial} km, in-track ${u.in_track} km, cross-track ${u.cross_track} km</div>
    <div class="dim">History: ${o.history.length} element sets${o.manoeuvres.length ? `, ${o.manoeuvres.length} manoeuvre(s) detected` : ""}</div></div>`;
  refresh();
}
async function conj() {
  const d = await (await fetch("/api/conjunctions?limit=15")).json();
  if (!d.run) { document.getElementById("conjmeta").textContent = "No screening run yet."; return; }
  document.getElementById("conjmeta").textContent = `${d.conjunctions.toLocaleString()} under ${d.run.threshold_km} km in ${d.run.days} day(s) from ${String(d.run.window_start).slice(0, 16)} · ${d.run.objects.toLocaleString()} objects, ${d.run.runtime_s.toFixed(0)} s`;
  document.getElementById("conj").innerHTML = "<tr><th>Objects</th><th>Miss</th><th>Pc</th></tr>" + d.items.map((c) =>
    `<tr class="click" data-a="${c.a}" data-b="${c.b}" data-lat="${c.lat}" data-lon="${c.lon}" data-alt="${c.alt_m}"><td>${(c.name_a || c.a)}<br><span class="dim">${(c.name_b || c.b)}</span></td>
     <td>${(c.miss_km * 1000).toFixed(0)} m<br><span class="dim">${c.tca.slice(5, 16)}</span></td><td>${c.pc ? c.pc.toExponential(1) : "-"}</td></tr>`).join("");
  document.querySelectorAll("#conj tr.click").forEach((r) => r.onclick = async () => {
    const t = await (await fetch(`/api/conjunction/${r.dataset.a}/${r.dataset.b}`)).json();
    tracks.entities.removeAll();
    const cols = ["#5ea2ff", "#ff6b6b"];
    Object.values(t.tracks).forEach((tr, k) => tracks.entities.add({ polyline: { positions: C.Cartesian3.fromDegreesArrayHeights(tr.flatMap(([x, y, h]) => [x, y, h * 1000])), width: 3, material: C.Color.fromCssColorString(cols[k]), arcType: C.ArcType.NONE } }));
    viewer.camera.flyTo({ destination: C.Cartesian3.fromDegrees(+r.dataset.lon, +r.dataset.lat, +r.dataset.alt + 400000), duration: 1.5 });
  });
}
async function geo() {
  const d = await (await fetch("/api/geo")).json();
  document.getElementById("geo").innerHTML = "<tr><th>Satellite</th><th>Lon</th><th>Nearest (24 h)</th></tr>" + d.canadian
    .map((c) => { const b = d.belt.find((x) => x.norad === c.norad) || {}; const nb = c.neighbours[0] || {};
      return `<tr class="click" data-n="${c.norad}"><td>${c.name}<br><span class="dim">${b.state || ""}</span></td><td>${b.lon ?? ""}</td><td>${nb.name || ""}<br><span class="dim">${nb.min_km_24h ?? ""} km</span></td></tr>`; }).join("");
  document.querySelectorAll("#geo tr.click").forEach((r) => r.onclick = () => select(+r.dataset.n));
}
document.querySelectorAll("[data-c]").forEach((b) => b.onclick = () => {
  document.querySelectorAll("[data-c]").forEach((x) => x.classList.remove("on")); b.classList.add("on");
  S.mode = b.dataset.c; legend(); refresh();
});
let qt;
document.getElementById("q").oninput = (e) => { clearTimeout(qt); qt = setTimeout(async () => {
  const r = await (await fetch(`/api/search?q=${encodeURIComponent(e.target.value)}`)).json();
  document.getElementById("results").innerHTML = r.map((x) => `<div class="click" data-n="${x.norad}" style="cursor:pointer">${x.name} <span class="dim">#${x.norad}</span></div>`).join("");
  document.querySelectorAll("#results .click").forEach((d) => d.onclick = () => select(+d.dataset.n)); }, 250); };
new C.ScreenSpaceEventHandler(scene.canvas).setInputAction((c) => { const o = scene.pick(c.position); if (o && o.id && o.id.norad) select(o.id.norad); }, C.ScreenSpaceEventType.LEFT_CLICK);
legend(); refresh(); conj(); geo(); setInterval(refresh, 4000);
