// CI / benchmarks page (website/pages/ci.md): reads the JSON-lines results that
// docs/benchmarks/ci/run_ci.py appends, one record per (commit, dataset), and draws
// Latest / History / Table / Quantification views. No dependencies: plots are plain SVG.
// `?src=<url>` previews another results file (a fork's benchmarks branch, a local file).
(function () {
  const METRICS = [
    ["psms", "PSMs (1% FDR)"], ["peptides", "peptides"], ["protein_groups", "protein groups"],
    ["true_fdp", "true FDP", pct], ["wall_s", "runtime", dur], ["peak_rss_gb", "peak RSS", gb],
    ["cpu_s", "CPU time", dur], ["throughput", "throughput", (v) => num(v) + " sp/s"],
  ];
  const HISTORY = [["psms", "PSMs (1% FDR)"], ["wall_s", "runtime (s)"], ["peak_rss_gb", "peak RSS (GB)"]];
  const COLORS = { HUMAN: "#0072B2", YEAST: "#E69F00", ECOLI: "#8E6BE8" };

  function num(v) { return v == null ? "—" : Math.round(v).toLocaleString("en-US"); }
  function pct(v) { return v == null ? "—" : (100 * v).toFixed(2) + "%"; }
  function gb(v) { return v == null ? "—" : v.toFixed(2) + " GB"; }
  function f3(v) { return v == null ? "—" : v.toFixed(3); }
  function dur(s) {
    if (s == null) return "—";
    return s < 60 ? s.toFixed(1) + "s" : Math.floor(s / 60) + "m " + Math.round(s % 60) + "s";
  }
  function esc(s) { return String(s ?? "").replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]); }
  function when(r) { return new Date(r.date || r.run_at); }
  function sha(r) {
    return `<a class="ci-sha" href="https://github.com/bigbio/andes/commit/${esc(r.commit)}">${esc((r.commit || "").slice(0, 7))}</a>`;
  }
  function status(r) {
    return r.status === "ok" ? '<span class="ci-ok">OK</span>'
      : `<span class="ci-bad" title="${esc(r.error)}">${esc((r.status || "?").toUpperCase())}</span>`;
  }

  function latest(rows, sets) {
    return '<div class="ci-cards">' + sets.map((d) => {
      const r = rows.filter((x) => x.dataset === d).at(-1);
      const cells = METRICS.filter(([k]) => r[k] != null)
        .map(([k, label, fmt]) => `<div><small>${label}</small><b>${(fmt || num)(r[k])}</b></div>`).join("");
      return `<div class="ci-card"><h4>${esc(d)} <span>· ${esc(r.description)}</span></h4>
        <p class="ci-msg">${esc(r.message)}</p>${sha(r)}
        <div class="ci-grid">${cells}<div><small>status</small><b>${status(r)}</b></div></div></div>`;
    }).join("") + "</div>";
  }

  // One line chart: x = commit date, y from 0, a point per run.
  function chart(pts, label) {
    const W = 320, H = 170, L = 52, R = 10, T = 26, B = 22;
    if (!pts.length) return `<div class="ci-chart"><small>${label}</small><p>no data</p></div>`;
    const xs = pts.map((p) => p.x), ymax = Math.max(...pts.map((p) => p.y)) * 1.15 || 1;
    const x0 = Math.min(...xs), x1 = Math.max(...xs), span = x1 - x0 || 1;
    const X = (x) => L + (pts.length === 1 ? (W - L - R) / 2 : ((x - x0) / span) * (W - L - R));
    const Y = (y) => H - B - (y / ymax) * (H - T - B);
    let g = "";
    for (let i = 0; i <= 4; i++) {
      const v = (ymax * i) / 4;
      g += `<line x1="${L}" x2="${W - R}" y1="${Y(v)}" y2="${Y(v)}" class="ci-grid-line"/>` +
        `<text x="${L - 4}" y="${Y(v) + 3}" text-anchor="end">${v >= 100 ? num(v) : +v.toPrecision(2)}</text>`;
    }
    const fmtd = (x) => new Date(x).toLocaleDateString("en-US", { month: "short", day: "numeric" });
    g += `<text x="${L}" y="${H - 6}">${fmtd(x0)}</text><text x="${W - R}" y="${H - 6}" text-anchor="end">${fmtd(x1)}</text>`;
    const line = pts.map((p) => `${X(p.x)},${Y(p.y)}`).join(" ");
    const dots = pts.map((p) => {
      const tipText = `${label}: ${p.y.toLocaleString("en-US")}\n${p.t}`;
      return `<a href="https://github.com/bigbio/andes/commit/${esc(p.c)}" aria-label="${esc(tipText.replace("\n", ", "))}" data-tip="${esc(tipText)}"><circle cx="${X(p.x)}" cy="${Y(p.y)}" r="4"/></a>`;
    }).join("");
    return `<div class="ci-chart"><svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${label}">
      <text x="${L - 40}" y="14" class="ci-label">${label}</text>${g}<polyline points="${line}"/>${dots}</svg></div>`;
  }

  function historyView(rows, sets) {
    return sets.map((d) => {
      const rs = rows.filter((r) => r.dataset === d && r.status === "ok");
      const charts = HISTORY.map(([k, label]) => chart(rs.filter((r) => r[k] != null)
        .map((r) => ({ x: +when(r), y: r[k], c: r.commit, t: `${(r.commit || "").slice(0, 7)} · ${when(r).toLocaleDateString("en-CA")} · ${r.message || ""}` })), label)).join("");
      return `<h3>${esc(d)} <span class="ci-dim">· ${esc(rs.at(-1)?.description)}</span></h3><div class="ci-charts">${charts}</div>`;
    }).join("");
  }

  function table(rows, sets, el) {
    const pick = el.dataset.table || "all";
    const opts = ["all", ...sets].map((d) => `<option ${d === pick ? "selected" : ""}>${esc(d)}</option>`).join("");
    const body = rows.filter((r) => pick === "all" || r.dataset === pick).slice().reverse().map((r) =>
      `<tr><td>${when(r).toLocaleString("en-CA")}</td><td>${sha(r)}</td>
       <td>${esc(r.dataset)}<br><small>${esc(r.description)}</small></td>
       <td>${num(r.psms)}</td><td>${num(r.peptides)}</td><td>${num(r.protein_groups)}</td>
       <td>${dur(r.wall_s)}</td><td>${gb(r.peak_rss_gb)}</td><td><small>${esc(r.host)}</small></td><td>${status(r)}</td></tr>`).join("");
    return `<label>dataset <select data-act="table">${opts}</select></label>
      <table class="ci-table"><thead><tr><th>when</th><th>commit</th><th>dataset</th><th>PSMs</th>
      <th>peptides</th><th>protein groups</th><th>runtime</th><th>RSS</th><th>host</th><th>status</th></tr></thead>
      <tbody>${body}</tbody></table>`;
  }

  // Box plot of log2 ratios per species, diamond = expected ratio.
  function boxes(ratios) {
    const W = 760, rowH = 56, L = 150, R = 80, T = 10, n = ratios.length, H = T + n * rowH + 44;
    const vals = ratios.flatMap((g) => [g.lo, g.hi, g.expected ?? g.median]);
    const lo = Math.floor(Math.min(...vals) * 2) / 2 - 0.25, hi = Math.ceil(Math.max(...vals) * 2) / 2 + 0.25;
    const X = (v) => L + ((v - lo) / (hi - lo)) * (W - L - R);
    let s = "";
    for (let v = Math.ceil(lo * 2) / 2; v <= hi; v += 0.5) {
      s += `<line x1="${X(v)}" x2="${X(v)}" y1="${T}" y2="${T + n * rowH}" class="ci-grid-line${v === 0 ? " ci-zero" : ""}"/>` +
        `<text x="${X(v)}" y="${T + n * rowH + 16}" text-anchor="middle">${v.toFixed(1)}</text>`;
    }
    ratios.forEach((g, i) => {
      const y = T + i * rowH + rowH / 2, c = COLORS[g.group] || "#888";
      s += `<g stroke="${c}" fill="none" data-tip="${esc(g.group)} ${esc(g.comparison)}: median ${f3(g.median)}, middle 50% ${f3(g.q1)} to ${f3(g.q3)}, expected ${f3(g.expected)}, n=${num(g.n)}">
        <line x1="${X(g.lo)}" x2="${X(g.q1)}" y1="${y}" y2="${y}"/><line x1="${X(g.q3)}" x2="${X(g.hi)}" y1="${y}" y2="${y}"/>
        <line x1="${X(g.lo)}" x2="${X(g.lo)}" y1="${y - 7}" y2="${y + 7}"/><line x1="${X(g.hi)}" x2="${X(g.hi)}" y1="${y - 7}" y2="${y + 7}"/>
        <rect x="${X(g.q1)}" y="${y - 13}" width="${Math.max(1, X(g.q3) - X(g.q1))}" height="26" fill="${c}" fill-opacity="0.25"/>
        <line x1="${X(g.median)}" x2="${X(g.median)}" y1="${y - 13}" y2="${y + 13}" stroke-width="2"/></g>` +
        (g.expected == null ? "" : `<path d="M${X(g.expected)} ${y - 7}l7 7-7 7-7-7z" class="ci-expected"/>`) +
        `<text x="${L - 12}" y="${y + 4}" text-anchor="end" class="ci-label">${esc(g.group)} ${esc(g.comparison)}</text>` +
        `<text x="${W - R + 12}" y="${y + 4}">n=${num(g.n)}</text>`;
    });
    s += `<text x="${(L + W - R) / 2}" y="${H - 4}" text-anchor="middle" class="ci-label">log2 ratio (condition / reference)</text>`;
    return `<svg viewBox="0 0 ${W} ${H}" class="ci-boxes">${s}</svg>`;
  }

  function quant(rows, el) {
    // Keep a run without ratio estimates: its counts and CVs are still results.
    const q = rows.filter((r) => r.quant);
    if (!q.length) return "<p>No quantification runs yet.</p>";
    const sets = [...new Set(q.map((r) => r.dataset))];
    const set = sets.includes(el.dataset.qset) ? el.dataset.qset : sets[0];
    const runs = q.filter((r) => r.dataset === set).reverse();
    const id = (x) => `${x.commit}|${x.run_at}`;   // a commit can be benchmarked more than once
    const r = runs.find((x) => id(x) === el.dataset.qcommit) || runs[0], Q = r.quant;
    const ratios = Q.ratios || [], cv = Q.cv || [];
    const opt = (v, cur, label) => `<option value="${esc(v)}" ${v === cur ? "selected" : ""}>${esc(label)}</option>`;
    const ratioRows = ratios.map((g) => `<tr><td>${esc(g.group)}</td><td>${esc(g.comparison)}</td><td>${num(g.n)}</td>
      <td>${f3(g.median)}</td><td>${f3(g.q1)} to ${f3(g.q3)}</td><td>${f3(g.expected)}</td><td>${f3(g.bias)}</td><td>${f3(g.mae)}</td></tr>`).join("");
    const cvRows = cv.map((c) => `<tr><td>${esc(c.group)}</td><td>${esc(c.condition)}</td><td>${num(c.n)}</td><td>${(100 * c.median_cv).toFixed(1)}%</td></tr>`).join("");
    const legend = Object.entries(COLORS).map(([k, c]) => `<span><i style="background:${c}"></i>${k}</span>`).join("") +
      '<span><i class="ci-expected-key"></i>expected ratio</span>';
    return `<div class="ci-card ci-quant"><h4>Quantification accuracy and precision</h4>
      <div class="ci-controls"><label>Benchmark set <select data-act="qset">${sets.map((s) => opt(s, set, s)).join("")}</select></label>
      <label>Commit <select data-act="qcommit">${runs.map((x) => opt(id(x), id(r), x.commit.slice(0, 7) + " · " + when(x).toLocaleString("en-CA") + " · run " + x.run_at)).join("")}</select></label></div>
      <p class="ci-dim">${num(Q.quantified)} ${esc(Q.level)}s quantified (${esc(Q.rows)}). Ratios use mean log2 intensity
      per condition, at least 2 observed replicates per condition, no imputation.</p>
      ${ratios.length ? `<p class="ci-legend">${legend}</p>${boxes(ratios)}
      <p class="ci-dim">Boxes: middle 50%; center line: median. Whiskers: outermost observations within 1.5 × IQR. Hover a box for exact values.</p>
      <table class="ci-table"><thead><tr><th>Group</th><th>Comparison</th><th>n</th><th>Median log2</th><th>Middle 50%</th>
      <th>Expected</th><th>Bias</th><th>Median absolute error</th></tr></thead><tbody>${ratioRows}</tbody></table>`
      : "<p>No ratio estimates in this run: no precursor had at least 2 observed replicates in a condition and in the reference.</p>"}
      <h4>Replicate precision</h4><p class="ci-dim">Sample SD / mean on linear intensities; at least two observed replicates.</p>
      <table class="ci-table"><thead><tr><th>Group</th><th>Condition</th><th>n</th><th>Median CV</th></tr></thead><tbody>${cvRows}</tbody></table></div>`;
  }

  async function init() {
    const el = document.getElementById("andes-ci");
    if (!el || el.dataset.ready) return;
    el.dataset.ready = "1";
    const src = new URLSearchParams(location.search).get("src") || el.dataset.src;
    let rows;
    try {
      const res = await fetch(src, { cache: "no-cache" });
      if (res.status === 404) { el.innerHTML = "<p>No benchmark runs yet.</p>"; return; }
      if (!res.ok) throw new Error("HTTP " + res.status);
      const text = await res.text();
      rows = text.split("\n").filter((l) => l.trim()).map((l) => JSON.parse(l))
        .sort((a, b) => when(a) - when(b) || a.run_at.localeCompare(b.run_at));
    } catch (e) {
      el.innerHTML = `<p>Could not load <code>${esc(src)}</code>: ${esc(e.message)}</p>`;
      return;
    }
    if (!rows.length) { el.innerHTML = "<p>No benchmark runs yet.</p>"; return; }
    const sets = [...new Set(rows.map((r) => r.dataset))];
    const views = { Latest: () => latest(rows, sets), History: () => historyView(rows, sets),
      Table: () => table(rows, sets, el), Quantification: () => quant(rows, el) };
    let tab = Object.keys(views).find((t) => "#" + t.toLowerCase() === location.hash) || "Latest";
    const draw = () => {
      el.innerHTML = `<p class="ci-dim">${rows.length} runs · ${sets.length} datasets · ${new Set(rows.map((r) => r.commit)).size} commits</p>
        <div class="ci-tabs">${Object.keys(views).map((t) => `<button class="${t === tab ? "on" : ""}">${t}</button>`).join("")}</div>
        <div class="ci-view">${views[tab]()}</div>`;
    };
    el.addEventListener("click", (e) => { if (e.target.matches(".ci-tabs button")) { tab = e.target.textContent; history.replaceState(null, "", "#" + tab.toLowerCase()); draw(); } });
    el.addEventListener("change", (e) => {
      const act = e.target.dataset.act;
      if (act) { el.dataset[act] = e.target.value; if (act === "qset") delete el.dataset.qcommit; draw(); }
    });
    // Instant tooltip for chart points and quant boxes (anything with data-tip).
    const tip = document.querySelector(".ci-tip") || document.body.appendChild(document.createElement("div"));
    tip.className = "ci-tip";
    el.addEventListener("mousemove", (e) => {
      const t = e.target.closest("[data-tip]");
      tip.style.display = t ? "block" : "none";
      if (t) { tip.textContent = t.dataset.tip; tip.style.left = e.clientX + 12 + "px"; tip.style.top = e.clientY + 12 + "px"; }
    });
    el.addEventListener("mouseleave", () => { tip.style.display = "none"; });
    // The same details on keyboard focus (chart points are links, so Tab reaches them).
    el.addEventListener("focusin", (e) => {
      const t = e.target.closest("[data-tip]");
      if (!t) return;
      const b = t.getBoundingClientRect();
      tip.textContent = t.dataset.tip; tip.style.display = "block";
      tip.style.left = b.right + 8 + "px"; tip.style.top = b.bottom + 8 + "px";
    });
    el.addEventListener("focusout", () => { tip.style.display = "none"; });
    draw();
  }

  // Material's instant navigation swaps pages without reloading scripts.
  if (window.document$) window.document$.subscribe(init); else document.addEventListener("DOMContentLoaded", init);
})();
