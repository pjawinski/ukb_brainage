/* Brain age results browser — Jawinski et al. (2025) Nature Aging
   Static site: reads data/*.json produced by scripts/build_data.py
   Two result sets share one interface:
     pheno  phenome-wide associations (PHESANT) in all participants, women, men
     rg     genetic correlations with Neale lab UK Biobank GWAS (LDSC) */
(function () {
  "use strict";

  const MEASURES = { gwm: "Grey + white matter", gm: "Grey matter", wm: "White matter" };
  const MKEYS = ["gm", "wm", "gwm"];
  const SHORT = { gwm: "Grey + white", gm: "Grey", wm: "White" };
  const CAT_SHORT = {
    "Hospital Inpatient - Administration": "Hospital admin.", "Family history and early life factors": "Family & early life",
    "Maternity and sex-specific factors": "Maternity & sex-spec.", "Medical history and conditions": "Medical history",
    "Lifestyle and environment": "Lifestyle", "Diet by 24-hour recall": "Diet (24-h recall)",
  };
  const SAMPLES = { all: "All", female: "Women", male: "Men" };
  // one colour per category, in the alphabetical order of meta.categories
  const CAT_COLORS = [
    "#7d6b5d", "#c2417b", "#1c8fb0", "#7a8b22", "#b07a2a", "#3a6fd8", "#1f9e89", "#9b4fd9", "#6b7a8f", "#d4a017",
    "#e06c00", "#2e8b57", "#d63a5e", "#00838f", "#5c9ce6", "#a0761c", "#7b5ea7", "#e17ca4",
  ];
  const VIEWS = { pheno: ["manhattan", "volcano", "sex"], rg: ["manhattan", "volcano"], cmp: ["pg"], herit: [], loci: [], genes: [], rgsel: [], mr: [] };
  const OWN = { herit: () => drawHerit(), loci: () => drawLoci(), rgsel: () => drawRgSel(), mr: () => drawMr(), genes: () => drawGenes() }; // result sets with their own panels instead of plot + table
  const isOwn = (d) => d in OWN;
  const isGen = (d) => d === "rg" || d === "cmp"; // both use the genetic-correlation trait list
  const PAGE = 5;
  const $ = (s) => document.querySelector(s);
  const lt05 = (x) => x != null && x < 0.05;

  const state = { d: "pheno", v: "manhattan", m: "gwm", s: "all", t: null, q: "", cat: null, sig: "all", dir: "any", sort: "p", asc: true, page: 0, labels: 10 };
  const cache = {};
  let PMETA = null, RG = null, META = null;
  const P2RG = {}; // PheWAS trait index -> rg trait index
  let view = null;

  // ---------- helpers ----------
  const fmtInt = (n) => (n == null ? "–" : n.toLocaleString("en-US"));
  function fmtP(p, html = true) {
    if (p == null) return "–";
    if (p >= 0.001) return p.toPrecision(2);
    const [mant, ex] = p.toExponential(1).split("e");
    const e = String(+ex).replace("-", "−");
    return html ? `${mant}×10<sup>${e}</sup>` : `${mant}e${+ex}`;
  }
  const fmtR = (r) => (r == null ? "–" : (r < 0 ? "−" : "") + Math.abs(r).toFixed(3));
  const fmtB = (b, se) => (b == null ? "–" : `${b < 0 ? "−" : ""}${Math.abs(b).toPrecision(3)} (${se == null ? "–" : se.toPrecision(2)})`);
  const fmtF = (x, d = 3) => (x == null ? "–" : x.toFixed(d));
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
  function hexA(hex, a) {
    const n = parseInt(hex.slice(1), 16);
    return `rgba(${(n >> 16) & 255},${(n >> 8) & 255},${n & 255},${a})`;
  }
  const trunc = (s, n) => (s.length > n ? s.slice(0, n - 1) + "…" : s);
  const surfaceA = () => hexA(css("--surface").startsWith("#") ? css("--surface") : "#ffffff", 0.85);
  const reduced = () => matchMedia("(prefers-reduced-motion: reduce)").matches;

  function bh(ps) { // Benjamini–Hochberg on an array (nulls ignored)
    const idx = ps.map((p, i) => [p, i]).filter((x) => x[0] != null).sort((a, b) => a[0] - b[0]);
    const m = idx.length, q = new Array(ps.length).fill(null);
    let prev = 1;
    for (let k = m - 1; k >= 0; k--) {
      prev = Math.min(prev, (idx[k][0] * m) / (k + 1));
      q[idx[k][1]] = prev;
    }
    return q;
  }

  async function load(name) {
    // an offline copy can ship the data as a script (window.BAG_DATA), since browsers block fetch() from file://
    if (!cache[name] && window.BAG_DATA && window.BAG_DATA[name]) cache[name] = Promise.resolve(window.BAG_DATA[name]);
    if (!cache[name]) cache[name] = fetch(`data/${name}.json`).then((r) => {
      if (!r.ok) throw new Error(`Could not load data/${name}.json (${r.status})`);
      return r.json();
    });
    return cache[name];
  }

  // ---------- URL state ----------
  function readHash() {
    const h = new URLSearchParams(location.hash.slice(1));
    if (VIEWS[h.get("d")]) state.d = h.get("d");
    else if (!location.hash.slice(1)) state.d = "herit"; // the page opens on heritability; older links without d= stay on the phenotypic correlations
    if (VIEWS[state.d].includes(h.get("v"))) state.v = h.get("v");
    if (state.d === "cmp") state.v = "pg";
    if (MEASURES[h.get("m")]) state.m = h.get("m");
    if (state.d === "loci") LS.m = MEASURES[h.get("m")] ? h.get("m") : "all";
    if (state.d === "loci" && h.get("locus")) LS.sel = +h.get("locus"); // an open locus panel
    if (state.d === "mr" && h.get("trait")) mrSel = +h.get("trait"); // an open MR trait panel
    if (state.d === "genes" && h.get("gene")) GS.sel = h.get("gene");
    if (SAMPLES[h.get("s")] || h.get("s") === "sexdiff") state.s = h.get("s");
    if (h.get("t")) state.t = h.get("t");
    if (h.get("q")) state.q = h.get("q");
    if (h.get("sig")) state.sig = h.get("sig");
  }
  function writeHash() {
    if (isOwn(state.d)) {
      const o = new URLSearchParams({ d: state.d });
      if (state.d === "loci" && LS.m !== "all") o.set("m", LS.m);
      if (state.d === "loci" && LS.sel != null) o.set("locus", LS.sel);
      if (state.d === "mr" && mrSel != null) o.set("trait", mrSel);
      if (state.d === "genes") { o.set("m", state.m); if (GS.sel) o.set("gene", GS.sel); }
      return history.replaceState(null, "", "#" + o.toString());
    }
    const h = new URLSearchParams();
    h.set("d", state.d);
    if (state.v !== "manhattan" && state.d !== "cmp") h.set("v", state.v);
    h.set("m", state.m);
    if (state.d === "pheno") h.set("s", state.s);
    if (state.t) h.set("t", state.t);
    if (state.q) h.set("q", state.q);
    if (state.sig !== "all") h.set("sig", state.sig);
    history.replaceState(null, "", "#" + h.toString());
  }

  // ---------- build the current set of rows ----------
  async function buildView() {
    const m = state.m, rows = [];
    const counts = { gm: {}, wm: {}, gwm: {} };
    const tally = (mm, qs, catOf) => qs.forEach((q, i) => { if (lt05(q)) { const c = catOf(i); counts[mm][c] = (counts[mm][c] || 0) + 1; } });

    if (state.d === "cmp") {
      // trait pairs selected as in code/genetics/rgVSrp.R; "significant" here means FDR < 5% in both analyses
      const all = await load("all");
      const pairRows = (mm) => RG.pairs.ri.map((ri, k) => {
        const pi = RG.pairs.pi[k];
        if (RG[mm].p[ri] == null || all[mm].r[pi] == null) return null;
        return { i: ri, pidx: pi, p: RG[mm].p[ri], q: RG[mm].q[ri], r: RG[mm].r[ri], se: RG[mm].se[ri], n: RG.h2[ri],
          pr: all[mm].r[pi], pp: all[mm].p[pi], pq: all[mm].q[pi], b: all[mm].r[pi] };
      }).filter(Boolean);
      rows.push(...pairRows(m));
      MKEYS.forEach((mm) => pairRows(mm).forEach((r) => { if (lt05(r.q) && lt05(r.pq)) { const c = RG.cat[r.i]; counts[mm][c] = (counts[mm][c] || 0) + 1; } }));
      var pairs = rows;
    } else if (state.d === "rg") {
      const all = await load("all");
      const st = RG[m];
      for (let i = 0; i < RG.id.length; i++) {
        if (st.p[i] == null) continue;
        const pi = RG.phewas[i];
        rows.push({ i, p: st.p[i], q: st.q[i], r: st.r[i], se: st.se[i], b: st.se[i], n: RG.h2[i],
          pr: pi != null ? all[m].r[pi] : null, pq: pi != null ? all[m].q[pi] : null });
      }
      MKEYS.forEach((mm) => tally(mm, RG[mm].q, (i) => RG.cat[i]));
    } else if (state.s === "sexdiff") {
      const [sd, f, ma] = await Promise.all([load("sexdiff"), load("female"), load("male")]);
      const p = sd[m], q = bh(p);
      for (let i = 0; i < PMETA.id.length; i++) {
        if (p[i] == null) continue;
        const rf = f[m].r[i], rm = ma[m].r[i];
        rows.push({ i, p: p[i], q: q[i], r: rf != null && rm != null ? rf - rm : null, rf, rm, qf: f[m].q[i], qm: ma[m].q[i],
          b: null, se: null, n: (f.ntotal[i] || 0) + (ma.ntotal[i] || 0) });
      }
      MKEYS.forEach((mm) => tally(mm, bh(sd[mm]), (i) => PMETA.cat[i]));
    } else {
      const d = await load(state.s);
      const st = d[m];
      for (let i = 0; i < PMETA.id.length; i++) {
        if (st.p[i] == null) continue;
        rows.push({ i, p: st.p[i], q: st.q[i], r: st.r[i], b: st.b[i], se: st.se[i], n: d.ntotal[i] });
      }
      MKEYS.forEach((mm) => tally(mm, d[mm].q, (i) => PMETA.cat[i]));
    }
    // Manhattan x positions: categories present in this set, alphabetical, ppoints within category
    const byCat = {};
    rows.forEach((r) => (byCat[META.cat[r.i]] = byCat[META.cat[r.i]] || []).push(r));
    const cats = Object.keys(byCat).map(Number).sort((a, b) => a - b);
    cats.forEach((c, k) => byCat[c].forEach((r, j) => (r.x = k + (j + 0.5) / byCat[c].length)));
    const sigP = rows.filter((r) => lt05(r.q)).map((r) => r.p);
    view = { rows, byCat, cats, counts, pairs: typeof pairs !== "undefined" ? pairs : [],
      bonf: 0.05 / (state.d === "cmp" ? RG.id.length : rows.length), pbonf: 0.05 / PMETA.id.length, fdr: sigP.length ? Math.max(...sigP) : null };
  }

  function passes(r) {
    if (state.cat != null && META.cat[r.i] !== state.cat) return false;
    if (state.sig === "fdr" && !isSig(r)) return false;
    if (state.sig === "bonf" && !(r.p < view.bonf && (state.d !== "cmp" || r.pp < view.pbonf))) return false;
    if (state.dir === "pos" && !(r.r > 0)) return false;
    if (state.dir === "neg" && !(r.r < 0)) return false;
    if (state.q) {
      const t = state.q.toLowerCase(), f = META.field[r.i];
      if (!META.desc[r.i].toLowerCase().includes(t) && !META.id[r.i].toLowerCase().includes(t) && !(f != null && String(f).includes(t))) return false;
    }
    return true;
  }

  // ---------- main plot ----------
  let plotReady = false;
  const yl = (p) => -Math.log10(p);
  const isSig = (r) => lt05(r.q) && (state.d !== "cmp" || lt05(r.pq));
  function coords(r) {
    if (state.v === "volcano") return [r.r, yl(r.p)];
    if (state.v === "sex") return [r.rf, r.rm];
    if (state.v === "pg") return [r.r, r.pr];
    return [r.x, yl(r.p)];
  }
  function plotRows() {
    if (state.v === "sex") return view.rows.filter((r) => r.rf != null && r.rm != null);
    if (state.v === "pg") return view.pairs;
    if (state.v === "volcano") return view.rows.filter((r) => r.r != null);
    return view.rows;
  }
  // scatter views: colour = significant in at least one of the two compared analyses, outline = the highlighted contrast
  function scatterClass(r) {
    if (state.v === "sex") return { colour: lt05(r.qf) || lt05(r.qm), ring: isSig(r) };
    return { colour: lt05(r.q) || lt05(r.pq), ring: lt05(r.q) && lt05(r.pq) };
  }

  function plotLayout() {
    const narrow = window.innerWidth < 900;
    const ink = css("--ink"), muted = css("--muted"), rule = css("--rule"), band = css("--band");
    const v = state.v, rows = plotRows();
    const shapes = [], annotations = [];
    const sexdiff = state.d === "pheno" && state.s === "sexdiff";
    const zoomCat = v === "manhattan" && state.cat != null && view.cats.includes(state.cat);
    const grid = { gridcolor: rule, griddash: "dot", linecolor: rule, fixedrange: narrow };
    let xaxis, yaxis, x0, x1, y0, y1;

    if (v === "sex" || v === "pg") {
      const xs = rows.map((r) => coords(r)[0]), ys = rows.map((r) => coords(r)[1]);
      if (v === "sex") {
        const lim = Math.max(0.02, ...xs.map(Math.abs), ...ys.map(Math.abs)) * 1.12;
        [x0, x1, y0, y1] = [-lim, lim, -lim, lim];
        shapes.push({ type: "line", x0, x1, y0, y1, line: { color: muted, width: 1, dash: "dash" }, layer: "below" });
      } else {
        const lx = Math.max(0.1, ...xs.map(Math.abs)) * 1.08, ly = Math.max(0.02, ...ys.map(Math.abs)) * 1.12;
        [x0, x1, y0, y1] = [-lx, lx, -ly, ly];
        const ps = RG.pairStats[state.m]; // least-squares line, as geom_smooth(method = "lm") in the paper figure
        shapes.push({ type: "line", x0, x1, y0: ps.intercept + ps.slope * x0, y1: ps.intercept + ps.slope * x1, line: { color: ink, width: 1.4 } });
        annotations.push({ xref: "paper", yref: "paper", x: 0.5, y: 1, yanchor: "bottom", showarrow: false,
          text: `r = ${ps.r.toFixed(2)}  |  p = ${fmtP(ps.p)}  |  MAD = ${ps.mad.toFixed(3)}  |  ${ps.n} trait pairs`,
          font: { size: narrow ? 11.5 : 13, color: ink, family: "Source Sans 3, sans-serif" } });
      }
      xaxis = { ...grid, title: { text: v === "sex" ? "r in women" : "Genetic correlation r<sub>g</sub>", standoff: 6 }, range: [x0, x1], zeroline: true, zerolinecolor: muted };
      yaxis = { ...grid, title: { text: v === "sex" ? "r in men" : "Phenotypic correlation r", standoff: 6 }, range: [y0, y1], zeroline: true, zerolinecolor: muted };
    } else {
      y0 = 0;
      // Manhattan zooms into a selected category: its points spread over the full width and the y-axis fits them
      const yRows = zoomCat ? rows.filter((r) => META.cat[r.i] === state.cat) : rows;
      y1 = Math.max(zoomCat ? yl(view.bonf) * 1.15 : 8, ...yRows.map((r) => yl(r.p))) * (state.labels > 0 ? 1.2 : 1.06);
      shapes.push({ type: "line", xref: "paper", x0: 0, x1: 1, y0: yl(view.bonf), y1: yl(view.bonf), line: { color: ink, width: 1.1 } });
      if (view.fdr) shapes.push({ type: "line", xref: "paper", x0: 0, x1: 1, y0: yl(view.fdr), y1: yl(view.fdr), line: { color: ink, width: 1.1, dash: "dash" } });
      yaxis = { ...grid, title: { text: sexdiff ? "−log<sub>10</sub>(p) for sex difference" : "−log<sub>10</sub>(p)", standoff: 6 }, range: [y0, y1], zeroline: false };
      if (v === "volcano") {
        const lim = Math.max(0.02, ...rows.map((r) => Math.abs(r.r))) * 1.12;
        [x0, x1] = [-lim, lim];
        const xt = isGen(state.d) ? "Genetic correlation r<sub>g</sub>" : sexdiff ? "Δr (women − men)" : "r (correlation-scale effect size)";
        xaxis = { ...grid, title: { text: xt, standoff: 6 }, range: [x0, x1], zeroline: true, zerolinecolor: muted };
      } else {
        const cats = view.cats;
        if (zoomCat) {
          const k = cats.indexOf(state.cat);
          [x0, x1] = [k, k + 1];
          xaxis = { range: [x0, x1], fixedrange: true, showgrid: false, zeroline: false, tickvals: [k + 0.5],
            ticktext: [`${META.categories[state.cat]}: all ${fmtInt(view.byCat[state.cat].length)} traits`], tickangle: 0, ticks: "", linecolor: rule };
        } else {
          [x0, x1] = [0, cats.length];
          cats.forEach((c, k) => { if (k % 2 === 1) shapes.push({ type: "rect", xref: "x", yref: "paper", x0: k, x1: k + 1, y0: 0, y1: 1, fillcolor: band, line: { width: 0 }, layer: "below" }); });
          xaxis = { range: [x0, x1], fixedrange: true, showgrid: false, zeroline: false, tickvals: cats.map((_, k) => k + 0.5),
            ticktext: narrow ? cats.map(() => "") : cats.map((c) => CAT_SHORT[META.categories[c]] || META.categories[c]), tickangle: 40, ticks: "", linecolor: rule };
        }
      }
    }

    const margin = { l: 56, r: 12, t: v === "pg" ? 34 : 16, b: v === "manhattan" && !zoomCat ? (narrow ? 16 : 96) : v === "manhattan" ? 36 : 48 };
    if (state.labels > 0 || state.t != null) {
      // Label placement in pixel space: try several positions around each point and keep a label only where its
      // box stays inside the plot and does not overlap labels already placed.
      const el = $("#plot");
      const W = Math.max(100, el.clientWidth - margin.l - margin.r), H = Math.max(100, el.clientHeight - margin.t - margin.b);
      const fs = narrow ? 10.5 : 12, maxChars = narrow ? 26 : 40, lh = fs * 1.35 + 4;
      const px = (r) => ((coords(r)[0] - x0) / (x1 - x0)) * W, py = (r) => (1 - (coords(r)[1] - y0) / (y1 - y0)) * H;
      const sel = state.t != null ? rows.find((r) => META.id[r.i] === state.t) : null;
      const cand = rows.filter((r) => r !== sel && passes(r) && isSig(r)).sort((a, b) => a.p - b.p);
      const queue = (sel ? [sel] : []).concat(cand);
      const boxes = [], placed = [], seenTrait = new Set();
      const offsets = [[22, -20], [-22, -20], [22, 20], [-22, 20], [28, -42], [-28, -42], [28, 42], [-28, 42], [0, -46], [0, 46]];
      const overlap = (a, b) => a.x0 < b.x1 && b.x0 < a.x1 && a.y0 < b.y1 && b.y0 < a.y1;
      for (const r of queue) {
        if (placed.length >= state.labels + (sel ? 1 : 0)) break;
        if (seenTrait.has(r.i)) continue;
        const [x, y] = coords(r);
        if (x == null || y == null) continue;
        const text = trunc(META.desc[r.i], maxChars), w = text.length * fs * 0.53 + 8;
        const ax0 = px(r), ay0 = py(r);
        let hit = null;
        for (const [ox, oy] of offsets) {
          const tx = ax0 + ox, ty = ay0 + oy;
          const anchor = ox > 0 ? "left" : ox < 0 ? "right" : "center";
          const bx0 = anchor === "left" ? tx : anchor === "right" ? tx - w : tx - w / 2;
          const box = { x0: bx0, x1: bx0 + w, y0: ty - lh / 2, y1: ty + lh / 2 };
          if (box.x0 < 0 || box.x1 > W || box.y0 < 0 || box.y1 > H) continue;
          if (boxes.some((b) => overlap(b, box))) continue;
          hit = { ox, oy, anchor, box }; break;
        }
        if (!hit) continue;
        boxes.push(hit.box);
        boxes.push({ x0: ax0 - 5, x1: ax0 + 5, y0: ay0 - 5, y1: ay0 + 5 }); // keep later labels off this point
        seenTrait.add(r.i);
        placed.push(r);
        annotations.push({
          x, y, text: r === sel ? `<b>${esc(text)}</b>` : esc(text),
          showarrow: true, arrowhead: 0, arrowwidth: 0.8, arrowcolor: muted, ax: hit.ox, ay: hit.oy, xanchor: hit.anchor,
          font: { size: fs, color: ink, family: "Source Sans 3, sans-serif" }, bgcolor: surfaceA(),
        });
      }
    }
    return {
      margin,
      paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
      font: { family: "Source Sans 3, sans-serif", color: muted, size: 12.5 },
      xaxis, yaxis, shapes, annotations, showlegend: false, dragmode: narrow ? false : "zoom",
      hoverlabel: { bgcolor: css("--surface"), bordercolor: rule, font: { color: ink, family: "Source Sans 3, sans-serif", size: 13 }, align: "left" },
    };
  }

  function hoverText(r) {
    const head = `<b>${esc(trunc(META.desc[r.i], 70))}</b><br>${esc(META.categories[META.cat[r.i]])}` +
      (META.field[r.i] != null ? ` · field ${META.field[r.i]}` : "") + "<br>";
    if (isGen(state.d)) {
      const other = r.pidx != null ? ` (${esc(PMETA.id[r.pidx])})` : "";
      return head + `r<sub>g</sub> = ${fmtR(r.r)} (SE ${fmtF(r.se)}) · p = ${fmtP(r.p)} · FDR = ${fmtP(r.q)}` +
        (r.pr != null ? `<br>phenotypic r${other} = ${fmtR(r.pr)} · FDR = ${fmtP(r.pq)}` : "");
    }
    if (state.s === "sexdiff") return head + `r women = ${fmtR(r.rf)} · r men = ${fmtR(r.rm)}<br>p (difference) = ${fmtP(r.p)} · FDR = ${fmtP(r.q)}`;
    return head + (r.r != null ? `r = ${fmtR(r.r)} · ` : "") + `p = ${fmtP(r.p)}` + (r.q != null ? ` · FDR = ${fmtP(r.q)}` : "");
  }

  function plotTraces() {
    const traces = [];
    const v = state.v, rows = plotRows();
    const scatter = v === "sex" || v === "pg";
    const dirMarks = v === "manhattan" && !(state.d === "pheno" && state.s === "sexdiff");
    const grey = css("--faint"), ink = css("--ink");
    const filtering = state.cat != null || state.q || state.sig !== "all" || state.dir !== "any";
    const byCat = {};
    rows.forEach((r) => (byCat[META.cat[r.i]] = byCat[META.cat[r.i]] || []).push(r));
    // the selected category is drawn last so its points sit on top of the greyed-out rest
    Object.entries(byCat).sort((a, b) => (+a[0] === state.cat) - (+b[0] === state.cat)).forEach(([c, rs]) => {
      const col = CAT_COLORS[+c % CAT_COLORS.length];
      const x = [], y = [], sym = [], size = [], color = [], cd = [], text = [], lw = [];
      const rank = (r) => scatter ? (scatterClass(r).ring ? 2 : 0) + (scatterClass(r).colour ? 1 : 0) : isSig(r) ? 1 : 0;
      rs.slice().sort((a, b) => rank(a) - rank(b)).forEach((r) => { // significant points drawn last, on top
        const sig = isSig(r), on = !filtering || passes(r);
        const [px, py] = coords(r);
        x.push(px); y.push(py);
        sym.push(dirMarks && sig ? (r.r >= 0 ? "triangle-up" : "triangle-down") : "circle");
        if (scatter) {
          const k = scatterClass(r);
          size.push(k.ring ? 11 : k.colour ? 7 : filtering && on ? 5.5 : 4);
          color.push(!on ? hexA(grey, 0.12) : k.colour || k.ring ? col : filtering ? hexA(col, 0.55) : hexA(grey, 0.28));
          lw.push(k.ring && on ? 1.6 : 0);
        } else {
          // inside an active filter every point keeps its category colour; only points outside the filter turn grey
          size.push(sig ? (v === "manhattan" ? 9 : 8) : filtering && on ? 6 : (v === "manhattan" ? 5 : 4.5));
          color.push(on ? (sig ? col : hexA(col, filtering ? 0.6 : v === "manhattan" ? 0.38 : 0.3)) : hexA(grey, 0.16));
          lw.push(0);
        }
        cd.push(r.i);
        text.push(hoverText(r));
      });
      traces.push({
        type: "scattergl", mode: "markers", x, y, customdata: cd, text, hovertemplate: "%{text}<extra></extra>",
        marker: { symbol: sym, size, color, line: { width: lw, color: ink } },
      });
    });
    const sel = state.t != null ? rows.find((r) => META.id[r.i] === state.t) : null;
    const sc = sel ? coords(sel) : null;
    traces.push({
      type: "scatter", mode: "markers", x: sc ? [sc[0]] : [], y: sc ? [sc[1]] : [], hoverinfo: "skip",
      marker: { symbol: "circle-open", size: 18, color: ink, line: { width: 2 } },
    });
    return traces;
  }

  function drawPlot() {
    if (!view || isOwn(state.d)) return; // nothing to draw until a correlation view has been built
    const el = $("#plot");
    const cfg = { responsive: true, displaylogo: false, modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"],
      toImageButtonOptions: { filename: `brainage_${state.d}_${state.v}_${state.m}${state.d === "pheno" ? "_" + state.s : ""}`, scale: 3 } };
    Plotly.react(el, plotTraces(), plotLayout(), cfg).then(() => {
      $("#plot-loading").hidden = true;
      if (!plotReady) {
        plotReady = true;
        el.on("plotly_click", (ev) => {
          const pt = ev.points && ev.points[0];
          if (pt && pt.customdata != null) select(META.id[pt.customdata], true);
        });
      }
    });
  }

  // ---------- legend under the plot ----------
  function keyHtml() {
    const tri = (d) => `<svg width="12" height="12" viewBox="0 0 12 12"><path d="${d === "up" ? "M6 1 11 11H1z" : "M6 11 11 1H1z"}"/></svg>`;
    const line = (dash) => `<svg width="22" height="12" viewBox="0 0 22 12"><line x1="0" y1="6" x2="22" y2="6"${dash ? ' stroke-dasharray="4 3"' : ""}/></svg>`;
    const dot = (ring) => `<svg width="14" height="14" viewBox="0 0 14 14"><circle cx="7" cy="7" r="5" class="${ring ? "ring" : "fill"}"/></svg>`;
    const thresholds = `<span>${line(true)} FDR 5%</span><span>${line(false)} Bonferroni 5%</span>`;
    if (state.v === "sex") return `<span><svg width="16" height="16" viewBox="0 0 16 16"><line x1="1" y1="15" x2="15" y2="1" stroke-dasharray="3 2"/></svg> equal effect in women and men</span>` +
      `<span>${dot(false)} significant in women or men</span><span>${dot(true)} women and men differ (FDR 5%)</span>`;
    if (state.v === "pg") return `<span>${dot(false)} significant in one analysis (FDR 5%)</span><span>${dot(true)} significant in both</span>`;
    if (state.v === "volcano" || (state.d === "pheno" && state.s === "sexdiff")) return thresholds;
    const up = isGen(state.d) ? "positive genetic correlation" : "older-appearing brain with higher trait value";
    const down = isGen(state.d) ? "negative genetic correlation" : "older-appearing brain with lower trait value";
    return `<span>${tri("up")} ${up}</span><span>${tri("down")} ${down}</span>` + thresholds;
  }

  // ---------- category chips ----------
  function drawCats() {
    if (!view) return;
    const counts = view.counts[state.m];
    const box = $("#cats");
    box.classList.toggle("filtering", state.cat != null);
    box.innerHTML = view.cats.map((k) => {
      const c = META.categories[k];
      return `<button type="button" class="cat" style="--c:${CAT_COLORS[k]}" data-k="${k}" aria-pressed="${state.cat === k}" title="${counts[k] || 0} traits at FDR < 5%">` +
        `<i></i>${esc(c)}${counts[k] ? ` <b>${counts[k]}</b>` : ""}</button>`;
    }).join("");
    box.querySelectorAll(".cat").forEach((b) => b.addEventListener("click", () => setCat(+b.dataset.k)));
  }
  function setCat(k) {
    state.cat = state.cat === k ? null : k; state.page = 0;
    refresh(false);
  }

  // ---------- table ----------
  function sortedFiltered() {
    const rows = view.rows.filter(passes);
    const k = state.sort, dir = state.asc ? 1 : -1;
    const cmp = state.d === "cmp";
    const val = (r) => k === "desc" ? META.desc[r.i].toLowerCase() : k === "cat" ? META.categories[META.cat[r.i]] :
      k === "r" ? (r.r == null ? null : Math.abs(r.r)) : cmp && k === "b" ? Math.abs(r.pr) : cmp && k === "p" ? r.q : cmp && k === "q" ? r.pq : r[k];
    rows.sort((a, b) => {
      const va = val(a), vb = val(b);
      if (va == null) return 1; if (vb == null) return -1;
      return (va < vb ? -1 : va > vb ? 1 : 0) * dir || a.p - b.p;
    });
    return rows;
  }

  function drawTable() {
    if (!view) return;
    const rg = isGen(state.d), cmp = state.d === "cmp";
    const rows = sortedFiltered();
    $("#table").classList.toggle("cmp", cmp);
    const pages = Math.max(1, Math.ceil(rows.length / PAGE));
    state.page = Math.min(state.page, pages - 1);
    const slice = rows.slice(state.page * PAGE, state.page * PAGE + PAGE);
    const tb = $("#table tbody");
    tb.innerHTML = slice.length ? slice.map((r) => {
      const id = META.id[r.i], c = META.cat[r.i];
      let sub = rg ? (META.field[r.i] != null ? `field ${META.field[r.i]}` : "") : id;
      // a pair can join two different fields with the same name; then name the phenotypic one too
      if (r.pidx != null && PMETA.field[r.pidx] !== META.field[r.i]) sub += ` · phenotypic r: field ${PMETA.field[r.pidx]}${PMETA.qual[r.pidx] ? `, ${PMETA.qual[r.pidx]}` : ""}`;
      return `<tr data-id="${esc(id)}" class="${id === state.t ? "sel" : ""}" tabindex="0">` +
        `<td class="trait">${esc(META.desc[r.i])}<small>${esc(sub)}</small></td>` +
        `<td class="hide-s catcell"><span class="dot" style="background:${CAT_COLORS[c]}"></span>${esc(META.categories[c])}</td>` +
        `<td class="num hide-s">${rg ? fmtF(r.n) : fmtInt(r.n)}</td>` +
        `<td class="num${r.r < 0 ? " neg" : ""}">${fmtR(r.r)}</td>` +
        (cmp ? `<td class="num bcol${r.pr < 0 ? " neg" : ""}">${fmtR(r.pr)}</td>` +
          `<td class="num${lt05(r.q) ? " sig" : ""}">${fmtP(r.q)}</td><td class="num${lt05(r.pq) ? " sig" : ""}">${fmtP(r.pq)}</td></tr>`
        : `<td class="num hide-m">${rg ? fmtF(r.se) : fmtB(r.b, r.se)}</td>` +
          `<td class="num${isSig(r) ? " sig" : ""}">${fmtP(r.p)}</td>` +
          `<td class="num">${fmtP(r.q)}</td></tr>`);
    }).join("") : `<tr><td colspan="7">No traits match these filters. Clear the search or pick another category.</td></tr>`;
    tb.querySelectorAll("tr[data-id]").forEach((tr) => {
      const go = () => select(tr.dataset.id, false);
      tr.addEventListener("click", go);
      tr.addEventListener("keydown", (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); go(); } });
    });
    const nSig = rows.filter(isSig).length;
    $("#count").textContent = cmp
      ? `${fmtInt(rows.length)} of ${fmtInt(view.rows.length)} trait pairs shown, ${fmtInt(nSig)} significant in both analyses (FDR < 5%)`
      : `${fmtInt(rows.length)} of ${fmtInt(view.rows.length)} traits shown, ${fmtInt(nSig)} at FDR < 5%`;
    $("#pageinfo").textContent = `Page ${state.page + 1} of ${pages}`;
    $("#prev").disabled = state.page === 0;
    $("#next").disabled = state.page >= pages - 1;
    document.querySelectorAll("#table th").forEach((th) => {
      th.setAttribute("aria-sort", th.dataset.k === state.sort ? (state.asc ? "ascending" : "descending") : "none");
    });
    $("#table th[data-k='r']").innerHTML = rg ? "r<sub>g</sub>" : state.s === "sexdiff" ? "Δr" : "r";
    $("#table th[data-k='n']").innerHTML = rg ? "h²" : "N";
    $("#table th[data-k='b']").textContent = cmp ? "r" : rg ? "SE" : "β (SE)";
    $("#table th[data-k='p']").innerHTML = cmp ? "FDR r<sub>g</sub>" : "p";
    $("#table th[data-k='q']").innerHTML = cmp ? "FDR r" : "FDR";
    $("#table th[data-k='b']").classList.toggle("bcol", cmp);
    $("#table th[data-k='n']").title = rg ? "SNP heritability of the UK Biobank trait (LDSC, observed scale)" : "";
  }

  // ---------- side panel: category overview or one trait ----------
  function forest(el, series, xTitle, height) {
    const ink = css("--ink"), muted = css("--muted"), rule = css("--rule");
    const traces = series.map((s) => ({
      type: "scatter", mode: "markers", name: s.name, x: s.pts.map((p) => p.x), y: s.pts.map((p) => 2 - p.k + (s.off || 0)),
      text: s.pts.map((p) => p.text), hovertemplate: "%{text}<extra></extra>",
      error_x: { type: "data", symmetric: false, array: s.pts.map((p) => p.hi - p.x), arrayminus: s.pts.map((p) => p.x - p.lo), color: s.color, thickness: 1.4, width: 0 },
      marker: { color: s.color, size: 7, symbol: s.symbol || "circle" },
    }));
    // the x-axis always includes 0, so one can see whether a confidence interval excludes it
    const vals = series.flatMap((s) => s.pts.flatMap((p) => [p.lo, p.hi, p.x])).filter((v) => v != null && isFinite(v));
    let x0 = Math.min(0, ...vals), x1 = Math.max(0, ...vals);
    const pad = (x1 - x0 || 0.1) * 0.08; x0 -= pad; x1 += pad;
    Plotly.react(el, traces, {
      margin: { l: 82, r: 10, t: 6, b: series.length > 1 ? 54 : 36 }, height,
      paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
      font: { family: "Source Sans 3, sans-serif", color: muted, size: 12 },
      xaxis: { range: [x0, x1], zeroline: true, zerolinecolor: ink, zerolinewidth: 1, gridcolor: rule, griddash: "dot", fixedrange: true, title: { text: xTitle, standoff: 4 } },
      yaxis: { tickvals: [2, 1, 0], ticktext: MKEYS.map((m) => SHORT[m]), range: [-0.5, 2.5], fixedrange: true, showgrid: false, zeroline: false },
      showlegend: series.length > 1, legend: { orientation: "h", x: 0, y: -0.32, font: { size: 12.5, color: ink } },
      hoverlabel: { bgcolor: css("--surface"), bordercolor: rule, font: { color: ink, size: 12.5 } },
    }, { displayModeBar: false, responsive: true });
  }

  function drawOverview() {
    const box = $("#detail");
    box.classList.add("is-overview");
    const c = view.counts;
    const cats = view.cats.filter((k) => MKEYS.some((m) => c[m][k]));
    cats.sort((a, b) => (c[state.m][b] || 0) - (c[state.m][a] || 0) || MKEYS.reduce((s, m) => s + (c[m][b] || 0) - (c[m][a] || 0), 0));
    const what = state.d === "rg" ? "genetic correlations at FDR &lt; 5%" : state.d === "cmp" ? "trait pairs significant in both analyses (FDR &lt; 5%)"
      : state.s === "sexdiff" ? "sex differences at FDR &lt; 5%" : "associations at FDR &lt; 5%";
    box.innerHTML = `<h2>FDR hits by category</h2>` +
      `<p class="lede">Number of ${what} for each brain age model. Select a bar to filter the plot and table; select a trait for its details.</p>` +
      (cats.length ? `<div class="bar-key" aria-hidden="true">${MKEYS.map((m) => `<span class="${m === state.m ? "on" : ""}"><i style="background:var(--m-${m})"></i>${SHORT[m]}</span>`).join("")}</div><div id="catbars" class="catbars"></div>` : `<p class="lede">No category has results at FDR &lt; 5% in this view.</p>`) +
      (state.cat != null ? `<button type="button" class="btn" id="clear-cat">Show all categories</button>` : "");
    if (state.cat != null) $("#clear-cat").addEventListener("click", () => setCat(state.cat));
    if (!cats.length) return;
    const ink = css("--ink"), muted = css("--muted"), rule = css("--rule");
    const names = cats.map((k) => CAT_SHORT[META.categories[k]] || META.categories[k]);
    const traces = MKEYS.map((m) => ({
      type: "bar", orientation: "h", name: SHORT[m], y: names, x: cats.map((k) => c[m][k] || 0), customdata: cats,
      marker: { color: css(`--m-${m}`), opacity: cats.map((k) => (state.cat == null || state.cat === k ? 1 : 0.3)),
        line: { width: m === state.m ? 1.5 : 0, color: ink } },
      hovertemplate: `%{y}<br>${MEASURES[m]}: %{x}${state.d === "cmp" ? " pairs significant in both" : " at FDR < 5%"}<extra></extra>`,
    }));
    const el = $("#catbars");
    Plotly.react(el, traces, {
      height: 34 + cats.length * 34, barmode: "group", bargap: 0.28, bargroupgap: 0.08,
      margin: { l: 128, r: 10, t: 4, b: 30 },
      paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
      font: { family: "Source Sans 3, sans-serif", color: muted, size: 12 },
      xaxis: { gridcolor: rule, griddash: "dot", fixedrange: true, zeroline: false, rangemode: "tozero" },
      yaxis: { autorange: "reversed", fixedrange: true, ticks: "", tickfont: { color: ink } },
      showlegend: false,
      hoverlabel: { bgcolor: css("--surface"), bordercolor: rule, font: { color: ink, size: 12.5 } },
    }, { displayModeBar: false, responsive: true }).then(() => {
      el.on("plotly_click", (ev) => { const pt = ev.points && ev.points[0]; if (pt) setCat(pt.customdata); });
    });
  }


  async function drawDetail() {
    if (!view || isOwn(state.d)) return;
    // the selected trait opens in its own panel under the table; the category overview stays in the sidebar
    const box = $("#trait");
    const i = state.t != null ? META.id.indexOf(state.t) : -1;
    box.hidden = i < 0;
    if (i < 0) { box.innerHTML = ""; return; }
    const ink = css("--ink");
    const field = META.field[i];
    const header = (extra) =>
      `<button type="button" class="btn close" id="close-detail">Close</button>` + `<button type="button" class="btn close share" data-share>Copy link</button>` +
      `<h2>${esc(META.desc[i])}</h2><div class="meta">` +
      `<div><span class="dot" style="background:${CAT_COLORS[META.cat[i]]}"></span>${esc(META.categories[META.cat[i]])}</div>` + extra +
      `<div class="path">${esc(META.paths[META.path[i]])}</div></div>`;
    const showcase = (f, code) => f != null ? `<div><a href="https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=${f}">UK Biobank field ${f}</a>${code ? `, coding ${esc(code)}` : ""}</div>` : "";
    const all = await load("all");

    if (isGen(state.d)) {
      const rgRows = MKEYS.map((m) => `<tr${m === state.m ? ' class="sel"' : ""}><td>${SHORT[m]}</td><td class="num">${fmtR(RG[m].r[i])}</td>` +
        `<td class="num">${fmtF(RG[m].se[i])}</td><td class="num${lt05(RG[m].q[i]) ? " sig" : ""}">${fmtP(RG[m].p[i])}</td><td class="num">${fmtP(RG[m].q[i])}</td></tr>`).join("");
      box.innerHTML = header(showcase(field, null) +
        `<div>SNP heritability h² = ${fmtF(RG.h2[i])} (SE ${fmtF(RG.h2se[i])})</div>`) +
        `<div class="trait-body"><div><h3>Genetic correlation r<sub>g</sub> with 95% CI</h3><div class="forest" id="forest"></div></div><div>` +
        `<div class="mini-wrap"><table class="mini"><thead><tr><th>Model</th><th class="num">r<sub>g</sub></th><th class="num">SE</th><th class="num">p</th><th class="num">FDR</th></tr></thead><tbody>${rgRows}</tbody></table></div>` +
        `</div></div>`;
      forest("forest", [{ name: "rg", color: ink, symbol: "diamond", pts: MKEYS.map((m, k) => {
        const r = RG[m].r[i], se = RG[m].se[i];
        return r == null ? null : { k, x: r, lo: r - 1.96 * se, hi: r + 1.96 * se, text: `${MEASURES[m]}<br>r<sub>g</sub> = ${fmtR(r)} [${fmtR(r - 1.96 * se)}, ${fmtR(r + 1.96 * se)}]<br>p = ${fmtP(RG[m].p[i])}` };
      }).filter(Boolean) }], "r<sub>g</sub>", 190);
    } else {
      const [f, ma, sd] = await Promise.all(["female", "male", "sexdiff"].map(load));
      const data = { all, female: f, male: ma };
      const id = META.id[i];
      const code = id.includes("#") ? id.split("#")[1] : id.includes("-") ? id.split("-").slice(1).join("-") : null;
      const nstr = (all.n[i] || f.n[i] || ma.n[i] || "");
      const nm = /^(\d+)\/(\d+)\((\d+)\)$/.exec(nstr);
      const nText = nm ? `N = ${fmtInt(+nm[3])} (${fmtInt(+nm[1])} / ${fmtInt(+nm[2])} by outcome)` : all.ntotal[i] != null ? `N = ${fmtInt(all.ntotal[i])}` : "";
      const rowsHtml = MKEYS.map((m) => {
        const s = all[m];
        if (s.p[i] == null) return `<tr><td>${SHORT[m]}</td><td class="num" colspan="4">not tested in full sample</td></tr>`;
        return `<tr${m === state.m ? ' class="sel"' : ""}><td>${SHORT[m]}</td><td class="num">${fmtR(s.r[i])}</td>` +
          `<td class="num hide-s">${fmtB(s.b[i], s.se[i])}</td><td class="num${lt05(s.q[i]) ? " sig" : ""}">${fmtP(s.p[i])}</td><td class="num">${fmtP(s.q[i])}</td></tr>`;
      }).join("");
      const sexRows = MKEYS.map((m) => sd[m][i] == null ? "" :
        `<tr><td>${SHORT[m]}</td><td class="num">${fmtR(f[m].r[i])}</td><td class="num">${fmtR(ma[m].r[i])}</td><td class="num">${fmtP(sd[m][i])}</td></tr>`).join("");

      box.innerHTML = header(showcase(field, code) +
        `<div>${esc(META.resTypes[META.rt[i]].toLowerCase().replace("-", " "))} regression${nText ? `, ${nText}` : ""}</div>`) +
        `<div class="trait-body"><div><h3>Effect size r with 95% CI</h3><div class="forest" id="forest"></div></div><div>` +
        `<h3>All participants</h3><div class="mini-wrap"><table class="mini"><thead><tr><th>Model</th><th class="num">r</th><th class="num hide-s">β (SE)</th><th class="num">p</th><th class="num">FDR</th></tr></thead><tbody>${rowsHtml}</tbody></table></div>` +
        (sexRows ? `<h3>Women versus men</h3><div class="mini-wrap"><table class="mini"><thead><tr><th>Model</th><th class="num">r women</th><th class="num">r men</th><th class="num">p diff.</th></tr></thead><tbody>${sexRows}</tbody></table></div>` : "") +
        `<p class="note">Model rows: brain age gap from grey matter, white matter, or both combined.</p></div></div>`;

      const sColors = { all: ink, female: "#b03a7a", male: "#2f74c0" };
      const offs = { all: 0.22, female: 0, male: -0.22 };
      forest("forest", Object.keys(SAMPLES).map((s) => ({
        name: SAMPLES[s], color: sColors[s], symbol: s === "all" ? "diamond" : "circle", off: offs[s],
        pts: MKEYS.map((m, k) => {
          const d = data[s], r = d[m].r[i], n = d.ntotal[i];
          if (r == null || !n) return null;
          const zf = Math.atanh(r), se = 1 / Math.sqrt(n - 3);
          const lo = Math.tanh(zf - 1.96 * se), hi = Math.tanh(zf + 1.96 * se);
          return { k, x: r, lo, hi, text: `${SAMPLES[s]}, ${MEASURES[m].toLowerCase()}<br>r = ${fmtR(r)} [${fmtR(lo)}, ${fmtR(hi)}]<br>p = ${fmtP(d[m].p[i])}, N = ${fmtInt(n)}` };
        }).filter(Boolean),
      })), "r", 250);
    }
    $("#close-detail").addEventListener("click", () => select(null, false));
    box.querySelectorAll(".jump").forEach((b) => b.addEventListener("click", () => setMode(b.dataset.d)));
  }


  // ---------- heritability panels (values as reported in the result tables) ----------
  function modeKey(hollow) {
    return MKEYS.map((m) => `<span><i style="background:var(--m-${m})"></i>${SHORT[m]}</span>`).join("") +
      (hollow ? `<span><i class="hollow"></i>not significant (FDR ≥ 5%)</span>` : "");
  }
  function hLayout(extra) {
    const ink = css("--ink"), muted = css("--muted"), rule = css("--rule");
    return Object.assign({
      paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)", showlegend: false,
      font: { family: "Source Sans 3, sans-serif", color: muted, size: 12.5 },
      hoverlabel: { bgcolor: css("--surface"), bordercolor: rule, font: { color: ink, size: 12.5 }, align: "left" },
      modebar: { bgcolor: "rgba(0,0,0,0)", color: muted, activecolor: ink },
    }, extra);
  }
  const H_CFG = { displaylogo: false, responsive: true, modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d", "zoom2d", "pan2d", "zoomIn2d", "zoomOut2d"] };

  async function drawHerit() {
    const H = await load("herit");
    const ink = css("--ink"), muted = css("--muted"), rule = css("--rule");
    const narrow = window.innerWidth < 640;
    document.querySelectorAll("#herit .model-key").forEach((el, k) => (el.innerHTML = modeKey(k > 0)));
    const off = { gm: 0.24, wm: 0, gwm: -0.24 };

    // SNP heritability: one row per sample, three models per row
    const n = H.ldsc.length, labels = H.ldsc.map((r) => r.label);
    const h2traces = MKEYS.map((m) => ({
      type: "scatter", mode: "markers", name: SHORT[m],
      x: H.ldsc.map((r) => r[m].h2), y: H.ldsc.map((_, k) => n - 1 - k + off[m]),
      error_x: { type: "data", array: H.ldsc.map((r) => 1.96 * r[m].se), color: css(`--m-${m}`), thickness: 1.6, width: 0 },
      marker: { color: css(`--m-${m}`), size: 9, line: { width: 1, color: ink } },
      text: H.ldsc.map((r) => `<b>${r.label}</b>${r.n ? `, n = ${fmtInt(r.n)}` : ""}<br>${MEASURES[m]}<br>h² = ${r[m].h2.toFixed(3)} (SE ${r[m].se.toFixed(3)})<br>95% CI ${(r[m].h2 - 1.96 * r[m].se).toFixed(3)} to ${(r[m].h2 + 1.96 * r[m].se).toFixed(3)}<br>LDSC intercept = ${r[m].intercept.toFixed(3)} (SE ${r[m].intercept_se.toFixed(3)})` +
        (r.key === "female" || r.key === "male" ? "<br><i>UK Biobank only</i>" : "")),
      hovertemplate: "%{text}<extra></extra>",
    }));
    Plotly.react("h2plot", h2traces, hLayout({
      height: 64 + n * 52, margin: { l: 118, r: 12, t: 8, b: 40 },
      xaxis: { title: { text: "SNP heritability h²", standoff: 4 }, range: [0, 0.42], gridcolor: rule, griddash: "dot", zeroline: false, fixedrange: true },
      yaxis: { tickvals: labels.map((_, k) => n - 1 - k),
        ticktext: H.ldsc.map((r) => r.n ? `${r.label}<br><span style="font-size:11px;color:${muted}">n = ${fmtInt(r.n)}</span>` : r.label),
        range: [-0.6, n - 0.4], fixedrange: true, showgrid: false, zeroline: false, tickfont: { color: ink } },
      shapes: labels.slice(1).map((_, k) => ({ type: "line", xref: "paper", x0: 0, x1: 1, y0: n - 1.5 - k, y1: n - 1.5 - k, line: { color: rule, width: 1 } })),
    }), H_CFG);

    // GENESIS: number of SNPs with non-zero effect; brain age models in model colours, reference traits grey
    const G = H.genesis, gn = G.length;
    const gname = (r) => r.model ? MEASURES[r.model] : r.label;
    const gshort = (r) => r.model ? MEASURES[r.model] : r.label.replace(/ \(.*\)/, "") + " (ref.)";
    const gcol = (r) => r.model ? css(`--m-${r.model}`) : hexA(css("--faint").startsWith("#") ? css("--faint") : "#8a93a6", 0.55);
    Plotly.react("genplot", [{
      type: "bar", orientation: "h", x: G.map((r) => r.causal), y: G.map((_, k) => gn - 1 - k),
      error_x: { type: "data", array: G.map((r) => r.causal_se), color: ink, thickness: 1.2, width: 4 },
      marker: { color: G.map(gcol), line: { width: 1, color: ink } },
      text: G.map((r) => `<b>${esc(gname(r))}</b><br>SNPs with non-zero effect: ${fmtInt(r.causal)} (SE ${fmtInt(r.causal_se)})<br>in the large-effect component: ${fmtInt(r.causal_large)} (SE ${fmtInt(Math.round(r.causal_large_se))})` +
        `<br>N to explain 80% of h²: ${fmtN(r.reqsample)}<br>expected loci at that N: ${fmtInt(Math.round(r.reqloci))}`),
      hovertemplate: "%{text}<extra></extra>", textposition: "none",
    }], hLayout({
      height: 64 + gn * 40, margin: { l: 140, r: 16, t: 8, b: 40 }, bargap: 0.62,
      xaxis: { title: { text: "SNPs with non-zero effect", standoff: 4 }, gridcolor: rule, griddash: "dot", zeroline: false, fixedrange: true, rangemode: "tozero" },
      yaxis: { tickvals: G.map((_, k) => gn - 1 - k), ticktext: G.map(gshort),
        fixedrange: true, showgrid: false, zeroline: false, tickfont: { color: ink } },
      shapes: [{ type: "line", xref: "paper", x0: 0, x1: 1, y0: 1.5, y1: 1.5, line: { color: rule, width: 1, dash: "dot" } }],
    }), H_CFG);

    // partitioned heritability and cell-type groups: enrichment per annotation, filled = FDR < 5%
    const enrichPlot = (el, P) => {
      const k = P.annotation.length;
      const order = P.annotation.map((_, i) => i).sort((a, b) => P.gwm.enr[b] - P.gwm.enr[a]);
      const yOf = (i) => k - 1 - order.indexOf(i);
      const traces = MKEYS.map((m) => {
        const col = css(`--m-${m}`);
        return {
          type: "scatter", mode: "markers", name: SHORT[m],
          x: order.map((i) => P[m].enr[i]), y: order.map((i) => yOf(i) + off[m] * 0.75),
          marker: { size: 8, color: order.map((i) => lt05(P[m].q[i]) ? col : "rgba(0,0,0,0)"), line: { width: 1.6, color: order.map((i) => lt05(P[m].q[i]) ? ink : col) } },
          text: order.map((i) => `<b>${esc(P.annotation[i])}</b><br>${MEASURES[m]}<br>enrichment = ${fmtF(P[m].enr[i], 2)}<br>share of SNPs = ${(100 * P.propSnps[i]).toFixed(1)}%, share of h² = ${(100 * P[m].propH2[i]).toFixed(1)}%<br>p (one-sided) = ${fmtP(P[m].p[i])} · FDR = ${fmtP(P[m].q[i])}`),
          hovertemplate: "%{text}<extra></extra>",
        };
      });
      const xs = MKEYS.flatMap((m) => P[m].enr);
      Plotly.react(el, traces, hLayout({
        height: 66 + k * (narrow ? 30 : 26), margin: { l: narrow ? 165 : 230, r: 12, t: 22, b: 40 },
        xaxis: { title: { text: "Heritability enrichment", standoff: 4 }, range: [Math.min(0, ...xs) - 0.5, Math.max(...xs) * 1.06], gridcolor: rule, griddash: "dot", zeroline: false, fixedrange: true },
        yaxis: { tickvals: order.map(yOf), ticktext: order.map((i) => narrow ? trunc(P.annotation[i], 24) : P.annotation[i]), tickfont: { color: ink, size: narrow ? 11 : 12.5 }, range: [-0.7, k - 0.3], fixedrange: true, showgrid: false, zeroline: false },
        shapes: order.filter((_, j) => j % 2 === 0).map((i) => ({ type: "rect", xref: "paper", x0: 0, x1: 1, y0: yOf(i) - 0.5, y1: yOf(i) + 0.5,
            fillcolor: css("--band"), line: { width: 0 }, layer: "below" }))
          .concat([{ type: "line", yref: "paper", x0: 1, x1: 1, y0: 0, y1: 1, line: { color: muted, width: 1, dash: "dash" } }]),
        annotations: [{ x: 1, y: 1, xref: "x", yref: "paper", yanchor: "bottom", text: "no enrichment", showarrow: false, font: { size: 11, color: muted } }],
      }), H_CFG);
    };
    enrichPlot("baseplot", H.baseline);
    enrichPlot("ctgplot", H.celltype);
  }
  const fmtN = (x) => x >= 1e6 ? `${(x / 1e6).toFixed(2).replace(/\.?0+$/, "")} million` : `${fmtInt(Math.round(x / 1000))},000`;


  // ---------- GWAS loci (values and evidence strings as reported in the result tables) ----------
  // chromosome lengths in GRCh37, used only to sort loci by genome position
  const CHR_LEN = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747,
    135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560];
  const CHR_START = CHR_LEN.reduce((a, l, k) => (a.push(k ? a[k - 1] + CHR_LEN[k - 1] : 0), a), []);
  const chrIndex = (c) => (c === "X" || c === "XY" ? 22 : +c - 1); // XY = pseudoautosomal region, shown with X
  const gpos = (chr, bp) => CHR_START[chrIndex(chr)] + bp;
  const LS = { q: "", m: "all", nov: "all", sel: null, sort: "p", asc: true, more: false }; // loci open sorted by p, strongest first
  const LOCI_SHOW = 10; // rows shown before "Show all" // m: a brain age model, or "all" for the three together
  let LOCI = null, lociBound = false;

  const lociHit = (l) => (LS.m === "all" ? l.hits.reduce((a, b) => (b.p < a.p ? b : a)) : l.hits.find((h) => h.model === LS.m)); // the selected model's lead variant, or the strongest one
  function lociPasses(l) {
    if (LS.m !== "all" && !l.models.includes(LS.m)) return false;
    if (LS.nov === "novel" && !l.novel) return false;
    if (LS.nov === "known" && l.novel) return false;
    if (LS.q) {
      const t = LS.q.toLowerCase();
      const hay = [l.gene, l.cytoband, "chr" + l.chr, ...l.hits.flatMap((h) => [h.id, h.nearest || "", h.prioritized || ""])].join(" ").toLowerCase();
      if (!hay.includes(t)) return false;
    }
    return true;
  }
  const mchip = (m, best) => `<span class="mchip ${m}${best ? " best" : ""}" style="--mc:var(--m-${m})">${{ gm: "GM", wm: "WM", gwm: "GWM" }[m]}</span>`;

  async function drawLoci() {
    if (!LOCI) LOCI = (await load("loci")).loci;
    if (!lociBound) bindLoci();
    const all = LOCI, nNovel = all.filter((l) => l.novel).length;
    const all3 = LS.m === "all", label = all3 ? "" : MEASURES[LS.m].toLowerCase() + " ";
    const mine = all3 ? all : all.filter((l) => l.models.includes(LS.m)), mNovel = mine.filter((l) => l.novel).length;
    $("#loci-title").textContent = `${mine.length} genome-wide significant loci for ${label}brain age gap`;
    $("#loci-tally").textContent = all3
      ? `Combined European meta-analysis, n = 54,890. ${nNovel} novel and ${all.length - nNovel} previously reported loci across the grey matter, white matter and combined brain age models.`
      : `Combined European meta-analysis, n = 54,890. ${mNovel} of these loci are novel. Across the three brain age models, ${all.length} distinct loci, ${nNovel} of them novel.`;
    $("#loci-note").innerHTML = (all3
      ? "Manhattan plot showing, for each variant, the smallest p value across the three brain age models."
      : "Manhattan plot of the combined European meta-analysis for the selected brain age model.") +
      " Diamonds mark the lead variants of independent loci; the line marks p = 5×10<sup>−8</sup>, and the y-axis is truncated at −log<sub>10</sub>(p) = 40. " +
      (all3 ? "The table lists all loci; select a row for details." : "The table lists all loci found for this model; select a row for details.");
    const img = $("#loci-img");
    img.src = `img/manhattan.${LS.m}.png`;
    img.alt = `Manhattan plot for ${all3 ? "the three brain age models together" : label + "brain age gap"}`;
    if (LS.sel != null && !mine.some((l) => l.locus === LS.sel)) LS.sel = null;
    drawLociTable(); drawLocus();
  }

  function drawLociTable() {
    const rows = LOCI.filter(lociPasses).map((l) => ({ l, h: lociHit(l) }));
    const dir = LS.asc ? 1 : -1;
    const key = { pos: (r) => gpos(r.l.chr, r.l.bp), gene: (r) => r.l.gene.toLowerCase(), id: (r) => r.h.id, models: (r) => r.l.models.length, beta: (r) => Math.abs(r.h.beta), p: (r) => r.h.p }[LS.sort];
    rows.sort((a, b) => (key(a) < key(b) ? -1 : key(a) > key(b) ? 1 : 0) * dir);
    const tb = $("#loci-table tbody");
    const shown = LS.more ? rows : rows.slice(0, LOCI_SHOW);
    tb.innerHTML = shown.map(({ l, h }) =>
      `<tr data-locus="${l.locus}" class="${l.locus === LS.sel ? "sel" : ""}" tabindex="0">` +
      `<td class="trait">${esc(l.cytoband)}<small>chr${esc(l.chr)}:${fmtInt(h.bp)}</small></td>` +
      `<td class="gene"><b><i>${esc(l.gene)}</i></b>${l.novel ? '<span class="novel-tag">novel</span>' : ""}</td>` +
      `<td class="hide-s vid">${esc(h.id)}<small>${esc(h.a1)}/${esc(h.a2)}, freq. ${fmtF(h.freq, 2)}</small></td>` +
      `<td class="models hide-s">${MKEYS.filter((m) => l.models.includes(m)).map((m) => mchip(m, LS.m === "all" && m === h.model && l.models.length > 1)).join("")}</td>` +
      `<td class="num hide-m">${fmtB(h.beta, h.se)}</td>` +
      `<td class="num sig">${fmtP(h.p)}</td></tr>`).join("") ||
      `<tr><td colspan="6">No loci match these filters. Clear the search or choose another option.</td></tr>`;
    tb.querySelectorAll("tr[data-locus]").forEach((tr) => {
      const go = () => selectLocus(+tr.dataset.locus, true);
      tr.addEventListener("click", go);
      tr.addEventListener("keydown", (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); go(); } });
    });
    const nn = rows.filter((r) => r.l.novel).length;
    const total = LOCI.filter((l) => LS.m === "all" || l.models.includes(LS.m)).length;
    const btn = $("#loci-more");
    btn.hidden = rows.length <= LOCI_SHOW;
    btn.textContent = LS.more ? `Show first ${LOCI_SHOW} only` : `Show all ${rows.length} loci`;
    $("#loci-count").textContent = `${rows.length === total ? `${total} loci` : `${rows.length} of ${total} loci match`}, ${nn} novel. ` + (LS.m === "all"
      ? "Statistics refer to the model with the smallest p value per locus (outlined in “Models”)."
      : `Statistics refer to ${MEASURES[LS.m].toLowerCase()} brain age gap; “Models” lists every model for which the locus was found.`);
    document.querySelectorAll("#loci-table th").forEach((th) => th.setAttribute("aria-sort", th.dataset.k === LS.sort ? (LS.asc ? "ascending" : "descending") : "none"));
  }

  // column definitions as given in the paper's supplementary tables
  const LDEF = {
    "Locus": "Independent discovery count, each containing up to three co-inherited index variants derived from the three genome-wide association analyses of brain age gap.",
    "Cytoband": "Cytogenetic band that contains the index variant.",
    "Position": "Position of the index variant in base pairs according to human genome build hg19 (GRCh37).",
    "Variant": "Identifier of the index variant.",
    "A1/A2": "A1 is the allele for which effects were calculated; A2 is the other allele.",
    "Freq.": "Frequency of the effect allele (A1).",
    "β (SE)": "Beta weight of the association between index variant and phenotype, with its standard error (years of brain age gap per A1 allele).",
    "p": "P value of the association between index variant and phenotype.",
    "ηp²": "Partial eta squared; proportion of variance of the phenotype explained by the index variant (adjusted for sex, age, age², total intracranial volume, scanner site, type of array, and the first twenty genetic principal components).",
    "N": "Sample size.",
    "Nearest gene": "HGNC symbol of the nearest gene based on ANNOVAR annotations (hg19 RefSeq gene table updated September 29, 2019), with the most relevant functional category of the index variant and the distance to the gene in base pairs. ANNOVAR prioritizes the most deleterious annotation for variants located in genomic regions where multiple genes overlap. Description and biotype follow the RefSeq gene annotation file in GFF3 format updated November 5, 2019.",
    "Credible set size": "Size of the 95% credible set of variants derived from applying SBayesRC, susieR and FINEMAP.",
    "SBayesRC genes": "Genes nominated by SBayesRC credible variant analysis. Brackets include the cumulative posterior probability of variants that have been annotated with the corresponding gene.",
    "Nonsynonymous variants": "Nonsynonymous exonic variants from the 95% credible set. Cells show genes whose transcripts are affected by the respective exonic variants. Brackets contain the number of identified nonsynonymous variants, the top nonsynonymous variant, and its CADD deleteriousness score.",
    "SMR eQTL": "Genes whose expression levels putatively mediate the effect of a locus variant on brain age gap, identified using summary-data-based Mendelian randomization (SMR) and the BrainMeta v2 eQTL dataset. Brackets contain the locus variant and the SMR raw p value.",
    "SMR sQTL": "Genes whose RNA splicing putatively mediates the effect of a locus variant on brain age gap, identified using SMR and the BrainMeta v2 sQTL dataset. Brackets contain the locus variant and the SMR raw p value.",
    "GTEx single tissue": "Regulated genes identified by mapping GWAS results to single-tissue expression quantitative trait loci of the Genotype-Tissue Expression (GTEx) database. Brackets include the number of tissues with a significant eQTL and the minimum p value across tissues.",
    "GTEx multi-tissue": "Regulated genes identified by mapping GWAS results to multi-tissue expression quantitative trait loci of the GTEx database. Brackets include the number of tissues where the eQTL had a posterior probability ≥ 0.9, and the Han and Eskin RE2 p value.",
    "PoPS": "Genes implicated by the Polygenic Priority Score (PoPS) analysis. Brackets include the polygenic priority score.",
    "Prioritized gene": "Gene prioritized by aggregating the results of the seven gene nomination strategies.",
    "GWAS Catalog": "NHGRI-EBI GWAS Catalog results showing other complex traits previously associated with the index variant or any other genome-wide significant variant in strong linkage disequilibrium with it (r² > 0.8).",
    "Literature": "Studies that previously identified the locus, with the reported variant with the strongest p value in brackets.",
  };
  const tip = (k) => (LDEF[k] ? ` title="${esc(LDEF[k])}"` : "");

  function evidenceList(str) { // "GENE (x) | GENE (y) | ..." as a short list, the rest folded away
    if (!str) return '<span class="neg">–</span>';
    const items = str.split(" | ").map((t) => `<li>${esc(t)}</li>`);
    if (items.length <= 4) return `<ul>${items.join("")}</ul>`;
    return `<ul>${items.slice(0, 4).join("")}</ul><details><summary>${items.length - 4} more</summary><ul>${items.slice(4).join("")}</ul></details>`;
  }

  function drawLocus() {
    const box = $("#locus");
    const l = LS.sel != null ? LOCI.find((x) => x.locus === LS.sel) : null;
    box.hidden = !l;
    if (!l) { box.innerHTML = ""; return; }
    const hits = MKEYS.map((m) => l.hits.find((h) => h.model === m)).filter(Boolean);
    const statRows = hits.map((h) => `<tr><td>${SHORT[h.model]}</td><td>${esc(h.id)}<small>chr${esc(l.chr)}:${fmtInt(h.bp)}</small></td><td>${esc(h.a1)}/${esc(h.a2)}</td><td class="num hide-s">${fmtF(h.freq, 2)}</td>` +
      `<td class="num">${fmtB(h.beta, h.se)}</td><td class="num sig">${fmtP(h.p)}</td><td class="num hide-s">${h.eta2 != null ? String(+h.eta2.toPrecision(3)) : "–"}</td><td class="num hide-s">${fmtInt(h.n)}</td></tr>`).join("");
    const ev = (label, f) => `<tr><th scope="row"${tip(label)}>${label}</th>${hits.map((h) => `<td>${f(h)}</td>`).join("")}</tr>`;
    const cs = (h) => Object.entries(h.cs).map(([k, v]) => `${k}: ${v ? esc(v) : "–"}`).join("<br>");
    const evidRows = ev("Prioritized gene", (h) => (h.prioritized ? `<b><i>${esc(h.prioritized)}</i></b>` : "–")) +
      ev("Nearest gene", (h) => (h.nearest ? `<i>${esc(h.nearest)}</i> (${esc(h.region || "")}${h.distance ? `, ${fmtInt(h.distance)} bp` : ""})` +
        (h.nearest_desc ? `<small>${esc(h.nearest_desc)}${h.nearest_type ? `, ${esc(h.nearest_type.replace(/_/g, " "))}` : ""}</small>` : "") : "–")) +
      ev("Credible set size", cs) +
      Object.keys(hits[0].evidence).map((k) => ev(k, (h) => evidenceList(h.evidence[k]))).join("");
    const catalog = hits.map((h) => h.catalog).find(Boolean);
    const nCat = catalog ? catalog.split(" | ").length : 0;
    box.innerHTML =
      `<button type="button" class="btn close" id="close-locus">Close</button>` + `<button type="button" class="btn close share" data-share>Copy link</button>` +
      `<h2><i>${esc(l.gene)}</i> · ${esc(l.cytoband)}${l.novel ? '<span class="novel-tag">novel</span>' : ""}</h2>` +
      `<div class="locus-meta"><p>Chromosome ${esc(l.chr)}, found for ${MKEYS.filter((m) => l.models.includes(m)).map((m) => MEASURES[m].toLowerCase()).join(", ")} brain age gap.</p>` +
      (l.literature ? `<p>Previously reported: ${esc(l.literature.replace(/_/g, " "))}</p>` : `<p>Not reported in earlier GWAS of brain age gap.</p>`) + `</div>` +
      `<h3>Lead variants</h3><div class="mini-wrap"><table class="mini"><thead><tr><th>Model</th>${["Variant", "A1/A2"].map((k) => `<th${tip(k)}>${k}</th>`).join("")}<th class="num hide-s"${tip("Freq.")}>Freq.</th><th class="num"${tip("β (SE)")}>β (SE)</th><th class="num"${tip("p")}>p</th><th class="num hide-s"${tip("ηp²")}>η<sub>p</sub>²</th><th class="num hide-s"${tip("N")}>N</th></tr></thead><tbody>${statRows}</tbody></table></div>` +
      `<h3>Gene prioritization</h3><div class="mini-wrap"><table class="mini evid"><thead><tr><th></th>${hits.map((h) => `<th>${SHORT[h.model]}</th>`).join("")}</tr></thead><tbody>${evidRows}</tbody></table></div>` +
      (catalog ? `<details class="catalog"><summary>Associated with ${nCat} trait${nCat > 1 ? "s" : ""} in the GWAS Catalog</summary><p>${esc(catalog.split(" | ").join("; "))}</p></details>` : "") +
      `<details class="defs"><summary>What the columns mean</summary><dl>${["Variant", "A1/A2", "Freq.", "β (SE)", "p", "ηp²", "N", "Prioritized gene", "Nearest gene", "Credible set size", ...Object.keys(hits[0].evidence), "GWAS Catalog", "Literature"]
        .filter((k) => LDEF[k]).map((k) => `<dt>${k === "ηp²" ? "η<sub>p</sub>²" : esc(k)}</dt><dd>${esc(LDEF[k])}</dd>`).join("")}</dl></details>`;
    $("#close-locus").addEventListener("click", () => selectLocus(null, false));
  }

  function selectLocus(id, scroll) {
    LS.sel = id;
    writeHash(); drawLociTable(); drawLocus();
    if (id != null && scroll) $("#locus").scrollIntoView({ behavior: reduced() ? "auto" : "smooth", block: "start" });
  }

  function bindLoci() {
    lociBound = true;
    let t;
    $("#loci-q").addEventListener("input", (e) => { clearTimeout(t); t = setTimeout(() => { LS.q = e.target.value.trim(); drawLociTable(); }, 180); });
    $("#loci-nov").addEventListener("change", (e) => { LS.nov = e.target.value; drawLociTable(); });
    $("#loci-more").addEventListener("click", () => { LS.more = !LS.more; drawLociTable(); if (!LS.more) $("#loci-table").scrollIntoView({ behavior: reduced() ? "auto" : "smooth", block: "start" }); });
    document.querySelectorAll("#loci-table th").forEach((th) => th.addEventListener("click", () => {
      const k = th.dataset.k;
      if (LS.sort === k) LS.asc = !LS.asc; else { LS.sort = k; LS.asc = !["beta", "models"].includes(k); }
      drawLociTable();
    }));
  }

  // ---------- genetic correlations with 38 selected GWAS (values as in the result table) ----------
  const stars = (x) => (x.fdr < 0.05 ? "**" : x.p < 0.05 ? "*" : "");
  async function drawRgSel() {
    const R = await load("rgsel"), T = R.traits;
    const ink = css("--ink"), muted = css("--muted"), surf = css("--surface"), band = css("--band");
    const narrow = window.innerWidth < 640;
    const nf = MKEYS.map((m) => `${T.filter((t) => t[m].fdr < 0.05).length} for ${SHORT[m].toLowerCase()}`);
    $("#rgsel-tally").textContent = `Combined European meta-analysis, n = 54,890. Correlations passing FDR < 5%: ${nf.join(", ")}.`;
    const n = T.length, ys = T.map((_, k) => n - 1 - k);
    const xl = (m) => (narrow ? { gm: "GM", wm: "WM", gwm: "GWM" }[m] : SHORT[m]);
    const xs = MKEYS.map(xl);
    const text = T.map((t) => MKEYS.map((m) => `<b>${esc(t.trait)}</b> (${esc(t.ref)})<br>${MEASURES[m]} brain age gap<br>r<sub>g</sub> = ${fmtR(t[m].rg)} (SE ${fmtF(t[m].se)})<br>p = ${fmtP(t[m].p)}, FDR = ${fmtP(t[m].fdr)}<br><span style="color:${muted}">h² of trait = ${fmtF(t.h2)} (SE ${fmtF(t.h2_se)})</span>`));
    const sig = (t, m) => t[m].p < 0.05;
    // two layers: a faint fill for p ≥ 0.05, the r_g colour scale for p < 0.05 (as in the paper figure)
    const base = { type: "heatmap", x: xs, y: ys, text, hovertemplate: "%{text}<extra></extra>", xgap: 2, ygap: 2 };
    const faint = { ...base, z: T.map((t) => MKEYS.map((m) => (sig(t, m) ? null : 0))), colorscale: [[0, band], [1, band]], showscale: false };
    const blue = "#2166ac", red = "#b2182b";
    const H = 50 + n * (narrow ? 20 : 22) + 70;
    const col = { ...base, z: T.map((t) => MKEYS.map((m) => (sig(t, m) ? Math.max(-0.3, Math.min(0.3, t[m].rg)) : null))), zmin: -0.3, zmax: 0.3,
      colorscale: [[0, blue], [0.5, surf.startsWith("#") ? surf : "#ffffff"], [1, red]],
      colorbar: { orientation: "h", title: { text: "Genetic correlation r<sub>g</sub>", side: "top" }, thickness: 10, len: narrow ? 0.9 : 0.8, lenmode: "fraction", tickangle: 0,
        x: 0.5, xanchor: "center", y: -10 / H, yanchor: "top", tickvals: narrow ? [-0.3, 0, 0.3] : [-0.3, -0.15, 0, 0.15, 0.3], ticktext: narrow ? ["−0.3", "0", "0.3"] : ["−0.3", "−0.15", "0", "0.15", "0.3"], outlinewidth: 0, tickfont: { size: 11 } } };
    const annotations = [];
    T.forEach((t, k) => MKEYS.forEach((m, j) => { const st = stars(t[m]); if (st) annotations.push({ x: xs[j], y: ys[k], text: st, showarrow: false, font: { size: 13, color: Math.abs(t[m].rg) > 0.17 ? "#fff" : ink }, yshift: -3 }); }));
    // frame around the cells, domain dividers and labels
    const line = { color: muted, width: 1 };
    const shapes = [{ type: "rect", xref: "x", yref: "y", x0: -0.5, x1: 2.5, y0: -0.5, y1: n - 0.5, line, fillcolor: "rgba(0,0,0,0)" }];
    T.forEach((t, k) => {
      if (k && t.domain !== T[k - 1].domain) shapes.push({ type: "line", xref: "x", yref: "y", x0: -0.5, x1: 2.5, y0: ys[k] + 0.5, y1: ys[k] + 0.5, line });
      if (!k || t.domain !== T[k - 1].domain) {
        const last = T.reduce((a, u, j) => (u.domain === t.domain ? j : a), k);
        annotations.push({ xref: "paper", x: 1.03, y: (ys[k] + ys[last]) / 2, text: narrow ? { Psychiatric: "Psych.", "Substance use": "Subst.", Neurological: "Neuro.", Personality: "Person.", Sleep: "Sleep", Cognition: "Cogn.", Anthropometric: "Anthrop.", Cardiovascular: "Cardio." }[t.domain] : t.domain,
          showarrow: false, xanchor: "left", font: { size: narrow ? 10.5 : 12, color: muted } });
      }
    });
    Plotly.react("rgselplot", [faint, col], hLayout({
      height: H, margin: { l: narrow ? 168 : 200, r: narrow ? 58 : 120, t: 34, b: 78 },
      xaxis: { side: "top", tickangle: 0, fixedrange: true, showgrid: false, zeroline: false, ticks: "outside", ticklen: 4, tickcolor: "rgba(0,0,0,0)", tickfont: { color: ink, size: narrow ? 11.5 : 12.5 } },
      yaxis: { tickvals: ys, ticktext: T.map((t) => esc(t.trait)), ticks: "outside", ticklen: 6, tickcolor: "rgba(0,0,0,0)", fixedrange: true, showgrid: false, zeroline: false, tickfont: { color: ink, size: narrow ? 10.5 : 12 }, range: [-0.5, n - 0.5] },
      shapes, annotations,
    }), { ...H_CFG, toImageButtonOptions: { filename: "brainage_rg_selected_traits", scale: 3 } });
  }

  // ---------- Mendelian randomization (GSMR and sensitivity methods, values as in the result table) ----------
  const MR_DIRS = [["to", "Trait → brain age gap"], ["from", "Brain age gap → trait"]];
  let MRD = null, mrSel = null;
  async function drawMr() {
    MRD = MRD || (await load("mr"));
    const T = MRD.traits, n = T.length;
    const ink = css("--ink"), muted = css("--muted"), surf = css("--surface"), band = css("--band");
    const narrow = window.innerWidth < 640;
    const S = PMETA.summary.mr;
    $("#mr-tally").textContent = `Combined European meta-analysis, n = 54,890. ${S.to_any} of ${n} traits show an effect on brain age gap and ${S.from_any} an effect of brain age gap at FDR < 5% for at least one model.`;
    // six columns: two directions × three models
    const cols = MR_DIRS.flatMap(([dk]) => MKEYS.map((m) => ({ dk, m })));
    const xs = cols.map((_, j) => j + (j >= 3 ? 1 : 0)), xAll = [0, 1, 2, 3, 4, 5, 6]; // column 3 is an empty spacer between the directions
    const spread = (f) => T.map((t) => { const r = cols.map((c) => f(t, c)); r.splice(3, 0, null); return r; });
    const ys = T.map((_, k) => n - 1 - k);
    const short = (m) => (narrow ? { gm: "G", wm: "W", gwm: "GW" }[m] : { gm: "Grey", wm: "White", gwm: "G + W" }[m]);
    const cell = (t, c) => t[c.dk][c.m];
    const hover = (t, c) => {
      const e = cell(t, c), dir = c.dk === "to" ? `${t.trait} → ${MEASURES[c.m].toLowerCase()} brain age gap` : `${MEASURES[c.m]} brain age gap → ${t.trait.toLowerCase()}`;
      if (!e) return `<b>${esc(t.trait)}</b> (${esc(t.ref)})<br>${esc(dir)}<br>no estimate available`;
      return `<b>${esc(t.trait)}</b> (${esc(t.ref)})<br>${esc(dir)}<br>GSMR b = ${fmtB(e.b, e.se)}<br>p = ${fmtP(e.p)}, FDR = ${fmtP(e.fdr)}` +
        `<br>${e.nheidi} of ${e.nsnp} instruments kept after HEIDI<br>${e.n05} of 10 MR methods with p < 0.05`;
    };
    const z = (t, c) => { const e = cell(t, c); return e && e.p < 0.05 ? Math.max(-6, Math.min(6, e.b / e.se)) : null; };
    const base = { type: "heatmap", x: xAll, y: ys, text: spread(hover), hovertemplate: "%{text}<extra></extra>", hoverongaps: false, xgap: 2, ygap: 2 };
    const faint = { ...base, z: spread((t, c) => (z(t, c) == null ? 0 : null)), colorscale: [[0, band], [1, band]], showscale: false };
    const H = 86 + n * (narrow ? 24 : 28) + 70;
    const col = { ...base, z: spread(z), zmin: -6, zmax: 6,
      colorscale: [[0, "#2166ac"], [0.5, surf.startsWith("#") ? surf : "#ffffff"], [1, "#b2182b"]],
      colorbar: { orientation: "h", title: { text: "GSMR z-score (b / SE)", side: "top" }, thickness: 10, len: narrow ? 0.9 : 0.6, tickangle: 0,
        x: 0.5, xanchor: "center", y: -10 / H, yanchor: "top", tickvals: [-6, -3, 0, 3, 6], ticktext: ["≤ −6", "−3", "0", "3", "≥ 6"], outlinewidth: 0, tickfont: { size: 11 } } };
    const annotations = [];
    T.forEach((t, k) => cols.forEach((c, j) => { const e = cell(t, c); const st = e ? (e.fdr < 0.05 ? "**" : e.p < 0.05 ? "*" : "") : "–";
      if (st) annotations.push({ x: xs[j], y: ys[k], text: st, showarrow: false, font: { size: 13, color: e && Math.abs(e.b / e.se) > 3.5 ? "#fff" : e ? ink : muted }, yshift: e ? -3 : 0 }); }));
    MR_DIRS.forEach(([dk, label], g) => annotations.push({ x: (xs[g * 3] + xs[g * 3 + 2]) / 2, y: 1, yref: "paper", yanchor: "bottom", yshift: 26, text: `<b>${narrow ? label.replace(/brain age gap/i, (x) => (x[0] === "B" ? "BAG" : "BAG")) : label}</b>`, showarrow: false, font: { size: narrow ? 11 : 12.5, color: ink } }));
    const line = { color: muted, width: 1 };
    const shapes = [0, 1].map((g) => ({ type: "rect", xref: "x", yref: "y", x0: xs[g * 3] - 0.5, x1: xs[g * 3 + 2] + 0.5, y0: -0.5, y1: n - 0.5, line, fillcolor: "rgba(0,0,0,0)" }));
    if (mrSel != null) { const y = n - 1 - mrSel; shapes.push({ type: "rect", xref: "paper", yref: "y", x0: 0, x1: 1, y0: y - 0.5, y1: y + 0.5, line: { color: css("--accent"), width: 2 }, fillcolor: "rgba(0,0,0,0)" }); }
    Plotly.react("mrplot", [faint, col], hLayout({
      height: H, margin: { l: narrow ? 142 : 190, r: 6, t: 62, b: 78 },
      xaxis: { side: "top", range: [-0.5, 6.5], tickvals: xs, ticktext: cols.map((c) => short(c.m)), tickangle: 0, fixedrange: true, showgrid: false, zeroline: false, ticks: "outside", ticklen: 4, tickcolor: "rgba(0,0,0,0)", tickfont: { color: ink, size: narrow ? 11 : 12 } },
      yaxis: { tickvals: ys, ticktext: T.map((t) => esc(t.trait)), ticks: "outside", ticklen: 6, tickcolor: "rgba(0,0,0,0)", fixedrange: true, showgrid: false, zeroline: false, tickfont: { color: ink, size: narrow ? 10.5 : 12.5 }, range: [-0.5, n - 0.5] },
      shapes, annotations,
    }), { ...H_CFG, toImageButtonOptions: { filename: "brainage_mendelian_randomization", scale: 3 } })
      .then((el) => { if (!el._mrClick) { el._mrClick = true; el.on("plotly_click", (ev) => { const pt = ev.points && ev.points[0]; if (pt && pt.x !== 3) selectMr(n - 1 - pt.y); }); } });
    drawMrDetail();
  }
  function selectMr(k) {
    mrSel = k; writeHash(); drawMr();
    if (k != null) setTimeout(() => $("#mr-detail").scrollIntoView({ behavior: reduced() ? "auto" : "smooth", block: "start" }), 50);
  }
  // which of the ten methods agree (p < 0.05), as a compact dot matrix
  function drawMrDetail() {
    const box = $("#mr-detail"), t = mrSel != null ? MRD.traits[mrSel] : null;
    box.hidden = !t;
    if (!t) { box.innerHTML = ""; return; }
    const cols = MR_DIRS.flatMap(([dk, dl]) => MKEYS.map((m) => ({ dk, dl, m, e: t[dk][m] })));
    const where = (c) => `${c.dl}, ${MEASURES[c.m].toLowerCase()}`;
    const dot = (c, k, name) => {
      if (!c.e) return '<td class="dotc neg">–</td>';
      const p = k === "gsmr" ? c.e.p : c.e.methods[k];
      if (p == null) return '<td class="dotc neg">–</td>';
      const on = p < 0.05;
      return `<td class="dotc tap" tabindex="0" data-info="${esc(`<b>${name}</b> · ${where(c)}: p = ${fmtP(p)}`)}"><span class="mdot${on ? " on" : ""}${on && c.e.b < 0 ? " neg-b" : ""}"></span></td>`;
    };
    const head = `<tr><th></th>${MR_DIRS.map(([, l]) => `<th colspan="3" class="grp">${l}</th>`).join("")}</tr><tr><th></th>${cols.map((c) => `<th class="dotc">${{ gm: "Grey", wm: "White", gwm: "G + W" }[c.m]}</th>`).join("")}</tr>`;
    const rows = MRD.methods.map(([k, name]) => `<tr><th scope="row">${esc(name)}</th>${cols.map((c) => dot(c, k, name)).join("")}</tr>`).join("");
    const gsmr = `<tr class="est-row"><th scope="row">GSMR b</th>${cols.map((c) => {
      if (!c.e) return '<td class="dotc neg">–</td>';
      const st = c.e.fdr < 0.05 ? "**" : c.e.p < 0.05 ? "*" : "";
      return `<td class="dotc tap${c.e.fdr < 0.05 ? " sig" : ""}" tabindex="0" data-info="${esc(`<b>GSMR</b> · ${where(c)}: b = ${fmtB(c.e.b, c.e.se)}, p = ${fmtP(c.e.p)}, FDR = ${fmtP(c.e.fdr)}; ${c.e.nheidi} of ${c.e.nsnp} instruments kept after HEIDI`)}">${(c.e.b < 0 ? "−" : "") + Math.abs(c.e.b).toPrecision(2)}${st ? `<sup>${st}</sup>` : ""}</td>`;
    }).join("")}</tr>`;
    const tally = `<tr class="tally-row"><th scope="row">Methods with p &lt; 0.05</th>${cols.map((c) => `<td class="dotc">${c.e ? `${c.e.n05}/10` : "–"}</td>`).join("")}</tr>`;
    box.innerHTML = `<button type="button" class="btn close" id="close-mr">Close</button>` + `<button type="button" class="btn close share" data-share>Copy link</button>` +
      `<h2>${esc(t.trait)}</h2><div class="locus-meta"><p>${esc(t.ref)}</p></div>` +
      `<div class="mini-wrap"><table class="mini mrdots"><thead>${head}</thead><tbody>${gsmr}${rows}${tally}</tbody></table></div>` +
      `<p class="mr-info" id="mr-info" aria-live="polite">Tap or hover a value or dot for details.</p>` +
      `<p class="note"><span class="mdot on"></span> p &lt; 0.05 with a positive GSMR estimate, <span class="mdot on neg-b"></span> with a negative one, <span class="mdot"></span> p ≥ 0.05. GSMR b with * p &lt; 0.05 and ** FDR &lt; 0.05. GSMR uses the instruments kept after the HEIDI outlier test; the other methods use all instruments.</p>`;
    $("#close-mr").addEventListener("click", () => selectMr(null));
    // details for a value or dot: on hover with a mouse, on tap with a finger, on focus with the keyboard
    const info = $("#mr-info");
    box.querySelectorAll("td.tap").forEach((td) => {
      const show = () => { box.querySelectorAll("td.tap.on").forEach((x) => x.classList.remove("on")); td.classList.add("on"); info.innerHTML = td.dataset.info; };
      td.addEventListener("click", show); td.addEventListener("mouseenter", show); td.addEventListener("focus", show);
    });
  }

  // ---------- fastBAT gene-based tests (values as in the result table) ----------
  const GS = { q: "", sig: "bonf", sort: "p", asc: true, more: false, bound: false, sel: null };
  let GENES = null;
  async function drawGenes() {
    GENES = GENES || (await load("genes"));
    if (!GS.bound) bindGenes();
    const m = state.m, S = PMETA.summary.genes, sm = S[m], model = MEASURES[m].toLowerCase();
    $("#genes-title").textContent = `${fmtInt(sm.bonf)} genes associated with ${model} brain age gap`;
    $("#genes-tally").innerHTML = `fastBAT gene-based tests, combined European meta-analysis, n = 54,890. ${fmtInt(sm.bonf)} of ${fmtInt(S.nTested)} genes pass Bonferroni (p &lt; ${fmtP(S.bonf)}) in ${fmtInt(sm.loci_bonf)} independent loci; ${fmtInt(sm.fdr)} pass FDR &lt; 5%.`;
    const img = $("#genes-img"); img.src = `img/genes.${m}.png`; img.alt = `Gene-based Manhattan plot for ${model} brain age gap`;
    if (GS.sel && !GENES.genes.some((x) => x.gene === GS.sel)) GS.sel = null;
    drawGenesTable(); drawGene();
  }
  function drawGenesTable() {
    const m = state.m, bonf = GENES.bonf;
    const pass = (g, mm) => (GS.sig === "bonf" ? g[mm].p < bonf : g[mm].fdr < 0.05);
    const q = GS.q.toLowerCase();
    let rows = GENES.genes.filter((g) => pass(g, m) && (!q || `${g.gene} ${g.desc || ""} ${g.cyto}`.toLowerCase().includes(q)));
    const dir = GS.asc ? 1 : -1;
    const key = { gene: (g) => g.gene.toLowerCase(), pos: (g) => (g.chr === "X" ? 23 : +g.chr) * 1e9 + g.start, nsnp: (g) => g.nsnp,
      models: (g) => MKEYS.filter((mm) => pass(g, mm)).length, p: (g) => g[m].p, fdr: (g) => g[m].fdr, lead: (g) => (g[m].lead || "").toLowerCase() }[GS.sort];
    rows.sort((a, b) => (key(a) < key(b) ? -1 : key(a) > key(b) ? 1 : 0) * dir);
    const total = GENES.genes.filter((g) => pass(g, m)).length;
    const shown = GS.more ? rows : rows.slice(0, 10);
    $("#genes-table tbody").innerHTML = shown.map((g) => `<tr data-gene="${esc(g.gene)}" class="${g.gene === GS.sel ? "sel" : ""}" tabindex="0">` +
      `<td class="gene"><b><i>${esc(g.gene)}</i></b><small>${esc(g.desc || "")}</small></td>` +
      `<td class="hide-s">${esc(g.cyto)}<small>chr${esc(g.chr)}:${fmtInt(g.start)}–${fmtInt(g.end)}</small></td>` +
      `<td class="hide-m">${g[m].index ? '<span class="novel-tag lead-tag">lead gene</span>' : g[m].lead ? `<i>${esc(g[m].lead)}</i>` : "–"}</td>` +
      `<td class="num hide-m">${fmtInt(g.nsnp)}</td>` +
      `<td class="models hide-s">${MKEYS.filter((mm) => pass(g, mm)).map((mm) => mchip(mm)).join("")}</td>` +
      `<td class="num${g[m].p < bonf ? " sig" : ""}">${fmtP(g[m].p)}</td><td class="num hide-s">${fmtP(g[m].fdr)}</td></tr>`).join("") ||
      `<tr><td colspan="7">No genes match. Clear the search or choose another threshold.</td></tr>`;
    $("#genes-table tbody").querySelectorAll("tr[data-gene]").forEach((tr) => {
      const go = () => selectGene(tr.dataset.gene === GS.sel ? null : tr.dataset.gene, true);
      tr.addEventListener("click", go);
      tr.addEventListener("keydown", (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); go(); } });
    });
    const btn = $("#genes-more");
    btn.hidden = rows.length <= 10;
    btn.textContent = GS.more ? "Show first 10 only" : `Show all ${fmtInt(rows.length)} genes`;
    $("#genes-count").textContent = `${rows.length === total ? fmtInt(total) : `${fmtInt(rows.length)} of ${fmtInt(total)}`} genes pass ${GS.sig === "bonf" ? "Bonferroni" : "FDR < 5%"} for ${MEASURES[m].toLowerCase()} brain age gap; “Models” lists every model for which the gene passes.`;
    document.querySelectorAll("#genes-table th").forEach((th) => th.setAttribute("aria-sort", th.dataset.k === GS.sort ? (GS.asc ? "ascending" : "descending") : "none"));
  }
  function selectGene(name, scroll) {
    GS.sel = name; writeHash(); drawGenesTable(); drawGene();
    if (name && scroll) $("#gene-detail").scrollIntoView({ behavior: reduced() ? "auto" : "smooth", block: "start" });
  }
  function drawGene() {
    const box = $("#gene-detail"), g = GS.sel ? GENES.genes.find((x) => x.gene === GS.sel) : null;
    box.hidden = !g;
    if (!g) { box.innerHTML = ""; return; }
    const bonf = GENES.bonf;
    const rows = MKEYS.map((m) => { const e = g[m];
      return `<tr><td>${SHORT[m]}</td><td class="num${e.p < bonf ? " sig" : ""}">${fmtP(e.p)}</td><td class="num">${fmtP(e.fdr)}</td>` +
        `<td>${e.index ? "this gene" : e.lead ? `<a href="#" class="lead-jump" data-gene="${esc(e.lead)}"><i>${esc(e.lead)}</i></a>` : "–"}</td>` +
        `<td>${e.p < bonf ? "Bonferroni" : e.fdr < 0.05 ? "FDR &lt; 5%" : "–"}</td></tr>`; }).join("");
    const links = [`<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=${encodeURIComponent(g.gene)}" target="_blank" rel="noopener">GeneCards</a>`]
      .concat(g.entrez ? [`<a href="https://www.ncbi.nlm.nih.gov/gene/${encodeURIComponent(g.entrez)}" target="_blank" rel="noopener">NCBI Gene</a>`] : []);
    box.innerHTML = `<button type="button" class="btn close" id="close-gene">Close</button>` +
      `<button type="button" class="btn close share" data-share>Copy link</button>` +
      `<h2><i>${esc(g.gene)}</i> · ${esc(g.cyto)}</h2>` +
      `<div class="locus-meta"><p>${esc(g.desc || "")}</p><p>chr${esc(g.chr)}:${fmtInt(g.start)}–${fmtInt(g.end)} (hg19), ${fmtInt(g.nsnp)} SNPs tested · ${links.join(" · ")}</p></div>` +
      `<div class="mini-wrap"><table class="mini"><thead><tr><th>Model</th><th class="num">p</th><th class="num">FDR</th><th>Lead gene of locus</th><th>Passes</th></tr></thead><tbody>${rows}</tbody></table></div>` +
      `<p class="note">Genes within 3 Mb of a more significant gene form one locus; its lead gene is the most significant gene, separately for each brain age model. Select a lead gene to open it.</p>`;
    $("#close-gene").addEventListener("click", () => selectGene(null));
    box.querySelectorAll(".lead-jump").forEach((a) => a.addEventListener("click", (e) => {
      e.preventDefault();
      const name = a.dataset.gene;
      if (!GENES.genes.some((x) => x.gene === name)) return;
      GS.q = ""; $("#genes-q").value = ""; selectGene(name, true);
    }));
  }
  function bindGenes() {
    GS.bound = true;
    let t;
    $("#genes-q").addEventListener("input", (e) => { clearTimeout(t); t = setTimeout(() => { GS.q = e.target.value.trim(); drawGenesTable(); }, 180); });
    $("#genes-sig").addEventListener("change", (e) => { GS.sig = e.target.value; GS.more = false; drawGenesTable(); });
    $("#genes-more").addEventListener("click", () => { GS.more = !GS.more; drawGenesTable(); if (!GS.more) $("#genes-table").scrollIntoView({ behavior: reduced() ? "auto" : "smooth", block: "start" }); });
    document.querySelectorAll("#genes-table th").forEach((th) => th.addEventListener("click", () => {
      const k = th.dataset.k;
      if (GS.sort === k) GS.asc = !GS.asc; else { GS.sort = k; GS.asc = !["nsnp", "models"].includes(k); }
      drawGenesTable();
    }));
  }

  // ---------- one search across all results ----------
  async function searchIndex() {
    const [all, R, M, L, G] = await Promise.all([load("all"), load("rgsel"), load("mr"), load("loci"), load("genes")]);
    const sigP = (i) => MKEYS.some((m) => lt05(all[m].q[i]));
    const sigR = (i) => MKEYS.some((m) => lt05(RG[m].q[i]));
    const items = [];
    PMETA.desc.forEach((d, i) => items.push({ g: "Phenotypic correlations", text: d, sub: `${PMETA.id[i]}${sigP(i) ? " · FDR < 5%" : ""}`, sig: sigP(i), hay: `${d} ${PMETA.id[i]}`, go: { d: "pheno", t: PMETA.id[i] } }));
    RG.desc.forEach((d, i) => items.push({ g: "Genetic correlations, UK Biobank traits", text: d, sub: `field ${RG.field[i] ?? "–"}${sigR(i) ? " · FDR < 5%" : ""}`, sig: sigR(i), hay: d, go: { d: "rg", t: RG.id[i] } }));
    R.traits.forEach((t) => items.push({ g: "Genetic correlations, 38 selected traits", text: t.trait, sub: t.ref, sig: MKEYS.some((m) => t[m].fdr < 0.05), hay: `${t.trait} ${t.ref}`, go: { d: "rgsel" } }));
    M.traits.forEach((t, k) => items.push({ g: "Mendelian randomization", text: t.trait, sub: t.ref, sig: ["to", "from"].some((dk) => MKEYS.some((m) => t[dk][m] && t[dk][m].fdr < 0.05)), hay: `${t.trait} ${t.ref}`, go: { d: "mr", trait: k } }));
    L.loci.forEach((l) => items.push({ g: "GWAS loci", text: l.gene, sub: `${l.cytoband}${l.novel ? " · novel" : ""}`, sig: true, italic: true,
      hay: [l.gene, l.cytoband, ...l.hits.flatMap((h) => [h.id, h.nearest || "", h.prioritized || ""])].join(" "), go: { d: "loci", locus: l.locus } }));
    G.genes.forEach((g) => { const best = MKEYS.reduce((a, mm) => (g[mm].p < g[a].p ? mm : a), "gwm");
      items.push({ g: "Gene-based tests", text: g.gene, sub: `${g.cyto} · p = ${fmtP(g[best].p, false)}`, sig: MKEYS.some((mm) => g[mm].p < G.bonf), italic: true,
        hay: `${g.gene} ${g.desc || ""} ${g.cyto}`, go: { d: "genes", m: best, gq: g.gene } }); });
    items.forEach((x) => (x.hay = x.hay.toLowerCase()));
    return items;
  }
  const GROUPS = ["GWAS loci", "Gene-based tests", "Phenotypic correlations", "Genetic correlations, UK Biobank traits", "Genetic correlations, 38 selected traits", "Mendelian randomization"];
  function bindSearch() {
    const input = $("#gq"), list = $("#gq-list");
    let index = null, opts = [], active = -1, tmr;
    const close = () => { list.hidden = true; input.setAttribute("aria-expanded", "false"); active = -1; };
    const mark = (text, q) => { const i = text.toLowerCase().indexOf(q); return i < 0 ? esc(text) : `${esc(text.slice(0, i))}<mark>${esc(text.slice(i, i + q.length))}</mark>${esc(text.slice(i + q.length))}`; };
    const render = async () => {
      const q = input.value.trim().toLowerCase();
      if (q.length < 2) return close();
      index = index || (await searchIndex());
      opts = []; let html = "";
      GROUPS.forEach((g) => {
        const hits = index.filter((x) => x.g === g && x.hay.includes(q))
          .sort((a, b) => (b.text.toLowerCase().startsWith(q) - a.text.toLowerCase().startsWith(q)) || (b.sig - a.sig) || a.text.localeCompare(b.text));
        if (!hits.length) return;
        html += `<div class="gq-group">${esc(g)} <span>${hits.length}</span></div>`;
        hits.slice(0, 4).forEach((x) => { opts.push(x.go); html += `<div class="gq-opt" role="option" data-k="${opts.length - 1}"><span class="gq-t">${x.italic ? `<i>${mark(x.text, q)}</i>` : mark(x.text, q)}</span><span class="gq-s">${esc(x.sub)}</span></div>`; });
        if (hits.length > 4 && (g === "Phenotypic correlations" || g.startsWith("Genetic correlations, UK"))) {
          opts.push({ d: g === "Phenotypic correlations" ? "pheno" : "rg", q: input.value.trim() });
          html += `<div class="gq-opt gq-more" role="option" data-k="${opts.length - 1}">Show all ${fmtInt(hits.length)} matches</div>`;
        } else if (hits.length > 4) html += `<div class="gq-note">and ${hits.length - 4} more</div>`;
      });
      list.innerHTML = html || `<div class="gq-note">No results for “${esc(input.value.trim())}”.</div>`;
      list.hidden = false; input.setAttribute("aria-expanded", "true"); active = -1;
    };
    const pick = (k) => { const spec = opts[k]; if (!spec) return; close(); input.blur(); go(spec); };
    const highlight = () => list.querySelectorAll(".gq-opt").forEach((o) => o.classList.toggle("on", +o.dataset.k === active));
    input.addEventListener("input", () => { clearTimeout(tmr); tmr = setTimeout(render, 150); });
    input.addEventListener("focus", () => { if (input.value.trim().length >= 2) render(); else searchIndex().then((x) => (index = index || x)); });
    input.addEventListener("keydown", (e) => {
      if (e.key === "Escape") return close();
      if (list.hidden || !opts.length) return;
      if (e.key === "ArrowDown") { e.preventDefault(); active = Math.min(opts.length - 1, active + 1); highlight(); }
      else if (e.key === "ArrowUp") { e.preventDefault(); active = Math.max(0, active - 1); highlight(); }
      else if (e.key === "Enter") { e.preventDefault(); pick(active < 0 ? 0 : active); }
    });
    list.addEventListener("mousedown", (e) => { const o = e.target.closest(".gq-opt"); if (o) { e.preventDefault(); pick(+o.dataset.k); } });
    document.addEventListener("click", (e) => { if (!e.target.closest(".gsearch")) close(); });
  }

  // ---------- actions ----------
  function select(id, fromPlot) {
    state.t = id;
    writeHash();
    drawPlot();
    drawTable();
    drawDetail().then(() => {
      // clicking a point or a table row takes the reader to the trait panel under the table
      if (id) $("#trait").scrollIntoView({ behavior: reduced() ? "auto" : "smooth", block: "start" });
    });
  }

  const MODEL_TEXT = { gm: "grey matter", wm: "white matter", gwm: "grey + white matter" };
  // UK Biobank imaging discovery sample (results/mri/accuracy.sample.txt); traits have up to this many participants
  const PHEWAS_N = { all: "Up to 32,634 participants. ", female: "Up to 17,084 women. ", male: "Up to 15,550 men. ", sexdiff: "Up to 17,084 women and 15,550 men. " };
  const SAMPLE_TEXT = { all: "all participants", female: "women", male: "men", sexdiff: "women versus men" };
  function syncControls() {
    const cur = { d: state.d === "rgsel" ? "rg" : state.d === "genes" ? "loci" : state.d, g: state.d === "rgsel" ? "sel" : "ukb", lv: state.d === "genes" ? "genes" : "loci", v: state.v, m: state.d === "loci" ? LS.m : state.m, s: state.s, l: String(state.labels) };
    document.querySelectorAll(".seg").forEach((g) => g.querySelectorAll("button").forEach((b) => {
      const on = cur[g.dataset.key] === b.dataset.val;
      b.setAttribute("aria-checked", on); b.tabIndex = on ? 0 : -1;
      if (g.dataset.key === "v") b.hidden = !viewsFor().includes(b.dataset.val);
      if (g.dataset.key === "m" && b.dataset.val === "all") b.hidden = state.d !== "loci"; // "All together" exists only for the loci
    }));
    $("#ctl-sample").hidden = state.d !== "pheno";
    $("#ctl-plot").hidden = state.d === "cmp" || isOwn(state.d);
    $("#ctl-model").hidden = state.d === "herit" || state.d === "rgsel" || state.d === "mr"; // these panels show all three models side by side
    $("#ctl-rgset").hidden = state.d !== "rg" && state.d !== "rgsel";
    $("#ctl-level").hidden = state.d !== "loci" && state.d !== "genes";
    $("#genes").hidden = state.d !== "genes";
    $("#rgsel").hidden = state.d !== "rgsel"; $("#mr").hidden = state.d !== "mr";
    const own = isOwn(state.d);
    $("#explore-body").hidden = own; $("#split").hidden = own;
    $("#herit").hidden = state.d !== "herit"; $("#loci").hidden = state.d !== "loci";
    // "How to read this" shows only the explanations for the current result set
    document.querySelectorAll("#about [data-for]").forEach((el) => (el.hidden = !el.dataset.for.split(" ").includes(state.d)));
    updateCards(); updateNext();
    if (lastD !== state.d) { // a short fade when the result set changes
      const el = sectionEl();
      if (lastD && !reduced()) { el.classList.remove("enter"); void el.offsetWidth; el.classList.add("enter"); }
      lastD = state.d;
    }
  }
  let lastD = null;
  const SECTION = { herit: "#herit", loci: "#loci", genes: "#genes", rgsel: "#rgsel", mr: "#mr" };
  const sectionEl = () => $(SECTION[state.d] || "#explore-body");

  // ---------- result cards, jump chips and previous/next ----------
  function updateCards() {
    const S = PMETA && PMETA.summary;
    if (!S) return;
    const pct = (x) => Math.round(x * 100);
    const h = S.herit.ldsc_combined, hs = MKEYS.map((k) => h[k].h2);
    const set = (k, num, sub) => { $(`[data-num="${k}"]`).innerHTML = num; if (sub) $(`[data-sub="${k}"]`).innerHTML = sub; };
    // the same numbers whichever model is selected: correlations count traits at FDR < 5% for at least one model
    set("herit", `${pct(Math.min(...hs))}–${pct(Math.max(...hs))}%`, "explained by common genetic variants");
    set("loci", fmtInt(S.loci.n), `genome-wide significant, ${fmtInt(S.loci.novel)} novel`);
    set("pheno", fmtInt(S.samples.all.any), `of ${fmtInt(S.samples.all.nTested)} traits at FDR &lt; 5%`);
    set("rg", fmtInt(S.rg.any + S.rgsel.any), `of ${fmtInt(S.rg.nTested + S.rgsel.n)} traits at FDR &lt; 5%`); // UK Biobank traits plus the 38 selected traits
    set("cmp", `r vs. r<sub>g</sub>`, `side by side for ${fmtInt(S.rg.pairs.gwm.n)} traits`);
    set("mr", `${fmtInt(S.mr.to_any)} of ${fmtInt(S.mr.n)}`, "traits affect brain age gap at FDR &lt; 5%");
    $("#rc-note").textContent = "Across the three brain age models; correlation counts include traits with FDR < 5% for at least one model.";
  }
  const ORDER = ["herit", "loci", "genes", "pheno", "rg", "rgsel", "cmp", "mr"];
  const ORDER_NAME = { herit: "Heritability", loci: "Genomic loci", genes: "Gene-based tests", pheno: "Phenotypic correlations", rg: "Genetic correlations, UK Biobank traits",
    rgsel: "Genetic correlations, 38 selected traits", cmp: "Phenotypic vs. genetic correlations", mr: "Mendelian randomization" };
  function updateNext() {
    const k = ORDER.indexOf(state.d);
    [["#nn-prev", ORDER[k - 1]], ["#nn-next", ORDER[k + 1]]].forEach(([sel, d]) => {
      const b = $(sel); b.hidden = !d; b.dataset.d = d || ""; if (d) b.querySelector(".nn-name").textContent = ORDER_NAME[d];
    });
  }
  const scrollToResults = () => {
    const el = sectionEl();
    el.scrollIntoView({ behavior: reduced() ? "auto" : "smooth", block: "start" });
  };
  // opens a result set directly (used by previous/next), optionally with a locus or trait selected
  async function go(spec) {
    if (spec.m) state.m = spec.m;
    if (isOwn(spec.d)) {
      if (spec.d === "mr") mrSel = null;
      if (spec.d === "genes" && spec.gq) { GS.q = spec.gq; GS.sig = "fdr"; GS.more = false; $("#genes-q").value = GS.q; $("#genes-sig").value = "fdr"; }
      if (spec.d === "loci") { LS.m = "all"; LS.q = ""; LS.nov = "all"; LS.more = false; $("#loci-q").value = ""; $("#loci-nov").value = "all"; }
      if (!isOwn(state.d)) heritFrom = state.d;
      state.d = spec.d; state.t = null; LS.sel = null; GS.sel = null;
      syncControls(); writeHash();
      await OWN[spec.d]();
      const l = LOCI && (spec.locus != null ? LOCI.find((x) => x.locus === spec.locus) : spec.gene ? LOCI.find((x) => x.gene === spec.gene) : null);
      if (l) return selectLocus(l.locus, true);
      if (spec.d === "mr" && spec.trait != null) return selectMr(spec.trait);
      if (spec.d === "genes" && spec.gq) return selectGene(spec.gq, true);
      return scrollToResults();
    }
    state.d = spec.d; state.s = spec.s || "all"; state.v = state.s === "sexdiff" ? "sex" : "manhattan";
    state.q = spec.q || ""; $("#q").value = state.q; state.sig = "all"; $("#sig").value = "all"; state.cat = null; state.page = 0; state.t = null;
    META = isGen(state.d) ? RG : PMETA;
    syncControls();
    await refresh(true);
    await drawDetail();
    if (spec.t) select(spec.t); else scrollToResults();
  }
  function choose(key, val) {
    if (key === "l") { state.labels = +val; syncControls(); return drawPlot(); }
    if (key === "d") return setMode(val);
    if (key === "g") return setMode(val === "sel" ? "rgsel" : "rg");
    if (key === "lv") return setMode(val);
    if (key === "v") return setView(val);
    if (key === "m") {
      if (state.d === "genes") {
        if (val === state.m || val === "all") return;
        state.m = val; GS.more = false; syncControls(); writeHash(); return drawGenes();
      }
      if (state.d === "loci") {
        if (val === LS.m) return;
        LS.m = val; LS.more = false; if (val !== "all") state.m = val;
        syncControls(); writeHash(); return drawLoci();
      }
      if (val === state.m) return;
      state.m = val; syncControls();
      return refresh(true).then(drawDetail);
    }
    if (key === "s") {
      if (val === state.s) return;
      state.s = val; state.page = 0;
      // women vs. men opens on the scatter plot; the scatter plot exists only for women vs. men
      if (val === "sexdiff") state.v = "sex";
      else if (state.v === "sex") state.v = "manhattan";
      syncControls(); return refresh(true).then(drawDetail);
    }
  }
  function viewsFor() { // plot types offered for the current results and sample
    if (state.d !== "pheno") return VIEWS[state.d];
    return state.s === "sexdiff" ? ["sex", "manhattan", "volcano"] : ["manhattan", "volcano"];
  }
  function setView(v) {
    if (v === state.v || !viewsFor().includes(v)) return;
    state.v = v; state.page = 0;
    syncControls();
    refresh(true).then(drawDetail);
  }
  let heritFrom = "pheno"; // result set to return to when leaving the heritability panels
  function setMode(d) {
    if (d === state.d) return;
    if (isOwn(d)) {
      if (!isOwn(state.d)) heritFrom = state.d;
      state.d = d; state.t = null; LS.sel = null; mrSel = null; GS.sel = null;
      syncControls(); writeHash();
      return OWN[d]();
    }
    if (isOwn(state.d)) state.d = heritFrom;
    if (d === state.d) { syncControls(); return refresh(true).then(drawDetail); }
    state.t = null; // each result set opens on its overview, without a selected trait
    state.d = d;
    if (d === "cmp") state.v = "pg";
    else if (!viewsFor().includes(state.v)) state.v = d === "pheno" && state.s === "sexdiff" ? "sex" : "manhattan";
    META = isGen(d) ? RG : PMETA; state.page = 0;
    syncControls();
    refresh(true).then(drawDetail);
  }

  const VIEW_NOTES = {
    pheno: {
      manhattan: "Traits grouped by category, as in the paper figure. Height shows significance; triangles point in the direction of the association.",
      volcano: "Effect size against significance. Points far left or right are the strongest associations; the lines mark FDR and Bonferroni thresholds.",
      sex: "Each trait’s association in women against men. Coloured points are significant in at least one sex (FDR < 5%); outlined points differ between the sexes.",
    },
    rg: {
      manhattan: "Genetic correlations grouped by category. Height shows significance; triangles point in the direction of the genetic correlation.",
      volcano: "Genetic correlation against significance. Points far left or right share the most genetic signal with brain age gap.",
    },
    cmp: {
      pg: "Genetic correlation against phenotypic correlation (all participants), as in the paper. Trait pairs are matched by description; multinomial models are excluded and only traits with SNP heritability h²/SE > 1.96 are kept. Coloured points are significant in one analysis, outlined points in both.",
    },
  };

  function updateText() {
    const n = view.rows.length, nSig = view.rows.filter(isSig).length, nB = view.rows.filter((r) => r.p < view.bonf).length;
    const model = `<span class="pick">${MODEL_TEXT[state.m]}</span>`;
    const ps = RG.pairStats[state.m];
    $("#sentence").innerHTML = state.d === "cmp"
      ? `Genetic versus phenotypic correlations of ${model} brain age gap across <strong>${fmtInt(n)}</strong> trait pairs`
      : state.d === "rg"
      ? `Genetic correlations of ${model} brain age gap with <strong>${fmtInt(n)}</strong> UK Biobank traits`
      : `Associations of ${model} brain age gap with <strong>${fmtInt(n)}</strong> UK Biobank traits` +
        (state.s === "all" ? "" : ` in <span class="pick">${SAMPLE_TEXT[state.s]}</span>`); // sample size goes in the line below
    $("#tally").textContent = state.d === "cmp"
      ? `Across pairs, r = ${ps.r.toFixed(2)} between genetic and phenotypic correlations; ${fmtInt(nSig)} pairs are significant in both analyses (FDR < 5%).`
      : state.d === "rg"
      ? `${fmtInt(nSig)} genetic correlations pass FDR < 5%, ${fmtInt(nB)} pass Bonferroni correction.`
      : PHEWAS_N[state.s] + (state.s === "sexdiff"
        ? `${fmtInt(nSig)} traits differ between women and men at FDR < 5%, ${fmtInt(nB)} after Bonferroni correction.`
        : `${fmtInt(nSig)} associations pass FDR < 5%, ${fmtInt(nB)} pass Bonferroni correction.`);
    let note = VIEW_NOTES[state.d][state.v];
    $("#view-note").textContent = note;
    $("#keymarks").innerHTML = keyHtml();
    syncControls();
  }

  async function refresh(rebuild) {
    if (rebuild) {
      $("#plot-loading").hidden = false;
      await buildView();
      updateText();
    }
    writeHash();
    drawCats();
    drawPlot();
    drawTable();
    drawOverview();
  }

  function download() {
    const rows = sortedFiltered();
    const q = (v) => v == null ? "" : /[",\n]/.test(String(v)) ? `"${String(v).replace(/"/g, '""')}"` : v;
    let head, cells;
    if (state.d === "cmp") {
      head = ["ukb_field", "description", "category", "phewas_varName", "h2_obs", "rg", "rg_se", "rg_pvalue", "rg_fdr", "phenotypic_r", "phenotypic_pvalue", "phenotypic_fdr"];
      cells = (r) => [META.field[r.i], META.desc[r.i], META.categories[META.cat[r.i]], PMETA.id[r.pidx], r.n, r.r, r.se, r.p, r.q, r.pr, r.pp, r.pq];
    } else if (state.d === "rg") {
      head = ["ukb_field", "description", "category", "path", "h2_obs", "rg", "se", "pvalue", "fdr", "phenotypic_r", "phenotypic_fdr"];
      cells = (r) => [META.field[r.i], META.desc[r.i], META.categories[META.cat[r.i]], META.paths[META.path[r.i]], r.n, r.r, r.se, r.p, r.q, r.pr, r.pq];
    } else {
      head = ["varName", "ukb_field", "description", "category", "path", "regression", "n", "r", "beta", "se", "pvalue", "fdr"];
      cells = (r) => [META.id[r.i], META.field[r.i], META.desc[r.i], META.categories[META.cat[r.i]], META.paths[META.path[r.i]],
        META.resTypes[META.rt[r.i]], r.n, r.r, r.b, r.se, r.p, r.q];
    }
    const lines = [head.join(",")].concat(rows.map((r) => cells(r).map(q).join(",")));
    const blob = new Blob([lines.join("\n")], { type: "text/csv" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = state.d === "cmp" ? `brainage_rg_vs_phenotypic_${state.m}.csv` : state.d === "rg" ? `brainage_rg_${state.m}.csv` : `brainage_phewas_${state.m}_${state.s}.csv`;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  }

  const BIB = `@article{jawinski2025brainage,
  author  = {Jawinski, Philippe and Forstbach, Helena and Kirsten, Holger and Beyer, Frauke and Villringer, Arno and Witte, A. Veronica and Scholz, Markus and Ripke, Stephan and Markett, Sebastian},
  title   = {Genome-wide analysis of brain age identifies 59 associated loci and unveils relationships with mental and physical health},
  journal = {Nature Aging},
  volume  = {5},
  pages   = {2086--2103},
  year    = {2025},
  doi     = {10.1038/s43587-025-00962-7}
}`;

  function bind() {
    document.querySelectorAll(".seg").forEach((g) => {
      const key = g.dataset.key;
      g.querySelectorAll("button").forEach((btn) => btn.addEventListener("click", () => choose(key, btn.dataset.val)));
      g.addEventListener("keydown", (e) => { // arrow keys move within a group, skipping hidden buttons
        if (!["ArrowLeft", "ArrowRight"].includes(e.key)) return;
        const bs = [...g.querySelectorAll("button")].filter((x) => !x.hidden), k = bs.findIndex((x) => x.getAttribute("aria-checked") === "true");
        const nxt = bs[(k + (e.key === "ArrowRight" ? 1 : bs.length - 1)) % bs.length];
        e.preventDefault(); choose(key, nxt.dataset.val); nxt.focus();
      });
    });
    let tmr;
    $("#q").addEventListener("input", (e) => { clearTimeout(tmr); tmr = setTimeout(() => { state.q = e.target.value.trim(); state.page = 0; refresh(false); }, 180); });
    $("#sig").addEventListener("change", (e) => { state.sig = e.target.value; state.page = 0; refresh(false); });

    $("#prev").addEventListener("click", () => { state.page--; drawTable(); });
    $("#next").addEventListener("click", () => { state.page++; drawTable(); });
    $("#download").addEventListener("click", download);
    bindSearch();
    // "Copy link" in the detail panels copies the address of the current view, including the open trait or locus
    document.addEventListener("click", (e) => {
      const b = e.target.closest("[data-share]"); if (!b) return;
      const done = (t) => { b.textContent = t; setTimeout(() => (b.textContent = "Copy link"), 1800); };
      (navigator.clipboard ? navigator.clipboard.writeText(location.href) : Promise.reject()).then(() => done("Link copied"), () => done("Copy the address bar"));
    });
    document.querySelectorAll(".next-nav .nn").forEach((b) => b.addEventListener("click", async () => {
      const d = b.dataset.d; if (!d) return;
      if (isOwn(d)) await go({ d });
      else { setMode(d); await new Promise((r) => setTimeout(r, 60)); scrollToResults(); }
    }));
    $("#copy-bib").addEventListener("click", (e) => {
      navigator.clipboard.writeText(BIB).then(() => { e.target.textContent = "BibTeX copied"; setTimeout(() => (e.target.textContent = "Copy BibTeX"), 2000); },
        () => { e.target.textContent = "Copy failed, select the citation above instead"; });
    });
    document.querySelectorAll("#table th").forEach((th) => th.addEventListener("click", () => {
      const k = th.dataset.k;
      if (state.sort === k) state.asc = !state.asc; else { state.sort = k; state.asc = !["r", "n"].includes(k); }
      state.page = 0; drawTable();
    }));
    let rz, lastNarrow = window.innerWidth < 900;
    window.addEventListener("resize", () => { clearTimeout(rz); rz = setTimeout(() => {
      lastNarrow = window.innerWidth < 900; if (isOwn(state.d)) OWN[state.d](); else drawPlot(); }, 200); });
    matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => { if (isOwn(state.d)) OWN[state.d](); else { drawPlot(); drawDetail(); } });
  }

  async function init() {
    readHash();
    if (state.d === "pheno" && state.v === "sex") state.s = "sexdiff"; // the scatter plot exists only for women vs. men
    $("#q").value = state.q; $("#sig").value = state.sig;
    try {
      [PMETA, RG] = await Promise.all([load("meta"), load("rg")]);
      RG.categories = PMETA.categories;
      // fields that share a description (e.g. fluid intelligence at the assessment centre and online) get a short qualifier
      [PMETA, RG].forEach((M) => (M.qual || []).forEach((q, i) => { if (q) M.desc[i] += `, ${q}`; }));
      RG.phewas.forEach((pi, ri) => { if (pi != null) P2RG[pi] = ri; });
      bind();
      if (isOwn(state.d)) { // open straight on heritability or loci; other views build when first chosen
        heritFrom = "pheno"; META = PMETA; syncControls();
        await OWN[state.d]();
        const open = state.d === "loci" && LS.sel != null ? $("#locus") : state.d === "mr" && mrSel != null ? $("#mr-detail") : state.d === "genes" && GS.sel ? $("#gene-detail") : null;
        if (open && !open.hidden) open.scrollIntoView({ block: "start" });
        return;
      }
      META = isGen(state.d) ? RG : PMETA;
      syncControls();
      await refresh(true);
      await drawDetail();
      if (state.t && !$("#trait").hidden) $("#trait").scrollIntoView({ block: "start" }); // a shared link to one trait
    } catch (err) {
      $("#plot-loading").textContent = `${err.message}. If you opened index.html directly from disk, serve the folder instead (for example: python3 -m http.server).`;
    }
  }
  init();
})();
