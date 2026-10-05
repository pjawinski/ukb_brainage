/* Brain age results browser — Jawinski et al. (2025) Nature Aging
   Static site: reads data/*.json produced by scripts/build_data.py
   Two result sets share one interface:
     pheno  phenome-wide associations (PHESANT) in all participants, women, men
     rg     genetic correlations with Neale lab UK Biobank GWAS (LDSC) */
/* eslint-disable n/no-unsupported-features/node-builtins -- browser code, not Node.js */
/* eslint-disable security/detect-object-injection -- lookups use keys from the site's own data files and
   from controls validated against fixed lists, not arbitrary user input */
const measures = {
  gwm: "Grey + white matter",
  gm: "Grey matter",
  wm: "White matter",
};

const modelKeys = ["gm", "wm", "gwm"];

const shortModel = { gwm: "Grey + white", gm: "Grey", wm: "White" };

const catShort = {
  "Hospital Inpatient - Administration": "Hospital admin.",
  "Family history and early life factors": "Family & early life",
  "Maternity and sex-specific factors": "Maternity & sex-spec.",
  "Medical history and conditions": "Medical history",
  "Lifestyle and environment": "Lifestyle",
  "Diet by 24-hour recall": "Diet (24-h recall)",
};

const sampleNames = { all: "All", female: "Women", male: "Men" };

// one colour per category, in the alphabetical order of meta.categories
const catColors = [
  "#7d6b5d",
  "#c2417b",
  "#1c8fb0",
  "#7a8b22",
  "#b07a2a",
  "#3a6fd8",
  "#1f9e89",
  "#9b4fd9",
  "#6b7a8f",
  "#d4a017",
  "#e06c00",
  "#2e8b57",
  "#d63a5e",
  "#00838f",
  "#5c9ce6",
  "#a0761c",
  "#7b5ea7",
  "#e17ca4",
];

const viewsByMode = {
  pheno: ["manhattan", "volcano", "sex"],
  rg: ["manhattan", "volcano"],
  cmp: ["pg"],
  herit: [],
  loci: [],
  genes: [],
  rgsel: [],
  mr: [],
};

const pageSize = 5;

// first free position for a label around its point, or null when every candidate box collides
const labelOffsets = [
  [22, -20],
  [-22, -20],
  [22, 20],
  [-22, 20],
  [28, -42],
  [-28, -42],
  [28, 42],
  [-28, 42],
  [0, -46],
  [0, 46],
];

const heritPlotConfig = {
  displaylogo: false,
  responsive: true,
  modeBarButtonsToRemove: [
    "select2d",
    "lasso2d",
    "autoScale2d",
    "zoom2d",
    "pan2d",
    "zoomIn2d",
    "zoomOut2d",
  ],
};

// ---------- GWAS loci (values and evidence strings as reported in the result tables) ----------
// chromosome lengths in GRCh37, used only to sort loci by genome position
const chrLen = [
  249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
  146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
  102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895,
  51304566, 155270560,
];

const chrStart = chrLen.map((_, k) =>
  chrLen.slice(0, k).reduce((a, b) => a + b, 0),
);

const lociShow = 10; // rows shown before "Show all" // m: a brain age model, or "all" for the three together

// column definitions as given in the paper's supplementary tables
const lociColumnHelp = {
  Locus:
    "Independent discovery count, each containing up to three co-inherited index variants derived from the three genome-wide association analyses of brain age gap.",
  Cytoband: "Cytogenetic band that contains the index variant.",
  Position:
    "Position of the index variant in base pairs according to human genome build hg19 (GRCh37).",
  Variant: "Identifier of the index variant.",
  "A1/A2":
    "A1 is the allele for which effects were calculated; A2 is the other allele.",
  "Freq.": "Frequency of the effect allele (A1).",
  "β (SE)":
    "Beta weight of the association between index variant and phenotype, with its standard error (years of brain age gap per A1 allele).",
  p: "P value of the association between index variant and phenotype.",
  "ηp²":
    "Partial eta squared; proportion of variance of the phenotype explained by the index variant (adjusted for sex, age, age², total intracranial volume, scanner site, type of array, and the first twenty genetic principal components).",
  N: "Sample size.",
  "Nearest gene":
    "HGNC symbol of the nearest gene based on ANNOVAR annotations (hg19 RefSeq gene table updated September 29, 2019), with the most relevant functional category of the index variant and the distance to the gene in base pairs. ANNOVAR prioritizes the most deleterious annotation for variants located in genomic regions where multiple genes overlap. Description and biotype follow the RefSeq gene annotation file in GFF3 format updated November 5, 2019.",
  "Credible set size":
    "Size of the 95% credible set of variants derived from applying SBayesRC, susieR and FINEMAP.",
  "SBayesRC genes":
    "Genes nominated by SBayesRC credible variant analysis. Brackets include the cumulative posterior probability of variants that have been annotated with the corresponding gene.",
  "Nonsynonymous variants":
    "Nonsynonymous exonic variants from the 95% credible set. Cells show genes whose transcripts are affected by the respective exonic variants. Brackets contain the number of identified nonsynonymous variants, the top nonsynonymous variant, and its CADD deleteriousness score.",
  "SMR eQTL":
    "Genes whose expression levels putatively mediate the effect of a locus variant on brain age gap, identified using summary-data-based Mendelian randomization (SMR) and the BrainMeta v2 eQTL dataset. Brackets contain the locus variant and the SMR raw p value.",
  "SMR sQTL":
    "Genes whose RNA splicing putatively mediates the effect of a locus variant on brain age gap, identified using SMR and the BrainMeta v2 sQTL dataset. Brackets contain the locus variant and the SMR raw p value.",
  "GTEx single tissue":
    "Regulated genes identified by mapping GWAS results to single-tissue expression quantitative trait loci of the Genotype-Tissue Expression (GTEx) database. Brackets include the number of tissues with a significant eQTL and the minimum p value across tissues.",
  "GTEx multi-tissue":
    "Regulated genes identified by mapping GWAS results to multi-tissue expression quantitative trait loci of the GTEx database. Brackets include the number of tissues where the eQTL had a posterior probability ≥ 0.9, and the Han and Eskin RE2 p value.",
  PoPS: "Genes implicated by the Polygenic Priority Score (PoPS) analysis. Brackets include the polygenic priority score.",
  "Prioritized gene":
    "Gene prioritized by aggregating the results of the seven gene nomination strategies.",
  "GWAS Catalog":
    "NHGRI-EBI GWAS Catalog results showing other complex traits previously associated with the index variant or any other genome-wide significant variant in strong linkage disequilibrium with it (r² > 0.8).",
  Literature:
    "Studies that previously identified the locus, with the reported variant with the strongest p value in brackets.",
};

// ---------- Mendelian randomization (GSMR and sensitivity methods, values as in the result table) ----------
const mrDirs = [
  ["to", "Trait → brain age gap"],
  ["from", "Brain age gap → trait"],
];

const searchGroups = [
  "GWAS loci",
  "Gene-based tests",
  "Phenotypic correlations",
  "Genetic correlations, UK Biobank traits",
  "Genetic correlations, 38 selected traits",
  "Mendelian randomization",
];

const modelText = {
  gm: "grey matter",
  wm: "white matter",
  gwm: "grey + white matter",
};

// UK Biobank imaging discovery sample (results/mri/accuracy.sample.txt); traits have up to this many participants
const phewasN = {
  all: "Up to 32,634 participants. ",
  female: "Up to 17,084 women. ",
  male: "Up to 15,550 men. ",
  sexdiff: "Up to 17,084 women and 15,550 men. ",
};

const sampleText = {
  all: "all participants",
  female: "women",
  male: "men",
  sexdiff: "women versus men",
};

const sectionOf = {
  herit: "#herit",
  loci: "#loci",
  genes: "#genes",
  rgsel: "#rgsel",
  mr: "#mr",
};

const modeOrder = [
  "herit",
  "loci",
  "genes",
  "pheno",
  "rg",
  "rgsel",
  "cmp",
  "mr",
];

const modeName = {
  herit: "Heritability",
  loci: "Genomic loci",
  genes: "Gene-based tests",
  pheno: "Phenotypic correlations",
  rg: "Genetic correlations, UK Biobank traits",
  rgsel: "Genetic correlations, 38 selected traits",
  cmp: "Phenotypic vs. genetic correlations",
  mr: "Mendelian randomization",
};

const viewNotes = {
  pheno: {
    manhattan:
      "Traits grouped by category, as in the paper figure. Height shows significance; triangles point in the direction of the association.",
    volcano:
      "Effect size against significance. Points far left or right are the strongest associations; the lines mark FDR and Bonferroni thresholds.",
    sex: "Each trait’s association in women against men. Coloured points are significant in at least one sex (FDR < 5%); outlined points differ between the sexes.",
  },
  rg: {
    manhattan:
      "Genetic correlations grouped by category. Height shows significance; triangles point in the direction of the genetic correlation.",
    volcano:
      "Genetic correlation against significance. Points far left or right share the most genetic signal with brain age gap.",
  },
  cmp: {
    pg: "Genetic correlation against phenotypic correlation (all participants), as in the paper. Trait pairs are matched by description; multinomial models are excluded and only traits with SNP heritability h²/SE > 1.96 are kept. Coloured points are significant in one analysis, outlined points in both.",
  },
};

const bibtex = `@article{jawinski2025brainage,
  author  = {Jawinski, Philippe and Forstbach, Helena and Kirsten, Holger and Beyer, Frauke and Villringer, Arno and Witte, A. Veronica and Scholz, Markus and Ripke, Stephan and Markett, Sebastian},
  title   = {Genome-wide analysis of brain age identifies 59 associated loci and unveils relationships with mental and physical health},
  journal = {Nature Aging},
  volume  = {5},
  pages   = {2086--2103},
  year    = {2025},
  doi     = {10.1038/s43587-025-00962-7}
}`;

(function () {
  "use strict";

  const ownPanels = {
    herit: () => drawHerit(),
    loci: () => drawLoci(),
    rgsel: () => drawRgSel(),
    mr: () => drawMr(),
    genes: () => drawGenes(),
  }; // result sets with their own panels instead of plot + table
  const isOwn = (d) => d in ownPanels;
  const isGen = (d) => d === "rg" || d === "cmp"; // both use the genetic-correlation trait list
  const isNil = (x) => x === null || x === undefined;
  const qsel = (s) => document.querySelector(s);
  const lt05 = (x) => !isNil(x) && x < 0.05;

  const state = {
    d: "pheno",
    v: "manhattan",
    m: "gwm",
    s: "all",
    t: null,
    q: "",
    cat: null,
    sig: "all",
    dir: "any",
    sort: "p",
    asc: true,
    page: 0,
    labels: 10,
  };
  const cache = {};
  let pmeta = null,
    rgMeta = null,
    meta = null;
  const p2rg = {}; // PheWAS trait index -> rg trait index
  let view = null;

  // ---------- helpers ----------
  const fmtInt = (n) => (isNil(n) ? "–" : n.toLocaleString("en-US"));
  /**
   * Replace the content of an element with markup. All markup is built in this
   * file from the site's own data files; free text such as trait names or the
   * search input goes through esc() first. The markup is parsed in the context
   * of the element, so table rows stay table rows.
   */
  function setHtml(el, html) {
    const range = document.createRange();
    range.selectNodeContents(el);
    el.replaceChildren(range.createContextualFragment(html));
  }
  /** Format a p-value for display, as HTML (×10 superscript) or plain text. */
  function fmtP(p, html = true) {
    if (isNil(p)) return "–";
    if (p >= 0.001) return p.toPrecision(2);
    const [mant, ex] = p.toExponential(1).split("e");
    const e = String(+ex).replace("-", "−");
    return html ? `${mant}×10<sup>${e}</sup>` : `${mant}e${+ex}`;
  }
  const fmtR = (r) =>
    isNil(r) ? "–" : (r < 0 ? "−" : "") + Math.abs(r).toFixed(3);
  const fmtB = (b, se) =>
    isNil(b)
      ? "–"
      : `${b < 0 ? "−" : ""}${Math.abs(b).toPrecision(3)} (${isNil(se) ? "–" : se.toPrecision(2)})`;
  const fmtF = (x, d = 3) => (isNil(x) ? "–" : x.toFixed(d));
  const esc = (s) =>
    String(s).replace(
      /[&<>"]/g,
      (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c],
    );
  const css = (v) =>
    getComputedStyle(document.documentElement).getPropertyValue(v).trim();
  /** Convert a #rrggbb or #rgb colour to rgba() with the given opacity. */
  function hexA(hex, a) {
    let h = hex.slice(1);
    if (h.length === 3) h = h.replace(/./g, (c) => c + c);
    const n = parseInt(h, 16);
    return `rgba(${(n >> 16) & 255},${(n >> 8) & 255},${n & 255},${a})`;
  }
  const trunc = (s, n) => (s.length > n ? s.slice(0, n - 1) + "…" : s);
  const surfaceA = () =>
    hexA(css("--surface").startsWith("#") ? css("--surface") : "#ffffff", 0.85);
  const reduced = () => matchMedia("(prefers-reduced-motion: reduce)").matches;

  /** Benjamini–Hochberg on an array (nulls ignored). */
  function bh(ps) {
    const idx = ps
      .map((p, i) => [p, i])
      .filter((x) => !isNil(x[0]))
      .sort((a, b) => a[0] - b[0]);
    const m = idx.length,
      q = new Array(ps.length).fill(null);
    let prev = 1;
    for (let k = m - 1; k >= 0; k--) {
      prev = Math.min(prev, (idx[k][0] * m) / (k + 1));
      q[idx[k][1]] = prev;
    }
    return q;
  }

  /** Load a data file once and cache the promise. */
  async function load(name) {
    // an offline copy can ship the data as a script (window.BAG_DATA), since browsers block fetch() from file://
    if (!cache[name] && window.BAG_DATA && window.BAG_DATA[name])
      cache[name] = Promise.resolve(window.BAG_DATA[name]);
    if (!cache[name])
      cache[name] = window.fetch(`data/${name}.json`).then((r) => {
        if (!r.ok)
          throw new Error(`Could not load data/${name}.json (${r.status})`);
        return r.json();
      });
    return cache[name];
  }

  // ---------- URL state ----------
  /** Read the view state from the URL hash. */
  function readHash() {
    // eslint-disable-next-line compat/compat -- supported by every browser the site targets
    const h = new URLSearchParams(location.hash.slice(1));
    if (viewsByMode[h.get("d")]) state.d = h.get("d");
    else if (!location.hash.slice(1)) state.d = "herit"; // the page opens on heritability; older links without d= stay on the phenotypic correlations
    if (viewsByMode[state.d].includes(h.get("v"))) state.v = h.get("v");
    if (state.d === "cmp") state.v = "pg";
    if (measures[h.get("m")]) state.m = h.get("m");
    if (state.d === "loci")
      lociState.m = measures[h.get("m")] ? h.get("m") : "all";
    if (state.d === "loci" && h.get("locus")) lociState.sel = +h.get("locus"); // an open locus panel
    if (state.d === "mr" && h.get("trait")) mrSel = +h.get("trait"); // an open MR trait panel
    if (state.d === "genes" && h.get("gene")) geneState.sel = h.get("gene");
    if (sampleNames[h.get("s")] || h.get("s") === "sexdiff")
      state.s = h.get("s");
    if (h.get("t")) state.t = h.get("t");
    if (h.get("q")) state.q = h.get("q");
    if (h.get("sig")) state.sig = h.get("sig");
  }
  /** Write the current view state to the URL hash. */
  function writeHash() {
    if (isOwn(state.d)) {
      // eslint-disable-next-line compat/compat -- supported by every browser the site targets
      const o = new URLSearchParams({ d: state.d });
      if (state.d === "loci" && lociState.m !== "all") o.set("m", lociState.m);
      if (state.d === "loci" && !isNil(lociState.sel))
        o.set("locus", lociState.sel);
      if (state.d === "mr" && !isNil(mrSel)) o.set("trait", mrSel);
      if (state.d === "genes") {
        o.set("m", state.m);
        if (geneState.sel) o.set("gene", geneState.sel);
      }
      history.replaceState(null, "", "#" + o.toString());
      return;
    }
    // eslint-disable-next-line compat/compat -- supported by every browser the site targets
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
  /**
   * Build the rows, category order and counts for the current result set, model
   * and sample.
   */
  async function buildView() {
    const m = state.m,
      rows = [];
    const counts = { gm: {}, wm: {}, gwm: {} };
    const tally = (mm, qs, catOf) =>
      qs.forEach((q, i) => {
        if (lt05(q)) {
          const c = catOf(i);
          counts[mm][c] = (counts[mm][c] || 0) + 1;
        }
      });

    if (state.d === "cmp") {
      // trait pairs selected as in code/genetics/rgVSrp.R; "significant" here means FDR < 5% in both analyses
      const all = await load("all");
      const pairRows = (mm) =>
        rgMeta.pairs.ri
          .map((ri, k) => {
            const pi = rgMeta.pairs.pi[k];
            if (isNil(rgMeta[mm].p[ri]) || isNil(all[mm].r[pi])) return null;
            return {
              i: ri,
              pidx: pi,
              p: rgMeta[mm].p[ri],
              q: rgMeta[mm].q[ri],
              r: rgMeta[mm].r[ri],
              se: rgMeta[mm].se[ri],
              n: rgMeta.h2[ri],
              pr: all[mm].r[pi],
              pp: all[mm].p[pi],
              pq: all[mm].q[pi],
              b: all[mm].r[pi],
            };
          })
          .filter(Boolean);
      rows.push(...pairRows(m));
      modelKeys.forEach((mm) =>
        pairRows(mm).forEach((r) => {
          if (lt05(r.q) && lt05(r.pq)) {
            const c = rgMeta.cat[r.i];
            counts[mm][c] = (counts[mm][c] || 0) + 1;
          }
        }),
      );
      var pairs = rows;
    } else if (state.d === "rg") {
      const all = await load("all");
      const st = rgMeta[m];
      for (let i = 0; i < rgMeta.id.length; i++) {
        if (isNil(st.p[i])) continue;
        const pi = rgMeta.phewas[i];
        rows.push({
          i,
          p: st.p[i],
          q: st.q[i],
          r: st.r[i],
          se: st.se[i],
          b: st.se[i],
          n: rgMeta.h2[i],
          pr: !isNil(pi) ? all[m].r[pi] : null,
          pq: !isNil(pi) ? all[m].q[pi] : null,
        });
      }
      modelKeys.forEach((mm) => tally(mm, rgMeta[mm].q, (i) => rgMeta.cat[i]));
    } else if (state.s === "sexdiff") {
      const [sd, f, ma] = await Promise.all([
        load("sexdiff"),
        load("female"),
        load("male"),
      ]);
      const p = sd[m],
        q = bh(p);
      for (let i = 0; i < pmeta.id.length; i++) {
        if (isNil(p[i])) continue;
        const rf = f[m].r[i],
          rm = ma[m].r[i];
        rows.push({
          i,
          p: p[i],
          q: q[i],
          r: !isNil(rf) && !isNil(rm) ? rf - rm : null,
          rf,
          rm,
          qf: f[m].q[i],
          qm: ma[m].q[i],
          b: null,
          se: null,
          n: (f.ntotal[i] || 0) + (ma.ntotal[i] || 0),
        });
      }
      modelKeys.forEach((mm) => tally(mm, bh(sd[mm]), (i) => pmeta.cat[i]));
    } else {
      const d = await load(state.s);
      const st = d[m];
      for (let i = 0; i < pmeta.id.length; i++) {
        if (isNil(st.p[i])) continue;
        rows.push({
          i,
          p: st.p[i],
          q: st.q[i],
          r: st.r[i],
          b: st.b[i],
          se: st.se[i],
          n: d.ntotal[i],
        });
      }
      modelKeys.forEach((mm) => tally(mm, d[mm].q, (i) => pmeta.cat[i]));
    }
    // Manhattan x positions: categories present in this set, alphabetical, ppoints within category
    const byCat = {};
    rows.forEach((r) =>
      (byCat[meta.cat[r.i]] = byCat[meta.cat[r.i]] || []).push(r),
    );
    const cats = Object.keys(byCat)
      .map(Number)
      .sort((a, b) => a - b);
    cats.forEach((c, k) =>
      byCat[c].forEach((r, j) => (r.x = k + (j + 0.5) / byCat[c].length)),
    );
    const sigP = rows.filter((r) => lt05(r.q)).map((r) => r.p);
    view = {
      rows,
      byCat,
      cats,
      counts,
      pairs: typeof pairs !== "undefined" ? pairs : [],
      bonf: 0.05 / (state.d === "cmp" ? rgMeta.id.length : rows.length),
      pbonf: 0.05 / pmeta.id.length,
      fdr: sigP.length ? Math.max(...sigP) : null,
    };
  }

  /**
   * Whether a row passes the current category, search, significance and
   * direction filters.
   */
  function passes(r) {
    if (!isNil(state.cat) && meta.cat[r.i] !== state.cat) return false;
    if (state.sig === "fdr" && !isSig(r)) return false;
    if (
      state.sig === "bonf" &&
      !(r.p < view.bonf && (state.d !== "cmp" || r.pp < view.pbonf))
    )
      return false;
    if (state.dir === "pos" && !(r.r > 0)) return false;
    if (state.dir === "neg" && !(r.r < 0)) return false;
    if (state.q) {
      const t = state.q.toLowerCase(),
        f = meta.field[r.i];
      if (
        !meta.desc[r.i].toLowerCase().includes(t) &&
        !meta.id[r.i].toLowerCase().includes(t) &&
        !(!isNil(f) && String(f).includes(t))
      )
        return false;
    }
    return true;
  }

  // ---------- main plot ----------
  let plotReady = false;
  const yl = (p) => -Math.log10(p);
  const isSig = (r) => lt05(r.q) && (state.d !== "cmp" || lt05(r.pq));
  /** X and y position of a row in the current plot type. */
  function coords(r) {
    if (state.v === "volcano") return [r.r, yl(r.p)];
    if (state.v === "sex") return [r.rf, r.rm];
    if (state.v === "pg") return [r.r, r.pr];
    return [r.x, yl(r.p)];
  }
  /** Rows that can be drawn in the current plot type. */
  function plotRows() {
    if (state.v === "sex")
      return view.rows.filter((r) => !isNil(r.rf) && !isNil(r.rm));
    if (state.v === "pg") return view.pairs;
    if (state.v === "volcano") return view.rows.filter((r) => !isNil(r.r));
    return view.rows;
  }
  /**
   * Scatter views: colour = significant in at least one of the two compared
   * analyses, outline = the highlighted contrast.
   */
  function scatterClass(r) {
    if (state.v === "sex")
      return { colour: lt05(r.qf) || lt05(r.qm), ring: isSig(r) };
    return { colour: lt05(r.q) || lt05(r.pq), ring: lt05(r.q) && lt05(r.pq) };
  }

  /**
   * Axes and reference shapes for the two scatter views (women vs men,
   * phenotypic vs genetic).
   */
  function scatterAxes(c) {
    const { v, rows, grid, shapes, annotations, ink, muted, narrow } = c;
    const xs = rows.map((r) => coords(r)[0]),
      ys = rows.map((r) => coords(r)[1]);
    let x0, x1, y0, y1;
    if (v === "sex") {
      const lim =
        Math.max(0.02, ...xs.map(Math.abs), ...ys.map(Math.abs)) * 1.12;
      [x0, x1, y0, y1] = [-lim, lim, -lim, lim];
      shapes.push({
        type: "line",
        x0,
        x1,
        y0,
        y1,
        line: { color: muted, width: 1, dash: "dash" },
        layer: "below",
      });
    } else {
      const lx = Math.max(0.1, ...xs.map(Math.abs)) * 1.08,
        ly = Math.max(0.02, ...ys.map(Math.abs)) * 1.12;
      [x0, x1, y0, y1] = [-lx, lx, -ly, ly];
      const ps = rgMeta.pairStats[state.m]; // least-squares line, as geom_smooth(method = "lm") in the paper figure
      shapes.push({
        type: "line",
        x0,
        x1,
        y0: ps.intercept + ps.slope * x0,
        y1: ps.intercept + ps.slope * x1,
        line: { color: ink, width: 1.4 },
      });
      annotations.push({
        xref: "paper",
        yref: "paper",
        x: 0.5,
        y: 1,
        yanchor: "bottom",
        showarrow: false,
        text: `r = ${ps.r.toFixed(2)}  |  p = ${fmtP(ps.p)}  |  MAD = ${ps.mad.toFixed(3)}  |  ${ps.n} trait pairs`,
        font: {
          size: narrow ? 11.5 : 13,
          color: ink,
          family: "Source Sans 3, sans-serif",
        },
      });
    }
    const xt = v === "sex" ? "r in women" : "Genetic correlation r<sub>g</sub>";
    const yt = v === "sex" ? "r in men" : "Phenotypic correlation r";
    return {
      x0,
      x1,
      y0,
      y1,
      xaxis: {
        ...grid,
        title: { text: xt, standoff: 6 },
        range: [x0, x1],
        zeroline: true,
        zerolinecolor: muted,
      },
      yaxis: {
        ...grid,
        title: { text: yt, standoff: 6 },
        range: [y0, y1],
        zeroline: true,
        zerolinecolor: muted,
      },
    };
  }

  /**
   * X-axis of the Manhattan plot: one column per category, or a single zoomed-
   * in category.
   */
  function manhattanXaxis(c) {
    const { shapes, rule, band, narrow, zoomCat } = c;
    const cats = view.cats;
    if (zoomCat) {
      const k = cats.indexOf(state.cat);
      return {
        x0: k,
        x1: k + 1,
        xaxis: {
          range: [k, k + 1],
          fixedrange: true,
          showgrid: false,
          zeroline: false,
          tickvals: [k + 0.5],
          ticktext: [
            `${meta.categories[state.cat]}: all ${fmtInt(view.byCat[state.cat].length)} traits`,
          ],
          tickangle: 0,
          ticks: "",
          linecolor: rule,
        },
      };
    }
    cats.forEach((cat, k) => {
      if (k % 2 === 1)
        shapes.push({
          type: "rect",
          xref: "x",
          yref: "paper",
          x0: k,
          x1: k + 1,
          y0: 0,
          y1: 1,
          fillcolor: band,
          line: { width: 0 },
          layer: "below",
        });
    });
    const labels = narrow
      ? cats.map(() => "")
      : cats.map(
          (cat) => catShort[meta.categories[cat]] || meta.categories[cat],
        );
    return {
      x0: 0,
      x1: cats.length,
      xaxis: {
        range: [0, cats.length],
        fixedrange: true,
        showgrid: false,
        zeroline: false,
        tickvals: cats.map((_, k) => k + 0.5),
        ticktext: labels,
        tickangle: 40,
        ticks: "",
        linecolor: rule,
      },
    };
  }

  /**
   * Axes and threshold lines for the views with -log10(p) on the y-axis
   * (Manhattan, volcano).
   */
  function pAxes(c) {
    const { v, rows, grid, shapes, ink, muted, zoomCat } = c;
    const sexdiff = state.d === "pheno" && state.s === "sexdiff";
    // Manhattan zooms into a selected category: its points spread over the full width and the y-axis fits them
    const yRows = zoomCat
      ? rows.filter((r) => meta.cat[r.i] === state.cat)
      : rows;
    const y0 = 0,
      y1 =
        Math.max(
          zoomCat ? yl(view.bonf) * 1.15 : 8,
          ...yRows.map((r) => yl(r.p)),
        ) * (state.labels > 0 ? 1.2 : 1.06);
    shapes.push({
      type: "line",
      xref: "paper",
      x0: 0,
      x1: 1,
      y0: yl(view.bonf),
      y1: yl(view.bonf),
      line: { color: ink, width: 1.1 },
    });
    if (view.fdr)
      shapes.push({
        type: "line",
        xref: "paper",
        x0: 0,
        x1: 1,
        y0: yl(view.fdr),
        y1: yl(view.fdr),
        line: { color: ink, width: 1.1, dash: "dash" },
      });
    const yt = sexdiff
      ? "−log<sub>10</sub>(p) for sex difference"
      : "−log<sub>10</sub>(p)";
    const yaxis = {
      ...grid,
      title: { text: yt, standoff: 6 },
      range: [y0, y1],
      zeroline: false,
    };
    if (v !== "volcano") return { ...manhattanXaxis(c), y0, y1, yaxis };
    const lim = Math.max(0.02, ...rows.map((r) => Math.abs(r.r))) * 1.12;
    let xt = "r (correlation-scale effect size)";
    if (isGen(state.d)) xt = "Genetic correlation r<sub>g</sub>";
    else if (sexdiff) xt = "Δr (women − men)";
    return {
      x0: -lim,
      x1: lim,
      y0,
      y1,
      yaxis,
      xaxis: {
        ...grid,
        title: { text: xt, standoff: 6 },
        range: [-lim, lim],
        zeroline: true,
        zerolinecolor: muted,
      },
    };
  }

  /**
   * First free position for a label around its point, or null when every
   * candidate box collides.
   */
  function labelSpot(ax0, ay0, w, lh, W, H, boxes) {
    const overlap = (a, b) =>
      a.x0 < b.x1 && b.x0 < a.x1 && a.y0 < b.y1 && b.y0 < a.y1;
    for (const [ox, oy] of labelOffsets) {
      const tx = ax0 + ox,
        ty = ay0 + oy;
      let anchor = "center",
        bx0 = tx - w / 2;
      if (ox > 0) [anchor, bx0] = ["left", tx];
      else if (ox < 0) [anchor, bx0] = ["right", tx - w];
      const box = { x0: bx0, x1: bx0 + w, y0: ty - lh / 2, y1: ty + lh / 2 };
      if (box.x0 < 0 || box.x1 > W || box.y0 < 0 || box.y1 > H) continue;
      if (boxes.some((b) => overlap(b, box))) continue;
      return { ox, oy, anchor, box };
    }
    return null;
  }

  /**
   * Label placement in pixel space: try several positions around each point and
   * keep a label only where its box stays inside the plot and does not overlap
   * labels already placed.
   */
  function placeLabels(c, ax, margin) {
    const { rows, annotations, ink, muted, narrow } = c;
    const { x0, x1, y0, y1 } = ax;
    const el = qsel("#plot");
    const plotW = Math.max(100, el.clientWidth - margin.l - margin.r),
      plotH = Math.max(100, el.clientHeight - margin.t - margin.b);
    const fs = narrow ? 10.5 : 12,
      maxChars = narrow ? 26 : 40,
      lh = fs * 1.35 + 4;
    const sel = !isNil(state.t)
      ? rows.find((r) => meta.id[r.i] === state.t)
      : null;
    const cand = rows
      .filter((r) => r !== sel && passes(r) && isSig(r))
      .sort((a, b) => a.p - b.p);
    const queue = (sel ? [sel] : []).concat(cand);
    const boxes = [],
      seenTrait = new Set();
    const limit = state.labels + (sel ? 1 : 0);
    for (const r of queue) {
      if (seenTrait.size >= limit) break;
      const [x, y] = coords(r);
      if (seenTrait.has(r.i) || isNil(x) || isNil(y)) continue;
      const text = trunc(meta.desc[r.i], maxChars),
        w = text.length * fs * 0.53 + 8;
      const ax0 = ((x - x0) / (x1 - x0)) * plotW,
        ay0 = (1 - (y - y0) / (y1 - y0)) * plotH;
      const hit = labelSpot(ax0, ay0, w, lh, plotW, plotH, boxes);
      if (!hit) continue;
      boxes.push(hit.box);
      boxes.push({ x0: ax0 - 5, x1: ax0 + 5, y0: ay0 - 5, y1: ay0 + 5 }); // keep later labels off this point
      seenTrait.add(r.i);
      annotations.push({
        x,
        y,
        text: r === sel ? `<b>${esc(text)}</b>` : esc(text),
        showarrow: true,
        arrowhead: 0,
        arrowwidth: 0.8,
        arrowcolor: muted,
        ax: hit.ox,
        ay: hit.oy,
        xanchor: hit.anchor,
        font: { size: fs, color: ink, family: "Source Sans 3, sans-serif" },
        bgcolor: surfaceA(),
      });
    }
  }

  /** Plot margins for the current plot type and screen width. */
  function plotMargin(v, zoomCat, narrow) {
    let bottom = 48;
    if (v === "manhattan") bottom = zoomCat ? 36 : narrow ? 16 : 96;
    return { l: 56, r: 12, t: v === "pg" ? 34 : 16, b: bottom };
  }

  /** Plotly layout of the main plot: axes, threshold lines and labels. */
  function plotLayout() {
    const narrow = window.innerWidth < 900;
    const ink = css("--ink"),
      muted = css("--muted"),
      rule = css("--rule"),
      band = css("--band");
    const v = state.v;
    const c = {
      v,
      rows: plotRows(),
      shapes: [],
      annotations: [],
      ink,
      muted,
      rule,
      band,
      narrow,
      zoomCat:
        v === "manhattan" && !isNil(state.cat) && view.cats.includes(state.cat),
      grid: {
        gridcolor: rule,
        griddash: "dot",
        linecolor: rule,
        fixedrange: narrow,
      },
    };
    const ax = v === "sex" || v === "pg" ? scatterAxes(c) : pAxes(c);
    const margin = plotMargin(v, c.zoomCat, narrow);
    if (state.labels > 0 || !isNil(state.t)) placeLabels(c, ax, margin);
    return {
      margin,
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)",
      font: { family: "Source Sans 3, sans-serif", color: muted, size: 12.5 },
      xaxis: ax.xaxis,
      yaxis: ax.yaxis,
      shapes: c.shapes,
      annotations: c.annotations,
      showlegend: false,
      dragmode: narrow ? false : "zoom",
      hoverlabel: {
        bgcolor: css("--surface"),
        bordercolor: rule,
        font: { color: ink, family: "Source Sans 3, sans-serif", size: 13 },
        align: "left",
      },
    };
  }

  /** Hover text of one point in the main plot. */
  function hoverText(r) {
    const head =
      `<b>${esc(trunc(meta.desc[r.i], 70))}</b><br>${esc(meta.categories[meta.cat[r.i]])}` +
      (!isNil(meta.field[r.i]) ? ` · field ${meta.field[r.i]}` : "") +
      "<br>";
    if (isGen(state.d)) {
      const other = !isNil(r.pidx) ? ` (${esc(pmeta.id[r.pidx])})` : "";
      return (
        head +
        `r<sub>g</sub> = ${fmtR(r.r)} (SE ${fmtF(r.se)}) · p = ${fmtP(r.p)} · FDR = ${fmtP(r.q)}` +
        (!isNil(r.pr)
          ? `<br>phenotypic r${other} = ${fmtR(r.pr)} · FDR = ${fmtP(r.pq)}`
          : "")
      );
    }
    if (state.s === "sexdiff")
      return (
        head +
        `r women = ${fmtR(r.rf)} · r men = ${fmtR(r.rm)}<br>p (difference) = ${fmtP(r.p)} · FDR = ${fmtP(r.q)}`
      );
    return (
      head +
      (!isNil(r.r) ? `r = ${fmtR(r.r)} · ` : "") +
      `p = ${fmtP(r.p)}` +
      (!isNil(r.q) ? ` · FDR = ${fmtP(r.q)}` : "")
    );
  }

  /** Size and colour of a point in the two scatter views. */
  function scatterPointStyle(r, col, grey, filtering, on) {
    const k = scatterClass(r);
    let size = 4;
    if (k.ring) size = 11;
    else if (k.colour) size = 7;
    else if (filtering && on) size = 5.5;
    let color = hexA(grey, 0.28);
    if (!on) color = hexA(grey, 0.12);
    else if (k.colour || k.ring) color = col;
    else if (filtering) color = hexA(col, 0.55);
    return { size, color, lw: k.ring && on ? 1.6 : 0 };
  }
  /**
   * Size and colour of a point in the Manhattan and volcano views; inside an
   * active filter every point keeps its category colour, only points outside
   * the filter turn grey.
   */
  function pPointStyle(sig, col, c, on) {
    const man = c.v === "manhattan";
    let size = man ? 5 : 4.5;
    if (sig) size = man ? 9 : 8;
    else if (c.filtering && on) size = 6;
    let color = hexA(c.grey, 0.16);
    if (on && sig) color = col;
    else if (on) color = hexA(col, c.filtering ? 0.6 : man ? 0.38 : 0.3);
    return { size, color, lw: 0 };
  }
  /**
   * Marker symbol, size, colour and outline of one point in the correlation
   * plots.
   */
  function pointStyle(r, col, c) {
    const sig = isSig(r),
      on = !c.filtering || passes(r);
    const style = c.scatter
      ? scatterPointStyle(r, col, c.grey, c.filtering, on)
      : pPointStyle(sig, col, c, on);
    style.sym = "circle";
    if (c.dirMarks && sig)
      style.sym = r.r >= 0 ? "triangle-up" : "triangle-down";
    return style;
  }

  /** Plotly traces of the main plot, one per category. */
  function plotTraces() {
    const traces = [];
    const v = state.v,
      rows = plotRows();
    const scatter = v === "sex" || v === "pg";
    const dirMarks =
      v === "manhattan" && !(state.d === "pheno" && state.s === "sexdiff");
    const grey = css("--faint"),
      ink = css("--ink");
    const filtering =
      !isNil(state.cat) ||
      state.q ||
      state.sig !== "all" ||
      state.dir !== "any";
    const c0 = { v, scatter, dirMarks, grey, filtering };
    const byCat = {};
    rows.forEach((r) =>
      (byCat[meta.cat[r.i]] = byCat[meta.cat[r.i]] || []).push(r),
    );
    // the selected category is drawn last so its points sit on top of the greyed-out rest
    Object.entries(byCat)
      .sort((a, b) => (+a[0] === state.cat) - (+b[0] === state.cat))
      .forEach(([c, rs]) => {
        const col = catColors[+c % catColors.length];
        const x = [],
          y = [],
          sym = [],
          size = [],
          color = [],
          cd = [],
          text = [],
          lw = [];
        const rank = (r) =>
          scatter
            ? (scatterClass(r).ring ? 2 : 0) + (scatterClass(r).colour ? 1 : 0)
            : isSig(r)
              ? 1
              : 0;
        rs.slice()
          .sort((a, b) => rank(a) - rank(b))
          .forEach((r) => {
            // significant points drawn last, on top
            const [px, py] = coords(r);
            x.push(px);
            y.push(py);
            const st = pointStyle(r, col, c0);
            sym.push(st.sym);
            size.push(st.size);
            color.push(st.color);
            lw.push(st.lw);
            cd.push(r.i);
            text.push(hoverText(r));
          });
        traces.push({
          type: "scattergl",
          mode: "markers",
          x,
          y,
          customdata: cd,
          text,
          hovertemplate: "%{text}<extra></extra>",
          marker: { symbol: sym, size, color, line: { width: lw, color: ink } },
        });
      });
    const sel = !isNil(state.t)
      ? rows.find((r) => meta.id[r.i] === state.t)
      : null;
    const sc = sel ? coords(sel) : null;
    traces.push({
      type: "scatter",
      mode: "markers",
      x: sc ? [sc[0]] : [],
      y: sc ? [sc[1]] : [],
      hoverinfo: "skip",
      marker: {
        symbol: "circle-open",
        size: 18,
        color: ink,
        line: { width: 2 },
      },
    });
    return traces;
  }

  /** Draw or update the main plot. */
  function drawPlot() {
    if (!view || isOwn(state.d)) return; // nothing to draw until a correlation view has been built
    const el = qsel("#plot");
    const cfg = {
      responsive: true,
      displaylogo: false,
      modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"],
      toImageButtonOptions: {
        filename: `brainage_${state.d}_${state.v}_${state.m}${state.d === "pheno" ? "_" + state.s : ""}`,
        scale: 3,
      },
    };
    Plotly.react(el, plotTraces(), plotLayout(), cfg).then(() => {
      qsel("#plot-loading").hidden = true;
      if (!plotReady) {
        plotReady = true;
        el.on("plotly_click", (ev) => {
          const pt = ev.points && ev.points[0];
          if (pt && !isNil(pt.customdata)) select(meta.id[pt.customdata], true);
        });
      }
    });
  }

  // ---------- legend under the plot ----------
  /** Legend under the main plot for the current plot type. */
  function keyHtml() {
    const tri = (d) =>
      `<svg width="12" height="12" viewBox="0 0 12 12"><path d="${d === "up" ? "M6 1 11 11H1z" : "M6 11 11 1H1z"}"/></svg>`;
    const line = (dash) =>
      `<svg width="22" height="12" viewBox="0 0 22 12"><line x1="0" y1="6" x2="22" y2="6"${dash ? ' stroke-dasharray="4 3"' : ""}/></svg>`;
    const dot = (ring) =>
      `<svg width="14" height="14" viewBox="0 0 14 14"><circle cx="7" cy="7" r="5" class="${ring ? "ring" : "fill"}"/></svg>`;
    const thresholds = `<span>${line(true)} FDR 5%</span><span>${line(false)} Bonferroni 5%</span>`;
    if (state.v === "sex")
      return (
        `<span><svg width="16" height="16" viewBox="0 0 16 16"><line x1="1" y1="15" x2="15" y2="1" stroke-dasharray="3 2"/></svg> equal effect in women and men</span>` +
        `<span>${dot(false)} significant in women or men</span><span>${dot(true)} women and men differ (FDR 5%)</span>`
      );
    if (state.v === "pg")
      return `<span>${dot(false)} significant in one analysis (FDR 5%)</span><span>${dot(true)} significant in both</span>`;
    if (state.v === "volcano" || (state.d === "pheno" && state.s === "sexdiff"))
      return thresholds;
    const up = isGen(state.d)
      ? "positive genetic correlation"
      : "older-appearing brain with higher trait value";
    const down = isGen(state.d)
      ? "negative genetic correlation"
      : "older-appearing brain with lower trait value";
    return (
      `<span>${tri("up")} ${up}</span><span>${tri("down")} ${down}</span>` +
      thresholds
    );
  }

  // ---------- category chips ----------
  /** Draw the category chips above the table. */
  function drawCats() {
    if (!view) return;
    const counts = view.counts[state.m];
    const box = qsel("#cats");
    box.classList.toggle("filtering", !isNil(state.cat));
    setHtml(
      box,
      view.cats
        .map((k) => {
          const c = meta.categories[k];
          return (
            `<button type="button" class="cat" style="--c:${catColors[k]}" data-k="${k}" aria-pressed="${state.cat === k}" title="${counts[k] || 0} traits at FDR < 5%">` +
            `<i></i>${esc(c)}${counts[k] ? ` <b>${counts[k]}</b>` : ""}</button>`
          );
        })
        .join(""),
    );
    box
      .querySelectorAll(".cat")
      .forEach((b) => b.addEventListener("click", () => setCat(+b.dataset.k)));
  }
  /** Select a category, or clear it when it is already selected. */
  function setCat(k) {
    state.cat = state.cat === k ? null : k;
    state.page = 0;
    refresh(false);
  }

  // ---------- table ----------
  /** Table rows after filtering and sorting. */
  function sortedFiltered() {
    const rows = view.rows.filter(passes);
    const k = state.sort,
      dir = state.asc ? 1 : -1;
    const cmp = state.d === "cmp";
    const val = (r) =>
      k === "desc"
        ? meta.desc[r.i].toLowerCase()
        : k === "cat"
          ? meta.categories[meta.cat[r.i]]
          : k === "r"
            ? isNil(r.r)
              ? null
              : Math.abs(r.r)
            : cmp && k === "b"
              ? Math.abs(r.pr)
              : cmp && k === "p"
                ? r.q
                : cmp && k === "q"
                  ? r.pq
                  : r[k];
    rows.sort((a, b) => {
      const va = val(a),
        vb = val(b);
      if (isNil(va)) return 1;
      if (isNil(vb)) return -1;
      return (va < vb ? -1 : va > vb ? 1 : 0) * dir || a.p - b.p;
    });
    return rows;
  }

  /** Draw the paged results table. */
  function drawTable() {
    if (!view) return;
    const rg = isGen(state.d),
      cmp = state.d === "cmp";
    const rows = sortedFiltered();
    qsel("#table").classList.toggle("cmp", cmp);
    const pages = Math.max(1, Math.ceil(rows.length / pageSize));
    state.page = Math.min(state.page, pages - 1);
    const slice = rows.slice(
      state.page * pageSize,
      state.page * pageSize + pageSize,
    );
    const tb = qsel("#table tbody");
    setHtml(
      tb,
      slice.length
        ? slice
            .map((r) => {
              const id = meta.id[r.i],
                c = meta.cat[r.i];
              const own = rg
                ? !isNil(meta.field[r.i])
                  ? `field ${meta.field[r.i]}`
                  : ""
                : id;
              // a pair can join two different fields with the same name; then name the phenotypic one too
              const other =
                !isNil(r.pidx) && pmeta.field[r.pidx] !== meta.field[r.i]
                  ? ` · phenotypic r: field ${pmeta.field[r.pidx]}` +
                    (pmeta.qual[r.pidx] ? `, ${pmeta.qual[r.pidx]}` : "")
                  : "";
              const sub = own + other;
              return (
                `<tr data-id="${esc(id)}" class="${id === state.t ? "sel" : ""}" tabindex="0">` +
                `<td class="trait">${esc(meta.desc[r.i])}<small>${esc(sub)}</small></td>` +
                `<td class="hide-s catcell"><span class="dot" style="background:${catColors[c]}"></span>${esc(meta.categories[c])}</td>` +
                `<td class="num hide-s">${rg ? fmtF(r.n) : fmtInt(r.n)}</td>` +
                `<td class="num${r.r < 0 ? " neg" : ""}">${fmtR(r.r)}</td>` +
                (cmp
                  ? `<td class="num bcol${r.pr < 0 ? " neg" : ""}">${fmtR(r.pr)}</td>` +
                    `<td class="num${lt05(r.q) ? " sig" : ""}">${fmtP(r.q)}</td><td class="num${lt05(r.pq) ? " sig" : ""}">${fmtP(r.pq)}</td></tr>`
                  : `<td class="num hide-m">${rg ? fmtF(r.se) : fmtB(r.b, r.se)}</td>` +
                    `<td class="num${isSig(r) ? " sig" : ""}">${fmtP(r.p)}</td>` +
                    `<td class="num">${fmtP(r.q)}</td></tr>`)
              );
            })
            .join("")
        : `<tr><td colspan="7">No traits match these filters. Clear the search or pick another category.</td></tr>`,
    );
    tb.querySelectorAll("tr[data-id]").forEach((tr) => {
      const activate = () => select(tr.dataset.id, false);
      tr.addEventListener("click", activate);
      tr.addEventListener("keydown", (e) => {
        if (e.key === "Enter" || e.key === " ") {
          e.preventDefault();
          activate();
        }
      });
    });
    const nSig = rows.filter(isSig).length;
    qsel("#count").textContent = cmp
      ? `${fmtInt(rows.length)} of ${fmtInt(view.rows.length)} trait pairs shown, ${fmtInt(nSig)} significant in both analyses (FDR < 5%)`
      : `${fmtInt(rows.length)} of ${fmtInt(view.rows.length)} traits shown, ${fmtInt(nSig)} at FDR < 5%`;
    qsel("#pageinfo").textContent = `Page ${state.page + 1} of ${pages}`;
    qsel("#prev").disabled = state.page === 0;
    qsel("#next").disabled = state.page >= pages - 1;
    document.querySelectorAll("#table th").forEach((th) => {
      th.setAttribute(
        "aria-sort",
        th.dataset.k === state.sort
          ? state.asc
            ? "ascending"
            : "descending"
          : "none",
      );
    });
    setHtml(
      qsel("#table th[data-k='r']"),
      rg ? "r<sub>g</sub>" : state.s === "sexdiff" ? "Δr" : "r",
    );
    setHtml(qsel("#table th[data-k='n']"), rg ? "h²" : "N");
    qsel("#table th[data-k='b']").textContent = cmp
      ? "r"
      : rg
        ? "SE"
        : "β (SE)";
    setHtml(qsel("#table th[data-k='p']"), cmp ? "FDR r<sub>g</sub>" : "p");
    setHtml(qsel("#table th[data-k='q']"), cmp ? "FDR r" : "FDR");
    qsel("#table th[data-k='b']").classList.toggle("bcol", cmp);
    qsel("#table th[data-k='n']").title = rg
      ? "SNP heritability of the UK Biobank trait (LDSC, observed scale)"
      : "";
  }

  // ---------- side panel: category overview or one trait ----------
  /** Draw a small forest plot of estimates with 95% confidence intervals. */
  function forest(el, series, xTitle, height) {
    const ink = css("--ink"),
      muted = css("--muted"),
      rule = css("--rule");
    const traces = series.map((s) => ({
      type: "scatter",
      mode: "markers",
      name: s.name,
      x: s.pts.map((p) => p.x),
      y: s.pts.map((p) => 2 - p.k + (s.off || 0)),
      text: s.pts.map((p) => p.text),
      hovertemplate: "%{text}<extra></extra>",
      error_x: {
        type: "data",
        symmetric: false,
        array: s.pts.map((p) => p.hi - p.x),
        arrayminus: s.pts.map((p) => p.x - p.lo),
        color: s.color,
        thickness: 1.4,
        width: 0,
      },
      marker: { color: s.color, size: 7, symbol: s.symbol || "circle" },
    }));
    // the x-axis always includes 0, so one can see whether a confidence interval excludes it
    const vals = series
      .flatMap((s) => s.pts.flatMap((p) => [p.lo, p.hi, p.x]))
      .filter((v) => !isNil(v) && isFinite(v));
    let x0 = Math.min(0, ...vals),
      x1 = Math.max(0, ...vals);
    const pad = (x1 - x0 || 0.1) * 0.08;
    x0 -= pad;
    x1 += pad;
    Plotly.react(
      el,
      traces,
      {
        margin: { l: 82, r: 10, t: 6, b: series.length > 1 ? 54 : 36 },
        height,
        paper_bgcolor: "rgba(0,0,0,0)",
        plot_bgcolor: "rgba(0,0,0,0)",
        font: { family: "Source Sans 3, sans-serif", color: muted, size: 12 },
        xaxis: {
          range: [x0, x1],
          zeroline: true,
          zerolinecolor: ink,
          zerolinewidth: 1,
          gridcolor: rule,
          griddash: "dot",
          fixedrange: true,
          title: { text: xTitle, standoff: 4 },
        },
        yaxis: {
          tickvals: [2, 1, 0],
          ticktext: modelKeys.map((m) => shortModel[m]),
          range: [-0.5, 2.5],
          fixedrange: true,
          showgrid: false,
          zeroline: false,
        },
        showlegend: series.length > 1,
        legend: {
          orientation: "h",
          x: 0,
          y: -0.32,
          font: { size: 12.5, color: ink },
        },
        hoverlabel: {
          bgcolor: css("--surface"),
          bordercolor: rule,
          font: { color: ink, size: 12.5 },
        },
      },
      { displayModeBar: false, responsive: true },
    );
  }

  /** Side panel without a selected trait: significant traits per category. */
  function drawOverview() {
    const box = qsel("#detail");
    box.classList.add("is-overview");
    const c = view.counts;
    const cats = view.cats.filter((k) => modelKeys.some((m) => c[m][k]));
    cats.sort(
      (a, b) =>
        (c[state.m][b] || 0) - (c[state.m][a] || 0) ||
        modelKeys.reduce((s, m) => s + (c[m][b] || 0) - (c[m][a] || 0), 0),
    );
    const what =
      state.d === "rg"
        ? "genetic correlations at FDR &lt; 5%"
        : state.d === "cmp"
          ? "trait pairs significant in both analyses (FDR &lt; 5%)"
          : state.s === "sexdiff"
            ? "sex differences at FDR &lt; 5%"
            : "associations at FDR &lt; 5%";
    setHtml(
      box,
      `<h2>FDR hits by category</h2>` +
        `<p class="lede">Number of ${what} for each brain age model. Select a bar to filter the plot and table; select a trait for its details.</p>` +
        (cats.length
          ? `<div class="bar-key" aria-hidden="true">${modelKeys.map((m) => `<span class="${m === state.m ? "on" : ""}"><i style="background:var(--m-${m})"></i>${shortModel[m]}</span>`).join("")}</div><div id="catbars" class="catbars"></div>`
          : `<p class="lede">No category has results at FDR &lt; 5% in this view.</p>`) +
        (!isNil(state.cat)
          ? `<button type="button" class="btn" id="clear-cat">Show all categories</button>`
          : ""),
    );
    if (!isNil(state.cat))
      qsel("#clear-cat").addEventListener("click", () => setCat(state.cat));
    if (!cats.length) return;
    const ink = css("--ink"),
      muted = css("--muted"),
      rule = css("--rule");
    const names = cats.map(
      (k) => catShort[meta.categories[k]] || meta.categories[k],
    );
    const traces = modelKeys.map((m) => ({
      type: "bar",
      orientation: "h",
      name: shortModel[m],
      y: names,
      x: cats.map((k) => c[m][k] || 0),
      customdata: cats,
      marker: {
        color: css(`--m-${m}`),
        opacity: cats.map((k) =>
          isNil(state.cat) || state.cat === k ? 1 : 0.3,
        ),
        line: { width: m === state.m ? 1.5 : 0, color: ink },
      },
      hovertemplate: `%{y}<br>${measures[m]}: %{x}${state.d === "cmp" ? " pairs significant in both" : " at FDR < 5%"}<extra></extra>`,
    }));
    const el = qsel("#catbars");
    Plotly.react(
      el,
      traces,
      {
        height: 34 + cats.length * 34,
        barmode: "group",
        bargap: 0.28,
        bargroupgap: 0.08,
        margin: { l: 128, r: 10, t: 4, b: 30 },
        paper_bgcolor: "rgba(0,0,0,0)",
        plot_bgcolor: "rgba(0,0,0,0)",
        font: { family: "Source Sans 3, sans-serif", color: muted, size: 12 },
        xaxis: {
          gridcolor: rule,
          griddash: "dot",
          fixedrange: true,
          zeroline: false,
          rangemode: "tozero",
        },
        yaxis: {
          autorange: "reversed",
          fixedrange: true,
          ticks: "",
          tickfont: { color: ink },
        },
        showlegend: false,
        hoverlabel: {
          bgcolor: css("--surface"),
          bordercolor: rule,
          font: { color: ink, size: 12.5 },
        },
      },
      { displayModeBar: false, responsive: true },
    ).then(() => {
      el.on("plotly_click", (ev) => {
        const pt = ev.points && ev.points[0];
        if (pt) setCat(pt.customdata);
      });
    });
  }

  /**
   * Side panel for the selected trait, or the overview when none is selected.
   */
  async function drawDetail() {
    if (!view || isOwn(state.d)) return;
    // the selected trait opens in its own panel under the table; the category overview stays in the sidebar
    const box = qsel("#trait");
    const i = !isNil(state.t) ? meta.id.indexOf(state.t) : -1;
    box.hidden = i < 0;
    if (i < 0) {
      setHtml(box, "");
      return;
    }
    const ink = css("--ink");
    const field = meta.field[i];
    const header = (extra) =>
      `<button type="button" class="btn close" id="close-detail">Close</button>` +
      `<button type="button" class="btn close share" data-share>Copy link</button>` +
      `<h2>${esc(meta.desc[i])}</h2><div class="meta">` +
      `<div><span class="dot" style="background:${catColors[meta.cat[i]]}"></span>${esc(meta.categories[meta.cat[i]])}</div>` +
      extra +
      `<div class="path">${esc(meta.paths[meta.path[i]])}</div></div>`;
    const showcase = (f, code) =>
      !isNil(f)
        ? `<div><a href="https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=${f}">UK Biobank field ${f}</a>${code ? `, coding ${esc(code)}` : ""}</div>`
        : "";
    const all = await load("all");

    if (isGen(state.d)) {
      const rgRows = modelKeys
        .map(
          (m) =>
            `<tr${m === state.m ? ' class="sel"' : ""}><td>${shortModel[m]}</td><td class="num">${fmtR(rgMeta[m].r[i])}</td>` +
            `<td class="num">${fmtF(rgMeta[m].se[i])}</td><td class="num${lt05(rgMeta[m].q[i]) ? " sig" : ""}">${fmtP(rgMeta[m].p[i])}</td><td class="num">${fmtP(rgMeta[m].q[i])}</td></tr>`,
        )
        .join("");
      setHtml(
        box,
        header(
          showcase(field, null) +
            `<div>SNP heritability h² = ${fmtF(rgMeta.h2[i])} (SE ${fmtF(rgMeta.h2se[i])})</div>`,
        ) +
          `<div class="trait-body"><div><h3>Genetic correlation r<sub>g</sub> with 95% CI</h3><div class="forest" id="forest"></div></div><div>` +
          `<div class="mini-wrap"><table class="mini"><thead><tr><th>Model</th><th class="num">r<sub>g</sub></th><th class="num">SE</th><th class="num">p</th><th class="num">FDR</th></tr></thead><tbody>${rgRows}</tbody></table></div>` +
          `</div></div>`,
      );
      forest(
        "forest",
        [
          {
            name: "rg",
            color: ink,
            symbol: "diamond",
            pts: modelKeys
              .map((m, k) => {
                const r = rgMeta[m].r[i],
                  se = rgMeta[m].se[i];
                return isNil(r)
                  ? null
                  : {
                      k,
                      x: r,
                      lo: r - 1.96 * se,
                      hi: r + 1.96 * se,
                      text: `${measures[m]}<br>r<sub>g</sub> = ${fmtR(r)} [${fmtR(r - 1.96 * se)}, ${fmtR(r + 1.96 * se)}]<br>p = ${fmtP(rgMeta[m].p[i])}`,
                    };
              })
              .filter(Boolean),
          },
        ],
        "r<sub>g</sub>",
        190,
      );
    } else {
      const [f, ma, sd] = await Promise.all(
        ["female", "male", "sexdiff"].map(load),
      );
      const data = { all, female: f, male: ma };
      const id = meta.id[i];
      const code = id.includes("#")
        ? id.split("#")[1]
        : id.includes("-")
          ? id.split("-").slice(1).join("-")
          : null;
      const nstr = all.n[i] || f.n[i] || ma.n[i] || "";
      const nm = /^(\d+)\/(\d+)\((\d+)\)$/.exec(nstr);
      const nText = nm
        ? `N = ${fmtInt(+nm[3])} (${fmtInt(+nm[1])} / ${fmtInt(+nm[2])} by outcome)`
        : !isNil(all.ntotal[i])
          ? `N = ${fmtInt(all.ntotal[i])}`
          : "";
      const rowsHtml = modelKeys
        .map((m) => {
          const s = all[m];
          if (isNil(s.p[i]))
            return `<tr><td>${shortModel[m]}</td><td class="num" colspan="4">not tested in full sample</td></tr>`;
          return (
            `<tr${m === state.m ? ' class="sel"' : ""}><td>${shortModel[m]}</td><td class="num">${fmtR(s.r[i])}</td>` +
            `<td class="num hide-s">${fmtB(s.b[i], s.se[i])}</td><td class="num${lt05(s.q[i]) ? " sig" : ""}">${fmtP(s.p[i])}</td><td class="num">${fmtP(s.q[i])}</td></tr>`
          );
        })
        .join("");
      const sexRows = modelKeys
        .map((m) =>
          isNil(sd[m][i])
            ? ""
            : `<tr><td>${shortModel[m]}</td><td class="num">${fmtR(f[m].r[i])}</td><td class="num">${fmtR(ma[m].r[i])}</td><td class="num">${fmtP(sd[m][i])}</td></tr>`,
        )
        .join("");

      setHtml(
        box,
        header(
          showcase(field, code) +
            `<div>${esc(meta.resTypes[meta.rt[i]].toLowerCase().replace("-", " "))} regression${nText ? `, ${nText}` : ""}</div>`,
        ) +
          `<div class="trait-body"><div><h3>Effect size r with 95% CI</h3><div class="forest" id="forest"></div></div><div>` +
          `<h3>All participants</h3><div class="mini-wrap"><table class="mini"><thead><tr><th>Model</th><th class="num">r</th><th class="num hide-s">β (SE)</th><th class="num">p</th><th class="num">FDR</th></tr></thead><tbody>${rowsHtml}</tbody></table></div>` +
          (sexRows
            ? `<h3>Women versus men</h3><div class="mini-wrap"><table class="mini"><thead><tr><th>Model</th><th class="num">r women</th><th class="num">r men</th><th class="num">p diff.</th></tr></thead><tbody>${sexRows}</tbody></table></div>`
            : "") +
          `<p class="note">Model rows: brain age gap from grey matter, white matter, or both combined.</p></div></div>`,
      );

      const sColors = { all: ink, female: "#b03a7a", male: "#2f74c0" };
      const offs = { all: 0.22, female: 0, male: -0.22 };
      forest(
        "forest",
        Object.keys(sampleNames).map((s) => ({
          name: sampleNames[s],
          color: sColors[s],
          symbol: s === "all" ? "diamond" : "circle",
          off: offs[s],
          pts: modelKeys
            .map((m, k) => {
              const d = data[s],
                r = d[m].r[i],
                n = d.ntotal[i];
              if (isNil(r) || !n) return null;
              const zf = Math.atanh(r),
                se = 1 / Math.sqrt(n - 3);
              const lo = Math.tanh(zf - 1.96 * se),
                hi = Math.tanh(zf + 1.96 * se);
              return {
                k,
                x: r,
                lo,
                hi,
                text: `${sampleNames[s]}, ${measures[m].toLowerCase()}<br>r = ${fmtR(r)} [${fmtR(lo)}, ${fmtR(hi)}]<br>p = ${fmtP(d[m].p[i])}, N = ${fmtInt(n)}`,
              };
            })
            .filter(Boolean),
        })),
        "r",
        250,
      );
    }
    qsel("#close-detail").addEventListener("click", () => select(null, false));
    box
      .querySelectorAll(".jump")
      .forEach((b) => b.addEventListener("click", () => setMode(b.dataset.d)));
  }

  // ---------- heritability panels (values as reported in the result tables) ----------
  /** HTML key for the three brain age models. */
  function modeKey(hollow) {
    return (
      modelKeys
        .map(
          (m) =>
            `<span><i style="background:var(--m-${m})"></i>${shortModel[m]}</span>`,
        )
        .join("") +
      (hollow
        ? `<span><i class="hollow"></i>not significant (FDR ≥ 5%)</span>`
        : "")
    );
  }
  /** Shared Plotly layout of the heritability panels. */
  function hLayout(extra) {
    const ink = css("--ink"),
      muted = css("--muted"),
      rule = css("--rule");
    return Object.assign(
      {
        paper_bgcolor: "rgba(0,0,0,0)",
        plot_bgcolor: "rgba(0,0,0,0)",
        showlegend: false,
        font: { family: "Source Sans 3, sans-serif", color: muted, size: 12.5 },
        hoverlabel: {
          bgcolor: css("--surface"),
          bordercolor: rule,
          font: { color: ink, size: 12.5 },
          align: "left",
        },
        modebar: { bgcolor: "rgba(0,0,0,0)", color: muted, activecolor: ink },
      },
      extra,
    );
  }

  /**
   * Draw the heritability panels: SNP heritability, polygenicity and
   * partitioned heritability.
   */
  async function drawHerit() {
    const heritData = await load("herit");
    const ink = css("--ink"),
      muted = css("--muted"),
      rule = css("--rule");
    const narrow = window.innerWidth < 640;
    document
      .querySelectorAll("#herit .model-key")
      .forEach((el, k) => setHtml(el, modeKey(k > 0)));
    const off = { gm: 0.24, wm: 0, gwm: -0.24 };

    // SNP heritability: one row per sample, three models per row
    const n = heritData.ldsc.length,
      labels = heritData.ldsc.map((r) => r.label);
    const h2traces = modelKeys.map((m) => ({
      type: "scatter",
      mode: "markers",
      name: shortModel[m],
      x: heritData.ldsc.map((r) => r[m].h2),
      y: heritData.ldsc.map((_, k) => n - 1 - k + off[m]),
      error_x: {
        type: "data",
        array: heritData.ldsc.map((r) => 1.96 * r[m].se),
        color: css(`--m-${m}`),
        thickness: 1.6,
        width: 0,
      },
      marker: {
        color: css(`--m-${m}`),
        size: 9,
        line: { width: 1, color: ink },
      },
      text: heritData.ldsc.map(
        (r) =>
          `<b>${r.label}</b>${r.n ? `, n = ${fmtInt(r.n)}` : ""}<br>${measures[m]}<br>h² = ${r[m].h2.toFixed(3)} (SE ${r[m].se.toFixed(3)})<br>95% CI ${(r[m].h2 - 1.96 * r[m].se).toFixed(3)} to ${(r[m].h2 + 1.96 * r[m].se).toFixed(3)}<br>LDSC intercept = ${r[m].intercept.toFixed(3)} (SE ${r[m].intercept_se.toFixed(3)})` +
          (r.key === "female" || r.key === "male"
            ? "<br><i>UK Biobank only</i>"
            : ""),
      ),
      hovertemplate: "%{text}<extra></extra>",
    }));
    Plotly.react(
      "h2plot",
      h2traces,
      hLayout({
        height: 64 + n * 52,
        margin: { l: 118, r: 12, t: 8, b: 40 },
        xaxis: {
          title: { text: "SNP heritability h²", standoff: 4 },
          range: [0, 0.42],
          gridcolor: rule,
          griddash: "dot",
          zeroline: false,
          fixedrange: true,
        },
        yaxis: {
          tickvals: labels.map((_, k) => n - 1 - k),
          ticktext: heritData.ldsc.map((r) =>
            r.n
              ? `${r.label}<br><span style="font-size:11px;color:${muted}">n = ${fmtInt(r.n)}</span>`
              : r.label,
          ),
          range: [-0.6, n - 0.4],
          fixedrange: true,
          showgrid: false,
          zeroline: false,
          tickfont: { color: ink },
        },
        shapes: labels.slice(1).map((_, k) => ({
          type: "line",
          xref: "paper",
          x0: 0,
          x1: 1,
          y0: n - 1.5 - k,
          y1: n - 1.5 - k,
          line: { color: rule, width: 1 },
        })),
      }),
      heritPlotConfig,
    );

    // GENESIS: number of SNPs with non-zero effect; brain age models in model colours, reference traits grey
    const genesisData = heritData.genesis,
      gn = genesisData.length;
    const gname = (r) => (r.model ? measures[r.model] : r.label);
    const gshort = (r) =>
      r.model ? measures[r.model] : r.label.replace(/ \(.*\)/, "") + " (ref.)";
    const gcol = (r) =>
      r.model
        ? css(`--m-${r.model}`)
        : hexA(
            css("--faint").startsWith("#") ? css("--faint") : "#8a93a6",
            0.55,
          );
    Plotly.react(
      "genplot",
      [
        {
          type: "bar",
          orientation: "h",
          x: genesisData.map((r) => r.causal),
          y: genesisData.map((_, k) => gn - 1 - k),
          error_x: {
            type: "data",
            array: genesisData.map((r) => r.causal_se),
            color: ink,
            thickness: 1.2,
            width: 4,
          },
          marker: {
            color: genesisData.map(gcol),
            line: { width: 1, color: ink },
          },
          text: genesisData.map(
            (r) =>
              `<b>${esc(gname(r))}</b><br>SNPs with non-zero effect: ${fmtInt(r.causal)} (SE ${fmtInt(r.causal_se)})<br>in the large-effect component: ${fmtInt(r.causal_large)} (SE ${fmtInt(Math.round(r.causal_large_se))})` +
              `<br>N to explain 80% of h²: ${fmtN(r.reqsample)}<br>expected loci at that N: ${fmtInt(Math.round(r.reqloci))}`,
          ),
          hovertemplate: "%{text}<extra></extra>",
          textposition: "none",
        },
      ],
      hLayout({
        height: 64 + gn * 40,
        margin: { l: 140, r: 16, t: 8, b: 40 },
        bargap: 0.62,
        xaxis: {
          title: { text: "SNPs with non-zero effect", standoff: 4 },
          gridcolor: rule,
          griddash: "dot",
          zeroline: false,
          fixedrange: true,
          rangemode: "tozero",
        },
        yaxis: {
          tickvals: genesisData.map((_, k) => gn - 1 - k),
          ticktext: genesisData.map(gshort),
          fixedrange: true,
          showgrid: false,
          zeroline: false,
          tickfont: { color: ink },
        },
        shapes: [
          {
            type: "line",
            xref: "paper",
            x0: 0,
            x1: 1,
            y0: 1.5,
            y1: 1.5,
            line: { color: rule, width: 1, dash: "dot" },
          },
        ],
      }),
      heritPlotConfig,
    );

    // partitioned heritability and cell-type groups: enrichment per annotation, filled = FDR < 5%
    const enrichPlot = (el, panel) => {
      const k = panel.annotation.length;
      const order = panel.annotation
        .map((_, i) => i)
        .sort((a, b) => panel.gwm.enr[b] - panel.gwm.enr[a]);
      const yOf = (i) => k - 1 - order.indexOf(i);
      const traces = modelKeys.map((m) => {
        const col = css(`--m-${m}`);
        return {
          type: "scatter",
          mode: "markers",
          name: shortModel[m],
          x: order.map((i) => panel[m].enr[i]),
          y: order.map((i) => yOf(i) + off[m] * 0.75),
          marker: {
            size: 8,
            color: order.map((i) =>
              lt05(panel[m].q[i]) ? col : "rgba(0,0,0,0)",
            ),
            line: {
              width: 1.6,
              color: order.map((i) => (lt05(panel[m].q[i]) ? ink : col)),
            },
          },
          text: order.map(
            (i) =>
              `<b>${esc(panel.annotation[i])}</b><br>${measures[m]}<br>enrichment = ${fmtF(panel[m].enr[i], 2)}<br>share of SNPs = ${(100 * panel.propSnps[i]).toFixed(1)}%, share of h² = ${(100 * panel[m].propH2[i]).toFixed(1)}%<br>p (one-sided) = ${fmtP(panel[m].p[i])} · FDR = ${fmtP(panel[m].q[i])}`,
          ),
          hovertemplate: "%{text}<extra></extra>",
        };
      });
      const xs = modelKeys.flatMap((m) => panel[m].enr);
      Plotly.react(
        el,
        traces,
        hLayout({
          height: 66 + k * (narrow ? 30 : 26),
          margin: { l: narrow ? 165 : 230, r: 12, t: 22, b: 40 },
          xaxis: {
            title: { text: "Heritability enrichment", standoff: 4 },
            range: [Math.min(0, ...xs) - 0.5, Math.max(...xs) * 1.06],
            gridcolor: rule,
            griddash: "dot",
            zeroline: false,
            fixedrange: true,
          },
          yaxis: {
            tickvals: order.map(yOf),
            ticktext: order.map((i) =>
              narrow ? trunc(panel.annotation[i], 24) : panel.annotation[i],
            ),
            tickfont: { color: ink, size: narrow ? 11 : 12.5 },
            range: [-0.7, k - 0.3],
            fixedrange: true,
            showgrid: false,
            zeroline: false,
          },
          shapes: order
            .filter((_, j) => j % 2 === 0)
            .map((i) => ({
              type: "rect",
              xref: "paper",
              x0: 0,
              x1: 1,
              y0: yOf(i) - 0.5,
              y1: yOf(i) + 0.5,
              fillcolor: css("--band"),
              line: { width: 0 },
              layer: "below",
            }))
            .concat([
              {
                type: "line",
                yref: "paper",
                x0: 1,
                x1: 1,
                y0: 0,
                y1: 1,
                line: { color: muted, width: 1, dash: "dash" },
              },
            ]),
          annotations: [
            {
              x: 1,
              y: 1,
              xref: "x",
              yref: "paper",
              yanchor: "bottom",
              text: "no enrichment",
              showarrow: false,
              font: { size: 11, color: muted },
            },
          ],
        }),
        heritPlotConfig,
      );
    };
    enrichPlot("baseplot", heritData.baseline);
    enrichPlot("ctgplot", heritData.celltype);
  }
  const fmtN = (x) =>
    x >= 1000000
      ? `${(x / 1000000).toFixed(2).replace(/\.?0+$/, "")} million`
      : `${fmtInt(Math.round(x / 1000))},000`;

  const chrIndex = (c) => (c === "X" || c === "XY" ? 22 : +c - 1); // XY = pseudoautosomal region, shown with X
  const gpos = (chr, bp) => chrStart[chrIndex(chr)] + bp;
  const lociState = {
    q: "",
    m: "all",
    nov: "all",
    sel: null,
    sort: "p",
    asc: true,
    more: false,
  }; // loci open sorted by p, strongest first
  let lociData = null,
    lociBound = false;

  const lociHit = (l) =>
    lociState.m === "all"
      ? l.hits.reduce((a, b) => (b.p < a.p ? b : a))
      : l.hits.find((h) => h.model === lociState.m); // the selected model's lead variant, or the strongest one
  /** Whether a locus passes the current model, search and novelty filters. */
  function lociPasses(l) {
    if (lociState.m !== "all" && !l.models.includes(lociState.m)) return false;
    if (lociState.nov === "novel" && !l.novel) return false;
    if (lociState.nov === "known" && l.novel) return false;
    if (lociState.q) {
      const t = lociState.q.toLowerCase();
      const hay = [
        l.gene,
        l.cytoband,
        "chr" + l.chr,
        ...l.hits.flatMap((h) => [h.id, h.nearest || "", h.prioritized || ""]),
      ]
        .join(" ")
        .toLowerCase();
      if (!hay.includes(t)) return false;
    }
    return true;
  }
  const mchip = (m, best) =>
    `<span class="mchip ${m}${best ? " best" : ""}" style="--mc:var(--m-${m})">${{ gm: "GM", wm: "WM", gwm: "GWM" }[m]}</span>`;

  /** Draw the Manhattan plot image and the table of genomic loci. */
  async function drawLoci() {
    if (!lociData) lociData = (await load("loci")).loci;
    if (!lociBound) bindLoci();
    const all = lociData,
      nNovel = all.filter((l) => l.novel).length;
    const all3 = lociState.m === "all",
      label = all3 ? "" : measures[lociState.m].toLowerCase() + " ";
    const mine = all3 ? all : all.filter((l) => l.models.includes(lociState.m)),
      mNovel = mine.filter((l) => l.novel).length;
    qsel("#loci-title").textContent =
      `${mine.length} genome-wide significant loci for ${label}brain age gap`;
    qsel("#loci-tally").textContent = all3
      ? `Combined European meta-analysis, n = 54,890. ${nNovel} novel and ${all.length - nNovel} previously reported loci across the grey matter, white matter and combined brain age models.`
      : `Combined European meta-analysis, n = 54,890. ${mNovel} of these loci are novel. Across the three brain age models, ${all.length} distinct loci, ${nNovel} of them novel.`;
    setHtml(
      qsel("#loci-note"),
      (all3
        ? "Manhattan plot showing, for each variant, the smallest p value across the three brain age models."
        : "Manhattan plot of the combined European meta-analysis for the selected brain age model.") +
        " Diamonds mark the lead variants of independent loci; the line marks p = 5×10<sup>−8</sup>, and the y-axis is truncated at −log<sub>10</sub>(p) = 40. " +
        (all3
          ? "The table lists all loci; select a row for details."
          : "The table lists all loci found for this model; select a row for details."),
    );
    const img = qsel("#loci-img");
    img.src = `img/manhattan.${lociState.m}.png`;
    img.alt = `Manhattan plot for ${all3 ? "the three brain age models together" : label + "brain age gap"}`;
    if (!isNil(lociState.sel) && !mine.some((l) => l.locus === lociState.sel))
      lociState.sel = null;
    drawLociTable();
    drawLocus();
  }

  /** Draw the table of genomic loci. */
  function drawLociTable() {
    const rows = lociData.filter(lociPasses).map((l) => ({ l, h: lociHit(l) }));
    const dir = lociState.asc ? 1 : -1;
    const key = {
      pos: (r) => gpos(r.l.chr, r.l.bp),
      gene: (r) => r.l.gene.toLowerCase(),
      id: (r) => r.h.id,
      models: (r) => r.l.models.length,
      beta: (r) => Math.abs(r.h.beta),
      p: (r) => r.h.p,
    }[lociState.sort];
    rows.sort((a, b) => (key(a) < key(b) ? -1 : key(a) > key(b) ? 1 : 0) * dir);
    const tb = qsel("#loci-table tbody");
    const shown = lociState.more ? rows : rows.slice(0, lociShow);
    setHtml(
      tb,
      shown
        .map(
          ({ l, h }) =>
            `<tr data-locus="${l.locus}" class="${l.locus === lociState.sel ? "sel" : ""}" tabindex="0">` +
            `<td class="trait">${esc(l.cytoband)}<small>chr${esc(l.chr)}:${fmtInt(h.bp)}</small></td>` +
            `<td class="gene"><b><i>${esc(l.gene)}</i></b>${l.novel ? '<span class="novel-tag">novel</span>' : ""}</td>` +
            `<td class="hide-s vid">${esc(h.id)}<small>${esc(h.a1)}/${esc(h.a2)}, freq. ${fmtF(h.freq, 2)}</small></td>` +
            `<td class="models hide-s">${modelKeys
              .filter((m) => l.models.includes(m))
              .map((m) =>
                mchip(
                  m,
                  lociState.m === "all" && m === h.model && l.models.length > 1,
                ),
              )
              .join("")}</td>` +
            `<td class="num hide-m">${fmtB(h.beta, h.se)}</td>` +
            `<td class="num sig">${fmtP(h.p)}</td></tr>`,
        )
        .join("") ||
        `<tr><td colspan="6">No loci match these filters. Clear the search or choose another option.</td></tr>`,
    );
    tb.querySelectorAll("tr[data-locus]").forEach((tr) => {
      const activate = () => selectLocus(+tr.dataset.locus, true);
      tr.addEventListener("click", activate);
      tr.addEventListener("keydown", (e) => {
        if (e.key === "Enter" || e.key === " ") {
          e.preventDefault();
          activate();
        }
      });
    });
    const nn = rows.filter((r) => r.l.novel).length;
    const total = lociData.filter(
      (l) => lociState.m === "all" || l.models.includes(lociState.m),
    ).length;
    const btn = qsel("#loci-more");
    btn.hidden = rows.length <= lociShow;
    btn.textContent = lociState.more
      ? `Show first ${lociShow} only`
      : `Show all ${rows.length} loci`;
    qsel("#loci-count").textContent =
      `${rows.length === total ? `${total} loci` : `${rows.length} of ${total} loci match`}, ${nn} novel. ` +
      (lociState.m === "all"
        ? "Statistics refer to the model with the smallest p value per locus (outlined in “Models”)."
        : `Statistics refer to ${measures[lociState.m].toLowerCase()} brain age gap; “Models” lists every model for which the locus was found.`);
    document
      .querySelectorAll("#loci-table th")
      .forEach((th) =>
        th.setAttribute(
          "aria-sort",
          th.dataset.k === lociState.sort
            ? lociState.asc
              ? "ascending"
              : "descending"
            : "none",
        ),
      );
  }

  const tip = (k) =>
    lociColumnHelp[k] ? ` title="${esc(lociColumnHelp[k])}"` : "";

  /** "GENE (x) | GENE (y) | ..." as a short list, the rest folded away. */
  function evidenceList(str) {
    if (!str) return '<span class="neg">–</span>';
    const items = str.split(" | ").map((t) => `<li>${esc(t)}</li>`);
    if (items.length <= 4) return `<ul>${items.join("")}</ul>`;
    return `<ul>${items.slice(0, 4).join("")}</ul><details><summary>${items.length - 4} more</summary><ul>${items.slice(4).join("")}</ul></details>`;
  }

  /** Draw the detail panel of the selected locus. */
  function drawLocus() {
    const box = qsel("#locus");
    const l = !isNil(lociState.sel)
      ? lociData.find((x) => x.locus === lociState.sel)
      : null;
    box.hidden = !l;
    if (!l) {
      setHtml(box, "");
      return;
    }
    const hits = modelKeys
      .map((m) => l.hits.find((h) => h.model === m))
      .filter(Boolean);
    const statRows = hits
      .map(
        (h) =>
          `<tr><td>${shortModel[h.model]}</td><td>${esc(h.id)}<small>chr${esc(l.chr)}:${fmtInt(h.bp)}</small></td><td>${esc(h.a1)}/${esc(h.a2)}</td><td class="num hide-s">${fmtF(h.freq, 2)}</td>` +
          `<td class="num">${fmtB(h.beta, h.se)}</td><td class="num sig">${fmtP(h.p)}</td><td class="num hide-s">${!isNil(h.eta2) ? String(+h.eta2.toPrecision(3)) : "–"}</td><td class="num hide-s">${fmtInt(h.n)}</td></tr>`,
      )
      .join("");
    const ev = (label, f) =>
      `<tr><th scope="row"${tip(label)}>${label}</th>${hits.map((h) => `<td>${f(h)}</td>`).join("")}</tr>`;
    const cs = (h) =>
      Object.entries(h.cs)
        .map(([k, v]) => `${k}: ${v ? esc(v) : "–"}`)
        .join("<br>");
    const evidRows =
      ev("Prioritized gene", (h) =>
        h.prioritized ? `<b><i>${esc(h.prioritized)}</i></b>` : "–",
      ) +
      ev("Nearest gene", (h) =>
        h.nearest
          ? `<i>${esc(h.nearest)}</i> (${esc(h.region || "")}${h.distance ? `, ${fmtInt(h.distance)} bp` : ""})` +
            (h.nearest_desc
              ? `<small>${esc(h.nearest_desc)}${h.nearest_type ? `, ${esc(h.nearest_type.replace(/_/g, " "))}` : ""}</small>`
              : "")
          : "–",
      ) +
      ev("Credible set size", cs) +
      Object.keys(hits[0].evidence)
        .map((k) => ev(k, (h) => evidenceList(h.evidence[k])))
        .join("");
    const catalog = hits.map((h) => h.catalog).find(Boolean);
    const nCat = catalog ? catalog.split(" | ").length : 0;
    setHtml(
      box,
      `<button type="button" class="btn close" id="close-locus">Close</button>` +
        `<button type="button" class="btn close share" data-share>Copy link</button>` +
        `<h2><i>${esc(l.gene)}</i> · ${esc(l.cytoband)}${l.novel ? '<span class="novel-tag">novel</span>' : ""}</h2>` +
        `<div class="locus-meta"><p>Chromosome ${esc(l.chr)}, found for ${modelKeys
          .filter((m) => l.models.includes(m))
          .map((m) => measures[m].toLowerCase())
          .join(", ")} brain age gap.</p>` +
        (l.literature
          ? `<p>Previously reported: ${esc(l.literature.replace(/_/g, " "))}</p>`
          : `<p>Not reported in earlier GWAS of brain age gap.</p>`) +
        `</div>` +
        `<h3>Lead variants</h3><div class="mini-wrap"><table class="mini"><thead><tr><th>Model</th>${["Variant", "A1/A2"].map((k) => `<th${tip(k)}>${k}</th>`).join("")}<th class="num hide-s"${tip("Freq.")}>Freq.</th><th class="num"${tip("β (SE)")}>β (SE)</th><th class="num"${tip("p")}>p</th><th class="num hide-s"${tip("ηp²")}>η<sub>p</sub>²</th><th class="num hide-s"${tip("N")}>N</th></tr></thead><tbody>${statRows}</tbody></table></div>` +
        `<h3>Gene prioritization</h3><div class="mini-wrap"><table class="mini evid"><thead><tr><th></th>${hits.map((h) => `<th>${shortModel[h.model]}</th>`).join("")}</tr></thead><tbody>${evidRows}</tbody></table></div>` +
        (catalog
          ? `<details class="catalog"><summary>Associated with ${nCat} trait${nCat > 1 ? "s" : ""} in the GWAS Catalog</summary><p>${esc(catalog.split(" | ").join("; "))}</p></details>`
          : "") +
        `<details class="defs"><summary>What the columns mean</summary><dl>${[
          "Variant",
          "A1/A2",
          "Freq.",
          "β (SE)",
          "p",
          "ηp²",
          "N",
          "Prioritized gene",
          "Nearest gene",
          "Credible set size",
          ...Object.keys(hits[0].evidence),
          "GWAS Catalog",
          "Literature",
        ]
          .filter((k) => lociColumnHelp[k])
          .map(
            (k) =>
              `<dt>${k === "ηp²" ? "η<sub>p</sub>²" : esc(k)}</dt><dd>${esc(lociColumnHelp[k])}</dd>`,
          )
          .join("")}</dl></details>`,
    );
    qsel("#close-locus").addEventListener("click", () =>
      selectLocus(null, false),
    );
  }

  /**
   * Select a locus (or clear the selection) and optionally scroll to its panel.
   */
  function selectLocus(id, scroll) {
    lociState.sel = id;
    writeHash();
    drawLociTable();
    drawLocus();
    if (!isNil(id) && scroll)
      qsel("#locus").scrollIntoView({
        behavior: reduced() ? "auto" : "smooth",
        block: "start",
      });
  }

  /** Attach the event handlers of the loci controls. */
  function bindLoci() {
    lociBound = true;
    let t;
    qsel("#loci-q").addEventListener("input", (e) => {
      clearTimeout(t);
      t = setTimeout(() => {
        lociState.q = e.target.value.trim();
        drawLociTable();
      }, 180);
    });
    qsel("#loci-nov").addEventListener("change", (e) => {
      lociState.nov = e.target.value;
      drawLociTable();
    });
    qsel("#loci-more").addEventListener("click", () => {
      lociState.more = !lociState.more;
      drawLociTable();
      if (!lociState.more)
        qsel("#loci-table").scrollIntoView({
          behavior: reduced() ? "auto" : "smooth",
          block: "start",
        });
    });
    document.querySelectorAll("#loci-table th").forEach((th) =>
      th.addEventListener("click", () => {
        const k = th.dataset.k;
        if (lociState.sort === k) lociState.asc = !lociState.asc;
        else {
          lociState.sort = k;
          lociState.asc = !["beta", "models"].includes(k);
        }
        drawLociTable();
      }),
    );
  }

  // ---------- genetic correlations with 38 selected GWAS (values as in the result table) ----------
  const stars = (x) => (x.fdr < 0.05 ? "**" : x.p < 0.05 ? "*" : "");
  /** Draw the heatmap of genetic correlations with 38 selected traits. */
  async function drawRgSel() {
    const rgselData = await load("rgsel"),
      rgselTraits = rgselData.traits;
    const ink = css("--ink"),
      muted = css("--muted"),
      surf = css("--surface"),
      band = css("--band");
    const narrow = window.innerWidth < 640;
    const nf = modelKeys.map(
      (m) =>
        `${rgselTraits.filter((t) => t[m].fdr < 0.05).length} for ${shortModel[m].toLowerCase()}`,
    );
    qsel("#rgsel-tally").textContent =
      `Combined European meta-analysis, n = 54,890. Correlations passing FDR < 5%: ${nf.join(", ")}.`;
    const n = rgselTraits.length,
      ys = rgselTraits.map((_, k) => n - 1 - k);
    const xl = (m) =>
      narrow ? { gm: "GM", wm: "WM", gwm: "GWM" }[m] : shortModel[m];
    const xs = modelKeys.map(xl);
    const text = rgselTraits.map((t) =>
      modelKeys.map(
        (m) =>
          `<b>${esc(t.trait)}</b> (${esc(t.ref)})<br>${measures[m]} brain age gap<br>r<sub>g</sub> = ${fmtR(t[m].rg)} (SE ${fmtF(t[m].se)})<br>p = ${fmtP(t[m].p)}, FDR = ${fmtP(t[m].fdr)}<br><span style="color:${muted}">h² of trait = ${fmtF(t.h2)} (SE ${fmtF(t.h2_se)})</span>`,
      ),
    );
    const sig = (t, m) => t[m].p < 0.05;
    // two layers: a faint fill for p ≥ 0.05, the r_g colour scale for p < 0.05 (as in the paper figure)
    const base = {
      type: "heatmap",
      x: xs,
      y: ys,
      text,
      hovertemplate: "%{text}<extra></extra>",
      xgap: 2,
      ygap: 2,
    };
    const faint = {
      ...base,
      z: rgselTraits.map((t) => modelKeys.map((m) => (sig(t, m) ? null : 0))),
      colorscale: [
        [0, band],
        [1, band],
      ],
      showscale: false,
    };
    const blue = "#2166ac",
      red = "#b2182b";
    const plotHeight = 50 + n * (narrow ? 20 : 22) + 70;
    const col = {
      ...base,
      z: rgselTraits.map((t) =>
        modelKeys.map((m) =>
          sig(t, m) ? Math.max(-0.3, Math.min(0.3, t[m].rg)) : null,
        ),
      ),
      zmin: -0.3,
      zmax: 0.3,
      colorscale: [
        [0, blue],
        [0.5, surf.startsWith("#") ? surf : "#ffffff"],
        [1, red],
      ],
      colorbar: {
        orientation: "h",
        title: { text: "Genetic correlation r<sub>g</sub>", side: "top" },
        thickness: 10,
        len: narrow ? 0.9 : 0.8,
        lenmode: "fraction",
        tickangle: 0,
        x: 0.5,
        xanchor: "center",
        y: -10 / plotHeight,
        yanchor: "top",
        tickvals: narrow ? [-0.3, 0, 0.3] : [-0.3, -0.15, 0, 0.15, 0.3],
        ticktext: narrow
          ? ["−0.3", "0", "0.3"]
          : ["−0.3", "−0.15", "0", "0.15", "0.3"],
        outlinewidth: 0,
        tickfont: { size: 11 },
      },
    };
    const annotations = [];
    rgselTraits.forEach((t, k) =>
      modelKeys.forEach((m, j) => {
        const st = stars(t[m]);
        if (st)
          annotations.push({
            x: xs[j],
            y: ys[k],
            text: st,
            showarrow: false,
            font: { size: 13, color: Math.abs(t[m].rg) > 0.17 ? "#fff" : ink },
            yshift: -3,
          });
      }),
    );
    // frame around the cells, domain dividers and labels
    const line = { color: muted, width: 1 };
    const shapes = [
      {
        type: "rect",
        xref: "x",
        yref: "y",
        x0: -0.5,
        x1: 2.5,
        y0: -0.5,
        y1: n - 0.5,
        line,
        fillcolor: "rgba(0,0,0,0)",
      },
    ];
    rgselTraits.forEach((t, k) => {
      if (k && t.domain !== rgselTraits[k - 1].domain)
        shapes.push({
          type: "line",
          xref: "x",
          yref: "y",
          x0: -0.5,
          x1: 2.5,
          y0: ys[k] + 0.5,
          y1: ys[k] + 0.5,
          line,
        });
      if (!k || t.domain !== rgselTraits[k - 1].domain) {
        const last = rgselTraits.reduce(
          (a, u, j) => (u.domain === t.domain ? j : a),
          k,
        );
        annotations.push({
          xref: "paper",
          x: 1.03,
          y: (ys[k] + ys[last]) / 2,
          text: narrow
            ? {
                Psychiatric: "Psych.",
                "Substance use": "Subst.",
                Neurological: "Neuro.",
                Personality: "Person.",
                Sleep: "Sleep",
                Cognition: "Cogn.",
                Anthropometric: "Anthrop.",
                Cardiovascular: "Cardio.",
              }[t.domain]
            : t.domain,
          showarrow: false,
          xanchor: "left",
          font: { size: narrow ? 10.5 : 12, color: muted },
        });
      }
    });
    Plotly.react(
      "rgselplot",
      [faint, col],
      hLayout({
        height: plotHeight,
        margin: { l: narrow ? 168 : 200, r: narrow ? 58 : 120, t: 34, b: 78 },
        xaxis: {
          side: "top",
          tickangle: 0,
          fixedrange: true,
          showgrid: false,
          zeroline: false,
          ticks: "outside",
          ticklen: 4,
          tickcolor: "rgba(0,0,0,0)",
          tickfont: { color: ink, size: narrow ? 11.5 : 12.5 },
        },
        yaxis: {
          tickvals: ys,
          ticktext: rgselTraits.map((t) => esc(t.trait)),
          ticks: "outside",
          ticklen: 6,
          tickcolor: "rgba(0,0,0,0)",
          fixedrange: true,
          showgrid: false,
          zeroline: false,
          tickfont: { color: ink, size: narrow ? 10.5 : 12 },
          range: [-0.5, n - 0.5],
        },
        shapes,
        annotations,
      }),
      {
        ...heritPlotConfig,
        toImageButtonOptions: {
          filename: "brainage_rg_selected_traits",
          scale: 3,
        },
      },
    );
  }

  let mrData = null,
    mrSel = null;
  /** Draw the Mendelian randomization heatmap. */
  async function drawMr() {
    mrData = mrData || (await load("mr"));
    const mrTraits = mrData.traits,
      n = mrTraits.length;
    const ink = css("--ink"),
      muted = css("--muted"),
      surf = css("--surface"),
      band = css("--band");
    const narrow = window.innerWidth < 640;
    const mrSummary = pmeta.summary.mr;
    qsel("#mr-tally").textContent =
      `Combined European meta-analysis, n = 54,890. ${mrSummary.to_any} of ${n} traits show an effect on brain age gap and ${mrSummary.from_any} an effect of brain age gap at FDR < 5% for at least one model.`;
    // six columns: two directions × three models
    const cols = mrDirs.flatMap(([dk]) => modelKeys.map((m) => ({ dk, m })));
    const xs = cols.map((_, j) => j + (j >= 3 ? 1 : 0)),
      xAll = [0, 1, 2, 3, 4, 5, 6]; // column 3 is an empty spacer between the directions
    const spread = (f) =>
      mrTraits.map((t) => {
        const r = cols.map((c) => f(t, c));
        r.splice(3, 0, null);
        return r;
      });
    const ys = mrTraits.map((_, k) => n - 1 - k);
    const short = (m) =>
      narrow
        ? { gm: "G", wm: "W", gwm: "GW" }[m]
        : { gm: "Grey", wm: "White", gwm: "G + W" }[m];
    const cell = (t, c) => t[c.dk][c.m];
    const hover = (t, c) => {
      const e = cell(t, c),
        dir =
          c.dk === "to"
            ? `${t.trait} → ${measures[c.m].toLowerCase()} brain age gap`
            : `${measures[c.m]} brain age gap → ${t.trait.toLowerCase()}`;
      if (!e)
        return `<b>${esc(t.trait)}</b> (${esc(t.ref)})<br>${esc(dir)}<br>no estimate available`;
      return (
        `<b>${esc(t.trait)}</b> (${esc(t.ref)})<br>${esc(dir)}<br>GSMR b = ${fmtB(e.b, e.se)}<br>p = ${fmtP(e.p)}, FDR = ${fmtP(e.fdr)}` +
        `<br>${e.nheidi} of ${e.nsnp} instruments kept after HEIDI<br>${e.n05} of 10 MR methods with p < 0.05`
      );
    };
    const z = (t, c) => {
      const e = cell(t, c);
      return e && e.p < 0.05 ? Math.max(-6, Math.min(6, e.b / e.se)) : null;
    };
    const base = {
      type: "heatmap",
      x: xAll,
      y: ys,
      text: spread(hover),
      hovertemplate: "%{text}<extra></extra>",
      hoverongaps: false,
      xgap: 2,
      ygap: 2,
    };
    const faint = {
      ...base,
      z: spread((t, c) => (isNil(z(t, c)) ? 0 : null)),
      colorscale: [
        [0, band],
        [1, band],
      ],
      showscale: false,
    };
    const plotHeight = 86 + n * (narrow ? 24 : 28) + 70;
    const col = {
      ...base,
      z: spread(z),
      zmin: -6,
      zmax: 6,
      colorscale: [
        [0, "#2166ac"],
        [0.5, surf.startsWith("#") ? surf : "#ffffff"],
        [1, "#b2182b"],
      ],
      colorbar: {
        orientation: "h",
        title: { text: "GSMR z-score (b / SE)", side: "top" },
        thickness: 10,
        len: narrow ? 0.9 : 0.6,
        tickangle: 0,
        x: 0.5,
        xanchor: "center",
        y: -10 / plotHeight,
        yanchor: "top",
        tickvals: [-6, -3, 0, 3, 6],
        ticktext: ["≤ −6", "−3", "0", "3", "≥ 6"],
        outlinewidth: 0,
        tickfont: { size: 11 },
      },
    };
    const annotations = [];
    mrTraits.forEach((t, k) =>
      cols.forEach((c, j) => {
        const e = cell(t, c);
        const st = e ? (e.fdr < 0.05 ? "**" : e.p < 0.05 ? "*" : "") : "–";
        if (st)
          annotations.push({
            x: xs[j],
            y: ys[k],
            text: st,
            showarrow: false,
            font: {
              size: 13,
              color: e && Math.abs(e.b / e.se) > 3.5 ? "#fff" : e ? ink : muted,
            },
            yshift: e ? -3 : 0,
          });
      }),
    );
    mrDirs.forEach(([, label], g) =>
      annotations.push({
        x: (xs[g * 3] + xs[g * 3 + 2]) / 2,
        y: 1,
        yref: "paper",
        yanchor: "bottom",
        yshift: 26,
        text: `<b>${narrow ? label.replace(/brain age gap/i, "BAG") : label}</b>`,
        showarrow: false,
        font: { size: narrow ? 11 : 12.5, color: ink },
      }),
    );
    const line = { color: muted, width: 1 };
    const shapes = [0, 1].map((g) => ({
      type: "rect",
      xref: "x",
      yref: "y",
      x0: xs[g * 3] - 0.5,
      x1: xs[g * 3 + 2] + 0.5,
      y0: -0.5,
      y1: n - 0.5,
      line,
      fillcolor: "rgba(0,0,0,0)",
    }));
    if (!isNil(mrSel)) {
      const y = n - 1 - mrSel;
      shapes.push({
        type: "rect",
        xref: "paper",
        yref: "y",
        x0: 0,
        x1: 1,
        y0: y - 0.5,
        y1: y + 0.5,
        line: { color: css("--accent"), width: 2 },
        fillcolor: "rgba(0,0,0,0)",
      });
    }
    Plotly.react(
      "mrplot",
      [faint, col],
      hLayout({
        height: plotHeight,
        margin: { l: narrow ? 142 : 190, r: 6, t: 62, b: 78 },
        xaxis: {
          side: "top",
          range: [-0.5, 6.5],
          tickvals: xs,
          ticktext: cols.map((c) => short(c.m)),
          tickangle: 0,
          fixedrange: true,
          showgrid: false,
          zeroline: false,
          ticks: "outside",
          ticklen: 4,
          tickcolor: "rgba(0,0,0,0)",
          tickfont: { color: ink, size: narrow ? 11 : 12 },
        },
        yaxis: {
          tickvals: ys,
          ticktext: mrTraits.map((t) => esc(t.trait)),
          ticks: "outside",
          ticklen: 6,
          tickcolor: "rgba(0,0,0,0)",
          fixedrange: true,
          showgrid: false,
          zeroline: false,
          tickfont: { color: ink, size: narrow ? 10.5 : 12.5 },
          range: [-0.5, n - 0.5],
        },
        shapes,
        annotations,
      }),
      {
        ...heritPlotConfig,
        toImageButtonOptions: {
          filename: "brainage_mendelian_randomization",
          scale: 3,
        },
      },
    ).then((el) => {
      if (!el._mrClick) {
        el._mrClick = true;
        el.on("plotly_click", (ev) => {
          const pt = ev.points && ev.points[0];
          if (pt && pt.x !== 3) selectMr(n - 1 - pt.y);
        });
      }
    });
    drawMrDetail();
  }
  /** Select a trait in the Mendelian randomization heatmap. */
  function selectMr(k) {
    mrSel = k;
    writeHash();
    drawMr();
    if (!isNil(k))
      setTimeout(
        () =>
          qsel("#mr-detail").scrollIntoView({
            behavior: reduced() ? "auto" : "smooth",
            block: "start",
          }),
        50,
      );
  }
  /** Which of the ten methods agree (p < 0.05), as a compact dot matrix. */
  function drawMrDetail() {
    const box = qsel("#mr-detail"),
      t = !isNil(mrSel) ? mrData.traits[mrSel] : null;
    box.hidden = !t;
    if (!t) {
      setHtml(box, "");
      return;
    }
    const cols = mrDirs.flatMap(([dk, dl]) =>
      modelKeys.map((m) => ({ dk, dl, m, e: t[dk][m] })),
    );
    const where = (c) => `${c.dl}, ${measures[c.m].toLowerCase()}`;
    const dot = (c, k, name) => {
      if (!c.e) return '<td class="dotc neg">–</td>';
      const p = k === "gsmr" ? c.e.p : c.e.methods[k];
      if (isNil(p)) return '<td class="dotc neg">–</td>';
      const on = p < 0.05;
      return `<td class="dotc tap" tabindex="0" data-info="${esc(`<b>${name}</b> · ${where(c)}: p = ${fmtP(p)}`)}"><span class="mdot${on ? " on" : ""}${on && c.e.b < 0 ? " neg-b" : ""}"></span></td>`;
    };
    const head = `<tr><th></th>${mrDirs.map(([, l]) => `<th colspan="3" class="grp">${l}</th>`).join("")}</tr><tr><th></th>${cols.map((c) => `<th class="dotc">${{ gm: "Grey", wm: "White", gwm: "G + W" }[c.m]}</th>`).join("")}</tr>`;
    const rows = mrData.methods
      .map(
        ([k, name]) =>
          `<tr><th scope="row">${esc(name)}</th>${cols.map((c) => dot(c, k, name)).join("")}</tr>`,
      )
      .join("");
    const gsmr = `<tr class="est-row"><th scope="row">GSMR b</th>${cols
      .map((c) => {
        if (!c.e) return '<td class="dotc neg">–</td>';
        const st = c.e.fdr < 0.05 ? "**" : c.e.p < 0.05 ? "*" : "";
        return `<td class="dotc tap${c.e.fdr < 0.05 ? " sig" : ""}" tabindex="0" data-info="${esc(`<b>GSMR</b> · ${where(c)}: b = ${fmtB(c.e.b, c.e.se)}, p = ${fmtP(c.e.p)}, FDR = ${fmtP(c.e.fdr)}; ${c.e.nheidi} of ${c.e.nsnp} instruments kept after HEIDI`)}">${(c.e.b < 0 ? "−" : "") + Math.abs(c.e.b).toPrecision(2)}${st ? `<sup>${st}</sup>` : ""}</td>`;
      })
      .join("")}</tr>`;
    const tally = `<tr class="tally-row"><th scope="row">Methods with p &lt; 0.05</th>${cols.map((c) => `<td class="dotc">${c.e ? `${c.e.n05}/10` : "–"}</td>`).join("")}</tr>`;
    setHtml(
      box,
      `<button type="button" class="btn close" id="close-mr">Close</button>` +
        `<button type="button" class="btn close share" data-share>Copy link</button>` +
        `<h2>${esc(t.trait)}</h2><div class="locus-meta"><p>${esc(t.ref)}</p></div>` +
        `<div class="mini-wrap"><table class="mini mrdots"><thead>${head}</thead><tbody>${gsmr}${rows}${tally}</tbody></table></div>` +
        `<p class="mr-info" id="mr-info" aria-live="polite">Tap or hover a value or dot for details.</p>` +
        `<p class="note"><span class="mdot on"></span> p &lt; 0.05 with a positive GSMR estimate, <span class="mdot on neg-b"></span> with a negative one, <span class="mdot"></span> p ≥ 0.05. GSMR b with * p &lt; 0.05 and ** FDR &lt; 0.05. GSMR uses the instruments kept after the HEIDI outlier test; the other methods use all instruments.</p>`,
    );
    qsel("#close-mr").addEventListener("click", () => selectMr(null));
    // details for a value or dot: on hover with a mouse, on tap with a finger, on focus with the keyboard
    const info = qsel("#mr-info");
    box.querySelectorAll("td.tap").forEach((td) => {
      const show = () => {
        box
          .querySelectorAll("td.tap.on")
          .forEach((x) => x.classList.remove("on"));
        td.classList.add("on");
        setHtml(info, td.dataset.info);
      };
      td.addEventListener("click", show);
      td.addEventListener("mouseenter", show);
      td.addEventListener("focus", show);
    });
  }

  // ---------- fastBAT gene-based tests (values as in the result table) ----------
  const geneState = {
    q: "",
    sig: "bonf",
    sort: "p",
    asc: true,
    more: false,
    bound: false,
    sel: null,
  };
  let geneData = null;
  /** Draw the gene-based results: Manhattan image, tally and table. */
  async function drawGenes() {
    geneData = geneData || (await load("genes"));
    if (!geneState.bound) bindGenes();
    const m = state.m,
      geneSummary = pmeta.summary.genes,
      sm = geneSummary[m],
      model = measures[m].toLowerCase();
    qsel("#genes-title").textContent =
      `${fmtInt(sm.bonf)} genes associated with ${model} brain age gap`;
    setHtml(
      qsel("#genes-tally"),
      `fastBAT gene-based tests, combined European meta-analysis, n = 54,890. ${fmtInt(sm.bonf)} of ${fmtInt(geneSummary.nTested)} genes pass Bonferroni (p &lt; ${fmtP(geneSummary.bonf)}) in ${fmtInt(sm.loci_bonf)} independent loci; ${fmtInt(sm.fdr)} pass FDR &lt; 5%.`,
    );
    const img = qsel("#genes-img");
    img.src = `img/genes.${m}.png`;
    img.alt = `Gene-based Manhattan plot for ${model} brain age gap`;
    if (geneState.sel && !geneData.genes.some((x) => x.gene === geneState.sel))
      geneState.sel = null;
    drawGenesTable();
    drawGene();
  }
  /** Draw the table of gene-based results. */
  function drawGenesTable() {
    const m = state.m,
      bonf = geneData.bonf;
    const pass = (g, mm) =>
      geneState.sig === "bonf" ? g[mm].p < bonf : g[mm].fdr < 0.05;
    const q = geneState.q.toLowerCase();
    let rows = geneData.genes.filter(
      (g) =>
        pass(g, m) &&
        (!q || `${g.gene} ${g.desc || ""} ${g.cyto}`.toLowerCase().includes(q)),
    );
    const dir = geneState.asc ? 1 : -1;
    const key = {
      gene: (g) => g.gene.toLowerCase(),
      pos: (g) => (g.chr === "X" ? 23 : +g.chr) * 1e9 + g.start,
      nsnp: (g) => g.nsnp,
      models: (g) => modelKeys.filter((mm) => pass(g, mm)).length,
      p: (g) => g[m].p,
      fdr: (g) => g[m].fdr,
      lead: (g) => (g[m].lead || "").toLowerCase(),
    }[geneState.sort];
    rows.sort((a, b) => (key(a) < key(b) ? -1 : key(a) > key(b) ? 1 : 0) * dir);
    const total = geneData.genes.filter((g) => pass(g, m)).length;
    const shown = geneState.more ? rows : rows.slice(0, 10);
    setHtml(
      qsel("#genes-table tbody"),
      shown
        .map(
          (g) =>
            `<tr data-gene="${esc(g.gene)}" class="${g.gene === geneState.sel ? "sel" : ""}" tabindex="0">` +
            `<td class="gene"><b><i>${esc(g.gene)}</i></b><small>${esc(g.desc || "")}</small></td>` +
            `<td class="hide-s">${esc(g.cyto)}<small>chr${esc(g.chr)}:${fmtInt(g.start)}–${fmtInt(g.end)}</small></td>` +
            `<td class="hide-m">${g[m].index ? '<span class="novel-tag lead-tag">lead gene</span>' : g[m].lead ? `<i>${esc(g[m].lead)}</i>` : "–"}</td>` +
            `<td class="num hide-m">${fmtInt(g.nsnp)}</td>` +
            `<td class="models hide-s">${modelKeys
              .filter((mm) => pass(g, mm))
              .map((mm) => mchip(mm))
              .join("")}</td>` +
            `<td class="num${g[m].p < bonf ? " sig" : ""}">${fmtP(g[m].p)}</td><td class="num hide-s">${fmtP(g[m].fdr)}</td></tr>`,
        )
        .join("") ||
        `<tr><td colspan="7">No genes match. Clear the search or choose another threshold.</td></tr>`,
    );
    qsel("#genes-table tbody")
      .querySelectorAll("tr[data-gene]")
      .forEach((tr) => {
        const activate = () =>
          selectGene(
            tr.dataset.gene === geneState.sel ? null : tr.dataset.gene,
            true,
          );
        tr.addEventListener("click", activate);
        tr.addEventListener("keydown", (e) => {
          if (e.key === "Enter" || e.key === " ") {
            e.preventDefault();
            activate();
          }
        });
      });
    const btn = qsel("#genes-more");
    btn.hidden = rows.length <= 10;
    btn.textContent = geneState.more
      ? "Show first 10 only"
      : `Show all ${fmtInt(rows.length)} genes`;
    qsel("#genes-count").textContent =
      `${rows.length === total ? fmtInt(total) : `${fmtInt(rows.length)} of ${fmtInt(total)}`} genes pass ${geneState.sig === "bonf" ? "Bonferroni" : "FDR < 5%"} for ${measures[m].toLowerCase()} brain age gap; “Models” lists every model for which the gene passes.`;
    document
      .querySelectorAll("#genes-table th")
      .forEach((th) =>
        th.setAttribute(
          "aria-sort",
          th.dataset.k === geneState.sort
            ? geneState.asc
              ? "ascending"
              : "descending"
            : "none",
        ),
      );
  }
  /**
   * Select a gene (or clear the selection) and optionally scroll to its panel.
   */
  function selectGene(name, scroll) {
    geneState.sel = name;
    writeHash();
    drawGenesTable();
    drawGene();
    if (name && scroll)
      qsel("#gene-detail").scrollIntoView({
        behavior: reduced() ? "auto" : "smooth",
        block: "start",
      });
  }
  /** Draw the detail panel of the selected gene. */
  function drawGene() {
    const box = qsel("#gene-detail"),
      g = geneState.sel
        ? geneData.genes.find((x) => x.gene === geneState.sel)
        : null;
    box.hidden = !g;
    if (!g) {
      setHtml(box, "");
      return;
    }
    const bonf = geneData.bonf;
    const rows = modelKeys
      .map((m) => {
        const e = g[m];
        return (
          `<tr><td>${shortModel[m]}</td><td class="num${e.p < bonf ? " sig" : ""}">${fmtP(e.p)}</td><td class="num">${fmtP(e.fdr)}</td>` +
          `<td>${e.index ? "this gene" : e.lead ? `<a href="#" class="lead-jump" data-gene="${esc(e.lead)}"><i>${esc(e.lead)}</i></a>` : "–"}</td>` +
          `<td>${e.p < bonf ? "Bonferroni" : e.fdr < 0.05 ? "FDR &lt; 5%" : "–"}</td></tr>`
        );
      })
      .join("");
    const links = [
      `<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=${encodeURIComponent(g.gene)}" target="_blank" rel="noopener">GeneCards</a>`,
    ].concat(
      g.entrez
        ? [
            `<a href="https://www.ncbi.nlm.nih.gov/gene/${encodeURIComponent(g.entrez)}" target="_blank" rel="noopener">NCBI Gene</a>`,
          ]
        : [],
    );
    setHtml(
      box,
      `<button type="button" class="btn close" id="close-gene">Close</button>` +
        `<button type="button" class="btn close share" data-share>Copy link</button>` +
        `<h2><i>${esc(g.gene)}</i> · ${esc(g.cyto)}</h2>` +
        `<div class="locus-meta"><p>${esc(g.desc || "")}</p><p>chr${esc(g.chr)}:${fmtInt(g.start)}–${fmtInt(g.end)} (hg19), ${fmtInt(g.nsnp)} SNPs tested · ${links.join(" · ")}</p></div>` +
        `<div class="mini-wrap"><table class="mini"><thead><tr><th>Model</th><th class="num">p</th><th class="num">FDR</th><th>Lead gene of locus</th><th>Passes</th></tr></thead><tbody>${rows}</tbody></table></div>` +
        `<p class="note">Genes within 3 Mb of a more significant gene form one locus; its lead gene is the most significant gene, separately for each brain age model. Select a lead gene to open it.</p>`,
    );
    qsel("#close-gene").addEventListener("click", () => selectGene(null));
    box.querySelectorAll(".lead-jump").forEach((a) =>
      a.addEventListener("click", (e) => {
        e.preventDefault();
        const name = a.dataset.gene;
        if (!geneData.genes.some((x) => x.gene === name)) return;
        geneState.q = "";
        qsel("#genes-q").value = "";
        selectGene(name, true);
      }),
    );
  }
  /** Attach the event handlers of the gene-based controls. */
  function bindGenes() {
    geneState.bound = true;
    let t;
    qsel("#genes-q").addEventListener("input", (e) => {
      clearTimeout(t);
      t = setTimeout(() => {
        geneState.q = e.target.value.trim();
        drawGenesTable();
      }, 180);
    });
    qsel("#genes-sig").addEventListener("change", (e) => {
      geneState.sig = e.target.value;
      geneState.more = false;
      drawGenesTable();
    });
    qsel("#genes-more").addEventListener("click", () => {
      geneState.more = !geneState.more;
      drawGenesTable();
      if (!geneState.more)
        qsel("#genes-table").scrollIntoView({
          behavior: reduced() ? "auto" : "smooth",
          block: "start",
        });
    });
    document.querySelectorAll("#genes-table th").forEach((th) =>
      th.addEventListener("click", () => {
        const k = th.dataset.k;
        if (geneState.sort === k) geneState.asc = !geneState.asc;
        else {
          geneState.sort = k;
          geneState.asc = !["nsnp", "models"].includes(k);
        }
        drawGenesTable();
      }),
    );
  }

  // ---------- one search across all results ----------
  /** Build the index for the global search across all result sets. */
  async function searchIndex() {
    const [all, R, M, L, G] = await Promise.all([
      load("all"),
      load("rgsel"),
      load("mr"),
      load("loci"),
      load("genes"),
    ]);
    const sigP = (i) => modelKeys.some((m) => lt05(all[m].q[i]));
    const sigR = (i) => modelKeys.some((m) => lt05(rgMeta[m].q[i]));
    const items = [];
    pmeta.desc.forEach((d, i) =>
      items.push({
        g: "Phenotypic correlations",
        text: d,
        sub: `${pmeta.id[i]}${sigP(i) ? " · FDR < 5%" : ""}`,
        sig: sigP(i),
        hay: `${d} ${pmeta.id[i]}`,
        go: { d: "pheno", t: pmeta.id[i] },
      }),
    );
    rgMeta.desc.forEach((d, i) =>
      items.push({
        g: "Genetic correlations, UK Biobank traits",
        text: d,
        sub: `field ${rgMeta.field[i] ?? "–"}${sigR(i) ? " · FDR < 5%" : ""}`,
        sig: sigR(i),
        hay: d,
        go: { d: "rg", t: rgMeta.id[i] },
      }),
    );
    R.traits.forEach((t) =>
      items.push({
        g: "Genetic correlations, 38 selected traits",
        text: t.trait,
        sub: t.ref,
        sig: modelKeys.some((m) => t[m].fdr < 0.05),
        hay: `${t.trait} ${t.ref}`,
        go: { d: "rgsel" },
      }),
    );
    M.traits.forEach((t, k) =>
      items.push({
        g: "Mendelian randomization",
        text: t.trait,
        sub: t.ref,
        sig: ["to", "from"].some((dk) =>
          modelKeys.some((m) => t[dk][m] && t[dk][m].fdr < 0.05),
        ),
        hay: `${t.trait} ${t.ref}`,
        go: { d: "mr", trait: k },
      }),
    );
    L.loci.forEach((l) =>
      items.push({
        g: "GWAS loci",
        text: l.gene,
        sub: `${l.cytoband}${l.novel ? " · novel" : ""}`,
        sig: true,
        italic: true,
        hay: [
          l.gene,
          l.cytoband,
          ...l.hits.flatMap((h) => [
            h.id,
            h.nearest || "",
            h.prioritized || "",
          ]),
        ].join(" "),
        go: { d: "loci", locus: l.locus },
      }),
    );
    G.genes.forEach((g) => {
      const best = modelKeys.reduce(
        (a, mm) => (g[mm].p < g[a].p ? mm : a),
        "gwm",
      );
      items.push({
        g: "Gene-based tests",
        text: g.gene,
        sub: `${g.cyto} · p = ${fmtP(g[best].p, false)}`,
        sig: modelKeys.some((mm) => g[mm].p < G.bonf),
        italic: true,
        hay: `${g.gene} ${g.desc || ""} ${g.cyto}`,
        go: { d: "genes", m: best, gq: g.gene },
      });
    });
    items.forEach((x) => (x.hay = x.hay.toLowerCase()));
    return items;
  }
  /** Attach the global search box and its result list. */
  function bindSearch() {
    const input = qsel("#gq"),
      list = qsel("#gq-list");
    let index = null,
      opts = [],
      active = -1,
      tmr;
    const close = () => {
      list.hidden = true;
      input.setAttribute("aria-expanded", "false");
      active = -1;
    };
    const mark = (text, q) => {
      const i = text.toLowerCase().indexOf(q);
      return i < 0
        ? esc(text)
        : `${esc(text.slice(0, i))}<mark>${esc(text.slice(i, i + q.length))}</mark>${esc(text.slice(i + q.length))}`;
    };
    const render = async () => {
      const q = input.value.trim().toLowerCase();
      if (q.length < 2) {
        close();
        return;
      }
      index = index || (await searchIndex());
      opts = [];
      let html = "";
      searchGroups.forEach((g) => {
        const hits = index
          .filter((x) => x.g === g && x.hay.includes(q))
          .sort(
            (a, b) =>
              b.text.toLowerCase().startsWith(q) -
                a.text.toLowerCase().startsWith(q) ||
              b.sig - a.sig ||
              a.text.localeCompare(b.text),
          );
        if (!hits.length) return;
        html += `<div class="gq-group">${esc(g)} <span>${hits.length}</span></div>`;
        hits.slice(0, 4).forEach((x) => {
          opts.push(x.go);
          html += `<div class="gq-opt" role="option" data-k="${opts.length - 1}"><span class="gq-t">${x.italic ? `<i>${mark(x.text, q)}</i>` : mark(x.text, q)}</span><span class="gq-s">${esc(x.sub)}</span></div>`;
        });
        if (
          hits.length > 4 &&
          (g === "Phenotypic correlations" ||
            g.startsWith("Genetic correlations, UK"))
        ) {
          opts.push({
            d: g === "Phenotypic correlations" ? "pheno" : "rg",
            q: input.value.trim(),
          });
          html += `<div class="gq-opt gq-more" role="option" data-k="${opts.length - 1}">Show all ${fmtInt(hits.length)} matches</div>`;
        } else if (hits.length > 4)
          html += `<div class="gq-note">and ${hits.length - 4} more</div>`;
      });
      setHtml(
        list,
        html ||
          `<div class="gq-note">No results for “${esc(input.value.trim())}”.</div>`,
      );
      list.hidden = false;
      input.setAttribute("aria-expanded", "true");
      active = -1;
    };
    const pick = (k) => {
      const spec = opts[k];
      if (!spec) return;
      close();
      input.blur();
      go(spec);
    };
    const highlight = () =>
      list
        .querySelectorAll(".gq-opt")
        .forEach((o) => o.classList.toggle("on", +o.dataset.k === active));
    input.addEventListener("input", () => {
      clearTimeout(tmr);
      tmr = setTimeout(render, 150);
    });
    input.addEventListener("focus", () => {
      if (input.value.trim().length >= 2) render();
      else searchIndex().then((x) => (index = index || x));
    });
    input.addEventListener("keydown", (e) => {
      if (e.key === "Escape") {
        close();
        return;
      }
      if (list.hidden || !opts.length) return;
      if (e.key === "ArrowDown") {
        e.preventDefault();
        active = Math.min(opts.length - 1, active + 1);
        highlight();
      } else if (e.key === "ArrowUp") {
        e.preventDefault();
        active = Math.max(0, active - 1);
        highlight();
      } else if (e.key === "Enter") {
        e.preventDefault();
        pick(active < 0 ? 0 : active);
      }
    });
    list.addEventListener("mousedown", (e) => {
      const o = e.target.closest(".gq-opt");
      if (o) {
        e.preventDefault();
        pick(+o.dataset.k);
      }
    });
    document.addEventListener("click", (e) => {
      if (!e.target.closest(".gsearch")) close();
    });
  }

  // ---------- actions ----------
  /** Select a trait in the correlation views and show its side panel. */
  function select(id) {
    state.t = id;
    writeHash();
    drawPlot();
    drawTable();
    drawDetail().then(() => {
      // clicking a point or a table row takes the reader to the trait panel under the table
      if (id)
        qsel("#trait").scrollIntoView({
          behavior: reduced() ? "auto" : "smooth",
          block: "start",
        });
    });
  }

  /** Update all controls and visible panels to the current state. */
  function syncControls() {
    const cur = {
      d: state.d === "rgsel" ? "rg" : state.d === "genes" ? "loci" : state.d,
      g: state.d === "rgsel" ? "sel" : "ukb",
      lv: state.d === "genes" ? "genes" : "loci",
      v: state.v,
      m: state.d === "loci" ? lociState.m : state.m,
      s: state.s,
      l: String(state.labels),
    };
    document.querySelectorAll(".seg").forEach((g) =>
      g.querySelectorAll("button").forEach((b) => {
        const on = cur[g.dataset.key] === b.dataset.val;
        b.setAttribute("aria-checked", on);
        b.tabIndex = on ? 0 : -1;
        if (g.dataset.key === "v")
          b.hidden = !viewsFor().includes(b.dataset.val);
        if (g.dataset.key === "m" && b.dataset.val === "all")
          b.hidden = state.d !== "loci"; // "All together" exists only for the loci
      }),
    );
    qsel("#ctl-sample").hidden = state.d !== "pheno";
    qsel("#ctl-plot").hidden = state.d === "cmp" || isOwn(state.d);
    qsel("#ctl-model").hidden =
      state.d === "herit" || state.d === "rgsel" || state.d === "mr"; // these panels show all three models side by side
    qsel("#ctl-rgset").hidden = state.d !== "rg" && state.d !== "rgsel";
    qsel("#ctl-level").hidden = state.d !== "loci" && state.d !== "genes";
    qsel("#genes").hidden = state.d !== "genes";
    qsel("#rgsel").hidden = state.d !== "rgsel";
    qsel("#mr").hidden = state.d !== "mr";
    const own = isOwn(state.d);
    qsel("#explore-body").hidden = own;
    qsel("#split").hidden = own;
    qsel("#herit").hidden = state.d !== "herit";
    qsel("#loci").hidden = state.d !== "loci";
    // "How to read this" shows only the explanations for the current result set
    document
      .querySelectorAll("#about [data-for]")
      .forEach(
        (el) => (el.hidden = !el.dataset.for.split(" ").includes(state.d)),
      );
    updateCards();
    updateNext();
    if (lastD !== state.d) {
      // a short fade when the result set changes
      const el = sectionEl();
      if (lastD && !reduced()) {
        el.classList.remove("enter");
        void el.offsetWidth;
        el.classList.add("enter");
      }
      lastD = state.d;
    }
  }
  let lastD = null;
  const sectionEl = () => qsel(sectionOf[state.d] || "#explore-body");

  // ---------- result cards, jump chips and previous/next ----------
  /** Fill the result cards with the fixed headline numbers. */
  function updateCards() {
    const cardSummary = pmeta && pmeta.summary;
    if (!cardSummary) return;
    const pct = (x) => Math.round(x * 100);
    const h = cardSummary.herit.ldsc_combined,
      hs = modelKeys.map((k) => h[k].h2);
    const set = (k, num, sub) => {
      setHtml(qsel(`[data-num="${k}"]`), num);
      if (sub) setHtml(qsel(`[data-sub="${k}"]`), sub);
    };
    // the same numbers whichever model is selected: correlations count traits at FDR < 5% for at least one model
    set(
      "herit",
      `${pct(Math.min(...hs))}–${pct(Math.max(...hs))}%`,
      "explained by common genetic variants",
    );
    set(
      "loci",
      fmtInt(cardSummary.loci.n),
      `genome-wide significant, ${fmtInt(cardSummary.loci.novel)} novel`,
    );
    set(
      "pheno",
      fmtInt(cardSummary.samples.all.any),
      `of ${fmtInt(cardSummary.samples.all.nTested)} traits at FDR &lt; 5%`,
    );
    set(
      "rg",
      fmtInt(cardSummary.rg.any + cardSummary.rgsel.any),
      `of ${fmtInt(cardSummary.rg.nTested + cardSummary.rgsel.n)} traits at FDR &lt; 5%`,
    ); // UK Biobank traits plus the 38 selected traits
    set(
      "cmp",
      `r vs. r<sub>g</sub>`,
      `side by side for ${fmtInt(cardSummary.rg.pairs.gwm.n)} traits`,
    );
    set(
      "mr",
      `${fmtInt(cardSummary.mr.to_any)} of ${fmtInt(cardSummary.mr.n)}`,
      "traits affect brain age gap at FDR &lt; 5%",
    );
    qsel("#rc-note").textContent =
      "Across the three brain age models; correlation counts include traits with FDR < 5% for at least one model.";
  }
  /** Update the previous/next buttons below the results. */
  function updateNext() {
    const k = modeOrder.indexOf(state.d);
    [
      ["#nn-prev", modeOrder[k - 1]],
      ["#nn-next", modeOrder[k + 1]],
    ].forEach(([sel, d]) => {
      const b = qsel(sel);
      b.hidden = !d;
      b.dataset.d = d || "";
      if (d) b.querySelector(".nn-name").textContent = modeName[d];
    });
  }
  const scrollToResults = () => {
    const el = sectionEl();
    el.scrollIntoView({
      behavior: reduced() ? "auto" : "smooth",
      block: "start",
    });
  };
  /**
   * Opens a result set directly (used by previous/next), optionally with a
   * locus or trait selected reset the filters of a panel before a search result
   * opens in it.
   */
  function resetOwnFilters(spec) {
    if (spec.d === "mr") mrSel = null;
    if (spec.d === "genes" && spec.gq) {
      geneState.q = spec.gq;
      geneState.sig = "fdr";
      geneState.more = false;
      qsel("#genes-q").value = geneState.q;
      qsel("#genes-sig").value = "fdr";
    }
    if (spec.d === "loci") {
      lociState.m = "all";
      lociState.q = "";
      lociState.nov = "all";
      lociState.more = false;
      qsel("#loci-q").value = "";
      qsel("#loci-nov").value = "all";
    }
  }
  /** The locus a search result points to, by locus ID or lead gene. */
  function findLocus(spec) {
    if (!lociData) return null;
    if (!isNil(spec.locus)) return lociData.find((x) => x.locus === spec.locus);
    if (spec.gene) return lociData.find((x) => x.gene === spec.gene);
    return null;
  }
  /**
   * Open a panel-based result set (loci, genes, heritability, MR) from a search
   * result.
   */
  async function goOwn(spec) {
    resetOwnFilters(spec);
    if (!isOwn(state.d)) heritFrom = state.d;
    state.d = spec.d;
    state.t = null;
    lociState.sel = null;
    geneState.sel = null;
    syncControls();
    writeHash();
    await ownPanels[spec.d]();
    const l = findLocus(spec);
    if (l) return selectLocus(l.locus, true);
    if (spec.d === "mr" && !isNil(spec.trait)) return selectMr(spec.trait);
    if (spec.d === "genes" && spec.gq) return selectGene(spec.gq, true);
    return scrollToResults();
  }
  /**
   * Open a search result: switch result set and model, then select the trait,
   * locus or gene.
   */
  async function go(spec) {
    if (spec.m) state.m = spec.m;
    if (isOwn(spec.d)) {
      await goOwn(spec);
      return;
    }
    state.d = spec.d;
    state.s = spec.s || "all";
    state.v = state.s === "sexdiff" ? "sex" : "manhattan";
    state.q = spec.q || "";
    qsel("#q").value = state.q;
    state.sig = "all";
    qsel("#sig").value = "all";
    state.cat = null;
    state.page = 0;
    state.t = null;
    meta = isGen(state.d) ? rgMeta : pmeta;
    syncControls();
    await refresh(true);
    await drawDetail();
    if (spec.t) select(spec.t);
    else scrollToResults();
  }
  /** Handle a click on one of the segmented controls. */
  function choose(key, val) {
    if (key === "l") {
      state.labels = +val;
      syncControls();
      drawPlot();
      return;
    }
    if (key === "d") {
      setMode(val);
      return;
    }
    if (key === "g") {
      setMode(val === "sel" ? "rgsel" : "rg");
      return;
    }
    if (key === "lv") {
      setMode(val);
      return;
    }
    if (key === "v") {
      setView(val);
      return;
    }
    if (key === "m") {
      if (state.d === "genes") {
        if (val === state.m || val === "all") return;
        state.m = val;
        geneState.more = false;
        syncControls();
        writeHash();
        drawGenes();
        return;
      }
      if (state.d === "loci") {
        if (val === lociState.m) return;
        lociState.m = val;
        lociState.more = false;
        if (val !== "all") state.m = val;
        syncControls();
        writeHash();
        drawLoci();
        return;
      }
      if (val === state.m) return;
      state.m = val;
      syncControls();
      refresh(true).then(drawDetail);
      return;
    }
    if (key === "s") {
      if (val === state.s) return;
      state.s = val;
      state.page = 0;
      // women vs. men opens on the scatter plot; the scatter plot exists only for women vs. men
      if (val === "sexdiff") state.v = "sex";
      else if (state.v === "sex") state.v = "manhattan";
      syncControls();
      refresh(true).then(drawDetail);
      return;
    }
  }
  /** Plot types offered for the current results and sample. */
  function viewsFor() {
    if (state.d !== "pheno") return viewsByMode[state.d];
    return state.s === "sexdiff"
      ? ["sex", "manhattan", "volcano"]
      : ["manhattan", "volcano"];
  }
  /** Switch the plot type. */
  function setView(v) {
    if (v === state.v || !viewsFor().includes(v)) return;
    state.v = v;
    state.page = 0;
    syncControls();
    refresh(true).then(drawDetail);
  }
  let heritFrom = "pheno"; // result set to return to when leaving the heritability panels
  /** Switch to another result set. */
  function setMode(d) {
    if (d === state.d) return;
    if (isOwn(d)) {
      if (!isOwn(state.d)) heritFrom = state.d;
      state.d = d;
      state.t = null;
      lociState.sel = null;
      mrSel = null;
      geneState.sel = null;
      syncControls();
      writeHash();
      ownPanels[d]();
      return;
    }
    if (isOwn(state.d)) state.d = heritFrom;
    if (d === state.d) {
      syncControls();
      refresh(true).then(drawDetail);
      return;
    }
    state.t = null; // each result set opens on its overview, without a selected trait
    state.d = d;
    if (d === "cmp") state.v = "pg";
    else if (!viewsFor().includes(state.v))
      state.v = d === "pheno" && state.s === "sexdiff" ? "sex" : "manhattan";
    meta = isGen(d) ? rgMeta : pmeta;
    state.page = 0;
    syncControls();
    refresh(true).then(drawDetail);
  }

  /** Update the explanatory text for the current view. */
  function updateText() {
    const n = view.rows.length,
      nSig = view.rows.filter(isSig).length,
      nB = view.rows.filter((r) => r.p < view.bonf).length;
    const model = `<span class="pick">${modelText[state.m]}</span>`;
    const ps = rgMeta.pairStats[state.m];
    setHtml(
      qsel("#sentence"),
      state.d === "cmp"
        ? `Genetic versus phenotypic correlations of ${model} brain age gap across <strong>${fmtInt(n)}</strong> trait pairs`
        : state.d === "rg"
          ? `Genetic correlations of ${model} brain age gap with <strong>${fmtInt(n)}</strong> UK Biobank traits`
          : `Associations of ${model} brain age gap with <strong>${fmtInt(n)}</strong> UK Biobank traits` +
            (state.s === "all"
              ? ""
              : ` in <span class="pick">${sampleText[state.s]}</span>`),
    ); // sample size goes in the line below
    qsel("#tally").textContent =
      state.d === "cmp"
        ? `Across pairs, r = ${ps.r.toFixed(2)} between genetic and phenotypic correlations; ${fmtInt(nSig)} pairs are significant in both analyses (FDR < 5%).`
        : state.d === "rg"
          ? `${fmtInt(nSig)} genetic correlations pass FDR < 5%, ${fmtInt(nB)} pass Bonferroni correction.`
          : phewasN[state.s] +
            (state.s === "sexdiff"
              ? `${fmtInt(nSig)} traits differ between women and men at FDR < 5%, ${fmtInt(nB)} after Bonferroni correction.`
              : `${fmtInt(nSig)} associations pass FDR < 5%, ${fmtInt(nB)} pass Bonferroni correction.`);
    let note = viewNotes[state.d][state.v];
    qsel("#view-note").textContent = note;
    setHtml(qsel("#keymarks"), keyHtml());
    syncControls();
  }

  /** Rebuild (optionally) and redraw the plot, chips and table. */
  async function refresh(rebuild) {
    if (rebuild) {
      qsel("#plot-loading").hidden = false;
      await buildView();
      updateText();
    }
    writeHash();
    drawCats();
    drawPlot();
    drawTable();
    drawOverview();
  }

  /** Download the filtered table as CSV. */
  function download() {
    const rows = sortedFiltered();
    const q = (v) =>
      isNil(v)
        ? ""
        : /[",\n]/.test(String(v))
          ? `"${String(v).replace(/"/g, '""')}"`
          : v;
    let head, cells;
    if (state.d === "cmp") {
      head = [
        "ukb_field",
        "description",
        "category",
        "phewas_varName",
        "h2_obs",
        "rg",
        "rg_se",
        "rg_pvalue",
        "rg_fdr",
        "phenotypic_r",
        "phenotypic_pvalue",
        "phenotypic_fdr",
      ];
      cells = (r) => [
        meta.field[r.i],
        meta.desc[r.i],
        meta.categories[meta.cat[r.i]],
        pmeta.id[r.pidx],
        r.n,
        r.r,
        r.se,
        r.p,
        r.q,
        r.pr,
        r.pp,
        r.pq,
      ];
    } else if (state.d === "rg") {
      head = [
        "ukb_field",
        "description",
        "category",
        "path",
        "h2_obs",
        "rg",
        "se",
        "pvalue",
        "fdr",
        "phenotypic_r",
        "phenotypic_fdr",
      ];
      cells = (r) => [
        meta.field[r.i],
        meta.desc[r.i],
        meta.categories[meta.cat[r.i]],
        meta.paths[meta.path[r.i]],
        r.n,
        r.r,
        r.se,
        r.p,
        r.q,
        r.pr,
        r.pq,
      ];
    } else {
      head = [
        "varName",
        "ukb_field",
        "description",
        "category",
        "path",
        "regression",
        "n",
        "r",
        "beta",
        "se",
        "pvalue",
        "fdr",
      ];
      cells = (r) => [
        meta.id[r.i],
        meta.field[r.i],
        meta.desc[r.i],
        meta.categories[meta.cat[r.i]],
        meta.paths[meta.path[r.i]],
        meta.resTypes[meta.rt[r.i]],
        r.n,
        r.r,
        r.b,
        r.se,
        r.p,
        r.q,
      ];
    }
    const lines = [head.join(",")].concat(
      rows.map((r) => cells(r).map(q).join(",")),
    );
    const blob = new Blob([lines.join("\n")], { type: "text/csv" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download =
      state.d === "cmp"
        ? `brainage_rg_vs_phenotypic_${state.m}.csv`
        : state.d === "rg"
          ? `brainage_rg_${state.m}.csv`
          : `brainage_phewas_${state.m}_${state.s}.csv`;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  }

  /** Attach all event handlers. */
  function bind() {
    document.querySelectorAll(".seg").forEach((g) => {
      const key = g.dataset.key;
      g.querySelectorAll("button").forEach((btn) =>
        btn.addEventListener("click", () => choose(key, btn.dataset.val)),
      );
      g.addEventListener("keydown", (e) => {
        // arrow keys move within a group, skipping hidden buttons
        if (!["ArrowLeft", "ArrowRight"].includes(e.key)) return;
        const bs = [...g.querySelectorAll("button")].filter((x) => !x.hidden),
          k = bs.findIndex((x) => x.getAttribute("aria-checked") === "true");
        const nxt =
          bs[(k + (e.key === "ArrowRight" ? 1 : bs.length - 1)) % bs.length];
        e.preventDefault();
        choose(key, nxt.dataset.val);
        nxt.focus();
      });
    });
    let tmr;
    qsel("#q").addEventListener("input", (e) => {
      clearTimeout(tmr);
      tmr = setTimeout(() => {
        state.q = e.target.value.trim();
        state.page = 0;
        refresh(false);
      }, 180);
    });
    qsel("#sig").addEventListener("change", (e) => {
      state.sig = e.target.value;
      state.page = 0;
      refresh(false);
    });

    qsel("#prev").addEventListener("click", () => {
      state.page--;
      drawTable();
    });
    qsel("#next").addEventListener("click", () => {
      state.page++;
      drawTable();
    });
    qsel("#download").addEventListener("click", download);
    bindSearch();
    // "Copy link" in the detail panels copies the address of the current view, including the open trait or locus
    document.addEventListener("click", (e) => {
      const b = e.target.closest("[data-share]");
      if (!b) return;
      const done = (t) => {
        b.textContent = t;
        setTimeout(() => (b.textContent = "Copy link"), 1800);
      };
      (navigator.clipboard
        ? navigator.clipboard.writeText(location.href)
        : Promise.reject()
      ).then(
        () => done("Link copied"),
        () => done("Copy the address bar"),
      );
    });
    document.querySelectorAll(".next-nav .nn").forEach((b) =>
      b.addEventListener("click", async () => {
        const d = b.dataset.d;
        if (!d) return;
        if (isOwn(d)) await go({ d });
        else {
          setMode(d);
          await new Promise((r) => setTimeout(r, 60));
          scrollToResults();
        }
      }),
    );
    qsel("#copy-bib").addEventListener("click", (e) => {
      navigator.clipboard.writeText(bibtex).then(
        () => {
          e.target.textContent = "BibTeX copied";
          setTimeout(() => (e.target.textContent = "Copy BibTeX"), 2000);
        },
        () => {
          e.target.textContent =
            "Copy failed, select the citation above instead";
        },
      );
    });
    document.querySelectorAll("#table th").forEach((th) =>
      th.addEventListener("click", () => {
        const k = th.dataset.k;
        if (state.sort === k) state.asc = !state.asc;
        else {
          state.sort = k;
          state.asc = !["r", "n"].includes(k);
        }
        state.page = 0;
        drawTable();
      }),
    );
    let rz;
    window.addEventListener("resize", () => {
      clearTimeout(rz);
      rz = setTimeout(() => {
        if (isOwn(state.d)) ownPanels[state.d]();
        else drawPlot();
      }, 200);
    });
    matchMedia("(prefers-color-scheme: dark)").addEventListener(
      "change",
      () => {
        if (isOwn(state.d)) ownPanels[state.d]();
        else {
          drawPlot();
          drawDetail();
        }
      },
    );
  }

  /** Load the data, restore the state from the URL and draw the page. */
  async function init() {
    readHash();
    if (state.d === "pheno" && state.v === "sex") state.s = "sexdiff"; // the scatter plot exists only for women vs. men
    qsel("#q").value = state.q;
    qsel("#sig").value = state.sig;
    try {
      [pmeta, rgMeta] = await Promise.all([load("meta"), load("rg")]);
      rgMeta.categories = pmeta.categories;
      // fields that share a description (e.g. fluid intelligence at the assessment centre and online) get a short qualifier
      [pmeta, rgMeta].forEach((metaSet) =>
        (metaSet.qual || []).forEach((q, i) => {
          if (q) metaSet.desc[i] += `, ${q}`;
        }),
      );
      rgMeta.phewas.forEach((pi, ri) => {
        if (!isNil(pi)) p2rg[pi] = ri;
      });
      bind();
      if (isOwn(state.d)) {
        // open straight on heritability or loci; other views build when first chosen
        heritFrom = "pheno";
        meta = pmeta;
        syncControls();
        await ownPanels[state.d]();
        const open =
          state.d === "loci" && !isNil(lociState.sel)
            ? qsel("#locus")
            : state.d === "mr" && !isNil(mrSel)
              ? qsel("#mr-detail")
              : state.d === "genes" && geneState.sel
                ? qsel("#gene-detail")
                : null;
        if (open && !open.hidden) open.scrollIntoView({ block: "start" });
        return;
      }
      meta = isGen(state.d) ? rgMeta : pmeta;
      syncControls();
      await refresh(true);
      await drawDetail();
      if (state.t && !qsel("#trait").hidden)
        qsel("#trait").scrollIntoView({ block: "start" }); // a shared link to one trait
    } catch (err) {
      qsel("#plot-loading").textContent =
        `${err.message}. If you opened index.html directly from disk, serve the folder instead (for example: python3 -m http.server).`;
    }
  }
  init();
})();
