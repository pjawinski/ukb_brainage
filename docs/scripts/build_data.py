#!/usr/bin/env python3
"""Build the compact JSON files used by the interactive PheWAS browser (docs/)."""

# Input : results/combined/discovery*.phewas.txt (output of code/genetics/phesant.combine*.R)
# Output: docs/data/meta.json        trait annotations (shared by all samples)
#         docs/data/<sample>.json    per-sample statistics, column-oriented and aligned to meta.json
#         docs/data/rg.json          genetic correlations with Neale lab UK Biobank GWAS
#                                    (gwama.eur.rgNeale.txt)
#         docs/data/herit.json       SNP heritability (LDSC), polygenicity (GENESIS) and partitioned
#                                    heritability (LDSC),
#                                    copied from the result tables without recomputation
#         docs/data/loci.json        genome-wide significant loci with gene prioritisation evidence
#                                    (gwama.eur.snplevel/discoveries.main.txt and discoveries.suppl.txt,
#                                    copied as reported)
#         docs/data/rgsel.json       genetic correlations with 38 selected published GWAS (gwama.eur.rgSelection.txt),
#                                    in the order and domains of code/genetics/rg.plotSelection.R
#         docs/data/genes.json       fastBAT gene-based tests (gwama.eur.fastbat.txt),
#                                    genes with FDR < 5% in at least one model
#         docs/data/mr.json          Mendelian randomization in both directions (gwama.eur.gsmr.multi.labels.txt):
#                                    GSMR estimates and the p-values of the sensitivity methods, as reported
#
# Run from the repository root:
#     python3 docs/scripts/build_data.py

import json
import math
import os
import re
from urllib.parse import unquote

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RES = os.path.join(ROOT, "results", "combined")
OUT = os.path.join(ROOT, "docs", "data")

MEASURES = ["gm", "wm", "gwm"]
SAMPLES = {
    "all": "discovery.phewas.txt",
    "female": "discovery.sex.female.phewas.txt",
    "male": "discovery.sex.male.phewas.txt",
}
SEXDIFF = "discovery.sex.phewas.txt"
RGNEALE = "gwama.eur.rgNeale.txt"


def read(fname):
    """Read a tab-separated result table as strings, exactly as written."""
    return pd.read_csv(os.path.join(RES, fname), sep="\t", dtype=str, keep_default_na=False, quoting=3)


def num(x, sig=4):
    """Round to `sig` significant digits; keep JSON small. None for missing."""
    if x is None or x == "" or x == "NA":
        return None
    v = float(x)
    if math.isnan(v):
        return None
    if v == 0:
        return 0
    return float(f"{v:.{sig}g}")


def pval(x):
    """Round a p-value to 4 significant digits (p-values span hundreds of orders of magnitude)."""
    return num(x, 4)


frames = {k: read(v) for k, v in SAMPLES.items()}
sexdiff = read(SEXDIFF)
rgdf = read(RGNEALE)

# ---- union of traits, ordered like the paper figure: category (alphabetical), then varName ----
meta_rows = {}
for key in ["all", "female", "male"]:
    for _, r in frames[key].iterrows():
        if r.varName not in meta_rows:
            meta_rows[r.varName] = r
for _, r in sexdiff.iterrows():
    if r.varName not in meta_rows:
        r = r.copy()
        r["varType"], r["resType"] = r["a_varType"], r["a_resType"]
        meta_rows[r.varName] = r

ids = sorted(meta_rows, key=lambda v: (meta_rows[v].custom_category, v))
index = {v: i for i, v in enumerate(ids)}
# one category list for PheWAS and genetic correlations, so colours match across both
categories = sorted({meta_rows[v].custom_category for v in ids} | set(rgdf.custom_category))
paths = sorted({meta_rows[v].Path for v in ids})
vtypes = sorted({meta_rows[v].varType for v in ids})
rtypes = sorted({meta_rows[v].resType for v in ids})


def qualifiers(descs, paths_of, fields):
    """Return short labels that tell apart different fields sharing one description."""
    # e.g. the fluid intelligence score from the assessment centre and from the online follow-up;
    # the descriptions themselves stay unchanged
    short = {
        "UK Biobank Assessment Centre": "assessment centre",
        "Assessment Centre": "assessment centre",
        "Online follow-up": "online follow-up",
    }
    groups = {}
    for i, d in enumerate(descs):
        groups.setdefault(d, []).append(i)
    qual = [None] * len(descs)
    for idx in groups.values():
        if len({fields[i] for i in idx}) < 2:
            continue
        segs = [paths_of(i) for i in idx]
        k = next((j for j in range(min(map(len, segs))) if len({s[j] for s in segs}) > 1), None)
        for i, sg in zip(idx, segs):
            q = short.get(sg[k], sg[k][:1].lower() + sg[k][1:]) if k is not None else f"field {fields[i]}"
            qual[i] = q
    return qual


def field_id(var):
    """Return the UK Biobank field ID at the start of a PHESANT variable name."""
    m = re.match(r"^(\d+)", var)
    return int(m.group(1)) if m else None


meta = {
    "categories": categories,
    "paths": paths,
    "varTypes": vtypes,
    "resTypes": rtypes,
    "id": ids,
    "field": [field_id(v) for v in ids],
    "desc": [meta_rows[v].description.strip() for v in ids],
    "cat": [categories.index(meta_rows[v].custom_category) for v in ids],
    "path": [paths.index(meta_rows[v].Path) for v in ids],
    "vt": [vtypes.index(meta_rows[v].varType) for v in ids],
    "rt": [rtypes.index(meta_rows[v].resType) for v in ids],
}
meta["qual"] = qualifiers(meta["desc"], lambda i: paths[meta["path"][i]].split(" > "), meta["field"])

os.makedirs(OUT, exist_ok=True)
summary = {"nTraits": len(ids), "samples": {}}

for key, df in frames.items():
    n = len(ids)
    out = {"ntotal": [None] * n, "n": [None] * n}
    for m in MEASURES:
        out[m] = {k: [None] * n for k in ["b", "se", "r", "p", "q"]}
    for _, r in df.iterrows():
        i = index[r.varName]
        out["ntotal"][i] = int(float(r.ntotal)) if r.ntotal not in ("", "NA") else None
        out["n"][i] = r.n
        for m in MEASURES:
            b, se = num(r[f"gap_{m}_beta"]), num(r[f"gap_{m}_se"])
            if b == -999:  # PHESANT placeholder: multinomial models have no single coefficient
                b, se = None, None
            out[m]["b"][i] = b
            out[m]["se"][i] = se
            out[m]["r"][i] = num(r[f"gap_{m}_rho"])
            out[m]["p"][i] = pval(r[f"gap_{m}_pvalue"])
            out[m]["q"][i] = pval(r[f"gap_{m}_fdr"])
    out["nTested"] = len(df)
    with open(os.path.join(OUT, f"{key}.json"), "w") as f:
        json.dump(out, f, separators=(",", ":"))
    summary["samples"][key] = {
        "nTested": len(df),
        **{m: int((df[f"gap_{m}_fdr"].astype(float) < 0.05).sum()) for m in MEASURES},
        "any": int(
            np.logical_or.reduce([df[f"gap_{m}_fdr"].astype(float) < 0.05 for m in MEASURES]).sum()
        ),  # FDR < 5% for at least one model
    }

# ---- sex-difference p-values (z-test on female vs. male estimates) ----
n = len(ids)
sd = {m: [None] * n for m in MEASURES}
for _, r in sexdiff.iterrows():
    i = index[r.varName]
    for m in MEASURES:
        sd[m][i] = pval(r[f"gap_{m}_deltaP"])
sd["nTested"] = len(sexdiff)
with open(os.path.join(OUT, "sexdiff.json"), "w") as f:
    json.dump(sd, f, separators=(",", ":"))


# ---- genetic correlations with Neale lab UK Biobank GWAS ----
def norm(t):
    """Normalise a trait name for matching: lower case, alphanumerics only."""
    return re.sub(r"[^a-z0-9]+", " ", t.lower()).strip()


by_field = {}
for i in range(len(ids)):
    by_field.setdefault(meta["field"][i], []).append(i)


def match_phewas(field, desc):
    """Link an rg trait to a PheWAS trait only when the match is unambiguous."""
    cands = by_field.get(field, [])
    hit = [i for i in cands if norm(meta["desc"][i]) == norm(desc)]
    if len(hit) == 1:
        return hit[0]
    if len(cands) == 1 and ":" not in desc and ":" not in meta["desc"][cands[0]]:
        return cands[0]
    return None


rgdf["field"] = rgdf.showcase.str.extract(r"id=(\d+)")[0]
rgdf = rgdf.sort_values(["custom_category", "trait_description"]).reset_index(drop=True)
rpaths = sorted(set(rgdf.category))
rg = {"id": [], "field": [], "desc": [], "cat": [], "path": [], "phewas": [], "h2": [], "h2se": []}
for m in MEASURES:
    rg[m] = {k: [] for k in ["r", "se", "p", "q"]}
seen = {}
for _, r in rgdf.iterrows():
    f = int(r.field) if isinstance(r.field, str) else None
    k = seen[f] = seen.get(f, 0) + 1
    rg["id"].append(f"rg{f if f is not None else 'x'}.{k}")
    rg["field"].append(f)
    rg["desc"].append(r.trait_description.strip())
    rg["cat"].append(categories.index(r.custom_category))
    rg["path"].append(rpaths.index(r.category))
    rg["phewas"].append(match_phewas(f, r.trait_description) if f is not None else None)
    rg["h2"].append(num(r.h2_obs))
    rg["h2se"].append(num(r.h2_obs_se))
    for m in MEASURES:
        rg[m]["r"].append(num(r[f"gap_{m}_rg"]))
        rg[m]["se"].append(num(r[f"gap_{m}_se"]))
        rg[m]["p"].append(pval(r[f"gap_{m}_p"]))
        rg[m]["q"].append(pval(r[f"gap_{m}_FDR"]))
# ---- phenotypic vs genetic comparison, selected exactly as in code/genetics/rgVSrp.R ----
# join on trait description, drop multinomial models, keep traits with h2/SE > 1.96
ph_all = frames["all"].copy()
rgdf["_ri"] = range(len(rgdf))
pairs = rgdf.merge(ph_all, left_on="trait_description", right_on="description", how="inner")
pairs = pairs[
    (pairs.resType != "MULTINOMIAL-LOGISTIC") & (pairs.h2_obs.astype(float) / pairs.h2_obs_se.astype(float) > 1.96)
]
rg["pairs"] = {"ri": pairs._ri.astype(int).tolist(), "pi": [index[v] for v in pairs.varName]}
rg["pairStats"] = {}
for m in MEASURES:
    x = pairs[f"gap_{m}_rg"].astype(float).to_numpy()
    y = pairs[f"gap_{m}_rho"].astype(float).to_numpy()
    r, p = stats.pearsonr(x, y)
    slope, intercept = np.polyfit(x, y, 1)
    rg["pairStats"][m] = {
        "n": len(x),
        "r": float(r),
        "p": float(p),
        "mad": float(np.mean(np.abs(x - y))),
        "slope": float(slope),
        "intercept": float(intercept),
    }
rg["paths"] = rpaths
rg["qual"] = qualifiers(rg["desc"], lambda i: rpaths[rg["path"][i]].split(" - ")[::-1], rg["field"])
rg["nTested"] = len(rgdf)
with open(os.path.join(OUT, "rg.json"), "w") as f:
    json.dump(rg, f, separators=(",", ":"))
summary["rg"] = {
    "nTested": len(rgdf),
    "matched": sum(x is not None for x in rg["phewas"]),
    "pairs": rg["pairStats"],
    **{m: int((rgdf[f"gap_{m}_FDR"].astype(float) < 0.05).sum()) for m in MEASURES},
    "any": int(np.logical_or.reduce([rgdf[f"gap_{m}_FDR"].astype(float) < 0.05 for m in MEASURES]).sum()),
}


# ---- heritability: values taken as reported in the result tables ----
def ldsc_h2(fname):
    """Read LDSC SNP heritability per model from a result table."""
    t = read(fname).set_index("trait")
    return {
        m: {
            "h2": num(t.loc[f"gap_{m}", "h2"]),
            "se": num(t.loc[f"gap_{m}", "h2_se"]),
            "intercept": num(t.loc[f"gap_{m}", "intercept"]),
            "intercept_se": num(t.loc[f"gap_{m}", "intercept_se"]),
        }
        for m in MEASURES
    }


H2_SAMPLES = [  # (key, label, file); "combined" = European discovery + replication meta-analysis
    ("discovery", "Discovery", "discovery.ldsc.h2.txt"),
    ("replication", "Replication", "replicate.eur.ldsc.h2.txt"),
    ("combined", "Combined", "gwama.eur.ldsc.h2.txt"),
    ("female", "Combined, women*", "gwama.eur.sex.female.ldsc.h2.txt"),
    ("male", "Combined, men*", "gwama.eur.sex.male.ldsc.h2.txt"),
]
# GWAS sample sizes are not part of the LDSC tables. Discovery: results/mri/accuracy.sample.txt (UKB discovery);
# combined: maximum per-SNP N in results/combined/gwama.eur.snplevel/snplevel.txt; replication: their difference
# (confirmed by the authors). Women and men: provided by the authors; the sex-stratified analyses include UK Biobank
# only (27,789 + 25,268 = 53,057 = 54,890 minus the 1,833 LIFE-Adult participants in results/mri/accuracy.sample.txt).
H2_N = {"discovery": 32634, "replication": 22256, "combined": 54890, "female": 27789, "male": 25268}
herit = {"ldsc": [{"key": k, "label": lab, "file": f, "n": H2_N[k], **ldsc_h2(f)} for k, lab, f in H2_SAMPLES]}

g = read("gwama.eur.genesis.stats.txt")
GENESIS_KEYS = {
    "Grey matter brain age gap": "gm",
    "White matter brain age gap": "wm",
    "Grey and white matter brain age gap": "gwm",
}


def exact(x):
    """Keep GENESIS counts exactly as tabulated (no rounding)."""
    return float(x)


herit["genesis"] = [
    {
        "label": r.trait,
        "model": GENESIS_KEYS.get(r.trait),
        "nsnps": int(float(r.nsnps)),
        "causal": exact(r.sSnps),
        "causal_se": exact(r.sSnps_se),
        "causal_large": exact(r.sSnps_large),
        "causal_large_se": exact(r.sSnps_large_se),
        "reqsample": exact(r.reqsample),
        "reqloci": exact(r.reqsnps),
    }
    for _, r in g.iterrows()
]


def partitioned(fname):
    """Read partitioned heritability (LDSC) per annotation and model."""
    t = read(fname)
    out = {"annotation": t.annotation.str.strip().tolist(), "propSnps": [num(x) for x in t["Prop._SNPs"]]}
    for m in MEASURES:
        out[m] = {
            "propH2": [num(x) for x in t[f"gap_{m}_Prop._h2"]],
            "enr": [num(x) for x in t[f"gap_{m}_Enrichment"]],
            "p": [pval(x) for x in t[f"gap_{m}_Enrichment_p"]],
            "q": [pval(x) for x in t[f"gap_{m}_Enrichment_FDR"]],
        }
    return out


herit["baseline"] = partitioned("gwama.eur.ldsc.partitioned.txt")
herit["celltype"] = partitioned("gwama.eur.ldsc.ctg.txt")
with open(os.path.join(OUT, "herit.json"), "w") as f:
    json.dump(herit, f, separators=(",", ":"))
summary["herit"] = {"ldsc_combined": herit["ldsc"][2], "genesis": herit["genesis"][:3]}

# ---- GWAS loci: values and evidence strings copied as reported ----
LOCI_DIR = os.path.join(RES, "gwama.eur.snplevel")
lm = pd.read_csv(os.path.join(LOCI_DIR, "discoveries.main.txt"), sep="\t", dtype=str, keep_default_na=False)
ls = pd.read_csv(os.path.join(LOCI_DIR, "discoveries.suppl.txt"), sep="\t", dtype=str, keep_default_na=False)
TRAIT_KEY = {"grey matter": "gm", "white matter": "wm", "grey and white matter": "gwm"}
MAIN_KEY = {"GM": "gm", "WM": "wm", "GWM": "gwm"}


def dash(x):
    """Return None for empty, dash or NA cells, else the stripped text."""
    return None if x.strip() in ("", "-", "NA") else x.strip()


loci = []
for k, r in lm.iterrows():
    locus_id = k + 1  # discoveries.main.txt lists loci in LOCUS_COUNT order
    rows = ls[ls.LOCUS_COUNT.astype(int) == locus_id]
    if not len(rows) or rows.CHR.iloc[0] != r.CHR:
        raise ValueError(f"locus {locus_id} does not match")
    first = rows.iloc[0]
    loci.append(
        {
            "locus": locus_id,
            "chr": r.CHR,
            "bp": int(r.BPstr.replace(",", "").strip()),
            "cytoband": first.CYTOBAND,
            "gene": r.prioritized,
            "models": [MAIN_KEY[t] for t in r.TraitsPerLocus.split(",")],
            "novel": r.LITERATURE_SHORT.strip() == "*",
            "literature": dash(first.LITERATURE_LONG) if first.LITERATURE_LONG != "Novel" else None,
            "hits": [
                {
                    "model": TRAIT_KEY[h.TraitNamesSuppl],
                    "id": h.ID,
                    "bp": int(h.BP),
                    "a1": h.A1,
                    "a2": h.A2,
                    "freq": num(h.A1_FREQ),
                    "beta": num(h.BETA),
                    "se": num(h.SE),
                    "z": num(h.Z),
                    "p": pval(h.P),
                    "eta2": num(h.ETA2),
                    "n": int(float(h.N)),
                    "region": dash(h.REGION),
                    "nearest": dash(h.NEAREST_GENE),
                    "nearest_desc": dash(h.NEAREST_GENE_DESCRIPTION),
                    "nearest_type": dash(h.NEAREST_GENE_BIOTYPE),
                    "distance": int(float(h.DISTANCE)) if dash(h.DISTANCE) else None,
                    "cs": {
                        "SBayesRC": dash(h.sbayes_cssizes),
                        "susieR": dash(h.susieR_cssizes),
                        "FINEMAP": dash(h.finemap_cssize),
                    },
                    "evidence": {
                        "SBayesRC genes": dash(h.sbayes_genes),
                        "Nonsynonymous variants": dash(h.nonsyn_customPthresh),
                        "SMR eQTL": dash(h.SMR_eqtl),
                        "SMR sQTL": dash(h.SMR_sqtl),
                        "GTEx single tissue": dash(h.GTEx_singleTissue),
                        "GTEx multi-tissue": dash(h.GTEx_multiTissue),
                        "PoPS": dash(h.PoPS),
                    },
                    "prioritized": dash(h.prioritized),
                    "catalog": dash(h.GWAS_CATALOG),
                }
                for _, h in rows.iterrows()
            ],
        }
    )
with open(os.path.join(OUT, "loci.json"), "w") as f:
    json.dump({"loci": loci, "n": 54890}, f, separators=(",", ":"))
summary["loci"] = {"n": len(loci), "novel": sum(locus["novel"] for locus in loci), "hits": len(ls)}

# ---------- genetic correlations with 38 selected GWAS (as in rg.plotSelection.R) ----------
RGSEL_ORDER = [
    "01_mdd_howard_2019",
    "01_bip_mullins_2021",
    "01_scz_trubetskoy_2022",
    "01_adhd_demontis_2023",
    "01_asd_grove_2019",
    "01_ptsd_nievergelt_2019",
    "01_an_watson_2019",
    "02_smkInit_saunders_2022",
    "02_cpd_saunders_2022",
    "02_dpw_saunders_2022",
    "02_cof_zhong_2019",
    "02_cud_johnson_2020",
    "03_ad_wightman_2021",
    "03_stroke_mishra_2022",
    "03_epiGen_ilae_2022",
    "03_epiFoc_ilae_2022",
    "03_als_vanRheenen_2021",
    "04_well_baselmans_2019",
    "04_neur_baselmans_2019",
    "04_lone_day_2018",
    "04_risk_karlssonlinner_2019",
    "05_ins_jansen_2019",
    "05_sd_dashti_2019",
    "05_chron_jones_2019",
    "05_ds_wang_2019",
    "06_ea_okbay_2022",
    "06_int_savage_2018",
    "06_cp_lee_2018",
    "06_rt_davies_2018",
    "06_mem_davies_2016",
    "07_height_yengo_2022",
    "07_bmi_yengo_2018",
    "07_whr_pulit_2019",
    "08_cad_aragam_2022",
    "08_dbp_evangelou_2018",
    "08_sbp_evangelou_2018",
    "08_myocardial_hartiala_2021",
    "08_diabetes_xue_2018",
]
RGSEL_DOMAINS = {
    "01": "Psychiatric",
    "02": "Substance use",
    "03": "Neurological",
    "04": "Personality",
    "05": "Sleep",
    "06": "Cognition",
    "07": "Anthropometric",
    "08": "Cardiovascular",
}
rs = read("gwama.eur.rgSelection.txt")
rl = read("gwama.eur.rgSelection.labels.txt")
if not len(rs) == len(rl) == len(RGSEL_ORDER):
    raise ValueError("rgSelection tables and RGSEL_ORDER differ in length")
# the label file lists the same rows in the same order
if not (rs.gap_gwm_rg.values == rl.gap_gwm_rg.values).all():
    raise ValueError("rgSelection label file is not aligned with the result table")
rs["label"] = rl["label"].values
rs = rs.set_index("p2")
rgsel = []
for key in RGSEL_ORDER:
    r = rs.loc[key]
    m = re.match(r"^(.*) \((.*)\)$", r.label)
    rgsel.append(
        {
            "key": key,
            "trait": m.group(1),
            "ref": m.group(2),
            "domain": RGSEL_DOMAINS[key[:2]],
            "h2": num(r.h2_obs),
            "h2_se": num(r.h2_obs_se),
            **{
                t: {
                    "rg": num(r[f"gap_{t}_rg"]),
                    "se": num(r[f"gap_{t}_se"]),
                    "z": num(r[f"gap_{t}_z"]),
                    "p": pval(r[f"gap_{t}_p"]),
                    "fdr": pval(r[f"gap_{t}_FDR"]),
                }
                for t in MEASURES
            },
        }
    )
with open(os.path.join(OUT, "rgsel.json"), "w") as f:
    json.dump({"traits": rgsel, "n": 54890}, f, separators=(",", ":"))
summary["rgsel"] = {
    "n": len(rgsel),
    "fdr": {t: sum(1 for x in rgsel if x[t]["fdr"] < 0.05) for t in MEASURES},
    "any": sum(1 for x in rgsel if any(x[t]["fdr"] < 0.05 for t in MEASURES)),
}

# ---------- fastBAT gene-based tests ----------
fb = read("gwama.eur.fastbat.txt")
n_genes = len(fb)
bonf = 0.05 / n_genes
fbk = fb[fb.topFDR.astype(float) < 0.05]
# lead gene of each locus (per model): the gene with GENE_COUNT 1 in that locus
lead = {
    t: dict(
        zip(
            fb.loc[fb[f"gap_{t}_GENE_COUNT"].astype(float) == 1, f"gap_{t}_LOCUS_COUNT"].astype(float).astype(int),
            fb.loc[fb[f"gap_{t}_GENE_COUNT"].astype(float) == 1, "Gene"],
        )
    )
    for t in MEASURES
}
genes = []
for _, r in fbk.iterrows():
    # GFF gene descriptions are URL-encoded (e.g. %2C for a comma)
    genes.append(
        {
            "gene": r.Gene,
            "desc": unquote(dash(r.GENE_DESCRIPTION)) if dash(r.GENE_DESCRIPTION) else None,
            "cyto": r.CYTOBAND,
            "chr": r.Chr,
            "start": int(r.Start),
            "end": int(r.End),
            "nsnp": int(float(r["No.SNPs"])),
            "entrez": dash(r.ENTREZ_ID).replace("GeneID:", "") if dash(r.ENTREZ_ID) else None,
            **{
                t: {
                    "p": pval(r[f"gap_{t}_Pvalue"]),
                    "fdr": pval(r[f"gap_{t}_FDR"]),
                    "locus": int(float(r[f"gap_{t}_LOCUS_COUNT"])),
                    "lead": lead[t].get(int(float(r[f"gap_{t}_LOCUS_COUNT"]))),
                    "index": int(float(r[f"gap_{t}_GENE_COUNT"])) == 1,
                }
                for t in MEASURES
            },
        }
    )
with open(os.path.join(OUT, "genes.json"), "w") as f:
    json.dump({"genes": genes, "nTested": n_genes, "bonf": bonf}, f, separators=(",", ":"))
fbp = {t: fb[f"gap_{t}_Pvalue"].astype(float) for t in MEASURES}
summary["genes"] = {
    "nTested": n_genes,
    "bonf": bonf,
    **{
        t: {
            "bonf": int((fbp[t] < bonf).sum()),
            "fdr": int((fb[f"gap_{t}_FDR"].astype(float) < 0.05).sum()),
            "loci_bonf": int(((fbp[t] < bonf) & (fb[f"gap_{t}_GENE_COUNT"].astype(float) == 1)).sum()),
        }
        for t in MEASURES
    },
    "any_bonf": int(np.logical_or.reduce([fbp[t] < bonf for t in MEASURES]).sum()),
}

# ---------- Mendelian randomization (GSMR plus sensitivity methods, both directions) ----------
# in the table, "gap_<m>_outcome" columns have brain age gap as the outcome (trait -> brain age gap),
# "gap_<m>_exposure" columns have it as the exposure (brain age gap -> trait)
MR_METHODS = [
    ("gsmr", "GSMR"),
    ("ivw", "IVW"),
    ("divw", "Debiased IVW"),
    ("pivw", "Penalized IVW"),
    ("median", "Weighted median"),
    ("egger", "MR-Egger"),
    ("ml", "Maximum likelihood"),
    ("mbe", "Mode-based"),
    ("cm", "Contamination mixture"),
    ("mrpresso", "MR-PRESSO"),
]
mrdf = read("gwama.eur.gsmr.multi.labels.txt")
mr = []
for _, r in mrdf.iterrows():
    m_ = re.match(r"^(.*) \((.*)\)$", r.label)
    ref = m_.group(2).replace("Scptt", "Scott")  # typo in the label file
    entry = {"trait": m_.group(1), "ref": ref}
    for key, col in (("to", "outcome"), ("from", "exposure")):
        entry[key] = {
            t: (
                None
                if not r[f"gap_{t}_{col}_gsmr_beta"]
                else {
                    "nsnp": int(float(r[f"gap_{t}_{col}_n_snps_total"])),
                    "nheidi": int(float(r[f"gap_{t}_{col}_n_snps_HEIDI"])),
                    "b": num(r[f"gap_{t}_{col}_gsmr_beta"]),
                    "se": num(r[f"gap_{t}_{col}_gsmr_se"]),
                    "p": pval(r[f"gap_{t}_{col}_gsmr_p"]),
                    "fdr": pval(r[f"gap_{t}_{col}_gsmr_fdr"]),
                    "methods": {k: pval(r[f"gap_{t}_{col}_{k}_p"]) for k, _ in MR_METHODS[1:]},
                    "n05": int(float(r[f"gap_{t}_{col}_sum_p05"])),
                }
            )
            for t in MEASURES
        }
    mr.append(entry)
with open(os.path.join(OUT, "mr.json"), "w") as f:
    json.dump({"traits": mr, "methods": MR_METHODS, "n": 54890}, f, separators=(",", ":"))
summary["mr"] = {
    "n": len(mr),
    **{
        f"{key}_any": sum(1 for x in mr if any(x[key][t] and x[key][t]["fdr"] < 0.05 for t in MEASURES))
        for key in ("to", "from")
    },
}

meta["summary"] = summary
with open(os.path.join(OUT, "meta.json"), "w") as f:
    json.dump(meta, f, separators=(",", ":"))

print(json.dumps(summary, indent=1))
for fn in sorted(os.listdir(OUT)):
    print(fn, round(os.path.getsize(os.path.join(OUT, fn)) / 1e6, 2), "MB")
