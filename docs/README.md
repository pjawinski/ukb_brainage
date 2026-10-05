# Interactive results browser

Static website (GitHub Pages) for the results of Jawinski et al. (2025)
*Nature Aging*. Live at <https://pjawinski.github.io/ukb_brainage/>

```text
docs/
├── index.html          page
├── assets/app.js       plots, tables, filters, search, detail panels
├── assets/style.css
├── data/*.json         generated from the result tables (see below)
└── scripts/build_data.py
```

The data files cover SNP heritability, the GWAS loci, gene-based tests,
phenome-wide associations, genetic correlations (UK Biobank traits and
38 selected traits) and Mendelian randomization.

**Update the data** after changing the result tables:

```bash
python3 docs/scripts/build_data.py
```

**Preview locally:** `cd docs && python3 -m http.server`, then open
<http://localhost:8000>.

**Deploy:** repository Settings → Pages → Build and deployment →
Source: *Deploy from a branch*, Branch: `main`, folder: `/docs`.

Shareable links keep the view in the URL, e.g.
`#d=pheno&m=gwm&t=2443` after the site address opens diabetes for the
combined model.
