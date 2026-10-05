# Interactive results browser

Static website (GitHub Pages) for the phenome-wide association and genetic correlation results of
Jawinski et al. (2025) *Nature Aging*. Live at https://pjawinski.github.io/ukb_brainage/

```
docs/
├── index.html          page
├── assets/app.js       plot, table, filters, trait panel
├── assets/style.css
├── data/*.json         generated from results/combined (PheWAS, genetic correlations incl. 38 selected traits, heritability, GWAS loci, Mendelian randomization)
└── scripts/build_data.py
```

**Update the data** after changing the result tables:

```
python3 docs/scripts/build_data.py
```

**Preview locally:** `cd docs && python3 -m http.server`, then open http://localhost:8000.

**Deploy:** repository Settings → Pages → Build and deployment → Source: *Deploy from a branch*,
Branch: `main`, folder: `/docs`.

Shareable links keep the view in the URL, e.g.
`https://pjawinski.github.io/ukb_brainage/#m=gwm&s=all&t=2443` opens diabetes for the combined model.
