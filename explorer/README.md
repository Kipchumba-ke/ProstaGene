# ProstaGene Explorer

An interactive companion to the ProstaGene report, with three role-based views:

| View | For | What it shows |
|---|---|---|
| **Learn** | Students and the public | A six-step guided story in plain language: pairing, marker genes, programmes, validation, ancestry, limits |
| **Clinician** | Clinicians | Established markers and their clinical context, subtypes by ancestry, progression models, equity notes, where the evidence stops |
| **Researcher** | Researchers | Gene search across every analysis, a 20,000-point volcano plot, Hallmark comparison across seven models, data downloads |

Design choices: a calm sage, muted-teal and warm-neutral palette; generous spacing; an 18 px base size with
Atkinson Hyperlegible (designed for low-vision readers); body and secondary text at 7:1 or more (WCAG AAA) and links, badges and buttons at 5.5:1 or more (WCAG AA); teal/clay chart colours that
stay distinguishable with common colour-vision deficiencies; keyboard focus outlines; reduced-motion support.

## Data

The app reads JSON from `public/data/`, written by the pipeline:

```bash
snakemake --cores 4 explorer/public/data/summary.json   # runs workflow/scripts/12_export_explorer.R
```

Commit those JSON files: the GitHub Actions workflow builds the site from them. Only summary statistics and
de-identified open-access TCGA values are exported.

If `public/data/` is empty, `npm run mock` writes illustrative data with the same structure, and the app shows a
"Preview data" banner.

## Develop

```bash
cd explorer
npm install
npm run mock      # only if public/data is empty
npm run dev       # http://localhost:5173
npm run build     # type-check + production build to dist/
```

Requires Node 20 or newer.
