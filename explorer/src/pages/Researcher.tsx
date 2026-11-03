import { useMemo, useState } from "react";
import type { AppData, Pathway } from "../data";
import { fmt, fmtP, fold, num } from "../data";
import { Section, Metric } from "../components/Layout";
import { NesBars, PairedDots, Volcano } from "../components/Charts";

type Analysis = { id: string; label: string; nes: keyof Pathway; padj: keyof Pathway | null; note: string };
const ANALYSES: Analysis[] = [
  { id: "main", label: "Tumor vs normal", nes: "NES", padj: "padj", note: "Paired DESeq2, ~ patient + condition; fgsea on the Wald statistic." },
  { id: "adj", label: "Composition-adjusted", nes: "NES_adj", padj: "padj_adj", note: "Adds ESTIMATE stromal + immune, smooth-muscle and basal scores to the paired model." },
  { id: "geo", label: "GEO validation", nes: "NES_geo", padj: "padj_geo", note: "Same model in GSE104131 (13 high-grade pairs)." },
  { id: "race", label: "Black vs White tumors", nes: "NES_race", padj: "padj_race", note: "Tumor-only, self-reported race, adjusted for Gleason group and age." },
  { id: "afr", label: "Genetic ancestry (AFR)", nes: "NES_afr", padj: "padj_afr", note: "Tumor-only, per unit African-ancestry fraction; gene-set p-values overstate certainty for co-regulated sets." },
  { id: "ergpos", label: "ERG+ tumor effect", nes: "NES_ergpos", padj: null, note: "Tumor vs normal within ERG-positive pairs (descriptive)." },
  { id: "ergneg", label: "ERG− tumor effect", nes: "NES_ergneg", padj: null, note: "Tumor vs normal within ERG-negative pairs (descriptive)." },
];

export default function Researcher({ data }: { data: AppData }) {
  const G = data.genes, s = data.summary;
  const [query, setQuery] = useState("");
  const [gene, setGene] = useState(s.de.top[0] ?? G.symbol[0]);
  const [aid, setAid] = useState("main");
  const a = ANALYSES.find((x) => x.id === aid)!;

  const matches = useMemo(() => {
    const q = query.trim().toUpperCase(); if (q.length < 2) return [];
    const starts = G.symbol.filter((g) => g.toUpperCase().startsWith(q));
    return (starts.length ? starts : G.symbol.filter((g) => g.toUpperCase().includes(q))).slice(0, 8);
  }, [query, G.symbol]);
  const i = G.symbol.indexOf(gene);
  const pick = (g: string) => { setGene(g); setQuery(""); };

  const rows = useMemo(() => [...data.hallmark]
    .filter((p) => p[a.nes] !== null)
    .sort((x, y) => Math.abs(Number(y[a.nes])) - Math.abs(Number(x[a.nes])))
    .slice(0, 16)
    .map((p) => ({ label: p.label, value: p[a.nes] as number, sig: a.padj ? Number(p[a.padj]) < 0.05 : undefined })), [data.hallmark, a]);

  const stat = (label: string, lfc: number | null, padj: number | null, unit = "log2FC") => (
    <tr><td>{label}</td><td className="num">{fmt(lfc)} <span className="muted small">{unit}</span></td><td className="num">{fmtP(padj)}</td></tr>
  );

  return (
    <>
      <section className="hero" style={{ paddingBottom: 32 }}>
        <div className="wrap">
          <span className="eyebrow">Researcher view</span>
          <h1>Explore genes and pathways</h1>
          <p className="lead">{s.design.genes_tested.toLocaleString()} genes tested in {s.design.pairs} matched pairs; results across validation, composition, ancestry and subtype models.</p>
          <div className="grid g4" style={{ marginTop: 32 }}>
            <Metric value={s.de.up.toLocaleString()} label="up in tumor (FDR < 0.05, ≥ 2-fold)" />
            <Metric value={s.de.down.toLocaleString()} label="down in tumor" />
            <Metric value={fmt(num(s.validation, "spearman_lfc_all"))} label="fold-change ρ with GEO" />
            <Metric value={fmt(num(s.genetic_ancestry, "nes_r_afr_vs_race"))} label="pathway r, genetic vs self-reported ancestry" />
          </div>
        </div>
      </section>

      <Section eyebrow="Gene lookup" title="Search a gene">
        <div style={{ position: "relative", maxWidth: 520 }}>
          <label htmlFor="q" className="muted small">Gene symbol</label>
          <input id="q" className="input" placeholder="e.g. AMACR, ERG, MYC" value={query} autoComplete="off"
            onChange={(e) => setQuery(e.target.value)}
            onKeyDown={(e) => { if (e.key === "Enter" && matches[0]) pick(matches[0]); }}
            aria-describedby="q-help" />
          <span id="q-help" className="muted small">Type at least two letters, then choose a match or press Enter.</span>
          {matches.length > 0 && (
            <ul role="listbox" className="card" style={{ position: "absolute", zIndex: 5, left: 0, right: 0, listStyle: "none", padding: 8, margin: "6px 0 0" }}>
              {matches.map((m) => (
                <li key={m}><button className="chip" style={{ width: "100%", textAlign: "left", border: "none" }} onClick={() => pick(m)}>{m}</button></li>
              ))}
            </ul>
          )}
        </div>

        {i >= 0 && (
          <div className="grid g2" style={{ marginTop: 32, alignItems: "start" }}>
            <div className="card">
              <div style={{ display: "flex", alignItems: "center", gap: 14, flexWrap: "wrap" }}>
                <h3 style={{ margin: 0, fontSize: "1.8rem" }}>{gene}</h3>
                <span className={`pill ${G.dir[i]}`}>{G.dir[i] === "up" ? `${fold(G.lfc[i])} higher in tumor` : G.dir[i] === "down" ? `${fold(G.lfc[i])} lower in tumor` : "Not a DEG"}</span>
                {G.germline.includes(gene) && <span className="pill ns">Likely inherited (germline) signal</span>}
              </div>
              <p className="muted small" style={{ marginTop: 8 }}>{G.type[i]?.replace("_", " ")}</p>
              <table className="table" style={{ marginTop: 12 }}>
                <thead><tr><th>Analysis</th><th className="num">Effect</th><th className="num">FDR</th></tr></thead>
                <tbody>
                  {stat("TCGA tumor vs normal (shrunk)", G.lfc[i], G.padj[i])}
                  {stat("GEO GSE104131 validation", G.geo_lfc[i], G.geo_padj[i])}
                  {stat("Black vs White tumors", G.race_lfc[i], G.race_padj[i])}
                  {stat("Per unit African ancestry", G.afr_lfc[i], G.afr_padj[i])}
                </tbody>
              </table>
            </div>
            <div className="card">
              {data.expression.genes[gene]
                ? <PairedDots gene={gene} values={data.expression.genes[gene]} samples={data.expression.samples} />
                : <p className="muted">Per-sample expression is exported for the top {Object.keys(data.expression.genes).length} genes only. Statistics on the left cover all tested genes.</p>}
            </div>
          </div>
        )}
      </Section>

      <Section tone="sage" eyebrow="All genes" title="Volcano plot" intro="Tumor vs matched normal. Hover to identify a gene and click to open it above.">
        <div className="card"><Volcano genes={G} selected={gene} onSelect={(g) => { setGene(g); document.getElementById("q")?.scrollIntoView({ behavior: "smooth", block: "center" }); }} /></div>
      </Section>

      <Section eyebrow="Pathways" title="Compare Hallmark results across analyses">
        <div className="chips" role="group" aria-label="Choose analysis">
          {ANALYSES.map((x) => <button key={x.id} className="chip" aria-pressed={aid === x.id} onClick={() => setAid(x.id)}>{x.label}</button>)}
        </div>
        <p className="muted" style={{ marginTop: 16 }}>{a.note} Showing the 16 sets with the largest |NES|{a.padj ? "; faded bars have FDR ≥ 0.05" : ""}.</p>
        <div className="card" style={{ marginTop: 16 }}><NesBars rows={rows} /></div>
      </Section>

      <Section tone="sage" eyebrow="Data" title="Download the exported results">
        <div className="chips">
          {["summary", "genes", "hallmark", "survival", "subtypes", "composition"].map((f) => (
            <a key={f} className="chip" href={`${import.meta.env?.BASE_URL ?? "./"}data/${f}.json`} download>{f}.json</a>
          ))}
        </div>
        <p className="muted" style={{ marginTop: 16 }}>Full tables, figures and code are in the report and the GitHub repository.</p>
      </Section>
    </>
  );
}
