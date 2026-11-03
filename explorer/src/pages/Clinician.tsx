import { useState } from "react";
import type { AppData } from "../data";
import { fmt, fmtP, fold, num, prettyScore } from "../data";
import { Section, Metric } from "../components/Layout";
import { KMPlot, StackedBars } from "../components/Charts";

const MARKERS: { gene: string; context: string }[] = [
  { gene: "AMACR", context: "Immunohistochemistry marker (P504S) used to support a cancer diagnosis on biopsy." },
  { gene: "PCA3", context: "Measured in post-DRE urine assays." },
  { gene: "DLX1", context: "Measured, with HOXC6, in a urine mRNA assay." },
  { gene: "HPN", context: "Hepsin; consistently overexpressed in prostate adenocarcinoma." },
  { gene: "TDRD1", context: "ERG target; marks ERG-fusion-positive tumors." },
  { gene: "KRT5", context: "Basal-cell keratin; loss of basal cells supports malignancy on IHC." },
  { gene: "TP63", context: "Basal-cell marker (p63); absent in invasive carcinoma glands." },
];
const SUB_COLORS: Record<string, string> = { ERG: "#5e8c7a", ETS_other: "#8fb8a8", SPOP: "#c99a6b", FOXA1: "#e0c49c", Other: "#c9c4b8" };

export default function Clinician({ data }: { data: AppData }) {
  const s = data.summary, G = data.genes;
  const [score, setScore] = useState(data.survival.scores.includes("hm_E2F_TARGETS") ? "hm_E2F_TARGETS" : data.survival.scores[0]);
  const adj = data.survival.cox.filter((c) => c.model === "adjusted");
  const uni = data.survival.cox.filter((c) => c.model === "univariable" && c.padj < 0.05);
  const nAdjSig = adj.filter((c) => c.padj < 0.05).length;
  const ergLow = num(s.subtypes, "erg_pct_afr_low"), ergHigh = num(s.subtypes, "erg_pct_afr_high");

  return (
    <>
      <section className="hero" style={{ paddingBottom: 32 }}>
        <div className="wrap">
          <span className="eyebrow">Clinician view</span>
          <h1>Prostate adenocarcinoma: what the tissue shows</h1>
          <p className="lead">Findings from {s.design.pairs} matched tumor–normal pairs and {num(s.survival, "tumors_with_pfi")} primary tumors with follow-up in TCGA-PRAD.</p>
          <div className="notice" style={{ marginTop: 24 }} role="note">
            <strong>Research summary.</strong>
            <span>Population-level associations from public datasets. Not validated for diagnosis, prognosis or treatment decisions in individual patients.</span>
          </div>
        </div>
      </section>

      <Section eyebrow="Markers" title="Established markers behave as expected" intro="Fold change is the average tumor-vs-normal difference within patients; ‘confirmed’ means the same direction and FDR < 0.05 in the independent GEO cohort.">
        <div className="card" style={{ overflowX: "auto" }}>
          <table className="table">
            <thead><tr><th>Gene</th><th>In tumor</th><th className="num">FDR</th><th>Independent cohort</th><th>Clinical context</th></tr></thead>
            <tbody>
              {MARKERS.map(({ gene, context }) => {
                const i = G.symbol.indexOf(gene); if (i < 0) return null;
                const l = G.lfc[i] ?? 0, ok = (G.geo_padj[i] ?? 1) < 0.05 && Math.sign(G.geo_lfc[i] ?? 0) === Math.sign(l);
                return (
                  <tr key={gene}>
                    <td><strong>{gene}</strong></td>
                    <td><span className={`pill ${l > 0 ? "up" : "down"}`}>{fold(l)} {l > 0 ? "higher" : "lower"}</span></td>
                    <td className="num">{fmtP(G.padj[i])}</td>
                    <td>{ok ? <span className="pill ok">Confirmed</span> : <span className="pill ns">Not confirmed</span>}</td>
                    <td className="muted">{context}</td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
      </Section>

      <Section tone="sage" eyebrow="Molecular subtypes" title="ERG-fusion tumors and ancestry"
        intro={`TCGA's 2015 classes group tumors by mutually exclusive drivers. ERG-positive tumors were ${Number.isFinite(ergHigh) ? `${ergHigh}%` : "–"} among patients with majority African ancestry versus ${Number.isFinite(ergLow) ? `${ergLow}%` : "–"} in others.`}>
        <div className="grid g2">
          <div className="card">
            <StackedBars groups={["AFR < 0.5", "AFR >= 0.5"]} colors={SUB_COLORS}
              data={data.subtypes.by_ancestry.map((d) => ({ group: d.afr_group, key: d.subtype_group, n: d.N }))} />
            <p className="muted small" style={{ marginTop: 16 }}>AFR = genetic African-ancestry fraction (Carrot-Zhang et al. 2020).</p>
          </div>
          <div className="stack">
            {data.subtypes.cox.slice(0, 1).map((c) => (
              <Metric key={c.term} value={`HR ${fmt(c.HR)}`} label={<>ERG-positive vs negative for progression, adjusted for Gleason, T stage and age (95% CI {fmt(c.lower)}–{fmt(c.upper)}, p = {fmtP(c.p)})</>} />
            ))}
            <p className="muted">Subtype frequencies differ by ancestry; group sizes for SPOP and FOXA1 tumors, and for patients of African ancestry, are small.</p>
          </div>
        </div>
      </Section>

      <Section eyebrow="Progression" title="Does tumor biology add to Gleason grade and stage?"
        intro={`Progression-free interval (TCGA-CDR), ${num(s.survival, "pfi_events")} events, median follow-up ${num(s.survival, "median_followup_months")} months. Clinical model (Gleason group, T stage, age) C-index ${fmt(num(s.survival, "c_index_clinical"), 3)}.`}>
        <div className="grid g2">
          <div className="card">
            <label htmlFor="km" className="muted small">Score</label>
            <div style={{ margin: "8px 0 16px" }}>
              <select id="km" className="input" value={score} onChange={(e) => setScore(e.target.value)}>
                {data.survival.scores.map((sc) => <option key={sc} value={sc}>{prettyScore(sc)}</option>)}
              </select>
            </div>
            <KMPlot rows={data.survival.rows} score={score} />
          </div>
          <div className="stack">
            <div className="card sage">
              <h3>{nAdjSig === 0 ? "No score adds independent information" : `${nAdjSig} score(s) remain after adjustment`}</h3>
              <p style={{ margin: 0 }}>
                {uni.length > 0 && <>Unadjusted, {uni.map((u) => `${prettyScore(u.score)} (HR ${fmt(u.HR)} per SD)`).join(", ")} predicted faster progression. </>}
                {nAdjSig === 0 && <>Once Gleason group, T stage and age are included, hazard ratios fall to about 1: the signal is already captured by grade and stage.</>}
              </p>
            </div>
            <p className="muted small">Curves split each score at its median and are unadjusted; the Cox models are the formal test.</p>
          </div>
        </div>
      </Section>

      <Section tone="sage" eyebrow="Equity" title="Ancestry: what holds up">
        <div className="grid g3">
          <div className="card"><h3>Representation</h3><p style={{ margin: 0 }}>{num(s.ancestry, "tcga_tumors_black")} of {num(s.ancestry, "tcga_tumors_black") + num(s.ancestry, "tcga_tumors_white")} analysed tumors came from Black patients. Self-reported race and measured genetic ancestry disagreed for a few patients.</p></div>
          <div className="card"><h3>Interferon signal</h3><p style={{ margin: 0 }}>Higher interferon activity in tumors from Black patients follows immune-cell content; per patient the shift is modest (p = {fmtP(num(s.genetic_ancestry, "ifng_afr_p"))}).</p></div>
          <div className="card"><h3>Inherited signals</h3><p style={{ margin: 0 }}>{num(s.genetic_ancestry, "likely_germline_genes")} of {num(s.genetic_ancestry, "top_genes_tested")} top race-associated genes track genetic ancestry in normal tissue or within Black patients, pointing to inherited variation rather than tumor biology.</p></div>
        </div>
      </Section>

      <Section eyebrow="Limits" title="Where the evidence stops">
        <ul className="lead" style={{ fontSize: "1.1rem" }}>
          <li>Bulk RNA from surgical specimens; composition (stroma, immune cells) shapes many signals.</li>
          <li>Follow-up is short and outcomes are progression events, not prostate-cancer deaths.</li>
          <li>Few patients of African ancestry; results for them carry wide uncertainty.</li>
          <li>No external outcome cohort was used to validate prognostic scores.</li>
        </ul>
      </Section>
    </>
  );
}
