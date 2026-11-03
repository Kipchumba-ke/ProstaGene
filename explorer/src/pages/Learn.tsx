import { useState, type ReactNode } from "react";
import type { AppData } from "../data";
import { fold, num } from "../data";
import { NesBars, PairedDots } from "../components/Charts";

const MARKERS: { gene: string; note: string }[] = [
  { gene: "AMACR", note: "Pathologists stain for this protein to help recognise prostate cancer in biopsies." },
  { gene: "PCA3", note: "An RNA that urine tests for prostate cancer measure." },
  { gene: "HPN", note: "Hepsin, a cell-surface enzyme often raised in prostate cancer." },
  { gene: "DLX1", note: "A developmental gene also measured in a urine test for prostate cancer." },
  { gene: "TDRD1", note: "Switched on mainly in tumors that carry an ERG gene fusion." },
];

const FRIENDLY: Record<string, string> = {
  MYC_TARGETS_V1: "Growth driven by the MYC gene", E2F_TARGETS: "Copying DNA before cells divide",
  G2M_CHECKPOINT: "Cell division", ANDROGEN_RESPONSE: "Response to male hormones (androgens)",
  MYOGENESIS: "Muscle-like supporting tissue", EPITHELIAL_MESENCHYMAL_TRANSITION: "Connective-tissue programme",
  APICAL_JUNCTION: "Cell-to-cell contacts", HYPOXIA: "Response to low oxygen",
};

function Step({ n, title, children }: { n: number; title: string; children: ReactNode }) {
  return (
    <section className="section">
      <div className="wrap step">
        <span className="step-num" aria-hidden="true">{n}</span>
        <div><h2>{title}</h2>{children}</div>
      </div>
    </section>
  );
}

export default function Learn({ data }: { data: AppData }) {
  const s = data.summary;
  const available = MARKERS.filter((m) => data.expression.genes[m.gene]);
  const [gene, setGene] = useState(available[0]?.gene ?? "");
  const gi = data.genes.symbol.indexOf(gene);
  const paths = Object.keys(FRIENDLY)
    .map((k) => data.hallmark.find((p) => p.pathway === `HALLMARK_${k}`))
    .filter((p): p is NonNullable<typeof p> => !!p)
    .sort((a, b) => b.NES - a.NES);
  const smooth = data.composition.find((c) => c.score === "smooth_muscle");
  const black = num(s.ancestry, "tcga_tumors_black"), white = num(s.ancestry, "tcga_tumors_white");
  const ergLow = num(s.subtypes, "erg_pct_afr_low"), ergHigh = num(s.subtypes, "erg_pct_afr_high");

  return (
    <>
      <section className="hero" style={{ paddingBottom: 24 }}>
        <div className="wrap">
          <span className="eyebrow">Learn · guided story</span>
          <h1>How prostate cancer changes the tissue it grows in</h1>
          <p className="lead">Six short steps. No background needed: each idea is explained as it comes up.</p>
        </div>
      </section>

      <Step n={1} title="We compared each tumor with healthy tissue from the same person">
        <p>
          The Cancer Genome Atlas (TCGA) collected prostate tumors removed during surgery. For {s.design.pairs} men, it also kept a
          piece of nearby healthy prostate. Comparing the two samples from the <em>same</em> person removes differences between people,
          such as age or genetics, so what remains is mostly the effect of the cancer.
        </p>
        <p className="muted">We measured the activity of about {Math.round(s.design.genes_tested / 1000)},000 genes in all {s.design.samples} samples using RNA sequencing.</p>
      </Step>

      <Step n={2} title="Some genes switch on strongly in cancer">
        <p>{s.de.degs.toLocaleString()} genes changed at least two-fold. The strongest are well-known prostate cancer markers, a good sign the analysis works.</p>
        <div className="chips" role="group" aria-label="Choose a marker gene" style={{ margin: "20px 0" }}>
          {available.map((m) => (
            <button key={m.gene} className="chip" aria-pressed={gene === m.gene} onClick={() => setGene(m.gene)}>{m.gene}</button>
          ))}
        </div>
        {gene && (
          <div className="grid g2" style={{ alignItems: "center" }}>
            <div className="card"><PairedDots gene={gene} values={data.expression.genes[gene]} samples={data.expression.samples} /></div>
            <div className="stack">
              <h3>{gene}</h3>
              <p>{MARKERS.find((m) => m.gene === gene)?.note}</p>
              {gi >= 0 && <p className="lead">About <strong>{fold(data.genes.lfc[gi])}</strong> {(data.genes.lfc[gi] ?? 0) > 0 ? "higher" : "lower"} in tumors on average.</p>}
            </div>
          </div>
        )}
      </Step>

      <Step n={3} title="Cancer cells grow faster, and the supporting tissue shrinks">
        <p>Genes work in teams called <strong>programmes</strong>. Teal bars are programmes more active in tumors; brown bars are less active.</p>
        <div className="card" style={{ marginTop: 20 }}>
          <NesBars rows={paths.map((p) => ({ label: FRIENDLY[p.pathway.replace("HALLMARK_", "")], value: p.NES, sig: p.padj < 0.05 }))} />
        </div>
        <p className="muted" style={{ marginTop: 16 }}>
          Healthy prostate is rich in muscle-like support tissue{smooth ? ` — its score is lower in the tumor for ${smooth.pairs - smooth.tumor_higher_pairs} of ${smooth.pairs} people` : ""}.
          Part of the "less active" signal is simply that tumors contain less of this tissue.
        </p>
      </Step>

      <Step n={4} title="We checked the results in a different group of patients">
        <p>
          A separate study (GEO GSE104131) sequenced tumors and healthy tissue from {num(s.validation, "complete_pairs")} other men.
          <strong> {Math.round(100 * num(s.validation, "sign_agreement_tcga_DEGs"))}%</strong> of the changes we found pointed the same way there, and
          a score built from our top genes told tumor from healthy tissue almost perfectly.
        </p>
      </Step>

      <Step n={5} title="Who is in the data matters">
        <p>
          Men of African ancestry are more often diagnosed with prostate cancer and more often die from it, yet they are a small
          share of large genomic studies. Here, {black} tumors came from Black patients and {white} from White patients.
        </p>
        <div className="grid g2" style={{ marginTop: 16 }}>
          <div className="card sage">
            <h3>Subtypes and ancestry</h3>
            <p style={{ margin: 0 }}>Tumors carrying an <strong>ERG gene fusion</strong> were {Number.isFinite(ergHigh) ? `${ergHigh}%` : "less common"} among patients with mostly African ancestry, against {Number.isFinite(ergLow) ? `${ergLow}%` : "–"} in others{ergHigh < ergLow ? ", in line with studies from the US, Ghana and South Africa" : ". Studies from the US, Ghana and South Africa have reported lower rates in men of African ancestry"}.</p>
          </div>
          <div className="card sand">
            <h3>Immune signals, read carefully</h3>
            <p style={{ margin: 0 }}>Tumors from Black patients showed somewhat stronger interferon (immune) activity. That signal travels with the number of immune cells in the sample, and at the level of individual patients it is modest.</p>
          </div>
        </div>
      </Step>

      <Step n={6} title="What this does not tell you">
        <ul>
          <li>These are averages across groups of patients. They cannot predict what will happen to any one person.</li>
          <li>Race recorded in a hospital form is a social category; genetic ancestry was measured separately and used as a check.</li>
          <li>Bulk tissue mixes many cell types, so some changes reflect what the sample is made of rather than what cancer cells do.</li>
        </ul>
        <div className="card sage" style={{ marginTop: 24 }}>
          <h3>Words used here</h3>
          <dl className="grid g2" style={{ margin: 0 }}>
            {[
              ["Gene expression", "How actively a gene is being read to make RNA."],
              ["Fold change", "How many times higher (or lower) a gene is in tumor than in healthy tissue."],
              ["Programme (pathway)", "A team of genes that work together on one job, such as cell division."],
              ["Gene fusion", "Two genes broken and joined, often switching one of them on."],
            ].map(([t, d]) => <div key={t}><dt><strong>{t}</strong></dt><dd style={{ margin: 0 }} className="muted">{d}</dd></div>)}
          </dl>
        </div>
      </Step>
    </>
  );
}
