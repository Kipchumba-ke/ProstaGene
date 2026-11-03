import type { AppData } from "../data";
import { num } from "../data";
import { Metric } from "../components/Layout";

const icons = {
  learn: <svg width="26" height="26" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" aria-hidden="true"><path d="M2 7l10-4 10 4-10 4L2 7z" /><path d="M6 9v5c0 1.7 2.7 3 6 3s6-1.3 6-3V9" /></svg>,
  clinician: <svg width="26" height="26" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" aria-hidden="true"><path d="M6 3v6a4 4 0 008 0V3" /><path d="M10 13v3a5 5 0 0010 0v-2" /><circle cx="20" cy="12" r="2" /></svg>,
  researcher: <svg width="26" height="26" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" aria-hidden="true"><circle cx="11" cy="11" r="7" /><path d="M21 21l-4.3-4.3" /><path d="M8 13l2-3 2 2 2-4" /></svg>,
};

export default function Home({ data }: { data: AppData }) {
  const s = data.summary;
  const roles = [
    { id: "learn", title: "I'm learning about this", body: "A guided walk through what changes when prostate tissue becomes cancer, in plain language, step by step.", icon: icons.learn },
    { id: "clinician", title: "I'm a clinician", body: "Established markers, molecular subtypes, what predicts progression, and where the evidence stops.", icon: icons.clinician },
    { id: "researcher", title: "I'm a researcher", body: "Search any gene, explore the volcano plot, compare pathway results across analyses and models.", icon: icons.researcher },
  ];
  return (
    <>
      <section className="hero">
        <div className="wrap">
          <span className="eyebrow">Open prostate cancer genomics</span>
          <h1 style={{ maxWidth: "18ch" }}>What changes when prostate tissue becomes cancer?</h1>
          <p className="lead" style={{ maxWidth: "60ch" }}>
            ProstaGene compares tumors with healthy tissue from the same {s.design.pairs} patients in The Cancer Genome Atlas,
            checks the findings in an independent group, and asks how tissue make-up, molecular subtype and ancestry shape the picture.
          </p>
          <p className="muted" style={{ marginTop: 28 }}>Choose the view that fits you:</p>
          <div className="grid g3" style={{ marginTop: 12 }}>
            {roles.map((r) => (
              <a key={r.id} href={`#/${r.id}`} className="card role-card">
                <span className="icon-dot">{r.icon}</span>
                <h3 style={{ margin: 0 }}>{r.title}</h3>
                <p className="muted" style={{ margin: 0 }}>{r.body}</p>
                <span className="go">Open view →</span>
              </a>
            ))}
          </div>
        </div>
      </section>
      <section className="section" style={{ background: "var(--surface-sage)" }}>
        <div className="wrap">
          <span className="eyebrow">At a glance</span>
          <div className="grid g4">
            <Metric value={s.design.pairs} label="patients with tumor and healthy tissue" />
            <Metric value={s.de.degs.toLocaleString()} label="genes clearly changed in tumors" />
            <Metric value={`${s.pathways.hallmark_sig} / 50`} label="biological programmes altered" />
            <Metric value={`${Math.round(100 * num(s.validation, "sign_agreement_tcga_DEGs"))}%`} label="of changes confirmed in an independent cohort" />
          </div>
        </div>
      </section>
    </>
  );
}
