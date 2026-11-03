import type { ReactNode } from "react";

export type Route = "home" | "learn" | "clinician" | "researcher";
const ROLES: { id: Route; label: string }[] = [
  { id: "learn", label: "Learn" },
  { id: "clinician", label: "Clinician" },
  { id: "researcher", label: "Researcher" },
];

export function Header({ route }: { route: Route }) {
  return (
    <header className="topbar">
      <div className="wrap">
        <a className="brand" href="#/"><span className="brand-mark" aria-hidden="true" />ProstaGene Explorer</a>
        <nav className="roles" aria-label="Choose a view">
          {ROLES.map((r) => (
            <a key={r.id} className="role-link" href={`#/${r.id}`} aria-current={route === r.id ? "page" : undefined}>{r.label}</a>
          ))}
        </nav>
      </div>
    </header>
  );
}

export function MockBanner({ generated }: { generated: string }) {
  if (generated !== "mock") return null;
  return (
    <div className="wrap" style={{ paddingTop: 24 }}>
      <div className="notice mock" role="status">
        <strong>Preview data.</strong>
        <span>This build is using illustrative mock data. Run the pipeline's export step to show the real results.</span>
      </div>
    </div>
  );
}

export function Footer({ generated }: { generated: string }) {
  return (
    <footer className="footer">
      <div className="wrap grid g2">
        <div>
          <strong>ProstaGene</strong> · An open, reproducible analysis of TCGA prostate adenocarcinoma (TCGA-PRAD),
          validated in GEO GSE104131. Built by Kipchumba.
          <p className="small" style={{ marginTop: 10 }}>
            For research and education only. These results describe groups of patients in public datasets and must not be used to make decisions about any individual's care.
          </p>
        </div>
        <div className="links">
          <a href="../">Read the full report</a>
          <span className="small">Data export: {generated === "mock" ? "mock preview" : generated}</span>
        </div>
      </div>
    </footer>
  );
}

export function Section({ id, eyebrow, title, intro, children, tone }: {
  id?: string; eyebrow?: string; title: string; intro?: ReactNode; children?: ReactNode; tone?: "sage";
}) {
  return (
    <section id={id} className="section" style={tone === "sage" ? { background: "var(--surface-sage)" } : undefined}>
      <div className="wrap">
        {eyebrow && <span className="eyebrow">{eyebrow}</span>}
        <h2>{title}</h2>
        {intro && <p className="lead">{intro}</p>}
        <div style={{ marginTop: 32 }}>{children}</div>
      </div>
    </section>
  );
}

export function Metric({ value, label }: { value: ReactNode; label: ReactNode }) {
  return (
    <div className="card metric">
      <div className="value">{value}</div>
      <div className="label">{label}</div>
    </div>
  );
}
