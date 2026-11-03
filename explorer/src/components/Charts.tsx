import { useEffect, useMemo, useRef, useState, type MouseEvent } from "react";
import { scaleLinear, scaleBand } from "d3-scale";
import { line, curveStepAfter } from "d3-shape";
import type { Genes } from "../data";

const UP = "#3e7a85", DOWN = "#a8622f", NS = "#c9c4b8", NORMAL = "#9fb7ad", TUMOR = "#3f6b5b";

/* ---------- Paired tumor vs normal dots ---------- */
export function PairedDots({ values, samples, gene, height = 300 }: {
  values: number[]; samples: { patient: string; condition: string }[]; gene: string; height?: number;
}) {
  const W = 420, H = height, m = { t: 20, r: 24, b: 44, l: 52 };
  const byPatient = useMemo(() => {
    const map = new Map<string, { normal?: number; tumor?: number }>();
    samples.forEach((s, i) => {
      const e = map.get(s.patient) ?? {};
      e[s.condition as "normal" | "tumor"] = values[i];
      map.set(s.patient, e);
    });
    return [...map.values()].filter((e) => e.normal !== undefined && e.tumor !== undefined) as { normal: number; tumor: number }[];
  }, [values, samples]);
  const lo = Math.min(...values), hi = Math.max(...values);
  const y = scaleLinear().domain([lo - 0.3, hi + 0.3]).range([H - m.b, m.t]).nice();
  const x = scaleBand().domain(["normal", "tumor"]).range([m.l, W - m.r]).padding(0.5);
  const higher = byPatient.filter((p) => p.tumor > p.normal).length;
  return (
    <figure style={{ margin: 0 }}>
      <svg className="chart" viewBox={`0 0 ${W} ${H}`} width="100%" role="img"
        aria-label={`${gene} expression in matched normal and tumor samples; higher in tumor in ${higher} of ${byPatient.length} patients`}>
        {y.ticks(5).map((t) => (
          <g key={t}><line x1={m.l} x2={W - m.r} y1={y(t)} y2={y(t)} stroke="#eee9df" />
            <text x={m.l - 10} y={y(t) + 4} textAnchor="end">{t}</text></g>
        ))}
        {byPatient.map((p, i) => (
          <line key={i} x1={x("normal")! + x.bandwidth() / 2} x2={x("tumor")! + x.bandwidth() / 2}
            y1={y(p.normal)} y2={y(p.tumor)} stroke={p.tumor > p.normal ? "#b8cfc6" : "#e2cdbb"} strokeWidth={1.2} />
        ))}
        {byPatient.flatMap((p, i) => [
          <circle key={`n${i}`} cx={x("normal")! + x.bandwidth() / 2} cy={y(p.normal)} r={4} fill={NORMAL} />,
          <circle key={`t${i}`} cx={x("tumor")! + x.bandwidth() / 2} cy={y(p.tumor)} r={4} fill={TUMOR} />,
        ])}
        {(["normal", "tumor"] as const).map((c) => (
          <text key={c} x={x(c)! + x.bandwidth() / 2} y={H - 14} textAnchor="middle" style={{ fontSize: 15, fill: "#1e2b27", fontWeight: 700 }}>
            {c === "normal" ? "Healthy tissue" : "Tumor"}
          </text>
        ))}
        <text transform={`translate(16 ${(H - m.b + m.t) / 2}) rotate(-90)`} textAnchor="middle">Expression (VST)</text>
      </svg>
      <figcaption className="muted small">Each line is one person. Higher in the tumor for {higher} of {byPatient.length} people.</figcaption>
    </figure>
  );
}

/* ---------- Diverging NES bars ---------- */
export function NesBars({ rows }: { rows: { label: string; value: number | null; sig?: boolean }[] }) {
  const W = 720, rowH = 34, m = { l: 300, r: 40, t: 8, b: 30 };
  // headroom so value labels never collide with the pathway names
  const max = Math.max(2, Math.ceil(Math.max(...rows.map((r) => Math.abs(r.value ?? 0))) + 0.8));
  const H = m.t + m.b + rows.length * rowH;
  const x = scaleLinear().domain([-max, max]).range([m.l, W - m.r]);
  return (
    <svg className="chart" viewBox={`0 0 ${W} ${H}`} width="100%" role="img" aria-label="Pathway enrichment scores">
      <line x1={x(0)} x2={x(0)} y1={m.t} y2={H - m.b} stroke="#b9b2a4" />
      {rows.map((r, i) => {
        const v = r.value ?? 0, yy = m.t + i * rowH;
        return (
          <g key={r.label}>
            <text x={m.l - 14} y={yy + rowH / 2 + 5} textAnchor="end" style={{ fontSize: 14, fill: "#1e2b27" }}>{r.label}</text>
            <rect x={Math.min(x(0), x(v))} y={yy + 7} width={Math.abs(x(v) - x(0))} height={rowH - 14} rx={6}
              fill={v > 0 ? UP : DOWN} opacity={r.sig === false ? 0.35 : 1} />
            <text x={v >= 0 ? x(v) + 8 : x(v) - 8} y={yy + rowH / 2 + 5} textAnchor={v >= 0 ? "start" : "end"}>{r.value === null ? "–" : v.toFixed(2)}</text>
          </g>
        );
      })}
      {[-max, -max / 2, 0, max / 2, max].map((t) => (
        <text key={t} x={x(t)} y={H - 8} textAnchor="middle">{t}</text>
      ))}
    </svg>
  );
}

/* ---------- Volcano (canvas, 20k points) ---------- */
export function Volcano({ genes, selected, onSelect }: { genes: Genes; selected?: string; onSelect: (g: string) => void }) {
  const ref = useRef<HTMLCanvasElement>(null);
  const [hover, setHover] = useState<number | null>(null);
  const W = 760, H = 460, m = { t: 20, r: 20, b: 50, l: 60 };
  const pts = useMemo(() => genes.symbol.map((s, i) => ({
    s, x: genes.lfc[i] ?? 0, y: -Math.log10(Math.max(genes.padj[i] ?? 1, 1e-300)), d: genes.dir[i],
  })), [genes]);
  const xMax = Math.min(8, Math.max(...pts.map((p) => Math.abs(p.x)))) + 0.5;
  const yMax = Math.min(100, Math.max(...pts.map((p) => p.y))) + 2;
  const x = scaleLinear().domain([-xMax, xMax]).range([m.l, W - m.r]);
  const y = scaleLinear().domain([0, yMax]).range([H - m.b, m.t]);
  const px = (p: { x: number; y: number }) => [x(Math.max(-xMax, Math.min(xMax, p.x))), y(Math.min(yMax, p.y))];

  useEffect(() => {
    const c = ref.current; if (!c) return;
    const dpr = window.devicePixelRatio || 1;
    c.width = W * dpr; c.height = H * dpr;
    const ctx = c.getContext("2d")!; ctx.scale(dpr, dpr); ctx.clearRect(0, 0, W, H);
    ctx.strokeStyle = "#e7e1d6"; ctx.fillStyle = "#46524e"; ctx.font = "13px system-ui"; ctx.textAlign = "center";
    for (const t of x.ticks(8)) { ctx.beginPath(); ctx.moveTo(x(t), m.t); ctx.lineTo(x(t), H - m.b); ctx.stroke(); ctx.fillText(String(t), x(t), H - m.b + 20); }
    ctx.textAlign = "right";
    for (const t of y.ticks(6)) { ctx.beginPath(); ctx.moveTo(m.l, y(t)); ctx.lineTo(W - m.r, y(t)); ctx.stroke(); ctx.fillText(String(t), m.l - 8, y(t) + 4); }
    for (const pass of ["ns", "sig"]) {
      for (const p of pts) {
        if ((pass === "ns") !== (p.d === "ns")) continue;
        const [cx, cy] = px(p);
        ctx.fillStyle = p.d === "up" ? UP : p.d === "down" ? DOWN : NS;
        ctx.globalAlpha = p.d === "ns" ? 0.45 : 0.8;
        ctx.beginPath(); ctx.arc(cx, cy, p.d === "ns" ? 1.6 : 2.4, 0, Math.PI * 2); ctx.fill();
      }
    }
    ctx.globalAlpha = 1;
    const mark = (i: number | null, color: string) => {
      if (i === null) return; const p = pts[i]; const [cx, cy] = px(p);
      ctx.strokeStyle = color; ctx.lineWidth = 2.5; ctx.beginPath(); ctx.arc(cx, cy, 7, 0, Math.PI * 2); ctx.stroke();
      ctx.fillStyle = "#1e2b27"; ctx.font = "bold 14px system-ui"; ctx.textAlign = cx > W / 2 ? "right" : "left";
      ctx.fillText(p.s, cx + (cx > W / 2 ? -12 : 12), cy - 10);
    };
    mark(selected ? pts.findIndex((p) => p.s === selected) : null, "#1e2b27");
    mark(hover, "#2f5f8a");
    ctx.fillStyle = "#46524e"; ctx.font = "14px system-ui"; ctx.textAlign = "center";
    ctx.fillText("log2 fold change, tumor vs healthy tissue", (m.l + W - m.r) / 2, H - 10);
    ctx.save(); ctx.translate(16, (m.t + H - m.b) / 2); ctx.rotate(-Math.PI / 2); ctx.fillText("−log10 adjusted p", 0, 0); ctx.restore();
  }, [pts, selected, hover]); // eslint-disable-line react-hooks/exhaustive-deps

  const nearest = (e: MouseEvent<HTMLCanvasElement>) => {
    const r = e.currentTarget.getBoundingClientRect();
    const mx = ((e.clientX - r.left) / r.width) * W, my = ((e.clientY - r.top) / r.height) * H;
    let best = -1, bd = 144;
    pts.forEach((p, i) => { if (p.d === "ns") return; const [cx, cy] = px(p); const d = (cx - mx) ** 2 + (cy - my) ** 2; if (d < bd) { bd = d; best = i; } });
    return best >= 0 ? best : null;
  };
  return (
    <div>
      <canvas ref={ref} style={{ width: "100%", height: "auto", aspectRatio: `${W} / ${H}`, cursor: "crosshair" }}
        onMouseMove={(e) => setHover(nearest(e))} onMouseLeave={() => setHover(null)}
        onClick={(e) => { const i = nearest(e); if (i !== null) onSelect(pts[i].s); }}
        role="img" aria-label="Volcano plot of all tested genes. Use the gene search to select a gene." />
      <div className="legend" style={{ marginTop: 8 }}>
        <span><i className="swatch" style={{ background: UP }} />Higher in tumor</span>
        <span><i className="swatch" style={{ background: DOWN }} />Lower in tumor</span>
        <span><i className="swatch" style={{ background: NS }} />Not significant or &lt; 2-fold</span>
        <span>Click a coloured point to open that gene.</span>
      </div>
    </div>
  );
}

/* ---------- Kaplan–Meier, median split ---------- */
export function KMPlot({ rows, score }: { rows: Record<string, number | string>[]; score: string }) {
  const W = 640, H = 340, m = { t: 16, r: 20, b: 48, l: 56 };
  const curves = useMemo(() => {
    const vals = rows.map((r) => Number(r[score])).filter((v) => !Number.isNaN(v)).sort((a, b) => a - b);
    const med = vals[Math.floor(vals.length / 2)];
    const km = (sub: Record<string, number | string>[]) => {
      const d = sub.map((r) => ({ t: Number(r.pfi_months), e: Number(r.pfi) })).sort((a, b) => a.t - b.t);
      let atRisk = d.length, s = 1; const pts: [number, number][] = [[0, 1]];
      for (let i = 0; i < d.length;) {
        const t = d[i].t; let ev = 0, n = 0;
        while (i < d.length && d[i].t === t) { ev += d[i].e; n++; i++; }
        if (ev > 0) { s *= 1 - ev / atRisk; pts.push([t, s]); }
        atRisk -= n;
      }
      pts.push([d[d.length - 1]?.t ?? 0, s]);
      return { pts, n: d.length, events: d.reduce((a, b) => a + b.e, 0) };
    };
    return {
      high: km(rows.filter((r) => Number(r[score]) > med)),
      low: km(rows.filter((r) => Number(r[score]) <= med)),
    };
  }, [rows, score]);
  const tMax = Math.max(...rows.map((r) => Number(r.pfi_months)));
  const x = scaleLinear().domain([0, Math.min(tMax, 120)]).range([m.l, W - m.r]);
  const y = scaleLinear().domain([0.4, 1]).range([H - m.b, m.t]);
  const path = line<[number, number]>().x((d) => x(Math.min(d[0], 120))).y((d) => y(Math.max(0.4, d[1]))).curve(curveStepAfter);
  return (
    <figure style={{ margin: 0 }}>
      <svg className="chart" viewBox={`0 0 ${W} ${H}`} width="100%" role="img" aria-label="Kaplan–Meier curves of progression-free interval">
        {y.ticks(6).map((t) => <g key={t}><line x1={m.l} x2={W - m.r} y1={y(t)} y2={y(t)} stroke="#eee9df" /><text x={m.l - 8} y={y(t) + 4} textAnchor="end">{t.toFixed(1)}</text></g>)}
        {x.ticks(6).map((t) => <text key={t} x={x(t)} y={H - m.b + 20} textAnchor="middle">{t}</text>)}
        <path d={path(curves.low.pts) ?? ""} fill="none" stroke={UP} strokeWidth={2.6} />
        <path d={path(curves.high.pts) ?? ""} fill="none" stroke={DOWN} strokeWidth={2.6} />
        <text x={(m.l + W - m.r) / 2} y={H - 8} textAnchor="middle">Months after surgery</text>
        <text transform={`translate(16 ${(H - m.b + m.t) / 2}) rotate(-90)`} textAnchor="middle">Progression-free</text>
      </svg>
      <div className="legend">
        <span><i className="swatch" style={{ background: DOWN }} />Higher score ({curves.high.events} events / {curves.high.n})</span>
        <span><i className="swatch" style={{ background: UP }} />Lower score ({curves.low.events} events / {curves.low.n})</span>
      </div>
    </figure>
  );
}

/* ---------- 100% stacked bars ---------- */
export function StackedBars({ groups, data, colors }: {
  groups: string[]; data: { group: string; key: string; n: number }[]; colors: Record<string, string>;
}) {
  const keys = Object.keys(colors);
  return (
    <div className="stack">
      {groups.map((g) => {
        const rows = data.filter((d) => d.group === g); const tot = rows.reduce((a, b) => a + b.n, 0) || 1;
        return (
          <div key={g}>
            <div style={{ display: "flex", justifyContent: "space-between", marginBottom: 6 }}><strong>{g}</strong><span className="muted">n = {tot}</span></div>
            <div style={{ display: "flex", height: 28, borderRadius: 10, overflow: "hidden" }} role="img"
              aria-label={`${g}: ` + keys.map((k) => `${k} ${Math.round((100 * (rows.find((r) => r.key === k)?.n ?? 0)) / tot)}%`).join(", ")}>
              {keys.map((k) => { const n = rows.find((r) => r.key === k)?.n ?? 0;
                return n ? <div key={k} title={`${k}: ${n}`} style={{ width: `${(100 * n) / tot}%`, background: colors[k] }} /> : null; })}
            </div>
          </div>
        );
      })}
      <div className="legend">{keys.map((k) => <span key={k}><i className="swatch" style={{ background: colors[k] }} />{k.replace("_", " ")}</span>)}</div>
    </div>
  );
}
