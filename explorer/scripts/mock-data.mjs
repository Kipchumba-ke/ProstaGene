// Generates MOCK explorer data with the same shape as workflow/scripts/12_export_explorer.R.
// Used for local development and CI when real pipeline results are not available.
// The app shows a "mock data" banner whenever summary.generated === "mock".
import { mkdirSync, writeFileSync, existsSync } from "node:fs";
import { join, dirname } from "node:path";
import { fileURLToPath } from "node:url";

const out = join(dirname(fileURLToPath(import.meta.url)), "..", "public", "data");
if (existsSync(join(out, "summary.json")) && !process.argv.includes("--force")) {
  console.log("public/data already has data; use --force to overwrite with mock data.");
  process.exit(0);
}
mkdirSync(out, { recursive: true });

let seed = 42;
const rnd = () => ((seed = (seed * 16807) % 2147483647) / 2147483647);
const norm = () => Math.sqrt(-2 * Math.log(rnd() + 1e-12)) * Math.cos(2 * Math.PI * rnd());
const r3 = (x) => Math.round(x * 1000) / 1000;
const write = (name, obj) => writeFileSync(join(out, name), JSON.stringify(obj));

const markers = { DLX1: 6.35, ZIC2: 5.73, PCSEAT: 4.92, AMACR: 3.16, DLX2: 4.99, PCA3: 5.23, SLC45A2: 4.17, HPN: 2.71,
  TDRD1: 4.43, ERG: 1.8, KRT5: -3.2, TP63: -2.9, KRT15: -3.4, MYH11: -2.6, ACTG2: -2.8, KLK3: 0.3, AR: 0.4, MYC: 0.9 };

// ---- genes ----
const N = 20000;
const g = { symbol: [], type: [], baseMean: [], lfc: [], padj: [], dir: [], geo_lfc: [], geo_padj: [],
            race_lfc: [], race_padj: [], afr_lfc: [], afr_padj: [], germline: ["AC009268.2", "FP565260.1", "F8A1", "GATD3B"] };
const names = [...Object.keys(markers), ...Array.from({ length: N - Object.keys(markers).length }, (_, i) => `GENE${i + 1}`)];
names.forEach((s) => {
  const lfc = markers[s] ?? norm() * 0.9;
  const se = 0.08 + rnd() * 0.25;
  const z = lfc / se;
  const p = Math.min(1, Math.exp(-0.717 * Math.abs(z) - 0.416 * z * z) * 2);
  const padj = Math.min(1, p * 3);
  g.symbol.push(s); g.type.push(rnd() < 0.8 ? "protein_coding" : "lncRNA");
  g.baseMean.push(r3(1 + rnd() * 3.5)); g.lfc.push(r3(lfc)); g.padj.push(Number(padj.toPrecision(3)));
  g.dir.push(padj < 0.05 && Math.abs(lfc) >= 1 ? (lfc > 0 ? "up" : "down") : "ns");
  g.geo_lfc.push(r3(lfc * 0.7 + norm() * 0.5)); g.geo_padj.push(Number(Math.min(1, padj * 40).toPrecision(3)));
  g.race_lfc.push(r3(norm() * 0.3)); g.race_padj.push(Number((0.05 + rnd()).toPrecision(3)));
  g.afr_lfc.push(r3(norm() * 0.4)); g.afr_padj.push(Number((0.05 + rnd()).toPrecision(3)));
});
write("genes.json", g);

// ---- hallmark ----
const hm = [
  ["MYC_TARGETS_V1", 2.73, 1.94, 1.9, -1.05, -1.39], ["MYC_TARGETS_V2", 2.67, 2.14, 2.1, -0.8, -0.9], ["G2M_CHECKPOINT", 2.6, 1.97, 2.2, -0.9, -1.0],
  ["E2F_TARGETS", 2.44, 1.55, 2.0, -0.7, -0.8], ["UNFOLDED_PROTEIN_RESPONSE", 2.08, 1.29, 1.5, -1.48, -1.7], ["ANDROGEN_RESPONSE", 1.94, 2.46, 1.7, -1.65, -1.69],
  ["MITOTIC_SPINDLE", 1.6, 1.4, 1.3, -1.61, -1.68], ["MTORC1_SIGNALING", 1.52, 1.24, 1.2, -1.66, -1.87], ["INTERFERON_GAMMA_RESPONSE", -0.4, 0.6, -0.2, 1.84, 2.01],
  ["INTERFERON_ALPHA_RESPONSE", -0.3, 0.5, -0.1, 1.96, 2.18], ["ALLOGRAFT_REJECTION", -0.88, 1.79, -0.5, 1.61, 1.86], ["TNFA_SIGNALING_VIA_NFKB", -1.77, -1.19, -1.5, 0.75, 0.9],
  ["HYPOXIA", -1.96, -1.6, -1.6, -1.3, -1.23], ["APICAL_JUNCTION", -1.97, -1.46, -1.5, -1.55, -1.35], ["EPITHELIAL_MESENCHYMAL_TRANSITION", -2.2, -1.47, -1.8, -1.23, -1.1],
  ["MYOGENESIS", -2.66, -2.09, -2.3, -2.44, -2.2], ["P53_PATHWAY", -1.61, -1.55, -1.3, -0.6, -0.5], ["ADIPOGENESIS", -1.35, -1.35, -1.0, -1.92, -2.05],
  ["PROTEIN_SECRETION", -0.6, -0.4, -0.5, -2.47, -2.79], ["OXIDATIVE_PHOSPHORYLATION", 0.62, -1.47, 0.4, -1.13, -1.42],
];
const pv = (n) => Number(Math.min(1, Math.exp(-2.2 * Math.abs(n) * Math.abs(n) + 3)).toPrecision(3));
write("hallmark.json", hm.map(([p, NES, adj, geo, race, afr]) => ({
  pathway: `HALLMARK_${p}`, label: p.replaceAll("_", " "), NES, padj: pv(NES), size: 100 + Math.round(rnd() * 100),
  NES_adj: adj, padj_adj: pv(adj), composition: Math.abs(NES) > 1.3 && Math.abs(adj) > 1.3 ? "robust to composition" : Math.abs(NES) > 1.3 ? "composition-sensitive" : "neither",
  NES_geo: geo, padj_geo: pv(geo), NES_race: race, padj_race: pv(race), NES_afr: afr, padj_afr: pv(afr),
  NES_ergpos: r3(NES + norm() * 0.3), NES_ergneg: r3(NES + norm() * 0.3), p_erg_diff: Number((0.05 + rnd() * 0.9).toPrecision(2)),
})));

// ---- expression (52 pairs) ----
const samples = [];
for (let i = 0; i < 52; i++) for (const c of ["normal", "tumor"]) samples.push({ patient: `TCGA-XX-${1000 + i}`, condition: c });
const expr = { samples, genes: {} };
for (const s of Object.keys(markers)) {
  const base = 6 + rnd() * 4;
  expr.genes[s] = samples.map((x, i) => Math.round((base + (x.condition === "tumor" ? markers[s] * 0.6 : 0) + norm() * 0.6 + (i % 7) * 0.02) * 100) / 100);
}
write("expression.json", expr);

// ---- survival ----
const scores = ["tcga_tumor_signature", "hm_E2F_TARGETS", "hm_INTERFERON_GAMMA_RESPONSE", "estimate_stromal"];
const rows = Array.from({ length: 497 }, () => {
  const e2f = norm();
  const gg = rnd() < 0.09 ? "GS6" : rnd() < 0.55 ? "GS7" : "GS8_10";
  const hazard = 0.006 * Math.exp(0.35 * e2f + (gg === "GS8_10" ? 1.2 : gg === "GS7" ? 0.3 : -1.5));
  const t = -Math.log(rnd()) / hazard, cens = 5 + rnd() * 70;
  const row = { pfi_months: r3(Math.min(t, cens)), pfi: t < cens ? 1 : 0, gleason_group: gg };
  row.tcga_tumor_signature = r3(norm()); row.hm_E2F_TARGETS = r3(e2f);
  row.hm_INTERFERON_GAMMA_RESPONSE = r3(norm()); row.estimate_stromal = r3(norm());
  return row;
});
const cox = scores.flatMap((s, i) => [
  { score: s, model: "univariable", HR: [1.15, 1.43, 1.07, 1.1][i], lower: [0.91, 1.13, 0.85, 0.87][i], upper: [1.46, 1.81, 1.33, 1.4][i], padj: [0.527, 0.0386, 0.658, 0.527][i], c_index: [0.538, 0.616, 0.529, 0.532][i] },
  { score: s, model: "adjusted", HR: [1.07, 1.04, 0.95, 0.96][i], lower: [0.85, 0.81, 0.76, 0.77][i], upper: [1.35, 1.32, 1.19, 1.21][i], padj: 0.982, c_index: [0.731, 0.731, 0.714, 0.723][i] },
]);
write("survival.json", { rows, scores, cox });

// ---- subtypes ----
write("subtypes.json", {
  counts: [{ subtype_group: "ERG", N: 152 }, { subtype_group: "Other", N: 85 }, { subtype_group: "ETS_other", N: 41 }, { subtype_group: "SPOP", N: 37 }, { subtype_group: "FOXA1", N: 9 }],
  by_ancestry: [
    { afr_group: "AFR < 0.5", subtype_group: "ERG", N: 140 }, { afr_group: "AFR < 0.5", subtype_group: "Other", N: 70 }, { afr_group: "AFR < 0.5", subtype_group: "ETS_other", N: 38 },
    { afr_group: "AFR < 0.5", subtype_group: "SPOP", N: 33 }, { afr_group: "AFR < 0.5", subtype_group: "FOXA1", N: 8 },
    { afr_group: "AFR >= 0.5", subtype_group: "ERG", N: 9 }, { afr_group: "AFR >= 0.5", subtype_group: "Other", N: 13 }, { afr_group: "AFR >= 0.5", subtype_group: "ETS_other", N: 3 },
    { afr_group: "AFR >= 0.5", subtype_group: "SPOP", N: 4 }, { afr_group: "AFR >= 0.5", subtype_group: "FOXA1", N: 1 }],
  erg_ancestry: [{ afr_group: "AFR < 0.5", n: 391, erg_pos: 182, pct_erg: 46.5 }, { afr_group: "AFR >= 0.5", n: 53, erg_pos: 15, pct_erg: 28.3 }],
  hallmark: [], cox: [{ model: "ERG (expression) + Gleason + T + age", term: "erg_expr_positiveTRUE", HR: 0.92, lower: 0.6, upper: 1.4, p: 0.7, n: 461, events: 76 }],
});

write("composition.json", [
  { score: "smooth_muscle", median_normal: 1.06, median_tumor: 0.288, tumor_higher_pairs: 5, pairs: 52, padj: 4.9e-7 },
  { score: "basal_epithelium", median_normal: 1.161, median_tumor: 0.31, tumor_higher_pairs: 7, pairs: 52, padj: 2.23e-5 },
  { score: "mcp_Fibroblasts", median_normal: 6.49, median_tumor: 6.11, tumor_higher_pairs: 15, pairs: 52, padj: 8.25e-3 },
  { score: "estimate_immune", median_normal: -208.8, median_tumor: -281.6, tumor_higher_pairs: 22, pairs: 52, padj: 0.0807 },
]);

// ---- summary (illustrative values; flagged as mock so the app shows a banner) ----
write("summary.json", {
  generated: "mock",
  design: { pairs: 52, samples: 104, genes_tested: 20303 },
  de: { degs: 2233, up: 763, down: 1470, any_fc: 12752, top: ["DLX1", "ZIC2", "PCSEAT", "AMACR", "DLX2", "PCA3", "SLC45A2", "NKX2-3", "HPN", "TDRD1"] },
  sensitivity: { spearman_log2FC_all: 0.993, jaccard: 0.871 },
  pathways: { hallmark_sig: 27, up: 9, down: 18 },
  validation: { complete_pairs: 13, sign_agreement_tcga_DEGs: 0.919, signature_auc: 0.982, hallmark_replicated: 20, hallmark_nes_pearson: 0.726, spearman_lfc_all: 0.684 },
  ancestry: { tcga_tumors_black: 57, tcga_tumors_white: 384, tcga_mt_pseudogene_hits: 15 },
  composition: { cor_tumor_smooth_muscle: -0.553, hallmark_robust_full: 17, race_immune_diff_sd: 0.237 },
  survival: { tumors_with_pfi: 497, pfi_events: 93, median_followup_months: 25.9, c_index_clinical: 0.719, adjusted_scores_fdr05: 0, univariable_scores_fdr05: 1, events_black: 9, n_black: 57 },
  genetic_ancestry: { tumors_with_ancestry: 471, black_median_afr: 0.805, nes_r_afr_vs_race: 0.873, likely_germline_genes: 10, top_genes_tested: 30, ifng_afr_beta: 0.256, ifng_afr_p: 0.176 },
  subtypes: { tumors_with_paper_subtype: 324, erg_pos_expression: 197, erg_pct_afr_low: 46.5, erg_pct_afr_high: 28.3, n_afr_high: 53, erg_call_agreement_pct: 95 },
});
console.log("Mock data written to", out);
