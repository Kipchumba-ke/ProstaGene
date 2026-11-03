// Types and loader for the JSON written by workflow/scripts/12_export_explorer.R
import { useEffect, useState } from "react";

type Num = number | null;
type Metrics = Record<string, number | string | null>;

export interface Summary {
  generated: string;
  design: { pairs: number; samples: number; genes_tested: number };
  de: { degs: number; up: number; down: number; any_fc: number; top: string[] };
  sensitivity: Metrics;
  pathways: { hallmark_sig: number; up: number; down: number };
  validation: Metrics;
  ancestry: Metrics;
  composition: Metrics;
  survival: Metrics;
  genetic_ancestry: Metrics;
  subtypes: Metrics;
}

export interface Genes {
  symbol: string[]; type: string[]; baseMean: number[]; lfc: Num[]; padj: Num[]; dir: string[];
  geo_lfc: Num[]; geo_padj: Num[]; race_lfc: Num[]; race_padj: Num[]; afr_lfc: Num[]; afr_padj: Num[];
  germline: string[];
}

export interface Pathway {
  pathway: string; label: string; NES: number; padj: number; size: number;
  NES_adj: Num; padj_adj: Num; composition: string | null;
  NES_geo: Num; padj_geo: Num; NES_race: Num; padj_race: Num; NES_afr: Num; padj_afr: Num;
  NES_ergpos: Num; NES_ergneg: Num; p_erg_diff: Num;
}

export interface Expression {
  samples: { patient: string; condition: "normal" | "tumor" }[];
  genes: Record<string, number[]>;
}

export interface Survival {
  rows: Record<string, number | string>[];
  scores: string[];
  cox: { score: string; model: string; HR: number; lower: number; upper: number; padj: number; c_index: number }[];
}

export interface Subtypes {
  counts: { subtype_group: string; N: number }[];
  by_ancestry: { afr_group: string; subtype_group: string; N: number }[];
  erg_ancestry: { afr_group: string; n: number; erg_pos: number; pct_erg: number }[];
  hallmark: { contrast: string; pathway: string; NES: number; padj: number }[];
  cox: { model: string; term: string; HR: number; lower: number; upper: number; p: number; n: number; events: number }[];
}

export interface Composition {
  score: string; median_normal: number; median_tumor: number; tumor_higher_pairs: number; pairs: number; padj: number;
}

export interface AppData {
  summary: Summary; genes: Genes; hallmark: Pathway[]; expression: Expression;
  survival: Survival; subtypes: Subtypes; composition: Composition[];
}

const BASE: string = import.meta.env?.BASE_URL ?? "./";

async function getJson<T>(name: string): Promise<T> {
  const res = await fetch(`${BASE}data/${name}`);
  if (!res.ok) throw new Error(`Could not load ${name} (${res.status})`);
  return res.json() as Promise<T>;
}

export function useAppData() {
  const [data, setData] = useState<AppData | null>(null);
  const [error, setError] = useState<string | null>(null);
  useEffect(() => {
    Promise.all([
      getJson<Summary>("summary.json"), getJson<Genes>("genes.json"), getJson<Pathway[]>("hallmark.json"),
      getJson<Expression>("expression.json"), getJson<Survival>("survival.json"),
      getJson<Subtypes>("subtypes.json"), getJson<Composition[]>("composition.json"),
    ])
      .then(([summary, genes, hallmark, expression, survival, subtypes, composition]) =>
        setData({ summary, genes, hallmark, expression, survival, subtypes, composition }))
      .catch((e: Error) => setError(e.message));
  }, []);
  return { data, error };
}

// ---- formatting helpers ----
export const fmt = (x: Num | undefined, d = 2) => (x === null || x === undefined || Number.isNaN(x) ? "–" : x.toFixed(d));
export const fmtP = (p: Num | undefined) =>
  p === null || p === undefined ? "–" : p === 0 ? "< 1e-300" : p < 0.001 ? p.toExponential(1) : p.toPrecision(2);
export const fold = (lfc: Num) => (lfc === null ? "–" : `${(2 ** Math.abs(lfc)).toFixed(1)}×`);
export const num = (m: Metrics, k: string) => {
  const v = m?.[k];
  return typeof v === "number" ? v : v === null || v === undefined ? NaN : Number(v);
};
const SCORE_NAMES: Record<string, string> = {
  tcga_tumor_signature: "TCGA tumor signature", hm_E2F_TARGETS: "E2F targets", hm_G2M_CHECKPOINT: "G2M checkpoint",
  hm_ANDROGEN_RESPONSE: "Androgen response", hm_INTERFERON_GAMMA_RESPONSE: "Interferon-γ response",
  hm_INTERFERON_ALPHA_RESPONSE: "Interferon-α response", hm_MYOGENESIS: "Myogenesis", hm_MYC_TARGETS_V1: "MYC targets",
  hm_EPITHELIAL_MESENCHYMAL_TRANSITION: "EMT", estimate_stromal: "ESTIMATE stromal", estimate_immune: "ESTIMATE immune",
  smooth_muscle: "Smooth muscle", mcp_Fibroblasts: "Fibroblasts (MCP-counter)",
};
export const prettyScore = (s: string) => SCORE_NAMES[s] ?? s.replace(/^(hm|mcp|estimate)_/, "").replaceAll("_", " ");
