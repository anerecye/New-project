from __future__ import annotations

import argparse
import time
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
    import requests
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires pandas, numpy, and requests. Install them with: "
        "python -m pip install -r requirements.txt"
    ) from exc


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
SUPP_DIR = BASE_DIR / "supplementary_tables"

GENE_ORDER = [
    "KCNQ1",
    "KCNH2",
    "SCN5A",
    "KCNE1",
    "KCNE2",
    "RYR2",
    "CASQ2",
    "TRDN",
    "CALM1",
    "CALM2",
    "CALM3",
    "ANK2",
    "SCN4B",
    "KCNJ2",
    "HCN4",
    "CACNA1C",
    "CACNB2",
    "CACNA2D1",
    "AKAP9",
    "SNTA1",
]

PRIMARY_CHROMS = {str(value) for value in range(1, 23)} | {"X", "Y", "M", "MT"}
AF_ULTRA_RARE = 1e-5


def prefixed_name(prefix: str, name: str) -> str:
    return f"{prefix}_{name}" if prefix else name


def data_path(prefix: str, name: str) -> Path:
    return DATA_DIR / prefixed_name(prefix, name)


def supplementary_path(name: str) -> Path:
    return SUPP_DIR / name


def save_table(df: pd.DataFrame, output_path: Path, sep: str = ",") -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False, sep=sep)
    print(f"Saved {output_path} ({len(df)} rows)")


def normalize_chrom(value: object) -> str:
    chrom = str(value).strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return chrom


def select_primary_genomic_position(genomic_pos: object) -> dict[str, object]:
    positions = genomic_pos if isinstance(genomic_pos, list) else [genomic_pos]
    positions = [position for position in positions if isinstance(position, dict)]
    primary = [
        position
        for position in positions
        if normalize_chrom(position.get("chr", "")).upper() in PRIMARY_CHROMS
    ]
    candidates = primary or positions
    if not candidates:
        raise ValueError("No genomic_pos entries were returned.")
    return max(
        candidates,
        key=lambda position: int(position.get("end", 0)) - int(position.get("start", 0)),
    )


def fetch_gene_locus(gene: str, session: requests.Session) -> dict[str, object]:
    response = session.get(
        "https://mygene.info/v3/query",
        params={
            "q": f"symbol:{gene}",
            "fields": "symbol,entrezgene,genomic_pos",
            "species": "human",
            "size": 5,
        },
        timeout=30,
    )
    response.raise_for_status()
    hits = response.json().get("hits", [])
    hit = next((item for item in hits if str(item.get("symbol", "")).upper() == gene), None)
    if hit is None and hits:
        hit = hits[0]
    if not hit:
        raise ValueError(f"No MyGene.info hit for {gene}.")

    pos = select_primary_genomic_position(hit.get("genomic_pos"))
    chrom = normalize_chrom(pos["chr"])
    start = int(pos["start"])
    end = int(pos["end"])
    ucsc_chrom = chrom if chrom.lower().startswith("chr") else f"chr{chrom}"
    sequence_response = session.get(
        "https://api.genome.ucsc.edu/getData/sequence",
        params={"genome": "hg38", "chrom": ucsc_chrom, "start": start, "end": end},
        timeout=90,
    )
    sequence_response.raise_for_status()
    sequence = str(sequence_response.json().get("dna", "")).upper()
    gc_fraction = (
        (sequence.count("G") + sequence.count("C")) / len(sequence)
        if sequence
        else np.nan
    )
    return {
        "gene": gene,
        "entrezgene": hit.get("entrezgene", ""),
        "chrom": chrom,
        "grch38_start_0based": start,
        "grch38_end_0based_exclusive": end,
        "gene_length_bp": end - start,
        "gc_fraction": gc_fraction,
        "gc_percent": gc_fraction * 100 if pd.notna(gc_fraction) else np.nan,
        "locus_source": "MyGene.info genomic_pos GRCh38",
        "gc_source": "UCSC hg38 getData/sequence",
    }


def load_or_fetch_gene_loci(cache_path: Path, refresh: bool = False) -> pd.DataFrame:
    if cache_path.exists() and not refresh:
        return pd.read_csv(cache_path)

    rows: list[dict[str, object]] = []
    with requests.Session() as session:
        session.headers.update({"User-Agent": "arrhythmia-supplementary-tables/1.0"})
        for gene in GENE_ORDER:
            rows.append(fetch_gene_locus(gene, session))
            time.sleep(0.15)
    loci = pd.DataFrame(rows)
    save_table(loci, cache_path)
    return loci


def make_s2(prefix: str) -> pd.DataFrame:
    risk_path = data_path(prefix, "reclassification_risk.csv")
    if not risk_path.exists():
        raise FileNotFoundError(
            f"Missing {risk_path}. Run advanced_variant_analyses.py first."
        )
    risk = pd.read_csv(risk_path)
    preferred_columns = [
        "variant_key",
        "gene",
        "clinvar_id",
        "variation_id",
        "title",
        "clinsig",
        "review_status",
        "review_strength",
        "review_score",
        "submitter_count",
        "submitter_count_source",
        "variant_type",
        "functional_class",
        "global_af",
        "global_ac",
        "global_an",
        "popmax_af",
        "popmax_ac",
        "popmax_population",
        "max_frequency_signal",
        "max_frequency_source",
        "max_frequency_ac",
        "frequency_threshold",
        "frequency_signal_ac_filter",
        "frequency_signal_ac_ge_20",
        "weak_review_signal",
        "potential_misclassification_signal",
        "risk_level",
        "signal_reason",
    ]
    columns = [column for column in preferred_columns if column in risk.columns]
    columns.extend(column for column in risk.columns if column not in columns)
    return risk.loc[:, columns]


def make_s5(prefix: str, gene_loci: pd.DataFrame) -> pd.DataFrame:
    matched_path = data_path(prefix, "gnomad_matched.csv")
    fisher_path = data_path(prefix, "per_gene_fisher.csv")
    if not matched_path.exists():
        raise FileNotFoundError(f"Missing {matched_path}.")

    matched = pd.read_csv(matched_path)
    for column in ["gnomad_af", "gnomad_ac", "gnomad_an"]:
        if column not in matched.columns:
            matched[column] = np.nan
        matched[column] = pd.to_numeric(matched[column], errors="coerce")
    matched["match_category"] = matched["match_category"].fillna("").astype(str)
    matched["gene"] = matched["gene"].fillna("").astype(str)

    rows: list[dict[str, object]] = []
    for gene in GENE_ORDER:
        group = matched[matched["gene"].eq(gene)]
        total = len(group)
        exact = int(group["match_category"].eq("exact_match").sum())
        allele_discordance = int(group["match_category"].eq("allele_discordance").sum())
        no_record = int(group["match_category"].eq("no_gnomad_record").sum())
        query_error = int(group["match_category"].eq("query_error").sum())
        non_overlap = allele_discordance + no_record
        af_covered = group[group["match_category"].eq("exact_match") & group["gnomad_af"].notna()]
        outlier = int((af_covered["gnomad_af"] > AF_ULTRA_RARE).sum())
        rows.append(
            {
                "gene": gene,
                "clinvar_plp_variant_count": total,
                "exact_match_count": exact,
                "allele_discordance_count": allele_discordance,
                "no_gnomad_record_count": no_record,
                "query_error_count": query_error,
                "non_overlap_count": non_overlap,
                "non_overlap_percent": 100 * non_overlap / total if total else np.nan,
                "af_covered_exact_match_count": len(af_covered),
                "af_gt_1e_5_count": outlier,
                "af_gt_1e_5_percent_of_af_covered": (
                    100 * outlier / len(af_covered) if len(af_covered) else np.nan
                ),
            }
        )

    summary = pd.DataFrame(rows)
    gene_loci = gene_loci.copy()
    summary = summary.merge(gene_loci, on="gene", how="left")

    if fisher_path.exists():
        fisher = pd.read_csv(fisher_path)
        fisher = fisher.rename(
            columns={
                "n_af_covered": "fisher_n_af_covered",
                "n_outlier": "fisher_n_af_gt_1e_5",
                "pct_outlier": "fisher_pct_af_gt_1e_5",
                "fisher_p": "fisher_p_greater",
                "bh_q": "fisher_bh_q",
            }
        )
        keep = [
            "gene",
            "odds_ratio",
            "fisher_p_greater",
            "fisher_bh_q",
            "significant",
        ]
        summary = summary.merge(fisher.loc[:, [column for column in keep if column in fisher.columns]], on="gene", how="left")

    ordered_columns = [
        "gene",
        "chrom",
        "grch38_start_0based",
        "grch38_end_0based_exclusive",
        "gene_length_bp",
        "gc_percent",
        "clinvar_plp_variant_count",
        "exact_match_count",
        "allele_discordance_count",
        "no_gnomad_record_count",
        "query_error_count",
        "non_overlap_count",
        "non_overlap_percent",
        "af_covered_exact_match_count",
        "af_gt_1e_5_count",
        "af_gt_1e_5_percent_of_af_covered",
        "odds_ratio",
        "fisher_p_greater",
        "fisher_bh_q",
        "significant",
        "entrezgene",
        "locus_source",
        "gc_source",
    ]
    gene_order = {gene: index for index, gene in enumerate(GENE_ORDER)}
    summary["_gene_order"] = summary["gene"].map(gene_order).fillna(len(gene_order))
    summary = summary.sort_values(["_gene_order", "gene"]).drop(columns="_gene_order").reset_index(drop=True)
    return summary.loc[:, [column for column in ordered_columns if column in summary.columns]]


def make_s6(prefix: str) -> pd.DataFrame:
    vital_path = data_path(prefix, "vital_scores.csv")
    if not vital_path.exists():
        raise FileNotFoundError(
            f"Missing {vital_path}. Run advanced_variant_analyses.py first."
        )
    vital = pd.read_csv(vital_path)
    preferred_columns = [
        "variant_key",
        "gene",
        "clinvar_id",
        "variation_id",
        "title",
        "clinsig",
        "review_status",
        "review_strength",
        "submitter_count",
        "match_category",
        "frequency_evidence_status",
        "variant_type",
        "functional_class",
        "global_af",
        "global_ac",
        "popmax_af",
        "popmax_ac",
        "popmax_population",
        "max_frequency_signal",
        "max_frequency_source",
        "qualifying_frequency_ac",
        "standard_acmg_frequency_flag",
        "standard_acmg_high_frequency_flag",
        "frequency_signal_ac_ge_20",
        "weak_review_signal",
        "gene_frequency_constraint_proxy",
        "technical_detectability_index",
        "frequency_pressure_score",
        "ac_reliability_score",
        "popmax_enrichment_score",
        "variant_type_tension_score",
        "technical_detectability_score",
        "gene_constraint_score",
        "review_fragility_score",
        "vital_score",
        "vital_band",
        "vital_red_flag",
        "predicted_revision_direction",
        "vital_signal_reason",
    ]
    columns = [column for column in preferred_columns if column in vital.columns]
    columns.extend(column for column in vital.columns if column not in columns)
    return vital.loc[:, columns]


def load_prefixed_table(prefix: str, name: str) -> pd.DataFrame:
    path = data_path(prefix, name)
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}. Run advanced_variant_analyses.py first.")
    return pd.read_csv(path)


def load_optional_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def make_s14_cross_disease_summary() -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for label, path in [
        (
            "current_2026_independent_cross_disease",
            DATA_DIR / "vital_cross_disease_3000_cross_disease_validation_summary.csv",
        ),
        (
            "baseline_2023_independent_cross_disease",
            DATA_DIR / "vital_cross_disease_3000_2023_01_cross_disease_validation_summary.csv",
        ),
        (
            "historical_2023_to_2026_validation_metrics",
            DATA_DIR / "vital_cross_disease_3000_2023_01_to_current_vital_historical_validation.csv",
        ),
        (
            "historical_2023_to_2026_enrichment",
            DATA_DIR / "vital_cross_disease_3000_2023_01_to_current_vital_historical_enrichment.csv",
        ),
    ]:
        table = load_optional_table(path)
        if table.empty:
            continue
        table.insert(0, "analysis", label)
        table.insert(1, "source_file", path.name)
        frames.append(table)
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def make_s15_cross_disease_outliers() -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for label, path in [
        (
            "current_2026_high_score_or_red",
            DATA_DIR / "vital_cross_disease_3000_cross_disease_outliers.csv",
        ),
        (
            "baseline_2023_high_score_or_red",
            DATA_DIR / "vital_cross_disease_3000_2023_01_cross_disease_outliers.csv",
        ),
    ]:
        table = load_optional_table(path)
        if table.empty:
            continue
        table.insert(0, "analysis", label)
        table.insert(1, "source_file", path.name)
        frames.append(table)

    historical = load_optional_table(
        DATA_DIR / "vital_cross_disease_3000_2023_01_to_current_vital_historical_predictions.csv"
    )
    if not historical.empty and "vital_red_flag" in historical.columns:
        red = historical[
            historical["vital_red_flag"].fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})
        ].copy()
        red.insert(0, "analysis", "historical_2023_red_followup")
        red.insert(1, "source_file", "vital_cross_disease_3000_2023_01_to_current_vital_historical_predictions.csv")
        frames.append(red)

    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate machine-readable supplementary tables S2, S5, and VITAL supplements."
    )
    parser.add_argument("--output-prefix", default="arrhythmia")
    parser.add_argument("--refresh-gene-loci", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    prefix = args.output_prefix.strip()
    locus_cache = data_path(prefix, "gene_locus_gc_summary.csv")
    gene_loci = load_or_fetch_gene_loci(locus_cache, refresh=args.refresh_gene_loci)

    s2 = make_s2(prefix)
    s5 = make_s5(prefix, gene_loci)
    s6 = make_s6(prefix)
    s7 = load_prefixed_table(prefix, "vital_method_comparison.csv")
    s8 = load_prefixed_table(prefix, "review_fragility_summary.csv")
    s9 = load_prefixed_table(prefix, "vital_threshold_sweep.csv")
    s10 = load_prefixed_table(prefix, "vital_acmg_disagreement.csv")
    s11 = load_prefixed_table(prefix, "vital_top_suspicious.csv")
    s12 = load_prefixed_table(prefix, "vital_absence_detectability_bias.csv")
    s14 = make_s14_cross_disease_summary()
    s15 = make_s15_cross_disease_outliers()

    save_table(s2, supplementary_path("Supplementary_Table_S2_reclassification_candidates.tsv"), sep="\t")
    save_table(s5, supplementary_path("Supplementary_Table_S5_gene_non_overlap_summary.tsv"), sep="\t")
    save_table(s6, supplementary_path("Supplementary_Table_S6_VITAL_scores.tsv"), sep="\t")
    save_table(s7, supplementary_path("Supplementary_Table_S7_VITAL_method_comparison.tsv"), sep="\t")
    save_table(s8, supplementary_path("Supplementary_Table_S8_review_fragility_summary.tsv"), sep="\t")
    save_table(s9, supplementary_path("Supplementary_Table_S9_VITAL_threshold_sweep.tsv"), sep="\t")
    save_table(s10, supplementary_path("Supplementary_Table_S10_ACMG_VITAL_disagreement.tsv"), sep="\t")
    save_table(s11, supplementary_path("Supplementary_Table_S11_top_suspicious_variants.tsv"), sep="\t")
    save_table(s12, supplementary_path("Supplementary_Table_S12_absence_detectability_bias.tsv"), sep="\t")
    if not s14.empty:
        save_table(s14, supplementary_path("Supplementary_Table_S14_cross_disease_3000_validation_summary.tsv"), sep="\t")
    if not s15.empty:
        save_table(s15, supplementary_path("Supplementary_Table_S15_cross_disease_3000_outliers.tsv"), sep="\t")
    save_table(s5, data_path(prefix, "gene_non_overlap_summary.csv"))


if __name__ == "__main__":
    main()
