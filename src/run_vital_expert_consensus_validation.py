from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from run_vital_annotation_mvp import (
    infer_action,
    infer_context,
    infer_evaluability,
    infer_flag,
    infer_reason,
    infer_regime,
    infer_threshold,
)
from score_vital_from_variant_summary import score_preloaded_snapshot
from validate_vital_reclassification import ARRHYTHMIA_GENES, clinical_group


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
DEFAULT_VARIANT_SUMMARY = BASE_DIR / "data" / "raw" / "variant_summary_2026-04.txt.gz"
DEFAULT_PREFIX = "arrhythmia_expert_consensus_validation"


def collapse_first_nonempty(values: pd.Series) -> str:
    cleaned = [str(value) for value in values.dropna() if str(value) and str(value) != "-"]
    return cleaned[0] if cleaned else ""


def collapse_unique(values: pd.Series) -> str:
    cleaned = sorted({str(value) for value in values.dropna() if str(value) and str(value) != "-"})
    return "|".join(cleaned)


def make_vcv_accession(variation_id: object) -> str:
    try:
        return f"VCV{int(str(variation_id)):09d}"
    except ValueError:
        return f"VCV{variation_id}"


def review_bucket(series: pd.Series) -> pd.Series:
    text = series.fillna("").astype(str).str.lower()
    bucket = pd.Series("other", index=series.index, dtype="object")
    bucket.loc[text.str.contains("reviewed by expert panel", na=False)] = "expert_panel"
    bucket.loc[text.str.contains("practice guideline", na=False)] = "practice_guideline"
    return bucket


def load_expert_consensus_snapshot(variant_summary_path: Path) -> pd.DataFrame:
    usecols = [
        "Type",
        "Name",
        "GeneSymbol",
        "ClinicalSignificance",
        "ReviewStatus",
        "NumberSubmitters",
        "Assembly",
        "Chromosome",
        "PositionVCF",
        "ReferenceAlleleVCF",
        "AlternateAlleleVCF",
        "VariationID",
    ]
    chunks: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        variant_summary_path,
        sep="\t",
        compression="infer",
        usecols=lambda column: column in usecols,
        dtype=str,
        chunksize=250_000,
        low_memory=False,
    ):
        chunk = chunk[chunk["Assembly"].eq("GRCh38") & chunk["GeneSymbol"].isin(ARRHYTHMIA_GENES)].copy()
        if chunk.empty:
            continue
        chunk["clinical_group"] = chunk["ClinicalSignificance"].map(clinical_group)
        chunk["review_bucket"] = review_bucket(chunk["ReviewStatus"])
        chunk = chunk[
            chunk["review_bucket"].isin({"expert_panel", "practice_guideline"})
            & chunk["clinical_group"].isin({"P_LP", "B_LB"})
        ].copy()
        if chunk.empty:
            continue
        chunks.append(chunk)

    if not chunks:
        return pd.DataFrame()

    df = pd.concat(chunks, ignore_index=True)
    for column in ["Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF", "VariationID"]:
        df[column] = df[column].fillna("").astype(str)
    df = df[
        df["Chromosome"].ne("")
        & df["PositionVCF"].ne("")
        & df["ReferenceAlleleVCF"].ne("")
        & df["AlternateAlleleVCF"].ne("")
        & ~df["ReferenceAlleleVCF"].isin({"na", "-"})
        & ~df["AlternateAlleleVCF"].isin({"na", "-"})
    ].copy()
    df["variant_key"] = (
        df["Chromosome"]
        + ":"
        + df["PositionVCF"]
        + ":"
        + df["ReferenceAlleleVCF"]
        + ":"
        + df["AlternateAlleleVCF"]
    )
    grouped = (
        df.groupby(["VariationID", "variant_key"], dropna=False)
        .agg(
            gene=("GeneSymbol", collapse_first_nonempty),
            chrom=("Chromosome", collapse_first_nonempty),
            pos=("PositionVCF", collapse_first_nonempty),
            ref=("ReferenceAlleleVCF", collapse_first_nonempty),
            alt=("AlternateAlleleVCF", collapse_first_nonempty),
            clinsig=("ClinicalSignificance", collapse_unique),
            clinical_group=("clinical_group", collapse_unique),
            review_status=("ReviewStatus", collapse_unique),
            review_bucket=("review_bucket", collapse_unique),
            title=("Name", collapse_first_nonempty),
            variant_type=("Type", collapse_first_nonempty),
            submitter_count=("NumberSubmitters", "max"),
        )
        .reset_index()
        .rename(columns={"VariationID": "variation_id"})
    )
    grouped["clinvar_id"] = grouped["variation_id"].map(make_vcv_accession)
    grouped["gnomad_af"] = pd.NA
    grouped["gnomad_ac"] = pd.NA
    grouped["gnomad_an"] = pd.NA
    grouped["match_category"] = "exact_match"
    grouped["gnomad_status"] = ""
    grouped["gnomad_variants_at_pos_for_match"] = ""
    grouped["gnomad_region_query_status"] = ""
    grouped["freq_class"] = ""
    grouped["submitter_count"] = pd.to_numeric(grouped["submitter_count"], errors="coerce")
    return grouped.sort_values(["clinical_group", "gene", "variation_id"]).reset_index(drop=True)


def add_vital_regimes(scores: pd.DataFrame) -> pd.DataFrame:
    table = scores.copy()
    table["VITAL_evaluability"] = table.apply(infer_evaluability, axis=1)
    table["VITAL_context"] = table["gene"].map(infer_context)
    table["VITAL_regime"] = table.apply(infer_regime, axis=1)
    table["VITAL_flag"] = table.apply(infer_flag, axis=1)
    table["VITAL_action"] = table.apply(infer_action, axis=1)
    table["VITAL_reason"] = table.apply(infer_reason, axis=1)
    table["VITAL_threshold"] = table.apply(infer_threshold, axis=1)
    return table


def make_summary_table(scored: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    groups = [
        ("expert_curated_P_LP", scored[scored["clinical_group"].eq("P_LP")]),
        ("expert_curated_B_LB", scored[scored["clinical_group"].eq("B_LB")]),
        ("all_expert_curated_arrhythmia", scored),
    ]
    for label, sub in groups:
        tier1 = sub[sub["VITAL_evaluability"].eq("TIER_1")].copy()
        rows.append(
            {
                "group": label,
                "n_total": len(sub),
                "n_tier1": len(tier1),
                "n_tier2": int(sub["VITAL_evaluability"].eq("TIER_2").sum()),
                "n_unevaluable": int(sub["VITAL_evaluability"].eq("UNEVALUABLE").sum()),
                "ok_count": int(tier1["VITAL_flag"].eq("OK").sum()),
                "check_popmax_count": int(tier1["VITAL_flag"].eq("CHECK_POPMAX").sum()),
                "model_conflict_count": int(tier1["VITAL_flag"].eq("MODEL_CONFLICT").sum()),
                "ok_percent_of_tier1": pct(int(tier1["VITAL_flag"].eq("OK").sum()), len(tier1)),
                "check_popmax_percent_of_tier1": pct(int(tier1["VITAL_flag"].eq("CHECK_POPMAX").sum()), len(tier1)),
                "model_conflict_percent_of_tier1": pct(int(tier1["VITAL_flag"].eq("MODEL_CONFLICT").sum()), len(tier1)),
                "interpretation": interpretation_for_group(label),
            }
        )
    return pd.DataFrame(rows)


def interpretation_for_group(label: str) -> str:
    if label == "expert_curated_P_LP":
        return (
            "Positive-control comparator. MODEL_CONFLICT here means incompatibility with the tested "
            "unqualified dominant high-penetrance model, not automatic disagreement with expert curation."
        )
    if label == "expert_curated_B_LB":
        return (
            "Negative-control comparator. False positive concern is whether benign expert-curated records "
            "are spuriously routed into MODEL_CONFLICT in exact usable AF space."
        )
    return "Combined expert-curated arrhythmia subset."


def pct(numerator: int, denominator: int) -> float | None:
    if denominator == 0:
        return None
    return 100.0 * numerator / denominator


def make_case_table(scored: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "gene",
        "clinvar_id",
        "variation_id",
        "title",
        "clinical_group",
        "clinsig",
        "review_status",
        "review_strength",
        "frequency_evidence_status",
        "global_af",
        "global_ac",
        "popmax_af",
        "popmax_ac",
        "popmax_population",
        "standard_acmg_frequency_flag",
        "frequency_signal_ac_ge_20",
        "VITAL_evaluability",
        "VITAL_flag",
        "VITAL_regime",
        "VITAL_action",
        "VITAL_reason",
        "VITAL_threshold",
    ]
    present = [column for column in keep if column in scored.columns]
    cases = scored[
        scored["VITAL_evaluability"].eq("TIER_1")
        & scored["VITAL_flag"].isin(["CHECK_POPMAX", "MODEL_CONFLICT"])
    ].copy()
    return cases.loc[:, present].sort_values(
        ["clinical_group", "VITAL_flag", "gene", "variation_id"],
        ascending=[True, False, True, True],
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate VITAL regime calls against expert-curated arrhythmia ClinVar records."
    )
    parser.add_argument("--variant-summary", type=Path, default=DEFAULT_VARIANT_SUMMARY)
    parser.add_argument("--output-prefix", default=DEFAULT_PREFIX)
    parser.add_argument("--dataset", default="gnomad_r4")
    parser.add_argument("--gnomad-pause", type=float, default=0.1)
    parser.add_argument("--no-fetch-gnomad", action="store_true")
    parser.add_argument("--force-gnomad-fetch", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    snapshot = load_expert_consensus_snapshot(args.variant_summary)
    if snapshot.empty:
        raise SystemExit("No expert-curated GRCh38 arrhythmia variants found in the supplied ClinVar snapshot.")

    score_preloaded_snapshot(
        snapshot=snapshot,
        output_prefix=args.output_prefix,
        dataset=args.dataset,
        pause=args.gnomad_pause,
        fetch_gnomad=not args.no_fetch_gnomad,
        force_gnomad=args.force_gnomad_fetch,
    )

    scored_path = DATA_DIR / f"{args.output_prefix}_vital_scores.csv"
    scored = pd.read_csv(scored_path)
    scored["variation_id"] = scored["variation_id"].astype(str)
    scored["variant_key"] = scored["variant_key"].astype(str)
    snapshot["variation_id"] = snapshot["variation_id"].astype(str)
    snapshot["variant_key"] = snapshot["variant_key"].astype(str)
    scored = scored.merge(
        snapshot.loc[:, ["variation_id", "variant_key", "clinical_group", "review_bucket"]],
        on=["variation_id", "variant_key"],
        how="left",
    )
    scored = add_vital_regimes(scored)

    regime_path = DATA_DIR / f"{args.output_prefix}_regime_calls.csv"
    summary_path = DATA_DIR / f"{args.output_prefix}_summary.csv"
    cases_path = DATA_DIR / f"{args.output_prefix}_flagged_cases.csv"

    scored.to_csv(regime_path, index=False)
    make_summary_table(scored).to_csv(summary_path, index=False)
    make_case_table(scored).to_csv(cases_path, index=False)

    print(f"Saved {regime_path}")
    print(f"Saved {summary_path}")
    print(f"Saved {cases_path}")


if __name__ == "__main__":
    main()
