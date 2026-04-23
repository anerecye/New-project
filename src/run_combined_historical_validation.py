from __future__ import annotations

import argparse
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
except ModuleNotFoundError as exc:
    raise SystemExit("This script requires pandas and numpy.") from exc

from score_vital_from_variant_summary import score_snapshot
from validate_vital_reclassification import (
    ARRHYTHMIA_GENES,
    clinical_group,
    evaluate_reclassification,
    wilson_ci,
)


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"


PANEL_DEFINITIONS = {
    "cardiomyopathy": [
        "MYH7",
        "MYBPC3",
        "TNNT2",
        "TNNI3",
        "TPM1",
        "ACTC1",
        "MYL2",
        "MYL3",
        "LMNA",
        "FLNC",
        "PLN",
        "DSP",
        "PKP2",
        "DSG2",
        "DSC2",
    ],
    "epilepsy": [
        "SCN1A",
        "SCN2A",
        "SCN8A",
        "KCNQ2",
        "KCNQ3",
        "STXBP1",
        "CDKL5",
        "PCDH19",
        "SLC2A1",
        "GABRA1",
        "GABRG2",
        "DEPDC5",
    ],
    "hearing_loss": [
        "GJB2",
        "GJB6",
        "SLC26A4",
        "MYO7A",
        "OTOF",
        "TECTA",
        "TMC1",
        "CDH23",
        "USH2A",
        "MYO15A",
    ],
}


DOMAIN_CONFIGS = [
    {
        "domain": "arrhythmia",
        "baseline_prefix": "arrhythmia_2023_01",
        "historical_prefix": "arrhythmia_2023_01_to_current",
        "genes": ARRHYTHMIA_GENES,
        "sample_variants": None,
        "random_seed": 13,
    },
    {
        "domain": "cardiomyopathy",
        "baseline_prefix": "vital_domain_cardiomyopathy_2023_01",
        "historical_prefix": "vital_domain_cardiomyopathy_2023_01_to_current",
        "genes": PANEL_DEFINITIONS["cardiomyopathy"],
        "sample_variants": 300,
        "random_seed": 37,
    },
    {
        "domain": "epilepsy",
        "baseline_prefix": "vital_domain_epilepsy_2023_01",
        "historical_prefix": "vital_domain_epilepsy_2023_01_to_current",
        "genes": PANEL_DEFINITIONS["epilepsy"],
        "sample_variants": 300,
        "random_seed": 37,
    },
    {
        "domain": "hearing_loss",
        "baseline_prefix": "vital_domain_hearing_loss_2023_01",
        "historical_prefix": "vital_domain_hearing_loss_2023_01_to_current",
        "genes": PANEL_DEFINITIONS["hearing_loss"],
        "sample_variants": 300,
        "random_seed": 37,
    },
    {
        "domain": "random_clinvar_plp",
        "baseline_prefix": "vital_random_clinvar_plp_2023_01",
        "historical_prefix": "vital_random_clinvar_plp_2023_01_to_current",
        "genes": None,
        "sample_variants": 500,
        "random_seed": 37,
    },
]


def data_path(prefix: str, name: str) -> Path:
    return DATA_DIR / f"{prefix}_{name}"


def save_table(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"Saved {output_path} ({len(df)} rows)")


def bool_series(series: pd.Series) -> pd.Series:
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def safe_rate(numerator: int | float, denominator: int | float) -> float:
    return float(numerator / denominator) if denominator else np.nan


def ensure_score_table(
    config: dict[str, object],
    baseline_variant_summary: Path,
    dataset: str,
    gnomad_pause: float,
    fetch_gnomad: bool,
    force_score: bool,
    force_gnomad_fetch: bool,
) -> Path:
    prefix = str(config["baseline_prefix"])
    score_table = data_path(prefix, "vital_scores.csv")
    if score_table.exists() and not force_score:
        print(f"Reusing {score_table}")
        return score_table
    print(f"Scoring baseline snapshot for {config['domain']} -> {prefix}")
    score_snapshot(
        variant_summary_path=baseline_variant_summary,
        output_prefix=prefix,
        genes=config["genes"],  # type: ignore[arg-type]
        dataset=dataset,
        pause=gnomad_pause,
        fetch_gnomad=fetch_gnomad,
        force_gnomad=force_gnomad_fetch,
        max_variants=None,
        sample_variants=config["sample_variants"],  # type: ignore[arg-type]
        random_seed=int(config["random_seed"]),
    )
    return score_table


def run_domain_validation(
    config: dict[str, object],
    score_table: Path,
    baseline_variant_summary: Path,
    followup_variant_summary: Path,
) -> pd.DataFrame:
    historical_prefix = str(config["historical_prefix"])
    print(f"Validating {config['domain']} -> {historical_prefix}")
    _, predictions = evaluate_reclassification(
        score_table=score_table,
        baseline_summary=baseline_variant_summary,
        followup_summary=followup_variant_summary,
        output_prefix=historical_prefix,
        genes=config["genes"],  # type: ignore[arg-type]
    )
    predictions = predictions.copy()
    predictions.insert(0, "domain", config["domain"])
    predictions.insert(1, "historical_prefix", historical_prefix)
    return predictions


def metric_row(predictions: pd.DataFrame, domain: str, target_name: str, target_column: str) -> dict[str, object]:
    labels = bool_series(predictions[target_column])
    red = bool_series(predictions["vital_red_flag"])
    total_n = len(predictions)
    event_count = int(labels.sum())
    red_count = int(red.sum())
    red_events = int((red & labels).sum())
    nonred_count = total_n - red_count
    nonred_events = event_count - red_events
    event_rate = safe_rate(event_count, total_n)
    red_event_rate = safe_rate(red_events, red_count)
    nonred_event_rate = safe_rate(nonred_events, nonred_count)
    red_ci_low, red_ci_high = wilson_ci(red_events, red_count)
    overall_ci_low, overall_ci_high = wilson_ci(event_count, total_n)
    return {
        "domain": domain,
        "target": target_name,
        "baseline_variant_count": total_n,
        "event_count": event_count,
        "event_rate": event_rate,
        "event_rate_ci_low": overall_ci_low,
        "event_rate_ci_high": overall_ci_high,
        "vital_red_count": red_count,
        "vital_red_event_count": red_events,
        "vital_red_event_rate": red_event_rate,
        "vital_red_event_rate_ci_low": red_ci_low,
        "vital_red_event_rate_ci_high": red_ci_high,
        "vital_red_recall": safe_rate(red_events, event_count),
        "nonred_event_count": nonred_events,
        "nonred_event_rate": nonred_event_rate,
        "red_enrichment_vs_overall": safe_rate(red_event_rate, event_rate),
        "red_enrichment_vs_nonred": safe_rate(red_event_rate, nonred_event_rate),
        "interpretation": (
            "pooled_across_heterogeneous_domains"
            if domain == "combined_pooled"
            else "within_domain_descriptive"
        ),
    }


def stratified_bootstrap(
    predictions: pd.DataFrame,
    target_column: str,
    n_iter: int,
    random_seed: int,
) -> tuple[pd.DataFrame, dict[str, object]]:
    rng = np.random.default_rng(random_seed)
    domains = sorted(predictions["domain"].dropna().unique())
    rows: list[dict[str, object]] = []
    for iteration in range(n_iter):
        samples = []
        for domain in domains:
            sub = predictions[predictions["domain"].eq(domain)]
            if sub.empty:
                continue
            take = rng.integers(0, len(sub), size=len(sub))
            samples.append(sub.iloc[take])
        if not samples:
            continue
        boot = pd.concat(samples, ignore_index=True)
        red = bool_series(boot["vital_red_flag"])
        labels = bool_series(boot[target_column])
        red_count = int(red.sum())
        if red_count == 0:
            rows.append(
                {
                    "iteration": iteration,
                    "red_count": red_count,
                    "red_event_count": 0,
                    "red_event_rate": np.nan,
                }
            )
            continue
        red_events = int((red & labels).sum())
        rows.append(
            {
                "iteration": iteration,
                "red_count": red_count,
                "red_event_count": red_events,
                "red_event_rate": red_events / red_count,
            }
        )
    draws = pd.DataFrame(rows)
    usable = draws["red_event_rate"].dropna()
    usable_event_counts = draws.loc[draws["red_event_rate"].notna(), "red_event_count"]
    summary = {
        "target_column": target_column,
        "bootstrap_type": "stratified_by_domain_rows_resampled_with_replacement",
        "bootstrap_interpretation": "binary-noise sensitivity for sparse red sets, not a smooth stable effect-size estimator",
        "iterations_requested": n_iter,
        "iterations_with_red_variants": int(usable.count()),
        "iterations_without_red_variants": int(draws["red_event_rate"].isna().sum()) if not draws.empty else 0,
        "bootstrap_samples_with_zero_red_events": int((usable_event_counts == 0).sum()) if not usable_event_counts.empty else 0,
        "bootstrap_samples_with_one_red_event": int((usable_event_counts == 1).sum()) if not usable_event_counts.empty else 0,
        "bootstrap_fraction_with_zero_red_events": float((usable_event_counts == 0).mean()) if not usable_event_counts.empty else np.nan,
        "bootstrap_fraction_with_one_red_event": float((usable_event_counts == 1).mean()) if not usable_event_counts.empty else np.nan,
        "red_event_rate_bootstrap_mean": float(usable.mean()) if not usable.empty else np.nan,
        "red_event_rate_bootstrap_median": float(usable.median()) if not usable.empty else np.nan,
        "red_event_rate_bootstrap_ci_low": float(usable.quantile(0.025)) if not usable.empty else np.nan,
        "red_event_rate_bootstrap_ci_high": float(usable.quantile(0.975)) if not usable.empty else np.nan,
    }
    return draws, summary


def first_present_column(df: pd.DataFrame, candidates: list[str]) -> pd.Series:
    result = pd.Series("", index=df.index, dtype="object")
    for column in candidates:
        if column in df.columns:
            values = df[column].fillna("").astype(str)
            result = result.where(result.astype(str).ne(""), values)
    return result


def add_stable_identity_columns(df: pd.DataFrame) -> pd.DataFrame:
    result = df.copy()
    result["stable_variant_key"] = first_present_column(
        result, ["variant_key_baseline", "variant_key_score", "variant_key"]
    )
    result["stable_gene"] = first_present_column(result, ["gene_baseline", "gene_score", "gene"])
    return result


def deduplicate_for_analysis(predictions: pd.DataFrame, output_prefix: str) -> pd.DataFrame:
    table = add_stable_identity_columns(predictions)
    table["vital_score_numeric"] = pd.to_numeric(table.get("vital_score"), errors="coerce").fillna(-1)
    raw_rows = len(table)
    raw_red_count = int(bool_series(table["vital_red_flag"]).sum())
    table = table.sort_values(
        ["domain", "variation_id", "vital_score_numeric"],
        ascending=[True, True, False],
    )
    before_id = len(table)
    table = table.drop_duplicates(["domain", "variation_id"], keep="first").copy()
    dropped_by_id = before_id - len(table)

    keyed = table["stable_variant_key"].astype(str).ne("")
    table_keyed = table.loc[keyed].sort_values(
        ["domain", "stable_variant_key", "vital_score_numeric"],
        ascending=[True, True, False],
    )
    duplicate_key_drop_index = table_keyed[table_keyed.duplicated(["domain", "stable_variant_key"], keep="first")].index
    dropped_by_key = len(duplicate_key_drop_index)
    dropped = table.loc[duplicate_key_drop_index].copy()
    if not dropped.empty:
        save_table(dropped, data_path(output_prefix, "historical_dropped_within_domain_duplicate_rows.csv"))
    table = table.drop(index=duplicate_key_drop_index).copy()
    save_table(
        pd.DataFrame(
            [
                {
                    "dedupe_stage": "within_domain",
                    "raw_rows_before_within_domain_dedupe": raw_rows,
                    "rows_after_within_domain_dedupe": len(table),
                    "removed_duplicate_variation_id_rows": dropped_by_id,
                    "removed_duplicate_variant_key_rows": dropped_by_key,
                    "vital_red_before_dedupe": raw_red_count,
                    "vital_red_after_dedupe": int(bool_series(table["vital_red_flag"]).sum()),
                    "vital_red_removed_by_dedupe": raw_red_count - int(bool_series(table["vital_red_flag"]).sum()),
                    "kept_rule": "within each domain, keep the max score per variation_id, then the max score per variant_key",
                }
            ]
        ),
        data_path(output_prefix, "historical_within_domain_dedupe_summary.csv"),
    )
    table["deduplication_note"] = (
        f"within_domain_deduplicated; dropped_by_variation_id={dropped_by_id}; "
        f"dropped_by_variant_key={dropped_by_key}"
    )
    return table.drop(columns=["vital_score_numeric"], errors="ignore")


def deduplicate_for_pooled(predictions: pd.DataFrame, output_prefix: str) -> pd.DataFrame:
    priority = {
        "arrhythmia": 0,
        "cardiomyopathy": 1,
        "epilepsy": 2,
        "hearing_loss": 3,
        "random_clinvar_plp": 4,
    }
    table = add_stable_identity_columns(predictions)
    table["domain_priority"] = table["domain"].map(priority).fillna(99)
    table["vital_score_numeric"] = pd.to_numeric(table.get("vital_score"), errors="coerce").fillna(-1)
    red_before_pooled_dedupe = int(bool_series(table["vital_red_flag"]).sum())
    decision_rows: list[dict[str, object]] = []
    for key_type, key_column in [("variation_id", "variation_id"), ("variant_key", "stable_variant_key")]:
        candidates = table.loc[table[key_column].fillna("").astype(str).ne("")].copy()
        group_sizes = candidates.groupby(key_column)["domain"].nunique().reset_index(name="domain_count")
        overlap_keys = set(group_sizes.loc[group_sizes["domain_count"] > 1, key_column].astype(str))
        for key_value, group in candidates[candidates[key_column].astype(str).isin(overlap_keys)].groupby(key_column):
            ranked = group.sort_values(
                ["domain_priority", "vital_score_numeric"],
                ascending=[True, False],
            )
            kept = ranked.iloc[0]
            reason_kept = (
                "tie_break_rule_applied"
                if ranked["domain_priority"].nunique() == 1
                else "higher_priority_domain"
            )
            decision_rows.append(
                {
                    "dedupe_key_type": key_type,
                    "dedupe_key": key_value,
                    "variation_ids_present": "|".join(sorted(set(group["variation_id"].dropna().astype(str)))),
                    "variant_keys_present": "|".join(sorted(set(group["stable_variant_key"].dropna().astype(str)))),
                    "domains_present": "|".join(sorted(set(group["domain"].dropna().astype(str)))),
                    "kept_domain": kept["domain"],
                    "kept_variation_id": kept["variation_id"],
                    "kept_variant_key": kept["stable_variant_key"],
                    "kept_vital_score": kept.get("vital_score", np.nan),
                    "dropped_domains": "|".join(
                        sorted(set(ranked.iloc[1:]["domain"].dropna().astype(str)))
                    ),
                    "candidate_count": len(group),
                    "reason_kept": reason_kept,
                    "tie_break_rule": "max_vital_score_within_same_domain_priority",
                }
            )
    if decision_rows:
        save_table(
            pd.DataFrame(decision_rows),
            data_path(output_prefix, "historical_cross_domain_dedupe_decisions.csv"),
        )
    table = table.sort_values(
        ["variation_id", "domain_priority", "vital_score_numeric"],
        ascending=[True, True, False],
    )
    duplicate_id_drop_index = table[table.duplicated(["variation_id"], keep="first")].index
    dropped_by_id = table.loc[duplicate_id_drop_index].copy()
    table = table.drop(index=duplicate_id_drop_index).copy()

    keyed = table["stable_variant_key"].astype(str).ne("")
    table_keyed = table.loc[keyed].sort_values(
        ["stable_variant_key", "domain_priority", "vital_score_numeric"],
        ascending=[True, True, False],
    )
    duplicate_key_drop_index = table_keyed[table_keyed.duplicated(["stable_variant_key"], keep="first")].index
    dropped_by_key = table.loc[duplicate_key_drop_index].copy()
    dropped = pd.concat([dropped_by_id, dropped_by_key], ignore_index=True)
    if not dropped.empty:
        save_table(dropped, data_path(output_prefix, "historical_dropped_cross_domain_duplicate_rows.csv"))
    table = table.drop(index=duplicate_key_drop_index).copy()
    save_table(
        pd.DataFrame(
            [
                {
                    "dedupe_stage": "pooled_cross_domain",
                    "raw_rows_before_pooled_dedupe": len(predictions),
                    "rows_after_within_domain_input": len(add_stable_identity_columns(predictions)),
                    "rows_after_pooled_dedupe": len(table),
                    "removed_cross_domain_variation_id_rows": len(dropped_by_id),
                    "removed_cross_domain_variant_key_rows": len(dropped_by_key),
                    "vital_red_before_pooled_dedupe": red_before_pooled_dedupe,
                    "vital_red_after_pooled_dedupe": int(bool_series(table["vital_red_flag"]).sum()),
                    "vital_red_removed_by_pooled_dedupe": red_before_pooled_dedupe
                    - int(bool_series(table["vital_red_flag"]).sum()),
                    "kept_rule": (
                        "prefer domain priority arrhythmia > cardiomyopathy > epilepsy > hearing_loss > "
                        "random_clinvar_plp; within the same key/domain priority keep the max score"
                    ),
                }
            ]
        ),
        data_path(output_prefix, "historical_pooled_dedupe_summary.csv"),
    )
    table["deduplication_note"] = (
        f"pooled_deduplicated; dropped_cross_domain_variation_id={len(dropped_by_id)}; "
        f"dropped_cross_domain_variant_key={len(dropped_by_key)}"
    )
    return table.drop(columns=["domain_priority", "vital_score_numeric"], errors="ignore")


def load_current_key_lookup(followup_variant_summary: Path, variant_keys: set[str]) -> pd.DataFrame:
    if not variant_keys:
        return pd.DataFrame()
    usecols = [
        "VariationID",
        "GeneSymbol",
        "ClinicalSignificance",
        "ReviewStatus",
        "Assembly",
        "Chromosome",
        "PositionVCF",
        "ReferenceAlleleVCF",
        "AlternateAlleleVCF",
    ]
    rows: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        followup_variant_summary,
        sep="\t",
        compression="infer",
        usecols=lambda column: column in usecols,
        dtype=str,
        chunksize=250_000,
        low_memory=False,
    ):
        chunk = chunk[chunk["Assembly"].eq("GRCh38")].copy()
        if chunk.empty:
            continue
        for column in ["Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF"]:
            chunk[column] = chunk[column].fillna("").astype(str)
        chunk["stable_variant_key"] = (
            chunk["Chromosome"]
            + ":"
            + chunk["PositionVCF"]
            + ":"
            + chunk["ReferenceAlleleVCF"]
            + ":"
            + chunk["AlternateAlleleVCF"]
        )
        chunk = chunk[chunk["stable_variant_key"].isin(variant_keys)].copy()
        if chunk.empty:
            continue
        chunk["current_clinical_group_by_key"] = chunk["ClinicalSignificance"].map(clinical_group)
        rows.append(chunk)
    if not rows:
        return pd.DataFrame()
    current = pd.concat(rows, ignore_index=True)
    return (
        current.groupby("stable_variant_key", dropna=False)
        .agg(
            current_variation_ids_by_key=("VariationID", lambda values: "|".join(sorted(set(values.dropna().astype(str))))),
            current_genes_by_key=("GeneSymbol", lambda values: "|".join(sorted(set(values.dropna().astype(str))))),
            current_clinical_groups_by_key=("current_clinical_group_by_key", lambda values: "|".join(sorted(set(values.dropna().astype(str))))),
            current_clinsig_by_key=("ClinicalSignificance", lambda values: "|".join(sorted(set(values.dropna().astype(str))))),
            current_review_by_key=("ReviewStatus", lambda values: "|".join(sorted(set(values.dropna().astype(str))))),
        )
        .reset_index()
    )


def write_sanity_checks(predictions: pd.DataFrame, output_prefix: str, followup_variant_summary: Path) -> None:
    qa = add_stable_identity_columns(predictions)
    qa["is_red"] = bool_series(qa["vital_red_flag"])
    qa["is_missing_followup"] = qa["followup_clinical_group"].fillna("").astype(str).eq("missing")
    summary_rows = []
    for domain, sub in qa.groupby("domain", dropna=False):
        summary_rows.append(
            {
                "domain": domain,
                "baseline_rows_after_score_merge": len(sub),
                "unique_variation_ids": sub["variation_id"].nunique(dropna=True),
                "unique_variant_keys": sub["stable_variant_key"].replace("", pd.NA).nunique(dropna=True),
                "vital_red_count": int(sub["is_red"].sum()),
                "missing_followup_count": int(sub["is_missing_followup"].sum()),
                "duplicate_variation_id_rows": int(sub.duplicated(["variation_id"], keep=False).sum()),
                "duplicate_variant_key_rows": int(
                    sub.loc[sub["stable_variant_key"].astype(str).ne("")].duplicated(
                        ["stable_variant_key"], keep=False
                    ).sum()
                ),
            }
        )
    save_table(pd.DataFrame(summary_rows), data_path(output_prefix, "historical_sanity_domain_summary.csv"))

    duplicate_ids = qa[qa.duplicated(["domain", "variation_id"], keep=False)].copy()
    save_table(
        duplicate_ids.sort_values(["domain", "variation_id", "vital_score"], ascending=[True, True, False]),
        data_path(output_prefix, "historical_within_domain_duplicate_variation_ids.csv"),
    )
    duplicate_keys = qa[
        qa["stable_variant_key"].astype(str).ne("")
        & qa.duplicated(["domain", "stable_variant_key"], keep=False)
    ].copy()
    save_table(
        duplicate_keys.sort_values(["domain", "stable_variant_key", "vital_score"], ascending=[True, True, False]),
        data_path(output_prefix, "historical_within_domain_duplicate_variant_keys.csv"),
    )

    id_domain_counts = qa.groupby("variation_id")["domain"].nunique().reset_index(name="domain_count")
    overlapping_ids = id_domain_counts[id_domain_counts["domain_count"] > 1]
    save_table(
        qa.merge(overlapping_ids, on="variation_id", how="inner").sort_values(
            ["variation_id", "domain", "vital_score"], ascending=[True, True, False]
        ),
        data_path(output_prefix, "historical_cross_domain_variation_id_overlap.csv"),
    )

    key_domain_counts = (
        qa.loc[qa["stable_variant_key"].astype(str).ne("")]
        .groupby("stable_variant_key")["domain"]
        .nunique()
        .reset_index(name="domain_count")
    )
    overlapping_keys = key_domain_counts[key_domain_counts["domain_count"] > 1]
    save_table(
        qa.merge(overlapping_keys, on="stable_variant_key", how="inner").sort_values(
            ["stable_variant_key", "domain", "vital_score"], ascending=[True, True, False]
        ),
        data_path(output_prefix, "historical_cross_domain_variant_key_overlap.csv"),
    )

    missing_cases = qa.loc[qa["is_missing_followup"]].copy()
    current_key_lookup = load_current_key_lookup(
        followup_variant_summary,
        set(missing_cases["stable_variant_key"].dropna().astype(str)) - {""},
    )
    if not current_key_lookup.empty:
        missing_cases = missing_cases.merge(current_key_lookup, on="stable_variant_key", how="left")
    else:
        for column in [
            "current_variation_ids_by_key",
            "current_genes_by_key",
            "current_clinical_groups_by_key",
            "current_clinsig_by_key",
            "current_review_by_key",
        ]:
            missing_cases[column] = ""
    missing_cases["missing_alignment_category"] = np.select(
        [
            missing_cases["current_variation_ids_by_key"].fillna("").astype(str).ne(""),
        ],
        ["mapping_failure_or_variation_id_changed"],
        default="disappeared_completely_or_no_current_GRCh38_key",
    )
    missing_cases["manual_review_priority"] = np.select(
        [
            missing_cases["is_red"],
            pd.to_numeric(missing_cases.get("vital_score"), errors="coerce").fillna(0) >= 70,
            pd.to_numeric(missing_cases.get("vital_score"), errors="coerce").fillna(0) >= 60,
        ],
        ["red_missing", "score_ge_70_missing", "score_ge_60_missing"],
        default="other_missing",
    )
    keep = [
        "domain",
        "variation_id",
        "clinvar_id",
        "stable_gene",
        "stable_variant_key",
        "name",
        "clinical_group",
        "followup_clinical_group",
        "clinical_significance",
        "followup_clinical_significance",
        "review_status_baseline",
        "followup_review_status",
        "current_variation_ids_by_key",
        "current_clinical_groups_by_key",
        "current_clinsig_by_key",
        "missing_alignment_category",
        "vital_score",
        "vital_band",
        "vital_red_flag",
        "manual_review_priority",
    ]
    keep = [column for column in keep if column in missing_cases.columns]
    save_table(
        missing_cases.sort_values(["manual_review_priority", "vital_score"], ascending=[True, False]).loc[:, keep],
        data_path(output_prefix, "historical_missing_cases_for_manual_review.csv"),
    )
    alignment_rows = []
    for domain, sub in qa.groupby("domain", dropna=False):
        changed_class = sub.loc[bool_series(sub.get("broad_revised_or_destabilized", pd.Series(False, index=sub.index)))]
        alignment_rows.append(
            {
                "domain": domain,
                "baseline_rows_after_score_merge": len(sub),
                "changed_classification_or_missing_broad": len(changed_class),
                "missing_by_variation_id": int(sub["is_missing_followup"].sum()),
                "mapping_failure_or_variation_id_changed": int(
                    missing_cases.loc[
                        missing_cases["domain"].eq(domain)
                        & missing_cases["missing_alignment_category"].eq("mapping_failure_or_variation_id_changed")
                    ].shape[0]
                ),
                "disappeared_completely_or_no_current_GRCh38_key": int(
                    missing_cases.loc[
                        missing_cases["domain"].eq(domain)
                        & missing_cases["missing_alignment_category"].eq(
                            "disappeared_completely_or_no_current_GRCh38_key"
                        )
                    ].shape[0]
                ),
            }
        )
    save_table(pd.DataFrame(alignment_rows), data_path(output_prefix, "historical_alignment_breakdown.csv"))


def summarize_combined(
    predictions: pd.DataFrame,
    output_prefix: str,
    followup_variant_summary: Path,
    bootstrap_iterations: int,
    random_seed: int,
) -> None:
    if predictions.empty:
        raise ValueError("No combined historical predictions were produced.")
    save_table(predictions, data_path(output_prefix, "historical_predictions_raw.csv"))
    write_sanity_checks(predictions, output_prefix, followup_variant_summary)
    within_domain_predictions = deduplicate_for_analysis(predictions, output_prefix)
    pooled_predictions = deduplicate_for_pooled(within_domain_predictions, output_prefix)
    save_table(within_domain_predictions, data_path(output_prefix, "historical_predictions.csv"))
    save_table(pooled_predictions, data_path(output_prefix, "historical_predictions_pooled_deduplicated.csv"))
    targets = [
        ("strict_PLP_to_BLB_or_VUS", "strict_revised_to_b_or_vus"),
        ("broad_PLP_to_non_PLP_or_missing", "broad_revised_or_destabilized"),
        ("expanded_PLP_to_non_PLP_or_review_change", "expanded_revised_or_review_changed"),
    ]
    metrics_rows: list[dict[str, object]] = []
    bootstrap_summaries: list[dict[str, object]] = []
    for domain in sorted(within_domain_predictions["domain"].dropna().unique()):
        sub = within_domain_predictions[within_domain_predictions["domain"].eq(domain)].copy()
        for target_name, column in targets:
            if column not in predictions.columns:
                continue
            metrics_rows.append(metric_row(sub, domain, target_name, column))
    for target_name, column in targets:
        if column not in predictions.columns:
            continue
        pooled = metric_row(pooled_predictions, "combined_pooled", target_name, column)
        pooled["interpretation_note"] = (
            "Pooled estimate is shown alongside within-domain rows because disease architectures, "
            "baseline event rates, and error mechanisms differ by domain."
        )
        metrics_rows.append(pooled)
        draws, summary = stratified_bootstrap(
            pooled_predictions,
            target_column=column,
            n_iter=bootstrap_iterations,
            random_seed=random_seed,
        )
        draws.insert(0, "target", target_name)
        save_table(draws, data_path(output_prefix, f"historical_{target_name}_stratified_bootstrap.csv"))
        summary["target"] = target_name
        bootstrap_summaries.append(summary)
    save_table(pd.DataFrame(metrics_rows), data_path(output_prefix, "historical_endpoint_summary.csv"))
    save_table(pd.DataFrame(bootstrap_summaries), data_path(output_prefix, "historical_stratified_bootstrap_summary.csv"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run multi-domain historical validation with strict, broad, and expanded endpoints."
    )
    parser.add_argument("--baseline-variant-summary", type=Path, default=BASE_DIR / "data" / "raw" / "variant_summary_2023-01.txt.gz")
    parser.add_argument("--followup-variant-summary", type=Path, default=BASE_DIR / "data" / "variant_summary.txt.gz")
    parser.add_argument("--output-prefix", default="combined_2023_01_to_current_vital")
    parser.add_argument("--dataset", default="gnomad_r4")
    parser.add_argument("--gnomad-pause", type=float, default=0.2)
    parser.add_argument("--no-fetch-gnomad", action="store_true")
    parser.add_argument("--force-gnomad-fetch", action="store_true")
    parser.add_argument("--force-score", action="store_true")
    parser.add_argument("--skip-scoring", action="store_true")
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--random-seed", type=int, default=101)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    combined_predictions: list[pd.DataFrame] = []
    for config in DOMAIN_CONFIGS:
        score_table = data_path(str(config["baseline_prefix"]), "vital_scores.csv")
        if not args.skip_scoring:
            score_table = ensure_score_table(
                config,
                baseline_variant_summary=args.baseline_variant_summary,
                dataset=args.dataset,
                gnomad_pause=args.gnomad_pause,
                fetch_gnomad=not args.no_fetch_gnomad,
                force_score=args.force_score,
                force_gnomad_fetch=args.force_gnomad_fetch,
            )
        elif not score_table.exists():
            raise FileNotFoundError(f"Missing score table for {config['domain']}: {score_table}")
        combined_predictions.append(
            run_domain_validation(
                config,
                score_table=score_table,
                baseline_variant_summary=args.baseline_variant_summary,
                followup_variant_summary=args.followup_variant_summary,
            )
        )
    summarize_combined(
        pd.concat(combined_predictions, ignore_index=True),
        args.output_prefix,
        followup_variant_summary=args.followup_variant_summary,
        bootstrap_iterations=args.bootstrap_iterations,
        random_seed=args.random_seed,
    )


if __name__ == "__main__":
    main()
