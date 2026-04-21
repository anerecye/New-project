from __future__ import annotations

import argparse
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
except ModuleNotFoundError as exc:
    raise SystemExit("This script requires pandas and numpy.") from exc


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
AF_ULTRA_RARE = 1e-5
MIN_RECLASSIFICATION_AC = 20


TARGETS = {
    "strict_PLP_to_BLB_or_VUS": "strict_revised_to_b_or_vus",
    "broad_PLP_to_non_PLP_or_missing": "broad_revised_or_destabilized",
    "expanded_PLP_to_non_PLP_or_review_change": "expanded_revised_or_review_changed",
}


def data_path(prefix: str, name: str) -> Path:
    return DATA_DIR / f"{prefix}_{name}" if prefix else DATA_DIR / name


def save_table(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"Saved {output_path} ({len(df)} rows)")


def bool_series(series: pd.Series) -> pd.Series:
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def numeric_series(df: pd.DataFrame, column: str, default: float = np.nan) -> pd.Series:
    if column not in df.columns:
        return pd.Series(default, index=df.index)
    return pd.to_numeric(df[column], errors="coerce")


def derive_ac_supported(df: pd.DataFrame) -> tuple[pd.Series, str]:
    if "frequency_signal_ac_ge_20" in df.columns:
        return bool_series(df["frequency_signal_ac_ge_20"]), "frequency_signal_ac_ge_20"
    qualifying_ac = numeric_series(df, "qualifying_frequency_ac").fillna(0)
    if "max_frequency_signal" in df.columns:
        frequency_signal = numeric_series(df, "max_frequency_signal").fillna(0) > AF_ULTRA_RARE
        return (
            frequency_signal & qualifying_ac.ge(MIN_RECLASSIFICATION_AC),
            "max_frequency_signal_gt_1e-5_and_qualifying_frequency_ac_ge_20",
        )
    return qualifying_ac.ge(MIN_RECLASSIFICATION_AC), "qualifying_frequency_ac_ge_20"


def derive_weak_review(df: pd.DataFrame) -> tuple[pd.Series, str]:
    if "weak_review_signal" in df.columns:
        return bool_series(df["weak_review_signal"]), "weak_review_signal"
    review_score = numeric_series(df, "review_score")
    submitter_count = numeric_series(df, "submitter_count")
    if review_score.notna().any() or submitter_count.notna().any():
        weak_review = review_score.le(1).fillna(False) | submitter_count.le(1).fillna(False)
        return weak_review, "review_score_le_1_or_submitter_count_le_1"
    return pd.Series(False, index=df.index), "missing_weak_review_fields_default_false"


def safe_rate(numerator: int | float, denominator: int | float) -> float:
    return float(numerator / denominator) if denominator else np.nan


def calibration_rows(df: pd.DataFrame, domain: str, thresholds: list[int]) -> list[dict[str, object]]:
    score = pd.to_numeric(df.get("vital_score"), errors="coerce").fillna(-np.inf)
    ac_supported, ac_mapping_source = derive_ac_supported(df)
    weak_review, weak_review_mapping_source = derive_weak_review(df)
    rows: list[dict[str, object]] = []
    for target_name, target_column in TARGETS.items():
        if target_column not in df.columns:
            continue
        labels = bool_series(df[target_column])
        total_n = len(df)
        event_count = int(labels.sum())
        overall_rate = safe_rate(event_count, total_n)
        for calibration_mode in ["score_only", "red_gate_compatible"]:
            for threshold in thresholds:
                predicted = score.ge(threshold)
                if calibration_mode == "red_gate_compatible":
                    predicted = predicted & ac_supported & weak_review
                flagged_count = int(predicted.sum())
                flagged_events = int((predicted & labels).sum())
                flagged_event_rate = safe_rate(flagged_events, flagged_count)
                rows.append(
                    {
                        "domain": domain,
                        "target": target_name,
                        "calibration_mode": calibration_mode,
                        "score_threshold": threshold,
                        "baseline_variant_count": total_n,
                        "target_event_count": event_count,
                        "overall_event_rate": overall_rate,
                        "flagged_count": flagged_count,
                        "flagged_event_count": flagged_events,
                        "flagged_event_rate": flagged_event_rate,
                        "recall": safe_rate(flagged_events, event_count),
                        "enrichment_vs_overall": safe_rate(flagged_event_rate, overall_rate),
                        "is_prespecified_threshold_70": threshold == 70,
                        "ac_supported_flag_count": int(ac_supported.sum()),
                        "weak_review_flag_count": int(weak_review.sum()),
                        "ac_mapping_source": ac_mapping_source,
                        "weak_review_mapping_source": weak_review_mapping_source,
                        "field_mapping_note": (
                            "Red-gate calibration derives missing actionability fields from validator logic; "
                            "AC threshold is an actionability gate, not a scoring component."
                        ),
                    }
                )
    return rows


def field_mapping_audit_rows(df: pd.DataFrame, domain: str) -> list[dict[str, object]]:
    score = pd.to_numeric(df.get("vital_score"), errors="coerce").fillna(-np.inf)
    ac_supported, ac_mapping_source = derive_ac_supported(df)
    weak_review, weak_review_mapping_source = derive_weak_review(df)
    derived_red_at_70 = score.ge(70) & ac_supported & weak_review
    recorded_red = (
        bool_series(df["vital_red_flag"])
        if "vital_red_flag" in df.columns
        else pd.Series(False, index=df.index)
    )
    kcn_h2_mask = pd.Series(False, index=df.index)
    if "gene" in df.columns:
        kcn_h2_mask = df["gene"].fillna("").astype(str).eq("KCNH2")
    elif "gene_score" in df.columns:
        kcn_h2_mask = df["gene_score"].fillna("").astype(str).eq("KCNH2")
    return [
        {
            "domain": domain,
            "row_count": len(df),
            "frequency_signal_ac_ge_20_present": "frequency_signal_ac_ge_20" in df.columns,
            "weak_review_signal_present": "weak_review_signal" in df.columns,
            "qualifying_frequency_ac_present": "qualifying_frequency_ac" in df.columns,
            "review_score_present": "review_score" in df.columns,
            "submitter_count_present": "submitter_count" in df.columns,
            "ac_mapping_source": ac_mapping_source,
            "weak_review_mapping_source": weak_review_mapping_source,
            "ac_supported_count": int(ac_supported.sum()),
            "weak_review_count": int(weak_review.sum()),
            "recorded_vital_red_flag_count": int(recorded_red.sum()),
            "derived_red_at_threshold_70_count": int(derived_red_at_70.sum()),
            "derived_red_at_threshold_70_variants": variant_summary(df.loc[derived_red_at_70]),
            "kcn_h2_recorded_red_count": int((recorded_red & kcn_h2_mask).sum()),
            "kcn_h2_derived_red_at_70_count": int((derived_red_at_70 & kcn_h2_mask).sum()),
            "kcn_h2_derived_red_at_70_variants": variant_summary(df.loc[derived_red_at_70 & kcn_h2_mask]),
            "audit_note": (
                "Before this fix, absent weak_review_signal defaulted to False in threshold calibration, "
                "yielding zero red-gate-compatible flags despite recorded vital_red_flag rows."
            ),
        }
    ]


def variant_summary(df: pd.DataFrame) -> str:
    if df.empty:
        return ""
    labels = []
    for _, row in df.iterrows():
        gene = row.get("gene", row.get("gene_score", ""))
        clinvar_id = row.get("clinvar_id", row.get("variation_id", ""))
        score = pd.to_numeric(pd.Series([row.get("vital_score")]), errors="coerce").iloc[0]
        score_label = f"{score:.1f}" if pd.notna(score) else "NA"
        labels.append(f"{gene}:{clinvar_id}:VITAL={score_label}")
    return "; ".join(labels)


def summarize_best(calibration: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    grouped = calibration.copy()
    grouped["eligible_for_best"] = grouped["flagged_count"].ge(1) & grouped["enrichment_vs_overall"].notna()
    for (domain, target, mode), sub in grouped.groupby(["domain", "target", "calibration_mode"], dropna=False):
        eligible = sub[sub["eligible_for_best"]].copy()
        threshold_70 = sub[sub["score_threshold"].eq(70)].iloc[0] if not sub[sub["score_threshold"].eq(70)].empty else None
        if eligible.empty:
            best = None
        else:
            best = eligible.sort_values(
                ["enrichment_vs_overall", "flagged_event_count", "flagged_count"],
                ascending=[False, False, True],
            ).iloc[0]
        rows.append(
            {
                "domain": domain,
                "target": target,
                "calibration_mode": mode,
                "best_threshold_by_enrichment": best["score_threshold"] if best is not None else np.nan,
                "best_threshold_flagged_count": best["flagged_count"] if best is not None else np.nan,
                "best_threshold_flagged_event_count": best["flagged_event_count"] if best is not None else np.nan,
                "best_threshold_enrichment_vs_overall": best["enrichment_vs_overall"] if best is not None else np.nan,
                "threshold_70_flagged_count": threshold_70["flagged_count"] if threshold_70 is not None else np.nan,
                "threshold_70_flagged_event_count": threshold_70["flagged_event_count"] if threshold_70 is not None else np.nan,
                "threshold_70_enrichment_vs_overall": threshold_70["enrichment_vs_overall"] if threshold_70 is not None else np.nan,
                "threshold_70_matches_best": (
                    bool(best is not None and threshold_70 is not None and best["score_threshold"] == 70)
                ),
                "calibration_note": (
                    "Exploratory historical threshold calibration; sparse red sets make maxima unstable."
                ),
            }
        )
    return pd.DataFrame(rows)


def summarize_for_mapping_compare(calibration: pd.DataFrame, label: str) -> pd.DataFrame:
    table = calibration[calibration["calibration_mode"].eq("red_gate_compatible")].copy()
    if table.empty:
        return pd.DataFrame()
    threshold_70 = table[table["score_threshold"].eq(70)].copy()
    threshold_70 = threshold_70.rename(
        columns={
            "flagged_count": f"{label}_threshold_70_flagged_count",
            "flagged_event_count": f"{label}_threshold_70_flagged_event_count",
            "enrichment_vs_overall": f"{label}_threshold_70_enrichment_vs_overall",
        }
    )[
        [
            "domain",
            "target",
            "calibration_mode",
            f"{label}_threshold_70_flagged_count",
            f"{label}_threshold_70_flagged_event_count",
            f"{label}_threshold_70_enrichment_vs_overall",
        ]
    ]
    maxima = (
        table.groupby(["domain", "target", "calibration_mode"], dropna=False)
        .agg(
            **{
                f"{label}_max_flagged_count": ("flagged_count", "max"),
                f"{label}_max_flagged_event_count": ("flagged_event_count", "max"),
            }
        )
        .reset_index()
    )
    return maxima.merge(threshold_70, on=["domain", "target", "calibration_mode"], how="left")


def save_mapping_fix_comparison(calibration: pd.DataFrame, output_prefix: str) -> None:
    before_path = data_path(output_prefix, "historical_threshold_calibration_field_mapping_inconsistent.csv")
    if not before_path.exists():
        return
    before = pd.read_csv(before_path)
    before_summary = summarize_for_mapping_compare(before, "before_fix")
    after_summary = summarize_for_mapping_compare(calibration, "after_fix")
    if before_summary.empty or after_summary.empty:
        return
    comparison = before_summary.merge(
        after_summary,
        on=["domain", "target", "calibration_mode"],
        how="outer",
    )
    comparison["mapping_fix"] = (
        "derive weak_review from review_score<=1 or submitter_count<=1; "
        "fallback AC support from qualifying_frequency_ac>=20 when frequency_signal_ac_ge_20 is absent"
    )
    comparison["interpretation"] = (
        "The field-mapping-inconsistent run is retained as a negative-control audit artifact; "
        "the corrected run follows validator actionability-gate logic."
    )
    save_table(comparison, data_path(output_prefix, "historical_threshold_calibration_mapping_fix_summary.csv"))


def run_calibration(prediction_table: Path, output_prefix: str) -> None:
    predictions = pd.read_csv(prediction_table)
    thresholds = list(range(40, 96, 5))
    rows: list[dict[str, object]] = []
    audit_rows: list[dict[str, object]] = []
    if "domain" in predictions.columns:
        for domain, sub in predictions.groupby("domain", dropna=False):
            rows.extend(calibration_rows(sub.copy(), str(domain), thresholds))
            audit_rows.extend(field_mapping_audit_rows(sub.copy(), str(domain)))
    rows.extend(calibration_rows(predictions, "combined_or_single_table", thresholds))
    audit_rows.extend(field_mapping_audit_rows(predictions, "combined_or_single_table"))
    calibration = pd.DataFrame(rows)
    save_table(calibration, data_path(output_prefix, "historical_threshold_calibration.csv"))
    save_table(summarize_best(calibration), data_path(output_prefix, "historical_threshold_calibration_best.csv"))
    save_table(
        pd.DataFrame(audit_rows),
        data_path(output_prefix, "historical_threshold_calibration_field_mapping_audit.csv"),
    )
    save_mapping_fix_comparison(calibration, output_prefix)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Explore historical VITAL threshold calibration.")
    parser.add_argument("--prediction-table", type=Path, required=True)
    parser.add_argument("--output-prefix", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_calibration(args.prediction_table, args.output_prefix)


if __name__ == "__main__":
    main()
