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
VITAL_ACTION_THRESHOLD = 70

COMPONENT_COLUMNS = {
    "frequency_pressure_score": 45.0,
    "ac_reliability_score": 20.0,
    "popmax_enrichment_score": 10.0,
    "variant_type_tension_score": 6.0,
    "technical_detectability_score": 8.0,
    "gene_constraint_score": 10.0,
    "review_fragility_score": 10.0,
}

WEIGHT_PROFILES = {
    "primary_expert_weights": {
        "frequency_pressure_score": 45.0,
        "ac_reliability_score": 20.0,
        "popmax_enrichment_score": 10.0,
        "variant_type_tension_score": 6.0,
        "technical_detectability_score": 8.0,
        "gene_constraint_score": 10.0,
        "review_fragility_score": 10.0,
    },
    "balanced_equal_components": {
        "frequency_pressure_score": 15.0,
        "ac_reliability_score": 15.0,
        "popmax_enrichment_score": 15.0,
        "variant_type_tension_score": 15.0,
        "technical_detectability_score": 15.0,
        "gene_constraint_score": 15.0,
        "review_fragility_score": 15.0,
    },
    "frequency_dominant": {
        "frequency_pressure_score": 55.0,
        "ac_reliability_score": 20.0,
        "popmax_enrichment_score": 15.0,
        "variant_type_tension_score": 4.0,
        "technical_detectability_score": 4.0,
        "gene_constraint_score": 5.0,
        "review_fragility_score": 6.0,
    },
    "review_fragility_dominant": {
        "frequency_pressure_score": 35.0,
        "ac_reliability_score": 20.0,
        "popmax_enrichment_score": 8.0,
        "variant_type_tension_score": 5.0,
        "technical_detectability_score": 6.0,
        "gene_constraint_score": 8.0,
        "review_fragility_score": 25.0,
    },
    "technical_detectability_dominant": {
        "frequency_pressure_score": 35.0,
        "ac_reliability_score": 20.0,
        "popmax_enrichment_score": 8.0,
        "variant_type_tension_score": 10.0,
        "technical_detectability_score": 20.0,
        "gene_constraint_score": 8.0,
        "review_fragility_score": 10.0,
    },
    "reduced_af_pressure": {
        "frequency_pressure_score": 30.0,
        "ac_reliability_score": 20.0,
        "popmax_enrichment_score": 10.0,
        "variant_type_tension_score": 8.0,
        "technical_detectability_score": 10.0,
        "gene_constraint_score": 12.0,
        "review_fragility_score": 15.0,
    },
}


def data_path(prefix: str, name: str) -> Path:
    return DATA_DIR / f"{prefix}_{name}" if prefix else DATA_DIR / name


def save_table(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"Saved {output_path} ({len(df)} rows)")


def bool_series(series: pd.Series) -> pd.Series:
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def rescore(table: pd.DataFrame, weights: dict[str, float]) -> pd.Series:
    total = pd.Series(0.0, index=table.index)
    for column, original_max in COMPONENT_COLUMNS.items():
        values = pd.to_numeric(table.get(column), errors="coerce").fillna(0)
        normalized = values / original_max if original_max else 0
        total = total + normalized.clip(lower=0, upper=1) * weights[column]
    return total.clip(lower=0, upper=100)


def spearman_correlation(left: pd.Series, right: pd.Series) -> float:
    pair = pd.DataFrame(
        {
            "left": pd.to_numeric(left, errors="coerce"),
            "right": pd.to_numeric(right, errors="coerce"),
        }
    ).dropna()
    if len(pair) < 2:
        return np.nan
    return float(pair["left"].rank(method="average").corr(pair["right"].rank(method="average")))


def variant_label(row: pd.Series) -> str:
    for column in ["clinvar_id", "variation_id", "variant_key"]:
        value = row.get(column)
        if pd.notna(value) and str(value):
            return str(value)
    return ""


def run_weight_sensitivity(score_table: Path, output_prefix: str) -> None:
    scores = pd.read_csv(score_table)
    scores["variant_label"] = scores.apply(variant_label, axis=1)
    primary_red = bool_series(scores.get("vital_red_flag", pd.Series(False, index=scores.index)))
    standard_flag = bool_series(scores.get("standard_acmg_frequency_flag", pd.Series(False, index=scores.index)))
    ac_supported = bool_series(scores.get("frequency_signal_ac_ge_20", pd.Series(False, index=scores.index)))
    weak_review = bool_series(scores.get("weak_review_signal", pd.Series(False, index=scores.index)))
    review_score = pd.to_numeric(scores.get("review_score"), errors="coerce").fillna(0)
    retained_proxy_negative = review_score.ge(2) & ~ac_supported
    primary_red_set = set(scores.loc[primary_red, "variant_label"])

    summary_rows: list[dict[str, object]] = []
    scored_rows: list[pd.DataFrame] = []
    profile_scores_by_name: dict[str, pd.Series] = {}
    profile_red_by_name: dict[str, pd.Series] = {}
    for profile_name, weights in WEIGHT_PROFILES.items():
        alt_score = rescore(scores, weights)
        alt_red = alt_score.ge(VITAL_ACTION_THRESHOLD) & ac_supported & weak_review
        profile_scores_by_name[profile_name] = alt_score
        profile_red_by_name[profile_name] = alt_red
        alt_red_set = set(scores.loc[alt_red, "variant_label"])
        shared = sorted(primary_red_set & alt_red_set)
        gained = sorted(alt_red_set - primary_red_set)
        lost = sorted(primary_red_set - alt_red_set)
        proxy_fp = int((alt_red & retained_proxy_negative).sum())
        summary_rows.append(
            {
                "weight_profile": profile_name,
                "profile_role": (
                    "primary" if profile_name == "primary_expert_weights" else "sensitivity_profile"
                ),
                "component_weight_total": sum(weights.values()),
                "vital_red_count": int(alt_red.sum()),
                "shared_with_primary_red_count": len(shared),
                "gained_vs_primary_count": len(gained),
                "lost_vs_primary_count": len(lost),
                "shared_with_primary_red_variants": "; ".join(shared),
                "gained_vs_primary_variants": "; ".join(gained),
                "lost_vs_primary_variants": "; ".join(lost),
                "naive_af_flag_count": int(standard_flag.sum()),
                "ac_supported_flag_count": int(ac_supported.sum()),
                "compression_vs_naive_percent": (
                    (1 - (int(alt_red.sum()) / int(standard_flag.sum()))) * 100
                    if int(standard_flag.sum())
                    else np.nan
                ),
                "proxy_false_positive_count": proxy_fp,
                "max_score": float(alt_score.max()) if alt_score.notna().any() else np.nan,
                "p95_score": float(alt_score.quantile(0.95)) if alt_score.notna().any() else np.nan,
                "spearman_rank_correlation_with_primary": spearman_correlation(scores["vital_score"], alt_score),
                "weighting_note": "expert_weight_sensitivity_not_empirical_refitting",
            }
        )
        profile_scores = scores.loc[
            alt_red,
            [
                "variant_label",
                "gene",
                "clinvar_id",
                "variation_id",
                "variant_key",
                "vital_score",
                "vital_band",
                "vital_signal_reason",
            ],
        ].copy()
        profile_scores.insert(0, "weight_profile", profile_name)
        profile_scores["sensitivity_vital_score"] = alt_score.loc[alt_red].to_numpy()
        profile_scores["relation_to_primary"] = np.where(
            profile_scores["variant_label"].isin(primary_red_set),
            "shared_with_primary",
            "gained_vs_primary",
        )
        scored_rows.append(profile_scores)

    save_table(pd.DataFrame(summary_rows), data_path(output_prefix, "vital_weight_sensitivity_summary.csv"))
    save_table(
        pd.concat(scored_rows, ignore_index=True) if scored_rows else pd.DataFrame(),
        data_path(output_prefix, "vital_weight_sensitivity_red_sets.csv"),
    )

    alternative_profiles = [name for name in WEIGHT_PROFILES if name != "primary_expert_weights"]
    matrix_rows: list[dict[str, object]] = []
    variant_profile_rows: list[dict[str, object]] = []
    for index, row in scores.loc[primary_red].iterrows():
        alt_red_count = int(sum(bool(profile_red_by_name[name].loc[index]) for name in alternative_profiles))
        if alt_red_count == len(alternative_profiles):
            stability_class = "anchor"
        elif alt_red_count >= len(alternative_profiles) - 1:
            stability_class = "near-stable"
        else:
            stability_class = "borderline"
        matrix_row: dict[str, object] = {
            "variant_label": row["variant_label"],
            "gene": row.get("gene", ""),
            "clinvar_id": row.get("clinvar_id", ""),
            "variation_id": row.get("variation_id", ""),
            "variant_key": row.get("variant_key", ""),
            "primary_vital_score": row.get("vital_score", np.nan),
            "alternative_profiles_tested": len(alternative_profiles),
            "alternative_profiles_retained_red": alt_red_count,
            "alternative_profiles_lost_red": len(alternative_profiles) - alt_red_count,
            "stability_class": stability_class,
            "stability_rule": "anchor=retained_in_all_alternative_profiles; near-stable=retained_in_4_of_5; borderline=retained_in_0_to_3_of_5",
        }
        for profile_name in WEIGHT_PROFILES:
            matrix_row[f"{profile_name}_score"] = profile_scores_by_name[profile_name].loc[index]
            matrix_row[f"{profile_name}_red"] = bool(profile_red_by_name[profile_name].loc[index])
        matrix_rows.append(matrix_row)
        for profile_name in alternative_profiles:
            weights = WEIGHT_PROFILES[profile_name]
            retained = bool(profile_red_by_name[profile_name].loc[index])
            variant_profile_rows.append(
                {
                    "variant_label": row["variant_label"],
                    "gene": row.get("gene", ""),
                    "clinvar_id": row.get("clinvar_id", ""),
                    "variation_id": row.get("variation_id", ""),
                    "primary_vital_score": row.get("vital_score", np.nan),
                    "weight_profile": profile_name,
                    "sensitivity_vital_score": profile_scores_by_name[profile_name].loc[index],
                    "retained_red_under_profile": retained,
                    "red_set_change": "retained" if retained else "lost_vs_primary",
                    "stability_class": stability_class,
                    "frequency_pressure_weight": weights["frequency_pressure_score"],
                    "ac_reliability_weight": weights["ac_reliability_score"],
                    "popmax_enrichment_weight": weights["popmax_enrichment_score"],
                    "variant_type_tension_weight": weights["variant_type_tension_score"],
                    "technical_detectability_weight": weights["technical_detectability_score"],
                    "gene_constraint_weight": weights["gene_constraint_score"],
                    "review_fragility_weight": weights["review_fragility_score"],
                    "weighting_note": "five_alternative_profiles_no_outcome_refitting",
                }
            )
    matrix = pd.DataFrame(matrix_rows)
    variant_profile_table = pd.DataFrame(variant_profile_rows)
    save_table(matrix, data_path(output_prefix, "vital_weight_sensitivity_variant_matrix.csv"))
    save_table(
        variant_profile_table,
        data_path(output_prefix, "vital_weight_sensitivity_variant_profile_table.csv"),
    )

    supplementary_dir = BASE_DIR / "supplementary_tables"
    supplementary_dir.mkdir(parents=True, exist_ok=True)
    variant_profile_table.to_csv(
        supplementary_dir / "Supplementary_Table_S13_VITAL_weight_sensitivity.tsv",
        sep="\t",
        index=False,
    )
    print(
        f"Saved {supplementary_dir / 'Supplementary_Table_S13_VITAL_weight_sensitivity.tsv'} "
        f"({len(variant_profile_table)} rows)"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run VITAL expert-weight sensitivity profiles.")
    parser.add_argument("--score-table", type=Path, default=data_path("arrhythmia", "vital_scores.csv"))
    parser.add_argument("--output-prefix", default="arrhythmia")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_weight_sensitivity(args.score_table, args.output_prefix)


if __name__ == "__main__":
    main()
