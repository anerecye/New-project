from __future__ import annotations

import argparse
import re
from pathlib import Path

try:
    import pandas as pd
except ModuleNotFoundError as exc:
    raise SystemExit("This script requires pandas.") from exc


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"


def data_path(prefix: str, name: str) -> Path:
    return DATA_DIR / f"{prefix}_{name}" if prefix else DATA_DIR / name


def save_table(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"Saved {output_path} ({len(df)} rows)")


def bool_series(series: pd.Series) -> pd.Series:
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def tension_band(score: float) -> str:
    if score < 20:
        return "00_20_low_tension"
    if score < 40:
        return "20_40_moderate_tension"
    if score < 60:
        return "40_60_elevated_tension"
    if score < 80:
        return "60_80_high_tension"
    return "80_100_extreme_tension"


def tension_zone(row: pd.Series) -> str:
    if bool(row.get("vital_red_flag_bool", False)):
        return f"{tension_band(row.get('vital_score_numeric', 0))}_red_actionable"
    return tension_band(row.get("vital_score_numeric", 0))


def functional_discordance_class(row: pd.Series) -> str:
    functional_class = str(row.get("functional_class", "")).lower()
    frequency_signal = bool(row.get("standard_acmg_frequency_flag_bool", False))
    if not frequency_signal:
        return "not_high_tension"
    if functional_class == "lof":
        return "frequency_LOF_discordant"
    if functional_class == "splice_or_intronic":
        return "frequency_splice_region_discordant"
    if functional_class == "other_or_unresolved":
        return "frequency_unresolved_or_composite_discordant"
    if functional_class == "missense":
        return "frequency_missense_discordant"
    return "frequency_function_discordant_other"


def classify_lof_subtype(row: pd.Series) -> str:
    if str(row.get("functional_class", "")).lower() != "lof":
        return "non_LOF"
    title = str(row.get("title", "")).lower()
    if re.search(r"\bp\.[^)\s;]*fs\b", title) or "frameshift" in title:
        return "frameshift"
    if re.search(r"\bc\.[^()\s;]*[+-][12][acgt]>", title) or "splice donor" in title or "splice acceptor" in title:
        return "canonical_splice"
    if re.search(r"\bp\.[^)\s;]*(?:ter|\*)", title) or "stop gained" in title or "nonsense" in title:
        return "stop_gained"
    return "other_LOF"


def format_variant_list(df: pd.DataFrame) -> str:
    if df.empty:
        return ""
    labels = []
    for _, row in df.sort_values("vital_score_numeric", ascending=False).iterrows():
        labels.append(
            f"{row.get('gene', '')}:{row.get('clinvar_id', '')}:"
            f"{row.get('lof_subtype', row.get('functional_class', ''))}:"
            f"score={row.get('vital_score_numeric', 0):.1f}"
        )
    return "; ".join(labels)


def summarize_bands(scored: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    ordered_bands = [
        "00_20_low_tension",
        "20_40_moderate_tension",
        "40_60_elevated_tension",
        "60_80_high_tension",
        "80_100_extreme_tension",
    ]
    for band in ordered_bands:
        subset = scored.loc[scored["vital_score_band_20pt"].eq(band)].copy()
        if subset.empty:
            rows.append(
                {
                    "vital_score_band_20pt": band,
                    "af_observed_variant_count": len(scored),
                    "band_variant_count": 0,
                    "band_variant_percent": 0,
                    "weak_or_single_review_percent": 0,
                    "popmax_enriched_percent": 0,
                    "ac_supported_frequency_percent": 0,
                    "indel_or_duplication_percent": 0,
                    "canonical_or_atypical_function_percent": 0,
                    "median_max_frequency_signal": pd.NA,
                    "median_qualifying_ac": pd.NA,
                    "median_gene_frequency_constraint_proxy": pd.NA,
                    "band_variants": "",
                    "interpretation_note": "Empty band in this dataset.",
                }
            )
            continue
        weak_or_single = subset["review_strength"].isin(["single_submitter", "weak_or_no_assertion"])
        popmax_ratio_ge_3 = subset["popmax_global_ratio"].fillna(0).ge(3)
        popmax_ratio_ge_10 = subset["popmax_global_ratio"].fillna(0).ge(10)
        popmax_high_frequency = subset["popmax_af"].fillna(0).gt(1e-4)
        ac_supported = subset["frequency_signal_ac_ge_20_bool"]
        indel_or_duplication = subset["variant_type"].isin(["deletion", "insertion", "duplication"])
        canonical_or_atypical = subset["functional_class"].isin(
            ["LOF", "splice_or_intronic", "other_or_unresolved"]
        )
        frameshift = subset["lof_subtype"].eq("frameshift")
        stop_gained = subset["lof_subtype"].eq("stop_gained")
        canonical_splice = subset["lof_subtype"].eq("canonical_splice")
        other_lof = subset["lof_subtype"].eq("other_LOF")
        rows.append(
            {
                "vital_score_band_20pt": band,
                "af_observed_variant_count": len(scored),
                "band_variant_count": len(subset),
                "band_variant_percent": 100 * len(subset) / len(scored) if len(scored) else 0,
                "weak_or_single_review_count": int(weak_or_single.sum()),
                "weak_or_single_review_percent": 100 * weak_or_single.mean(),
                "popmax_ratio_ge_3_count": int(popmax_ratio_ge_3.sum()),
                "popmax_ratio_ge_3_percent": 100 * popmax_ratio_ge_3.mean(),
                "popmax_ratio_ge_10_count": int(popmax_ratio_ge_10.sum()),
                "popmax_ratio_ge_10_percent": 100 * popmax_ratio_ge_10.mean(),
                "popmax_af_gt_1e_4_count": int(popmax_high_frequency.sum()),
                "popmax_af_gt_1e_4_percent": 100 * popmax_high_frequency.mean(),
                "ac_supported_frequency_count": int(ac_supported.sum()),
                "ac_supported_frequency_percent": 100 * ac_supported.mean(),
                "indel_or_duplication_count": int(indel_or_duplication.sum()),
                "indel_or_duplication_percent": 100 * indel_or_duplication.mean(),
                "lof_count": int(subset["functional_class"].eq("LOF").sum()),
                "frameshift_count": int(frameshift.sum()),
                "frameshift_percent": 100 * frameshift.mean(),
                "stop_gained_count": int(stop_gained.sum()),
                "stop_gained_percent": 100 * stop_gained.mean(),
                "canonical_splice_count": int(canonical_splice.sum()),
                "canonical_splice_percent": 100 * canonical_splice.mean(),
                "other_lof_count": int(other_lof.sum()),
                "other_lof_percent": 100 * other_lof.mean(),
                "splice_or_intronic_count": int(
                    subset["functional_class"].eq("splice_or_intronic").sum()
                ),
                "other_or_unresolved_count": int(
                    subset["functional_class"].eq("other_or_unresolved").sum()
                ),
                "missense_count": int(subset["functional_class"].eq("missense").sum()),
                "canonical_or_atypical_function_count": int(
                    canonical_or_atypical.sum()
                ),
                "canonical_or_atypical_function_percent": 100 * canonical_or_atypical.mean(),
                "median_max_frequency_signal": subset["max_frequency_signal"].median(),
                "median_popmax_global_ratio": subset["popmax_global_ratio"].median(),
                "median_qualifying_ac": subset["qualifying_frequency_ac"].median(),
                "median_gene_frequency_constraint_proxy": subset[
                    "gene_frequency_constraint_proxy"
                ].median(),
                "band_variants": format_variant_list(subset),
                "interpretation_note": (
                    "Bands represent a pathogenicity tension continuum, not binary classification."
                ),
            }
        )
    return pd.DataFrame(rows)


def summarize_signal_reorganization(scored: pd.DataFrame) -> pd.DataFrame:
    masks = {
        "all_frequency_observed": pd.Series(True, index=scored.index),
        "naive_AF_gt_1e-5": scored["standard_acmg_frequency_flag_bool"],
        "AC_supported_AF_gt_1e-5": scored["frequency_signal_ac_ge_20_bool"],
        "VITAL_60_80_high_tension": scored["vital_score_band_20pt"].eq("60_80_high_tension"),
        "VITAL_80_100_extreme_tension": scored["vital_score_band_20pt"].eq("80_100_extreme_tension"),
        "VITAL_red_actionable": scored["vital_red_flag_bool"],
    }
    rows: list[dict[str, object]] = []
    for signal_layer, mask in masks.items():
        subset = scored.loc[mask].copy()
        if subset.empty:
            rows.append({"signal_layer": signal_layer, "variant_count": 0})
            continue
        weak_or_single = subset["review_strength"].isin(["single_submitter", "weak_or_no_assertion"])
        popmax_ratio_ge_3 = subset["popmax_global_ratio"].fillna(0).ge(3)
        popmax_ratio_ge_10 = subset["popmax_global_ratio"].fillna(0).ge(10)
        popmax_high_frequency = subset["popmax_af"].fillna(0).gt(1e-4)
        indel_or_duplication = subset["variant_type"].isin(["deletion", "insertion", "duplication"])
        canonical_or_atypical = subset["functional_class"].isin(
            ["LOF", "splice_or_intronic", "other_or_unresolved"]
        )
        frameshift = subset["lof_subtype"].eq("frameshift")
        stop_gained = subset["lof_subtype"].eq("stop_gained")
        canonical_splice = subset["lof_subtype"].eq("canonical_splice")
        rows.append(
            {
                "signal_layer": signal_layer,
                "variant_count": len(subset),
                "median_vital_score": subset["vital_score_numeric"].median(),
                "weak_or_single_review_percent": 100 * weak_or_single.mean(),
                "popmax_ratio_ge_3_percent": 100 * popmax_ratio_ge_3.mean(),
                "popmax_ratio_ge_10_percent": 100 * popmax_ratio_ge_10.mean(),
                "popmax_af_gt_1e_4_percent": 100 * popmax_high_frequency.mean(),
                "ac_supported_frequency_percent": 100 * subset["frequency_signal_ac_ge_20_bool"].mean(),
                "indel_or_duplication_percent": 100 * indel_or_duplication.mean(),
                "canonical_or_atypical_function_percent": 100 * canonical_or_atypical.mean(),
                "frameshift_percent": 100 * frameshift.mean(),
                "stop_gained_percent": 100 * stop_gained.mean(),
                "canonical_splice_percent": 100 * canonical_splice.mean(),
                "median_max_frequency_signal": subset["max_frequency_signal"].median(),
                "median_qualifying_ac": subset["qualifying_frequency_ac"].median(),
                "interpretation_note": (
                    "Naive AF, AC support, and score bands reorganize the same frequency signal "
                    "into increasingly review-fragile and function-frequency discordant subsets."
                ),
            }
        )
    return pd.DataFrame(rows)


def summarize_lof_subtypes(scored: pd.DataFrame) -> pd.DataFrame:
    lof = scored[scored["functional_class"].eq("LOF")].copy()
    ordered_subtypes = ["frameshift", "stop_gained", "canonical_splice", "other_LOF"]
    rows: list[dict[str, object]] = []
    for subtype in ordered_subtypes:
        subset = lof[lof["lof_subtype"].eq(subtype)].copy()
        if subset.empty:
            rows.append({"lof_subtype": subtype, "lof_variant_count": 0})
            continue
        naive = subset["standard_acmg_frequency_flag_bool"]
        ac_supported = subset["frequency_signal_ac_ge_20_bool"]
        high_band = subset["vital_score_band_20pt"].isin(["60_80_high_tension", "80_100_extreme_tension"])
        rows.append(
            {
                "lof_subtype": subtype,
                "lof_variant_count": len(subset),
                "naive_af_flag_count": int(naive.sum()),
                "naive_af_flag_percent": 100 * naive.mean(),
                "ac_supported_frequency_count": int(ac_supported.sum()),
                "ac_supported_frequency_percent": 100 * ac_supported.mean(),
                "vital_60_100_count": int(high_band.sum()),
                "vital_60_100_percent": 100 * high_band.mean(),
                "vital_red_count": int(subset["vital_red_flag_bool"].sum()),
                "max_vital_score": subset["vital_score_numeric"].max(),
                "median_vital_score": subset["vital_score_numeric"].median(),
                "median_max_frequency_signal": subset["max_frequency_signal"].median(),
                "max_frequency_signal": subset["max_frequency_signal"].max(),
                "median_qualifying_ac": subset["qualifying_frequency_ac"].median(),
                "max_qualifying_ac": subset["qualifying_frequency_ac"].max(),
                "top_discordant_variants": format_variant_list(
                    subset.sort_values("vital_score_numeric", ascending=False).head(8)
                ),
                "interpretation_note": (
                    "LOF subtypes do not contribute equally; subtype-level tension helps test whether "
                    "supposedly severe variants behave like canonical high-penetrance alleles."
                ),
            }
        )
    return pd.DataFrame(rows)


def run_continuum(score_table: Path, output_prefix: str) -> None:
    scores = pd.read_csv(score_table)
    table = scores[scores["frequency_evidence_status"].eq("frequency_observed")].copy()
    for column in [
        "vital_score",
        "max_frequency_signal",
        "qualifying_frequency_ac",
        "global_ac",
        "popmax_ac",
        "gene_frequency_constraint_proxy",
        "technical_detectability_index",
    ]:
        table[column] = pd.to_numeric(table.get(column), errors="coerce")
    table["vital_score_numeric"] = table["vital_score"].fillna(0)
    table["vital_score_band_20pt"] = table["vital_score_numeric"].map(tension_band)
    table["lof_subtype"] = table.apply(classify_lof_subtype, axis=1)
    table["vital_red_flag_bool"] = bool_series(table.get("vital_red_flag", pd.Series(False, index=table.index)))
    table["standard_acmg_frequency_flag_bool"] = bool_series(
        table.get("standard_acmg_frequency_flag", pd.Series(False, index=table.index))
    )
    table["frequency_signal_ac_ge_20_bool"] = bool_series(
        table.get("frequency_signal_ac_ge_20", pd.Series(False, index=table.index))
    )
    table["pathogenicity_tension_zone"] = table.apply(tension_zone, axis=1)
    table["frequency_function_discordance_class"] = table.apply(functional_discordance_class, axis=1)
    table["continuum_note"] = (
        "The score positions ClinVar P/LP assertions on a pathogenicity tension continuum; "
        "high-scoring variants are candidates for functional overcall, reduced penetrance, "
        "hypomorphic effect, transcript/domain context, or review re-evaluation."
    )

    keep_columns = [
        "gene",
        "clinvar_id",
        "variation_id",
        "title",
        "functional_class",
        "lof_subtype",
        "variant_type",
        "review_strength",
        "global_af",
        "global_ac",
        "popmax_af",
        "popmax_ac",
        "popmax_population",
        "popmax_global_ratio",
        "max_frequency_signal",
        "qualifying_frequency_ac",
        "gene_frequency_constraint_proxy",
        "technical_detectability_index",
        "vital_score",
        "vital_score_band_20pt",
        "vital_band",
        "vital_red_flag",
        "pathogenicity_tension_zone",
        "frequency_function_discordance_class",
        "continuum_note",
    ]
    save_table(table[keep_columns], data_path(output_prefix, "vital_pathogenicity_tension_continuum.csv"))
    save_table(
        summarize_bands(table),
        data_path(output_prefix, "vital_frequency_function_discordance_summary.csv"),
    )
    save_table(
        summarize_signal_reorganization(table),
        data_path(output_prefix, "vital_signal_reorganization_summary.csv"),
    )
    save_table(
        summarize_lof_subtypes(table),
        data_path(output_prefix, "vital_lof_subtype_discordance_summary.csv"),
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize the pathogenicity tension continuum.")
    parser.add_argument("--score-table", type=Path, default=data_path("arrhythmia", "vital_scores.csv"))
    parser.add_argument("--output-prefix", default="arrhythmia")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_continuum(args.score_table, args.output_prefix)


if __name__ == "__main__":
    main()
