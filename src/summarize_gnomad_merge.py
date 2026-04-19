from __future__ import annotations

from pathlib import Path

try:
    import matplotlib.pyplot as plt
    import pandas as pd
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires pandas and matplotlib. Install them with: "
        "python -m pip install pandas matplotlib"
    ) from exc


BASE_DIR = Path(__file__).resolve().parents[1]
INPUT_PATH = BASE_DIR / "data" / "processed" / "clinvar_gnomad_merged.csv"
MATCH_SUMMARY_PATH = BASE_DIR / "data" / "processed" / "gnomad_match_summary.csv"
AF_BINS_OVERALL_PATH = BASE_DIR / "data" / "processed" / "gnomad_af_bins_overall.csv"
AF_BINS_REVIEW_PATH = (
    BASE_DIR / "data" / "processed" / "gnomad_af_bins_by_review_status.csv"
)
AF_BINS_CONFIDENCE_PATH = (
    BASE_DIR / "data" / "processed" / "gnomad_af_bins_by_confidence.csv"
)
FIGURE_PATH = BASE_DIR / "figures" / "gnomad_af_distribution.png"

REQUIRED_COLUMNS = (
    "variant_key",
    "best_review_status",
    "is_high_confidence",
    "AF",
    "gnomad_match",
)

AF_BIN_ORDER = [
    "AF == 0 or missing",
    "0 < AF < 1e-5",
    "1e-5 <= AF < 1e-4",
    "1e-4 <= AF < 1e-3",
    "AF >= 1e-3",
]


def load_merged_variants(input_path: Path) -> pd.DataFrame:
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    return pd.read_csv(input_path, dtype=str, keep_default_na=False)


def validate_required_columns(df: pd.DataFrame) -> None:
    missing = [column for column in REQUIRED_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")


def parse_boolean_series(series: pd.Series) -> pd.Series:
    normalized = series.fillna("").astype(str).str.strip().str.lower()
    return normalized.isin({"true", "1", "yes"})


def prepare_merged_variants(df: pd.DataFrame) -> pd.DataFrame:
    prepared = df.copy()
    prepared["gnomad_match"] = parse_boolean_series(prepared["gnomad_match"])
    prepared["is_high_confidence"] = parse_boolean_series(prepared["is_high_confidence"])
    prepared["AF"] = pd.to_numeric(prepared["AF"], errors="coerce")
    prepared["best_review_status"] = (
        prepared["best_review_status"].fillna("").astype(str).str.strip().replace("", "Missing")
    )
    prepared["af_bin"] = assign_af_bin(prepared["AF"])
    return prepared


def assign_af_bin(af: pd.Series) -> pd.Series:
    bins = pd.Series("AF == 0 or missing", index=af.index, dtype="object")
    bins = bins.mask((af > 0) & (af < 1e-5), "0 < AF < 1e-5")
    bins = bins.mask((af >= 1e-5) & (af < 1e-4), "1e-5 <= AF < 1e-4")
    bins = bins.mask((af >= 1e-4) & (af < 1e-3), "1e-4 <= AF < 1e-3")
    bins = bins.mask(af >= 1e-3, "AF >= 1e-3")
    return pd.Categorical(bins, categories=AF_BIN_ORDER, ordered=True)


def summarize_match_status(df: pd.DataFrame) -> pd.DataFrame:
    total_variants = len(df)
    matched_mask = df["gnomad_match"]
    matched_count = int(matched_mask.sum())
    unmatched_count = total_variants - matched_count
    matched_missing_af = int((matched_mask & df["AF"].isna()).sum())

    summary = pd.DataFrame(
        [
            {
                "metric": "total_clinvar_variants",
                "count": total_variants,
                "proportion": 1.0 if total_variants else 0.0,
            },
            {
                "metric": "matched_to_gnomad",
                "count": matched_count,
                "proportion": matched_count / total_variants if total_variants else 0.0,
            },
            {
                "metric": "unmatched_to_gnomad",
                "count": unmatched_count,
                "proportion": unmatched_count / total_variants if total_variants else 0.0,
            },
            {
                "metric": "matched_with_missing_af",
                "count": matched_missing_af,
                "proportion": matched_missing_af / matched_count if matched_count else 0.0,
            },
        ]
    )
    return summary


def make_af_bin_table(df: pd.DataFrame, group_column: str | None = None) -> pd.DataFrame:
    matched = df.loc[df["gnomad_match"]].copy()

    if group_column is None:
        counts = (
            matched["af_bin"]
            .value_counts(sort=False, dropna=False)
            .reindex(AF_BIN_ORDER, fill_value=0)
            .rename_axis("af_bin")
            .reset_index(name="count")
        )
        total = len(matched)
        counts["proportion"] = counts["count"] / total if total else 0.0
        return counts

    group_values = df[group_column].drop_duplicates()
    if group_column == "is_high_confidence":
        group_values = group_values.map({True: "True", False: "False"})
        matched = matched.copy()
        matched[group_column] = matched[group_column].map({True: "True", False: "False"})

    full_index = pd.MultiIndex.from_product(
        [sorted(group_values.tolist()), AF_BIN_ORDER],
        names=[group_column, "af_bin"],
    )

    if matched.empty:
        counts = full_index.to_frame(index=False)
        counts["count"] = 0
        counts["proportion"] = 0.0
        return counts

    counts = (
        matched.groupby([group_column, "af_bin"], dropna=False, observed=False)
        .size()
        .reindex(full_index, fill_value=0)
        .reset_index(name="count")
    )

    totals = counts.groupby(group_column)["count"].transform("sum")
    counts["proportion"] = counts["count"] / totals.where(totals > 0, other=1)
    counts.loc[totals == 0, "proportion"] = 0.0
    return counts.reset_index(drop=True)


def save_table(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)


def create_af_distribution_figure(overall_af_bins: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(overall_af_bins["af_bin"], overall_af_bins["count"], color="#2f6db2")
    ax.set_title("gnomAD AF Distribution for Matched ClinVar Variants")
    ax.set_xlabel("AF bin")
    ax.set_ylabel("Variant count")
    ax.tick_params(axis="x", rotation=25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle="--", alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def print_summary(summary_df: pd.DataFrame) -> None:
    summary_lookup = summary_df.set_index("metric")
    total = int(summary_lookup.loc["total_clinvar_variants", "count"])
    matched = int(summary_lookup.loc["matched_to_gnomad", "count"])
    unmatched = int(summary_lookup.loc["unmatched_to_gnomad", "count"])
    matched_prop = float(summary_lookup.loc["matched_to_gnomad", "proportion"])
    unmatched_prop = float(summary_lookup.loc["unmatched_to_gnomad", "proportion"])
    missing_af = int(summary_lookup.loc["matched_with_missing_af", "count"])

    print(f"Total ClinVar variants: {total}")
    print(f"Matched to gnomAD: {matched} ({matched_prop:.1%})")
    print(f"Unmatched to gnomAD: {unmatched} ({unmatched_prop:.1%})")
    print(f"Matched variants with missing AF: {missing_af}")


def main() -> None:
    merged = load_merged_variants(INPUT_PATH)
    validate_required_columns(merged)
    merged = prepare_merged_variants(merged)

    match_summary = summarize_match_status(merged)
    af_bins_overall = make_af_bin_table(merged)
    af_bins_by_review = make_af_bin_table(merged, group_column="best_review_status")
    af_bins_by_confidence = make_af_bin_table(merged, group_column="is_high_confidence")

    save_table(match_summary, MATCH_SUMMARY_PATH)
    save_table(af_bins_overall, AF_BINS_OVERALL_PATH)
    save_table(af_bins_by_review, AF_BINS_REVIEW_PATH)
    save_table(af_bins_by_confidence, AF_BINS_CONFIDENCE_PATH)
    create_af_distribution_figure(af_bins_overall, FIGURE_PATH)

    print_summary(match_summary)


if __name__ == "__main__":
    main()
