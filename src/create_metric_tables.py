from __future__ import annotations

from pathlib import Path

try:
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires pandas, numpy, and matplotlib. Install them with: "
        "python -m pip install -r requirements.txt"
    ) from exc


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
FIGURE_DIR = BASE_DIR / "figures"

CLINVAR_GNOMAD_PATH = DATA_DIR / "clinvar_gnomad_merged.csv"
GNOMAD_SUBSET_PATH = DATA_DIR / "gnomad_subset_with_header.tsv"
GENOMES_OUTLIER_PATH = DATA_DIR / "high_af_variants_genomes.csv"

MAIN_METRICS_PATH = DATA_DIR / "main_metrics.csv"
AF_CATEGORY_COUNTS_PATH = DATA_DIR / "af_category_counts.csv"
CLINVAR_VS_GNOMAD_PATH = DATA_DIR / "clinvar_vs_gnomad_af_summary.csv"
OUTLIERS_PATH = DATA_DIR / "outlier_variants.csv"
GENE_STATS_PATH = DATA_DIR / "gene_af_stats.csv"
OUTLIERS_BY_GENE_PATH = DATA_DIR / "outlier_counts_by_gene.csv"
TOP_VARIANTS_PATH = DATA_DIR / "top_variants.csv"
EXOMES_GENOMES_PATH = DATA_DIR / "exomes_genomes_af_comparison.csv"

AF_DISTRIBUTION_FIGURE = FIGURE_DIR / "af_distribution_updated.png"
AF_CATEGORIES_FIGURE = FIGURE_DIR / "af_categories.png"

ULTRA_RARE_THRESHOLD = 1e-5
COMMON_AF_THRESHOLD = 1e-4
OUTLIER_THRESHOLD = ULTRA_RARE_THRESHOLD
AF_CATEGORY_ORDER = ["ultra_rare", "rare", "common"]


def load_clinvar_gnomad(path: Path = CLINVAR_GNOMAD_PATH) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    df = pd.read_csv(path)
    df["AF"] = pd.to_numeric(df["AF"], errors="coerce")
    df["gnomad_match"] = parse_boolean_series(df["gnomad_match"])
    return df


def parse_boolean_series(series: pd.Series) -> pd.Series:
    normalized = series.fillna("").astype(str).str.strip().str.lower()
    return normalized.isin({"true", "1", "yes"})


def assign_af_categories(af: pd.Series) -> pd.Categorical:
    categories = pd.cut(
        af,
        bins=[0, ULTRA_RARE_THRESHOLD, COMMON_AF_THRESHOLD, np.inf],
        labels=AF_CATEGORY_ORDER,
        include_lowest=True,
    )
    return categories


def explode_numeric_af(series: pd.Series) -> pd.Series:
    af = (
        series.fillna("")
        .astype(str)
        .str.split(",")
        .explode()
        .str.strip()
        .replace("", pd.NA)
    )
    return pd.to_numeric(af, errors="coerce").dropna()


def load_gnomad_af(path: Path = GNOMAD_SUBSET_PATH) -> pd.Series:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    gnomad = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if "AF" not in gnomad.columns:
        raise ValueError(f"Missing AF column in {path}")
    return explode_numeric_af(gnomad["AF"])


def make_main_metrics(df: pd.DataFrame) -> pd.DataFrame:
    af_df = df.dropna(subset=["AF"]).copy()
    metrics = [
        ("total_clinvar_variants", len(df), "count"),
        ("variants_with_af", len(af_df), "count"),
        ("matched_to_gnomad", int(df["gnomad_match"].sum()), "count"),
        (
            "af_gt_1e_5_count",
            int((af_df["AF"] > ULTRA_RARE_THRESHOLD).sum()),
            "count",
        ),
        (
            "af_gt_1e_5_fraction",
            fraction(af_df["AF"] > ULTRA_RARE_THRESHOLD),
            "fraction",
        ),
        (
            "af_gt_1e_4_count",
            int((af_df["AF"] > COMMON_AF_THRESHOLD).sum()),
            "count",
        ),
        (
            "af_gt_1e_4_fraction",
            fraction(af_df["AF"] > COMMON_AF_THRESHOLD),
            "fraction",
        ),
        ("max_af", float(af_df["AF"].max()) if not af_df.empty else np.nan, "af"),
        ("median_af", float(af_df["AF"].median()) if not af_df.empty else np.nan, "af"),
    ]
    return pd.DataFrame(
        [
            {"metric": metric, "value": format_metric_value(value, value_type)}
            for metric, value, value_type in metrics
        ]
    )


def format_metric_value(value: float | int, value_type: str) -> str:
    if pd.isna(value):
        return ""
    if value_type == "count":
        return str(int(value))
    if value_type == "fraction":
        return f"{float(value):.6f}"
    return f"{float(value):.6g}"


def fraction(mask: pd.Series) -> float:
    if len(mask) == 0:
        return 0.0
    return float(mask.sum() / len(mask))


def make_af_category_counts(df: pd.DataFrame) -> pd.DataFrame:
    af_df = df.dropna(subset=["AF"]).copy()
    af_df["AF_category"] = assign_af_categories(af_df["AF"])

    counts = (
        af_df["AF_category"]
        .value_counts(sort=False)
        .reindex(AF_CATEGORY_ORDER, fill_value=0)
        .rename_axis("af_category")
        .reset_index(name="count")
    )
    total = counts["count"].sum()
    counts["fraction"] = counts["count"] / total if total else 0.0
    return counts


def describe_af(label: str, af: pd.Series) -> dict[str, float | str]:
    clean_af = pd.to_numeric(af, errors="coerce").dropna()
    description = clean_af.describe(percentiles=[0.25, 0.5, 0.75])
    row: dict[str, float | str] = {"dataset": label}
    for metric, value in description.items():
        row[str(metric)] = float(value)
    row["af_gt_1e_5_count"] = int((clean_af > ULTRA_RARE_THRESHOLD).sum())
    row["af_gt_1e_4_count"] = int((clean_af > COMMON_AF_THRESHOLD).sum())
    row["af_gt_1e_5_fraction"] = fraction(clean_af > ULTRA_RARE_THRESHOLD)
    row["af_gt_1e_4_fraction"] = fraction(clean_af > COMMON_AF_THRESHOLD)
    return row


def make_clinvar_vs_gnomad_summary(df: pd.DataFrame, gnomad_af: pd.Series) -> pd.DataFrame:
    clinvar_af = df.dropna(subset=["AF"])["AF"]
    rows = [
        describe_af("clinvar_matched_variants", clinvar_af),
        describe_af("gnomad_population_subset", gnomad_af),
    ]
    return pd.DataFrame(rows)


def make_outlier_table(df: pd.DataFrame) -> pd.DataFrame:
    af_df = df.dropna(subset=["AF"]).copy()
    af_df["AF_category"] = assign_af_categories(af_df["AF"])
    return (
        af_df.loc[af_df["AF"] > OUTLIER_THRESHOLD]
        .sort_values("AF", ascending=False)
        .reset_index(drop=True)
    )


def make_gene_stats(df: pd.DataFrame) -> pd.DataFrame:
    af_df = df.dropna(subset=["AF"]).copy()
    stats = af_df.groupby("GeneSymbol")["AF"].describe().reset_index()
    return stats.sort_values(["max", "count"], ascending=[False, False]).reset_index(drop=True)


def make_outlier_counts_by_gene(df: pd.DataFrame, outliers: pd.DataFrame) -> pd.DataFrame:
    gene_totals = (
        df.dropna(subset=["AF"])
        .groupby("GeneSymbol")
        .size()
        .rename("variants_with_af")
        .reset_index()
    )
    outlier_counts = (
        outliers.groupby("GeneSymbol")
        .size()
        .rename("outlier_count_af_gt_1e_5")
        .reset_index()
    )
    counts = gene_totals.merge(outlier_counts, on="GeneSymbol", how="left")
    counts["outlier_count_af_gt_1e_5"] = (
        counts["outlier_count_af_gt_1e_5"].fillna(0).astype(int)
    )
    counts["within_gene_fraction"] = (
        counts["outlier_count_af_gt_1e_5"] / counts["variants_with_af"]
    )
    total_outliers = counts["outlier_count_af_gt_1e_5"].sum()
    counts["fraction_of_all_outliers"] = (
        counts["outlier_count_af_gt_1e_5"] / total_outliers if total_outliers else 0.0
    )
    return counts.sort_values(
        ["outlier_count_af_gt_1e_5", "within_gene_fraction", "GeneSymbol"],
        ascending=[False, False, True],
    ).reset_index(drop=True)


def make_top_variants(df: pd.DataFrame, n: int = 10) -> pd.DataFrame:
    af_df = df.dropna(subset=["AF"]).copy()
    af_df["AF_category"] = assign_af_categories(af_df["AF"])
    preferred_columns = [
        "variant_key",
        "GeneSymbol",
        "AF",
        "AF_category",
        "AC",
        "AN",
        "ClinicalSignificance_values",
        "best_review_status",
        "is_high_confidence",
    ]
    columns = [column for column in preferred_columns if column in af_df.columns]
    return af_df.sort_values("AF", ascending=False).head(n).loc[:, columns]


def make_exomes_genomes_comparison(
    exomes_path: Path = GNOMAD_SUBSET_PATH,
    genomes_path: Path = GENOMES_OUTLIER_PATH,
) -> pd.DataFrame:
    columns = ["variant_key", "AF_exomes", "AF_genomes", "AF_diff"]
    if not exomes_path.exists() or not genomes_path.exists():
        return pd.DataFrame(columns=columns)

    exomes = pd.read_csv(exomes_path, sep="\t", dtype=str, keep_default_na=False)
    genomes = pd.read_csv(genomes_path, dtype=str, keep_default_na=False)
    required = ["CHROM", "POS", "REF", "ALT", "AF"]
    if any(column not in exomes.columns for column in required) or any(
        column not in genomes.columns for column in required
    ):
        return pd.DataFrame(columns=columns)

    for table in (exomes, genomes):
        table["variant_key"] = table[["CHROM", "POS", "REF", "ALT"]].agg(":".join, axis=1)

    comparison = exomes.loc[:, ["variant_key", "AF"]].merge(
        genomes.loc[:, ["variant_key", "AF"]],
        on="variant_key",
        how="inner",
        suffixes=("_exomes", "_genomes"),
    )
    if comparison.empty:
        return pd.DataFrame(columns=columns)

    comparison["AF_exomes"] = pd.to_numeric(comparison["AF_exomes"], errors="coerce")
    comparison["AF_genomes"] = pd.to_numeric(comparison["AF_genomes"], errors="coerce")
    comparison["AF_diff"] = comparison["AF_exomes"] - comparison["AF_genomes"]
    return comparison.loc[:, columns].sort_values("AF_diff", ascending=False)


def create_af_distribution_figure(df: pd.DataFrame, output_path: Path) -> None:
    af = df.dropna(subset=["AF"])["AF"]
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(np.log10(af + 1e-10), bins=20, color="#4c78a8", edgecolor="white")
    ax.axvline(
        np.log10(ULTRA_RARE_THRESHOLD),
        linestyle="--",
        color="#f58518",
        label="1e-5 threshold",
    )
    ax.axvline(
        np.log10(COMMON_AF_THRESHOLD),
        linestyle="--",
        color="#e45756",
        label="1e-4 threshold",
    )
    ax.set_xlabel("log10(AF)")
    ax.set_ylabel("Count")
    ax.set_title("ClinVar-matched variant allele frequencies")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def create_af_category_barplot(counts: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    bars = ax.bar(
        counts["af_category"], counts["count"], color=["#4c78a8", "#f58518", "#e45756"]
    )
    for bar in bars:
        height = bar.get_height()
        ax.annotate(
            f"{int(height)}",
            xy=(bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            va="bottom",
        )
    ax.set_ylabel("Number of variants")
    ax.set_xlabel("AF category")
    ax.set_title("ClinVar-matched variants by AF category")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def save_table(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)


def main() -> None:
    df = load_clinvar_gnomad()
    gnomad_af = load_gnomad_af()

    main_metrics = make_main_metrics(df)
    category_counts = make_af_category_counts(df)
    clinvar_vs_gnomad = make_clinvar_vs_gnomad_summary(df, gnomad_af)
    outliers = make_outlier_table(df)
    gene_stats = make_gene_stats(df)
    outliers_by_gene = make_outlier_counts_by_gene(df, outliers)
    top_variants = make_top_variants(df)
    exomes_genomes = make_exomes_genomes_comparison()

    save_table(main_metrics, MAIN_METRICS_PATH)
    save_table(category_counts, AF_CATEGORY_COUNTS_PATH)
    save_table(clinvar_vs_gnomad, CLINVAR_VS_GNOMAD_PATH)
    save_table(outliers, OUTLIERS_PATH)
    save_table(gene_stats, GENE_STATS_PATH)
    save_table(outliers_by_gene, OUTLIERS_BY_GENE_PATH)
    save_table(top_variants, TOP_VARIANTS_PATH)
    save_table(exomes_genomes, EXOMES_GENOMES_PATH)

    create_af_distribution_figure(df, AF_DISTRIBUTION_FIGURE)
    create_af_category_barplot(category_counts, AF_CATEGORIES_FIGURE)

    print("Main metrics")
    print(main_metrics.to_string(index=False))
    print("\nAF category counts")
    print(category_counts.to_string(index=False))
    print(f"\nTotal outliers above {OUTLIER_THRESHOLD:g}: {len(outliers)}")
    if exomes_genomes.empty:
        print("Exomes/genomes comparison: no overlapping variant keys available.")
    else:
        print("\nExomes/genomes AF differences")
        print(exomes_genomes.head().to_string(index=False))


if __name__ == "__main__":
    main()
