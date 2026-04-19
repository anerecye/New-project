from __future__ import annotations

from pathlib import Path

try:
    import numpy as np
    import pandas as pd
    from scipy.stats import fisher_exact
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires pandas, numpy, and scipy. Install them with: "
        "python -m pip install -r requirements.txt"
    ) from exc


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"

CLINVAR_GNOMAD_PATH = DATA_DIR / "clinvar_gnomad_merged.csv"
GNOMAD_SUBSET_PATH = DATA_DIR / "gnomad_subset_with_header.tsv"

MATCH_QC_SUMMARY_PATH = DATA_DIR / "match_qc_summary.csv"
MATCH_QC_BY_GENE_PATH = DATA_DIR / "match_qc_by_gene.csv"
GENE_OUTLIER_TESTS_PATH = DATA_DIR / "gene_outlier_enrichment_tests.csv"

OUTLIER_THRESHOLD = 1e-5


def parse_boolean_series(series: pd.Series) -> pd.Series:
    normalized = series.fillna("").astype(str).str.strip().str.lower()
    return normalized.isin({"true", "1", "yes"})


def standardize_chromosome(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.replace(r"^chr", "", regex=True, case=False)


def build_variant_key(
    chromosome: pd.Series,
    position: pd.Series,
    reference: pd.Series,
    alternate: pd.Series,
) -> pd.Series:
    key_df = pd.DataFrame(
        {
            "chromosome": standardize_chromosome(chromosome),
            "position": position.astype(str).str.strip(),
            "reference": reference.astype(str).str.strip(),
            "alternate": alternate.astype(str).str.strip(),
        }
    )
    return key_df.agg(":".join, axis=1)


def build_position_key(chromosome: pd.Series, position: pd.Series) -> pd.Series:
    key_df = pd.DataFrame(
        {
            "chromosome": standardize_chromosome(chromosome),
            "position": position.astype(str).str.strip(),
        }
    )
    return key_df.agg(":".join, axis=1)


def load_clinvar_gnomad() -> pd.DataFrame:
    df = pd.read_csv(CLINVAR_GNOMAD_PATH, dtype=str, keep_default_na=False)
    df["AF"] = pd.to_numeric(df["AF"].replace("", pd.NA), errors="coerce")
    df["gnomad_match"] = parse_boolean_series(df["gnomad_match"])
    df["_variant_key"] = build_variant_key(
        df["Chromosome"],
        df["Start"],
        df["ReferenceAllele"],
        df["AlternateAllele"],
    )
    df["_position_key"] = build_position_key(df["Chromosome"], df["Start"])
    return df


def load_gnomad_subset() -> pd.DataFrame:
    df = pd.read_csv(GNOMAD_SUBSET_PATH, sep="\t", dtype=str, keep_default_na=False)
    df["_variant_key"] = build_variant_key(df["CHROM"], df["POS"], df["REF"], df["ALT"])
    df["_position_key"] = build_position_key(df["CHROM"], df["POS"])
    return df


def classify_match_status(clinvar_df: pd.DataFrame, gnomad_df: pd.DataFrame) -> pd.Series:
    gnomad_positions = set(gnomad_df["_position_key"])

    status = pd.Series("no_gnomad_record_at_queried_position", index=clinvar_df.index)
    status = status.mask(
        clinvar_df["_position_key"].isin(gnomad_positions),
        "position_present_no_exact_allele_match",
    )
    status = status.mask(clinvar_df["gnomad_match"], "exact_variant_match")
    return status


def make_match_qc_tables(
    clinvar_df: pd.DataFrame, gnomad_df: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    qc = clinvar_df.copy()
    qc["match_status"] = classify_match_status(qc, gnomad_df)

    summary = (
        qc["match_status"]
        .value_counts()
        .rename_axis("match_status")
        .reset_index(name="count")
    )
    summary["fraction_of_clinvar_variants"] = summary["count"] / len(qc)

    by_gene = (
        qc.groupby(["GeneSymbol", "match_status"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    expected_columns = [
        "exact_variant_match",
        "position_present_no_exact_allele_match",
        "no_gnomad_record_at_queried_position",
    ]
    for column in expected_columns:
        if column not in by_gene.columns:
            by_gene[column] = 0

    by_gene["total_clinvar_variants"] = by_gene[expected_columns].sum(axis=1)
    by_gene["exact_match_fraction"] = (
        by_gene["exact_variant_match"] / by_gene["total_clinvar_variants"]
    )
    by_gene = by_gene.loc[
        :,
        [
            "GeneSymbol",
            "total_clinvar_variants",
            "exact_variant_match",
            "position_present_no_exact_allele_match",
            "no_gnomad_record_at_queried_position",
            "exact_match_fraction",
        ],
    ].sort_values("total_clinvar_variants", ascending=False)

    return summary, by_gene.reset_index(drop=True)


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    p = p_values.to_numpy(dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.empty(n, dtype=float)
    cumulative = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        cumulative = min(cumulative, ranked[i] * n / rank)
        adjusted[order[i]] = cumulative
    return pd.Series(adjusted, index=p_values.index)


def make_gene_outlier_tests(clinvar_df: pd.DataFrame) -> pd.DataFrame:
    af_df = clinvar_df.dropna(subset=["AF"]).copy()
    af_df["is_outlier_af_gt_1e_5"] = af_df["AF"] > OUTLIER_THRESHOLD
    total_variants = len(af_df)
    total_outliers = int(af_df["is_outlier_af_gt_1e_5"].sum())
    total_non_outliers = total_variants - total_outliers
    global_fraction = total_outliers / total_variants if total_variants else 0.0

    rows: list[dict[str, float | int | str]] = []
    for gene, group in af_df.groupby("GeneSymbol", sort=True):
        gene_total = len(group)
        gene_outliers = int(group["is_outlier_af_gt_1e_5"].sum())
        gene_non_outliers = gene_total - gene_outliers
        rest_outliers = total_outliers - gene_outliers
        rest_non_outliers = total_non_outliers - gene_non_outliers
        contingency = [[gene_outliers, gene_non_outliers], [rest_outliers, rest_non_outliers]]

        odds_ratio_greater, p_greater = fisher_exact(contingency, alternative="greater")
        odds_ratio_less, p_less = fisher_exact(contingency, alternative="less")
        odds_ratio_two_sided, p_two_sided = fisher_exact(contingency, alternative="two-sided")

        rows.append(
            {
                "GeneSymbol": gene,
                "variants_with_af": gene_total,
                "outlier_count_af_gt_1e_5": gene_outliers,
                "non_outlier_count": gene_non_outliers,
                "within_gene_fraction": gene_outliers / gene_total if gene_total else 0.0,
                "global_outlier_fraction": global_fraction,
                "expected_outliers_under_global_rate": gene_total * global_fraction,
                "odds_ratio_greater": odds_ratio_greater,
                "p_value_greater": p_greater,
                "odds_ratio_less": odds_ratio_less,
                "p_value_less": p_less,
                "odds_ratio_two_sided": odds_ratio_two_sided,
                "p_value_two_sided": p_two_sided,
            }
        )

    tests = pd.DataFrame(rows)
    tests["q_value_greater_bh"] = benjamini_hochberg(tests["p_value_greater"])
    tests["q_value_two_sided_bh"] = benjamini_hochberg(tests["p_value_two_sided"])
    return tests.sort_values(
        ["within_gene_fraction", "outlier_count_af_gt_1e_5"],
        ascending=[False, False],
    ).reset_index(drop=True)


def save_table(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)


def main() -> None:
    clinvar_df = load_clinvar_gnomad()
    gnomad_df = load_gnomad_subset()

    match_summary, match_by_gene = make_match_qc_tables(clinvar_df, gnomad_df)
    gene_tests = make_gene_outlier_tests(clinvar_df)

    save_table(match_summary, MATCH_QC_SUMMARY_PATH)
    save_table(match_by_gene, MATCH_QC_BY_GENE_PATH)
    save_table(gene_tests, GENE_OUTLIER_TESTS_PATH)

    print("Match QC summary")
    print(match_summary.to_string(index=False))
    print("\nGene outlier enrichment tests")
    print(gene_tests.to_string(index=False))


if __name__ == "__main__":
    main()
