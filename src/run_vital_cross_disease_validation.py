from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires pandas, numpy, and matplotlib. Install them with: "
        "python -m pip install -r requirements.txt"
    ) from exc

sys.path.insert(0, str(Path(__file__).resolve().parent))

import advanced_variant_analyses as ava
from run_vital_external_panels import PANEL_DEFINITIONS
from score_vital_from_variant_summary import load_plp_snapshot, score_preloaded_snapshot
from validate_vital_reclassification import ARRHYTHMIA_GENES, evaluate_reclassification


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
FIGURE_DIR = BASE_DIR / "figures"

KNOWN_FOUNDER_OR_CARRIER_CONTEXT_GENES = {
    "APC",
    "BRCA1",
    "BRCA2",
    "CFTR",
    "GBA",
    "HBB",
    "HEXA",
    "HFE",
    "MLH1",
    "MSH2",
    "MSH6",
    "PAH",
    "PMS2",
    "SERPINA1",
    "TPP1",
}

PREVIOUS_EXTERNAL_PANEL_GENES = {
    gene for genes in PANEL_DEFINITIONS.values() for gene in genes
}


def data_path(prefix: str, name: str) -> Path:
    return DATA_DIR / ava.prefixed_name(prefix, name)


def save_table(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"Saved {output_path} ({len(df)} rows)")


def as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def wilson_ci(successes: int, total: int, z: float = 1.96) -> tuple[float, float]:
    if total <= 0:
        return np.nan, np.nan
    proportion = successes / total
    denominator = 1 + (z**2 / total)
    center = (proportion + (z**2 / (2 * total))) / denominator
    margin = (
        z
        * math.sqrt((proportion * (1 - proportion) / total) + (z**2 / (4 * total**2)))
        / denominator
    )
    return max(0.0, center - margin), min(1.0, center + margin)


def rate_summary(label: str, successes: int, total: int) -> dict[str, object]:
    low, high = wilson_ci(successes, total)
    rate = successes / total if total else np.nan
    return {
        "metric": label,
        "count": successes,
        "denominator": total,
        "rate": rate,
        "rate_ci_low": low,
        "rate_ci_high": high,
        "expected_per_10000": rate * 10_000 if not np.isnan(rate) else np.nan,
        "expected_per_10000_ci_low": low * 10_000 if not np.isnan(low) else np.nan,
        "expected_per_10000_ci_high": high * 10_000 if not np.isnan(high) else np.nan,
    }


def variant_id(df: pd.DataFrame) -> pd.Series:
    return df["variation_id"].fillna("").astype(str) + "|" + df["variant_key"].fillna("").astype(str)


def build_cross_disease_sample(
    variant_summary: Path,
    output_prefix: str,
    top_gene_count: int,
    top_gene_sample_size: int,
    random_sample_size: int,
    random_seed: int,
    excluded_genes: set[str],
) -> pd.DataFrame:
    snapshot = load_plp_snapshot(
        variant_summary,
        genes=None,
        max_variants=None,
        sample_variants=None,
        random_seed=random_seed,
    )
    if snapshot.empty:
        raise ValueError(f"No P/LP GRCh38 variants found in {variant_summary}.")

    snapshot["gene"] = snapshot["gene"].fillna("").astype(str)
    snapshot = snapshot.loc[~snapshot["gene"].isin(excluded_genes)].copy()
    snapshot = snapshot.loc[snapshot["gene"].ne("")].copy()
    snapshot["sample_variant_id"] = variant_id(snapshot)

    gene_counts = (
        snapshot.groupby("gene", dropna=False)
        .size()
        .rename("plp_variant_count")
        .reset_index()
        .sort_values(["plp_variant_count", "gene"], ascending=[False, True])
        .reset_index(drop=True)
    )
    gene_counts["gene_plp_count_rank"] = np.arange(1, len(gene_counts) + 1)
    top_genes = set(gene_counts.head(top_gene_count)["gene"])
    save_table(gene_counts, data_path(output_prefix, "gene_plp_counts.csv"))

    rng_seed = random_seed
    top_pool = snapshot.loc[snapshot["gene"].isin(top_genes)].copy()
    top_n = min(top_gene_sample_size, len(top_pool))
    top_sample = top_pool.sample(n=top_n, random_state=rng_seed).copy() if top_n else top_pool.head(0).copy()
    top_sample["validation_stratum"] = "top100_gene_enriched"

    selected_ids = set(top_sample["sample_variant_id"])
    random_pool = snapshot.loc[~snapshot["sample_variant_id"].isin(selected_ids)].copy()
    random_n = min(random_sample_size, len(random_pool))
    random_sample = (
        random_pool.sample(n=random_n, random_state=rng_seed + 1).copy()
        if random_n
        else random_pool.head(0).copy()
    )
    random_sample["validation_stratum"] = "fully_random_remaining"

    selected = pd.concat([top_sample, random_sample], ignore_index=True)
    selected = selected.drop_duplicates(subset="sample_variant_id", keep="first").reset_index(drop=True)
    selected = selected.merge(gene_counts, on="gene", how="left")
    selected["top_gene_membership"] = selected["gene"].isin(top_genes)
    selected["known_founder_or_carrier_context_gene"] = selected["gene"].isin(
        KNOWN_FOUNDER_OR_CARRIER_CONTEXT_GENES
    )
    selected["sample_design_note"] = (
        f"{top_n} variants sampled from the top {top_gene_count} genes by number of "
        f"ClinVar P/LP submissions after exclusions; "
        f"{random_n} variants sampled from the remaining non-arrhythmia, non-external-panel pool"
    )
    save_table(selected, data_path(output_prefix, "sample_design.csv"))

    arrhythmia_leaks = selected["gene"].isin(ARRHYTHMIA_GENES)
    previous_panel_leaks = selected["gene"].isin(PREVIOUS_EXTERNAL_PANEL_GENES)
    stratum_summary = (
        selected.groupby("validation_stratum", dropna=False)
        .agg(
            variant_count=("sample_variant_id", "size"),
            unique_gene_count=("gene", "nunique"),
            top100_gene_membership_count=("top_gene_membership", "sum"),
            known_founder_or_carrier_context_count=("known_founder_or_carrier_context_gene", "sum"),
        )
        .reset_index()
    )
    save_table(stratum_summary, data_path(output_prefix, "sample_design_summary.csv"))
    sanity = pd.DataFrame(
        [
            {
                "check": "arrhythmia_training_gene_leakage",
                "variant_count": int(arrhythmia_leaks.sum()),
                "unique_gene_count": int(selected.loc[arrhythmia_leaks, "gene"].nunique()),
                "genes": "|".join(sorted(selected.loc[arrhythmia_leaks, "gene"].unique())),
            },
            {
                "check": "previous_external_panel_gene_overlap",
                "variant_count": int(previous_panel_leaks.sum()),
                "unique_gene_count": int(selected.loc[previous_panel_leaks, "gene"].nunique()),
                "genes": "|".join(sorted(selected.loc[previous_panel_leaks, "gene"].unique())),
            },
            {
                "check": "duplicate_variant_id_after_sampling",
                "variant_count": int(selected.duplicated(subset=["variation_id", "variant_key"]).sum()),
                "unique_gene_count": np.nan,
                "genes": "",
            },
        ]
    )
    save_table(sanity, data_path(output_prefix, "sample_sanity_checks.csv"))
    return selected.drop(columns=["sample_variant_id"])


def attach_sample_metadata(scores: pd.DataFrame, sample: pd.DataFrame) -> pd.DataFrame:
    metadata_columns = [
        "variation_id",
        "variant_key",
        "validation_stratum",
        "gene_plp_count_rank",
        "plp_variant_count",
        "top_gene_membership",
        "known_founder_or_carrier_context_gene",
    ]
    metadata = sample.loc[:, [column for column in metadata_columns if column in sample.columns]].copy()
    metadata["variation_id"] = metadata["variation_id"].astype(str)
    scores = scores.copy()
    scores["variation_id"] = scores["variation_id"].astype(str)
    return scores.merge(metadata, on=["variation_id", "variant_key"], how="left")


def summarize_scores(scores: pd.DataFrame, prefix: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    table = scores.copy()
    for column in ["standard_acmg_frequency_flag", "frequency_signal_ac_ge_20", "vital_red_flag"]:
        if column not in table.columns:
            table[column] = False
        table[column] = as_bool(table[column])
    table["vital_score"] = pd.to_numeric(table.get("vital_score"), errors="coerce")
    table["frequency_evidence_status"] = table.get("frequency_evidence_status", "").fillna("").astype(str)
    observed = table["frequency_evidence_status"].eq("frequency_observed")

    rows: list[dict[str, object]] = []
    scopes = {"overall": pd.Series(True, index=table.index)}
    if "validation_stratum" in table.columns:
        for stratum in sorted(table["validation_stratum"].dropna().astype(str).unique()):
            scopes[f"stratum:{stratum}"] = table["validation_stratum"].astype(str).eq(stratum)
    if "known_founder_or_carrier_context_gene" in table.columns:
        scopes["context:known_founder_or_carrier"] = table["known_founder_or_carrier_context_gene"].fillna(False).astype(bool)
        scopes["context:other_genes"] = ~table["known_founder_or_carrier_context_gene"].fillna(False).astype(bool)

    for scope, mask in scopes.items():
        sub = table.loc[mask].copy()
        n = len(sub)
        if n == 0:
            continue
        observed_sub = sub["frequency_evidence_status"].eq("frequency_observed")
        naive = int(sub["standard_acmg_frequency_flag"].sum())
        ac_supported = int(sub["frequency_signal_ac_ge_20"].sum())
        red = int(sub["vital_red_flag"].sum())
        rows.extend(
            [
                {"scope": scope, **rate_summary("frequency_observed", int(observed_sub.sum()), n)},
                {"scope": scope, **rate_summary("naive_af_gt_1e_5", naive, n)},
                {"scope": scope, **rate_summary("ac_supported_frequency_signal", ac_supported, n)},
                {"scope": scope, **rate_summary("vital_red", red, n)},
            ]
        )
        compression = naive / red if red else np.inf if naive else np.nan
        rows.append(
            {
                "scope": scope,
                "metric": "naive_to_vital_red_compression_ratio",
                "count": compression,
                "denominator": np.nan,
                "rate": np.nan,
                "rate_ci_low": np.nan,
                "rate_ci_high": np.nan,
                "expected_per_10000": np.nan,
                "expected_per_10000_ci_low": np.nan,
                "expected_per_10000_ci_high": np.nan,
            }
        )
        rows.append(
            {
                "scope": scope,
                "metric": "median_vital_score_frequency_observed",
                "count": float(sub.loc[observed_sub, "vital_score"].median()) if observed_sub.any() else np.nan,
                "denominator": np.nan,
                "rate": np.nan,
                "rate_ci_low": np.nan,
                "rate_ci_high": np.nan,
                "expected_per_10000": np.nan,
                "expected_per_10000_ci_low": np.nan,
                "expected_per_10000_ci_high": np.nan,
            }
        )
        rows.append(
            {
                "scope": scope,
                "metric": "p95_vital_score_frequency_observed",
                "count": float(sub.loc[observed_sub, "vital_score"].quantile(0.95)) if observed_sub.any() else np.nan,
                "denominator": np.nan,
                "rate": np.nan,
                "rate_ci_low": np.nan,
                "rate_ci_high": np.nan,
                "expected_per_10000": np.nan,
                "expected_per_10000_ci_low": np.nan,
                "expected_per_10000_ci_high": np.nan,
            }
        )

    summary = pd.DataFrame(rows)

    band_distribution = (
        table.groupby(["validation_stratum", "vital_band"], dropna=False)
        .size()
        .rename("variant_count")
        .reset_index()
    )
    band_distribution["variant_fraction_within_stratum"] = band_distribution.groupby(
        "validation_stratum"
    )["variant_count"].transform(lambda values: values / values.sum())

    outlier_columns = [
        "gene",
        "clinvar_id",
        "variation_id",
        "title",
        "clinsig",
        "review_status",
        "submitter_count",
        "validation_stratum",
        "known_founder_or_carrier_context_gene",
        "frequency_evidence_status",
        "variant_type",
        "functional_class",
        "global_af",
        "global_ac",
        "popmax_af",
        "popmax_ac",
        "popmax_population",
        "max_frequency_signal",
        "qualifying_frequency_ac",
        "standard_acmg_frequency_flag",
        "frequency_signal_ac_ge_20",
        "weak_review_signal",
        "vital_score",
        "vital_band",
        "vital_red_flag",
        "vital_signal_reason",
    ]
    present_outlier_columns = [column for column in outlier_columns if column in table.columns]
    outliers = table.loc[
        table["vital_red_flag"] | (table["vital_score"].fillna(0) >= 60),
        present_outlier_columns,
    ].sort_values(["vital_red_flag", "vital_score", "max_frequency_signal"], ascending=[False, False, False])
    if outliers.empty:
        outliers = table.loc[:, present_outlier_columns].sort_values(
            ["vital_score", "max_frequency_signal"], ascending=[False, False]
        ).head(20)

    save_table(summary, DATA_DIR / f"{prefix}_cross_disease_validation_summary.csv")
    save_table(band_distribution, DATA_DIR / f"{prefix}_cross_disease_band_distribution.csv")
    save_table(outliers, DATA_DIR / f"{prefix}_cross_disease_outliers.csv")
    return summary, band_distribution, outliers


def plot_score_distribution(scores: pd.DataFrame, output_path: Path) -> None:
    plot_df = scores.copy()
    plot_df["vital_score"] = pd.to_numeric(plot_df.get("vital_score"), errors="coerce")
    plot_df = plot_df.loc[
        plot_df["frequency_evidence_status"].eq("frequency_observed")
        & plot_df["vital_score"].notna()
    ].copy()
    if plot_df.empty:
        return

    labels = ["top100_gene_enriched", "fully_random_remaining"]
    labels = [label for label in labels if label in set(plot_df["validation_stratum"].dropna().astype(str))]
    if not labels:
        labels = sorted(plot_df["validation_stratum"].dropna().astype(str).unique())

    values = [plot_df.loc[plot_df["validation_stratum"].eq(label), "vital_score"].to_numpy() for label in labels]
    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    box = ax.boxplot(
        values,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#111827", "linewidth": 1.4},
        boxprops={"linewidth": 1.0},
        whiskerprops={"linewidth": 1.0},
        capprops={"linewidth": 1.0},
    )
    colors = ["#2f6f9f", "#7a7f38"]
    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.36)
        patch.set_edgecolor("#1f2937")

    rng = np.random.default_rng(20260421)
    for index, value_set in enumerate(values, start=1):
        if len(value_set) == 0:
            continue
        take = min(len(value_set), 600)
        sample = rng.choice(value_set, size=take, replace=False) if len(value_set) > take else value_set
        jitter = rng.normal(0, 0.035, len(sample))
        ax.scatter(
            np.full(len(sample), index) + jitter,
            sample,
            s=14,
            alpha=0.35,
            color=colors[(index - 1) % len(colors)],
            edgecolors="none",
        )

    ax.axhline(70, color="#b3261e", linestyle="--", linewidth=1.25, label="VITAL red threshold")
    ax.axhline(60, color="#d18f2f", linestyle=":", linewidth=1.05, label="watchlist threshold")
    ax.set_ylabel("VITAL score among frequency-observed variants")
    ax.set_ylim(-3, 105)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_title("Cross-disease ClinVar P/LP portability stress test")
    ax.grid(axis="y", color="#d8dde6", linewidth=0.8, alpha=0.72)
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)
    print(f"Saved {output_path}")


def run_current_validation(args: argparse.Namespace) -> None:
    excluded_genes = set(ARRHYTHMIA_GENES) | set(PREVIOUS_EXTERNAL_PANEL_GENES)
    sample = build_cross_disease_sample(
        variant_summary=args.variant_summary,
        output_prefix=args.output_prefix,
        top_gene_count=args.top_gene_count,
        top_gene_sample_size=args.top_gene_sample_size,
        random_sample_size=args.random_sample_size,
        random_seed=args.random_seed,
        excluded_genes=excluded_genes,
    )
    if not args.skip_scoring:
        score_preloaded_snapshot(
            sample,
            output_prefix=args.output_prefix,
            dataset=args.dataset,
            pause=args.gnomad_pause,
            fetch_gnomad=not args.no_fetch_gnomad,
            force_gnomad=args.force_gnomad_fetch,
        )

    scores_path = data_path(args.output_prefix, "vital_scores.csv")
    scores = pd.read_csv(scores_path)
    scores = attach_sample_metadata(scores, sample)
    save_table(scores, scores_path)
    summarize_scores(scores, args.output_prefix)
    plot_score_distribution(scores, FIGURE_DIR / f"{args.output_prefix}_score_distribution.png")


def run_historical_validation(args: argparse.Namespace) -> None:
    if args.historical_variant_summary is None or args.followup_variant_summary is None:
        return
    historical_prefix = args.historical_output_prefix
    excluded_genes = set(ARRHYTHMIA_GENES) | set(PREVIOUS_EXTERNAL_PANEL_GENES)
    sample = build_cross_disease_sample(
        variant_summary=args.historical_variant_summary,
        output_prefix=historical_prefix,
        top_gene_count=args.top_gene_count,
        top_gene_sample_size=args.top_gene_sample_size,
        random_sample_size=args.random_sample_size,
        random_seed=args.random_seed,
        excluded_genes=excluded_genes,
    )
    if not args.skip_scoring:
        score_preloaded_snapshot(
            sample,
            output_prefix=historical_prefix,
            dataset=args.dataset,
            pause=args.gnomad_pause,
            fetch_gnomad=not args.no_fetch_gnomad,
            force_gnomad=args.force_gnomad_fetch,
        )

    scores_path = data_path(historical_prefix, "vital_scores.csv")
    scores = pd.read_csv(scores_path)
    scores = attach_sample_metadata(scores, sample)
    save_table(scores, scores_path)
    summarize_scores(scores, historical_prefix)
    plot_score_distribution(scores, FIGURE_DIR / f"{historical_prefix}_score_distribution.png")

    sampled_genes = sorted(set(sample["gene"].dropna().astype(str)) - set(ARRHYTHMIA_GENES))
    evaluate_reclassification(
        score_table=scores_path,
        baseline_summary=args.historical_variant_summary,
        followup_summary=args.followup_variant_summary,
        output_prefix=f"{historical_prefix}_to_current",
        genes=sampled_genes,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run a large independent cross-disease VITAL portability validation "
            "using random ClinVar P/LP variants outside the arrhythmia training panel."
        )
    )
    parser.add_argument("--variant-summary", type=Path, default=BASE_DIR / "data" / "variant_summary.txt.gz")
    parser.add_argument("--output-prefix", default="vital_cross_disease_3000")
    parser.add_argument("--top-gene-count", type=int, default=100)
    parser.add_argument("--top-gene-sample-size", type=int, default=2000)
    parser.add_argument("--random-sample-size", type=int, default=1000)
    parser.add_argument("--random-seed", type=int, default=20260421)
    parser.add_argument("--dataset", default=ava.GNOMAD_DATASET)
    parser.add_argument("--gnomad-pause", type=float, default=0.1)
    parser.add_argument("--no-fetch-gnomad", action="store_true")
    parser.add_argument("--force-gnomad-fetch", action="store_true")
    parser.add_argument("--skip-scoring", action="store_true")
    parser.add_argument("--historical-variant-summary", type=Path)
    parser.add_argument("--followup-variant-summary", type=Path)
    parser.add_argument("--historical-output-prefix", default="vital_cross_disease_3000_2023_01")
    parser.add_argument(
        "--historical-only",
        action="store_true",
        help="Run only the historical baseline sample and follow-up validation.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.historical_only:
        run_current_validation(args)
    run_historical_validation(args)


if __name__ == "__main__":
    main()
