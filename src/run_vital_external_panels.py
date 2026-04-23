from __future__ import annotations

import argparse
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

import advanced_variant_analyses as ava
from score_vital_from_variant_summary import score_snapshot


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
FIGURE_DIR = BASE_DIR / "figures"

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


def data_path(prefix: str, name: str) -> Path:
    return DATA_DIR / ava.prefixed_name(prefix, name)


def save_table(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"Saved {output_path} ({len(df)} rows)")


def read_metric(summary: pd.DataFrame, metric: str) -> float:
    row = summary.loc[summary["metric"].eq(metric), "value"]
    return float(row.iloc[0]) if not row.empty else np.nan


def summarize_prefix(prefix: str, domain: str, domain_type: str) -> tuple[dict[str, object], pd.DataFrame]:
    scores_path = data_path(prefix, "vital_scores.csv")
    summary_path = data_path(prefix, "vital_summary.csv")
    if not scores_path.exists() or not summary_path.exists():
        raise FileNotFoundError(f"Missing score outputs for {prefix}.")
    scores = pd.read_csv(scores_path)
    summary = pd.read_csv(summary_path)
    vital_score = pd.to_numeric(scores.get("vital_score"), errors="coerce")
    red = scores.get("vital_red_flag", pd.Series(False, index=scores.index)).fillna(False).astype(str).str.lower().isin(
        {"true", "1", "yes"}
    )
    frequency_status = scores.get(
        "frequency_evidence_status",
        pd.Series("", index=scores.index, dtype="object"),
    ).fillna("").astype(str)
    observed = frequency_status.eq("frequency_observed")
    band_counts = (
        scores.groupby("vital_band", dropna=False)
        .size()
        .rename("variant_count")
        .reset_index()
    )
    band_counts.insert(0, "domain_type", domain_type)
    band_counts.insert(0, "domain", domain)
    band_counts.insert(0, "output_prefix", prefix)
    band_counts["variant_fraction"] = band_counts["variant_count"] / len(scores) if len(scores) else np.nan

    row = {
        "output_prefix": prefix,
        "domain": domain,
        "domain_type": domain_type,
        "variant_count": len(scores),
        "frequency_observed_count": int(observed.sum()),
        "no_frequency_evidence_count": int((~observed).sum()),
        "naive_af_gt_1e_5_flags": read_metric(summary, "standard_popmax_or_global_af_gt_1e_5_flags"),
        "ac_supported_frequency_flags": read_metric(summary, "ac_supported_frequency_flags"),
        "vital_red_count": int(red.sum()),
        "score_ge_50_count": int((vital_score >= 50).sum()),
        "score_ge_60_count": int((vital_score >= 60).sum()),
        "score_ge_70_count": int((vital_score >= 70).sum()),
        "score_ge_80_count": int((vital_score >= 80).sum()),
        "max_vital_score": float(vital_score.max()) if vital_score.notna().any() else np.nan,
        "p95_vital_score": float(vital_score.quantile(0.95)) if vital_score.notna().any() else np.nan,
        "median_observed_vital_score": float(vital_score.loc[observed].median()) if observed.any() else np.nan,
    }
    return row, band_counts


def score_distribution_for_prefix(prefix: str, domain: str, domain_type: str) -> pd.DataFrame:
    scores_path = data_path(prefix, "vital_scores.csv")
    if not scores_path.exists():
        raise FileNotFoundError(f"Missing score table for {prefix}.")
    scores = pd.read_csv(scores_path)
    keep_columns = [
        "variant_key",
        "gene",
        "clinvar_id",
        "variation_id",
        "title",
        "clinsig",
        "review_status",
        "review_strength",
        "frequency_evidence_status",
        "variant_type",
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
    present_columns = [column for column in keep_columns if column in scores.columns]
    distribution = scores[present_columns].copy()
    distribution.insert(0, "domain_type", domain_type)
    distribution.insert(0, "domain", domain)
    distribution.insert(0, "output_prefix", prefix)
    return distribution


def plot_external_score_distribution(score_distribution: pd.DataFrame, output_path: Path) -> None:
    plot_df = score_distribution.copy()
    plot_df["vital_score"] = pd.to_numeric(plot_df["vital_score"], errors="coerce")
    plot_df = plot_df.loc[
        plot_df["frequency_evidence_status"].eq("frequency_observed") & plot_df["vital_score"].notna()
    ].copy()
    domain_order = [
        "cardiomyopathy",
        "epilepsy",
        "hearing_loss",
        "random_clinvar_plp",
        "BRCA_MMR_APC_control",
    ]
    domain_order = [domain for domain in domain_order if domain in set(plot_df["domain"])]
    if not domain_order:
        return

    values = [plot_df.loc[plot_df["domain"].eq(domain), "vital_score"].to_numpy() for domain in domain_order]
    fig, ax = plt.subplots(figsize=(11, 5.6))
    box = ax.boxplot(
        values,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#101820", "linewidth": 1.5},
        boxprops={"linewidth": 1.1},
        whiskerprops={"linewidth": 1.0},
        capprops={"linewidth": 1.0},
    )
    palette = ["#3b7ddd", "#7b6fd6", "#d18f2f", "#60866a", "#7b8794"]
    for patch, color in zip(box["boxes"], palette):
        patch.set_facecolor(color)
        patch.set_alpha(0.38)
        patch.set_edgecolor("#263238")

    rng = np.random.default_rng(17)
    for index, scores in enumerate(values, start=1):
        if len(scores) == 0:
            continue
        jitter = rng.normal(0, 0.035, len(scores))
        ax.scatter(
            np.full(len(scores), index) + jitter,
            scores,
            s=18,
            alpha=0.42,
            color=palette[(index - 1) % len(palette)],
            edgecolors="none",
        )

    ax.axhline(70, color="#b3261e", linestyle="--", linewidth=1.3, label="Urgent-review cutoff")
    ax.axhline(60, color="#d18f2f", linestyle=":", linewidth=1.1, label="watchlist cutoff")
    ax.set_ylabel("Score among frequency-observed variants")
    ax.set_xlabel("")
    ax.set_ylim(-3, 105)
    ax.set_title("External-domain score distribution")
    ax.set_xticks(range(1, len(domain_order) + 1))
    ax.set_xticklabels(domain_order)
    ax.grid(axis="y", color="#d8dde6", linewidth=0.8, alpha=0.7)
    ax.legend(frameon=False, loc="upper right")
    ax.tick_params(axis="x", rotation=22)
    for label in ax.get_xticklabels():
        label.set_horizontalalignment("right")
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)
    print(f"Saved {output_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run small external-domain panels and summarize score distributions."
    )
    parser.add_argument("--variant-summary", type=Path, default=BASE_DIR / "data" / "variant_summary.txt.gz")
    parser.add_argument("--panel-sample-size", type=int, default=300)
    parser.add_argument("--random-sample-size", type=int, default=500)
    parser.add_argument("--random-seed", type=int, default=37)
    parser.add_argument("--dataset", default=ava.GNOMAD_DATASET)
    parser.add_argument("--gnomad-pause", type=float, default=0.2)
    parser.add_argument("--no-fetch-gnomad", action="store_true")
    parser.add_argument("--force-gnomad-fetch", action="store_true")
    parser.add_argument(
        "--skip-scoring",
        action="store_true",
        help="Only regenerate combined summaries from existing per-prefix outputs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    panel_specs: list[tuple[str, str, str, list[str] | None, int]] = [
        (
            "vital_domain_cardiomyopathy",
            "cardiomyopathy",
            "external_disease_panel",
            PANEL_DEFINITIONS["cardiomyopathy"],
            args.panel_sample_size,
        ),
        (
            "vital_domain_epilepsy",
            "epilepsy",
            "external_disease_panel",
            PANEL_DEFINITIONS["epilepsy"],
            args.panel_sample_size,
        ),
        (
            "vital_domain_hearing_loss",
            "hearing_loss",
            "external_disease_panel",
            PANEL_DEFINITIONS["hearing_loss"],
            args.panel_sample_size,
        ),
        (
            "vital_random_clinvar_plp",
            "random_clinvar_plp",
            "random_clinvar_sample",
            None,
            args.random_sample_size,
        ),
    ]

    if not args.skip_scoring:
        for prefix, domain, _, genes, sample_size in panel_specs:
            print(f"Scoring {domain} ({sample_size} sampled variants) -> {prefix}")
            score_snapshot(
                variant_summary_path=args.variant_summary,
                output_prefix=prefix,
                genes=genes,
                dataset=args.dataset,
                pause=args.gnomad_pause,
                fetch_gnomad=not args.no_fetch_gnomad,
                force_gnomad=args.force_gnomad_fetch,
                max_variants=None,
                sample_variants=sample_size,
                random_seed=args.random_seed,
            )

    summary_rows: list[dict[str, object]] = []
    band_frames: list[pd.DataFrame] = []
    score_frames: list[pd.DataFrame] = []
    for prefix, domain, domain_type, _, _ in panel_specs:
        row, bands = summarize_prefix(prefix, domain, domain_type)
        summary_rows.append(row)
        band_frames.append(bands)
        score_frames.append(score_distribution_for_prefix(prefix, domain, domain_type))
    if data_path("control", "vital_scores.csv").exists():
        row, bands = summarize_prefix("control", "BRCA_MMR_APC_control", "negative_control")
        summary_rows.append(row)
        band_frames.append(bands)
        score_frames.append(score_distribution_for_prefix("control", "BRCA_MMR_APC_control", "negative_control"))

    save_table(pd.DataFrame(summary_rows), DATA_DIR / "vital_external_panel_summary.csv")
    save_table(pd.concat(band_frames, ignore_index=True), DATA_DIR / "vital_external_panel_band_distribution.csv")
    score_distribution = pd.concat(score_frames, ignore_index=True)
    save_table(score_distribution, DATA_DIR / "vital_external_panel_score_distribution.csv")
    plot_external_score_distribution(score_distribution, FIGURE_DIR / "vital_external_panel_score_distribution.png")


if __name__ == "__main__":
    main()
