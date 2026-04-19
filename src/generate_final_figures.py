"""
Generate publication-ready figures from processed ClinVar/gnomAD outputs.

Outputs:
  figures/figure1_af_overlay_kde.png
  figures/figure1b_af_categories.png
  figures/figure2_per_gene_outlier_burden.png
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "processed"
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

GNOMAD_AF_PATH = DATA / "gnomad_subset_with_header.tsv"
CLINVAR_AF_PATH = DATA / "clinvar_gnomad_merged.csv"
AF_CATEGORY_PATH = DATA / "af_category_counts.csv"
OUTLIER_GENE_PATH = DATA / "outlier_counts_by_gene.csv"

GNOMAD_COLOR = "#4C72B0"
CLINVAR_COLOR = "#DD8452"
COMMON_COLOR = "#C44E52"
UNDERPOWERED_COLOR = "#B0B0B0"
ALPHA_FILL = 0.35

THRESH_1 = 1e-5
THRESH_2 = 1e-4
LOG_FLOOR = 1e-10


plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "figure.dpi": 150,
    }
)


def read_af_column(path: Path, sep: str = ",") -> pd.Series:
    df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
    if "AF" not in df.columns:
        raise ValueError(f"Missing AF column in {path}")

    af = (
        df["AF"]
        .astype(str)
        .str.split(",")
        .explode()
        .str.strip()
        .replace("", pd.NA)
    )
    return pd.to_numeric(af, errors="coerce").dropna()


def log_af(series: pd.Series) -> np.ndarray:
    return np.log10(series.clip(lower=LOG_FLOOR).to_numpy())


def kde_curve(values: np.ndarray, x_grid: np.ndarray) -> np.ndarray:
    kde = gaussian_kde(values, bw_method="scott")
    return kde(x_grid)


def find_column(columns: pd.Index, *tokens: str, exclude: tuple[str, ...] = ()) -> str:
    normalized = [str(column).lower().strip() for column in columns]
    for original, lowered in zip(columns, normalized):
        if all(token in lowered for token in tokens) and not any(
            token in lowered for token in exclude
        ):
            return str(original)
    token_text = ", ".join(tokens)
    raise ValueError(f"Could not detect column containing: {token_text}")


def figure1a() -> None:
    gnomad_af = log_af(read_af_column(GNOMAD_AF_PATH, sep="\t"))
    clinvar_af = log_af(read_af_column(CLINVAR_AF_PATH))

    x = np.linspace(
        min(gnomad_af.min(), clinvar_af.min()) - 0.2,
        max(gnomad_af.max(), clinvar_af.max()) + 0.2,
        500,
    )

    kde_gnomad = kde_curve(gnomad_af, x)
    kde_clinvar = kde_curve(clinvar_af, x)

    fig, ax = plt.subplots(figsize=(6.5, 3.8))

    ax.fill_between(x, kde_gnomad, alpha=ALPHA_FILL, color=GNOMAD_COLOR)
    ax.fill_between(x, kde_clinvar, alpha=ALPHA_FILL, color=CLINVAR_COLOR)
    ax.plot(
        x,
        kde_gnomad,
        color=GNOMAD_COLOR,
        lw=1.6,
        label=f"gnomAD exomes queried (n={len(gnomad_af):,})",
    )
    ax.plot(
        x,
        kde_clinvar,
        color=CLINVAR_COLOR,
        lw=1.6,
        label=f"ClinVar-matched (n={len(clinvar_af):,})",
    )

    y_top = max(kde_gnomad.max(), kde_clinvar.max())
    for threshold, label in [(THRESH_1, "10^-5"), (THRESH_2, "10^-4")]:
        ax.axvline(np.log10(threshold), color="black", lw=0.9, ls="--", alpha=0.7)
        ax.text(
            np.log10(threshold) + 0.03,
            y_top * 0.92,
            f"AF = {label}",
            fontsize=8,
            va="top",
            color="black",
            rotation=90,
            bbox={"boxstyle": "round,pad=0.15", "facecolor": "white", "edgecolor": "none", "alpha": 0.8},
        )

    ax.set_xlabel("log10(Allele frequency)", labelpad=6)
    ax.set_ylabel("Density", labelpad=6)
    ax.set_title(
        "Allele frequency distribution: queried gnomAD exomes vs ClinVar-matched variants",
        fontsize=10,
        pad=8,
    )
    ax.legend(frameon=False, fontsize=9, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()

    out = FIG / "figure1_af_overlay_kde.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out.relative_to(ROOT)}")


def figure1b() -> None:
    cats = pd.read_csv(AF_CATEGORY_PATH)
    cats.columns = cats.columns.str.lower().str.strip()

    cat_col = find_column(cats.columns, "cat")
    count_col = find_column(cats.columns, "count")

    order = ["ultra_rare", "rare", "common"]
    labels = ["Ultra-rare\n(AF < 10^-5)", "Rare\n(10^-5 - 10^-4)", "Common\n(AF > 10^-4)"]
    colors = [GNOMAD_COLOR, CLINVAR_COLOR, COMMON_COLOR]

    cats = cats.set_index(cat_col).reindex(order)
    counts = cats[count_col].fillna(0).astype(int)

    fig, ax = plt.subplots(figsize=(4.5, 3.5))
    bars = ax.bar(labels, counts, color=colors, width=0.55, edgecolor="white", linewidth=0.5)

    for bar, value in zip(bars, counts):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + counts.max() * 0.02,
            f"{value:,}",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    ax.set_ylabel("Number of variants", labelpad=6)
    ax.set_title("ClinVar-matched variants by AF category", fontsize=10, pad=8)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda value, _: f"{int(value):,}"))
    ax.set_ylim(0, counts.max() * 1.18)
    fig.tight_layout()

    out = FIG / "figure1b_af_categories.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out.relative_to(ROOT)}")


def figure2() -> None:
    gene_df = pd.read_csv(OUTLIER_GENE_PATH)
    gene_df.columns = gene_df.columns.str.lower().str.strip().str.replace(" ", "_")

    gene_col = find_column(gene_df.columns, "gene")
    fraction_col = "within_gene_fraction"
    denom_col = "variants_with_af"
    outlier_col = find_column(gene_df.columns, "outlier", exclude=("fraction", "all"))

    required = {fraction_col, denom_col}
    missing = required.difference(gene_df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")

    gene_df = gene_df.sort_values(fraction_col, ascending=False)

    genes = gene_df[gene_col].tolist()
    fractions = gene_df[fraction_col].tolist()
    outliers = gene_df[outlier_col].astype(int).tolist()
    denominators = gene_df[denom_col].astype(int).tolist()

    bar_colors = [
        UNDERPOWERED_COLOR if str(gene).upper() == "CALM2" else GNOMAD_COLOR
        for gene in genes
    ]

    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    bars = ax.bar(
        genes,
        fractions,
        color=bar_colors,
        width=0.6,
        edgecolor="white",
        linewidth=0.5,
    )

    y_padding = max(fractions) * 0.02 if max(fractions) else 0.01
    for bar, outlier_count, denominator in zip(bars, outliers, denominators):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + y_padding,
            f"{outlier_count}/{denominator}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    ax.axhline(0, color="black", lw=0.6)
    ax.set_ylabel("Within-gene outlier fraction\n(AF > 10^-5)", labelpad=6)
    ax.set_title("Per-gene outlier burden among AF-covered ClinVar variants", fontsize=10, pad=8)
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
    ax.set_ylim(0, max(fractions) * 1.25 if max(fractions) else 1)

    fig.text(
        0.5,
        -0.04,
        "Grey bar (CALM2): n = 1 AF-covered variant; underpowered.",
        ha="center",
        fontsize=8,
        color="#666666",
    )
    fig.tight_layout()

    out = FIG / "figure2_per_gene_outlier_burden.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out.relative_to(ROOT)}")


if __name__ == "__main__":
    figure1a()
    figure1b()
    figure2()
    print("All figures generated.")
