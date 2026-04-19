from __future__ import annotations

import argparse
import math
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
DEFAULT_GNOMAD = BASE_DIR / "data" / "processed" / "gnomad_subset_with_header.tsv"
DEFAULT_CLINVAR = BASE_DIR / "data" / "processed" / "clinvar_gnomad_merged.csv"
DEFAULT_OUTPUT = BASE_DIR / "figures" / "figure_final_publication.png"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create the publication figure comparing gnomAD and ClinVar AF distributions."
    )
    parser.add_argument("--gnomad", type=Path, default=DEFAULT_GNOMAD)
    parser.add_argument("--clinvar", type=Path, default=DEFAULT_CLINVAR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def load_af_values(path: Path, separator: str) -> pd.Series:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    df = pd.read_csv(path, sep=separator, dtype=str, keep_default_na=False)
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
    numeric_af = pd.to_numeric(af, errors="coerce")
    return numeric_af.loc[numeric_af > 0].dropna()


def log10_series(values: pd.Series) -> pd.Series:
    return values.map(math.log10)


def make_figure(gnomad_af: pd.Series, clinvar_af: pd.Series, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    gnomad_log = log10_series(gnomad_af)
    clinvar_log = log10_series(clinvar_af)
    lower = -7.0
    upper = -3.0
    bins = np.linspace(lower, upper, 55)

    fig, ax = plt.subplots(figsize=(12, 7))
    ax.hist(
        gnomad_log,
        bins=bins,
        density=True,
        alpha=0.45,
        color="#78a9cf",
        label="gnomAD (population)",
    )
    ax.hist(
        clinvar_log,
        bins=bins,
        density=True,
        alpha=0.65,
        color="#f28e2b",
        label="ClinVar (clinical)",
    )
    ax.axvline(math.log10(1e-5), color="#1f77b4", linestyle="--", linewidth=2)
    ax.axvline(math.log10(1e-4), color="#1f77b4", linestyle="--", linewidth=2)
    ax.plot([], [], color="#1f77b4", linestyle="--", label="Rare disease threshold (~1e-5)")
    ax.plot(
        [],
        [],
        color="#1f77b4",
        linestyle="--",
        label="Upper bound for Mendelian pathogenicity (~1e-4)",
    )

    ax.set_title(
        "Clinical databases capture a frequency-constrained subset of population variation",
        pad=12,
    )
    ax.set_xlabel("log10(Allele Frequency)")
    ax.set_ylabel("Density")
    ax.set_xlim(lower, upper)
    ax.legend(frameon=False, loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    gnomad_af = load_af_values(args.gnomad, separator="\t")
    clinvar_af = load_af_values(args.clinvar, separator=",")
    make_figure(gnomad_af, clinvar_af, args.output)
    print(f"Saved publication figure to: {args.output}")


if __name__ == "__main__":
    main()
