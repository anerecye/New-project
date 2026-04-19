from __future__ import annotations

import argparse
from pathlib import Path

try:
    import pandas as pd
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires pandas. Install it with: python -m pip install pandas"
    ) from exc


BASE_DIR = Path(__file__).resolve().parents[1]
CLINVAR_PATH = BASE_DIR / "data" / "processed" / "clinvar_variant_level.csv"
OUTPUT_PATH = BASE_DIR / "data" / "processed" / "clinvar_gnomad_merged.csv"

CLINVAR_REQUIRED_COLUMNS = (
    "variant_key",
    "Chromosome",
    "Start",
    "ReferenceAllele",
    "AlternateAllele",
)

GNOMAD_COLUMN_CANDIDATES = {
    "Chromosome": ("Chromosome", "chromosome", "#CHROM", "CHROM"),
    "Position": ("Position", "position", "POS", "Start"),
    "ReferenceAllele": ("ReferenceAllele", "referenceallele", "REF", "Ref"),
    "AlternateAllele": ("AlternateAllele", "alternateallele", "ALT", "Alt"),
    "AF": ("AF", "af"),
    "AC": ("AC", "ac"),
    "AN": ("AN", "an"),
}

MISSING_TOKENS = {"", ".", "na", "nan", "none", "null", "n/a"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge ClinVar variant-level data with a preprocessed gnomAD table."
    )
    parser.add_argument(
        "--gnomad",
        type=Path,
        required=True,
        help="Path to the preprocessed gnomAD CSV or TSV file.",
    )
    return parser.parse_args()


def infer_separator(path: Path) -> str | None:
    suffixes = [suffix.lower() for suffix in path.suffixes]
    if ".tsv" in suffixes or ".txt" in suffixes:
        return "\t"
    if ".csv" in suffixes:
        return ","
    return None


def load_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    separator = infer_separator(path)
    return pd.read_csv(
        path,
        sep=separator,
        engine="python" if separator is None else None,
        compression="infer",
        dtype=str,
        keep_default_na=False,
    )


def validate_required_columns(df: pd.DataFrame, required_columns: tuple[str, ...], label: str) -> None:
    missing = [column for column in required_columns if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in {label}: {', '.join(missing)}")


def resolve_gnomad_columns(df: pd.DataFrame) -> dict[str, str]:
    resolved: dict[str, str] = {}
    for canonical_name, candidates in GNOMAD_COLUMN_CANDIDATES.items():
        for candidate in candidates:
            if candidate in df.columns:
                resolved[canonical_name] = candidate
                break
        else:
            candidate_str = ", ".join(candidates)
            raise ValueError(
                f"Missing required gnomAD column for {canonical_name}. "
                f"Expected one of: {candidate_str}"
            )
    return resolved


def standardize_chromosome(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.strip()
        .str.replace(r"^chr", "", regex=True, case=False)
    )


def normalize_field(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip()


def normalize_field_with_missing(series: pd.Series) -> pd.Series:
    normalized = normalize_field(series)
    missing_mask = normalized.str.lower().isin(MISSING_TOKENS)
    return normalized.mask(missing_mask)


def build_match_key(
    chromosome: pd.Series,
    position: pd.Series,
    reference: pd.Series,
    alternate: pd.Series,
) -> pd.Series:
    key_df = pd.DataFrame(
        {
            "Chromosome": standardize_chromosome(chromosome),
            "Position": normalize_field_with_missing(position),
            "ReferenceAllele": normalize_field_with_missing(reference),
            "AlternateAllele": normalize_field_with_missing(alternate),
        }
    )
    has_missing = key_df.isna().any(axis=1)
    keys = pd.Series(pd.NA, index=key_df.index, dtype="object")

    if (~has_missing).any():
        valid_rows = key_df.loc[~has_missing].astype(str)
        keys.loc[~has_missing] = valid_rows.apply(lambda row: ":".join(row.tolist()), axis=1)

    return keys


def prepare_clinvar(df: pd.DataFrame) -> pd.DataFrame:
    validate_required_columns(df, CLINVAR_REQUIRED_COLUMNS, "ClinVar")
    prepared = df.copy()
    prepared["Chromosome"] = standardize_chromosome(prepared["Chromosome"])
    prepared["_match_key"] = build_match_key(
        prepared["Chromosome"],
        prepared["Start"],
        prepared["ReferenceAllele"],
        prepared["AlternateAllele"],
    )
    return prepared


def prepare_gnomad(df: pd.DataFrame) -> pd.DataFrame:
    resolved_columns = resolve_gnomad_columns(df)
    prepared = df.rename(columns={actual: canonical for canonical, actual in resolved_columns.items()})
    prepared = prepared.loc[
        :, ["Chromosome", "Position", "ReferenceAllele", "AlternateAllele", "AF", "AC", "AN"]
    ].copy()
    prepared["Chromosome"] = standardize_chromosome(prepared["Chromosome"])
    prepared["_match_key"] = build_match_key(
        prepared["Chromosome"],
        prepared["Position"],
        prepared["ReferenceAllele"],
        prepared["AlternateAllele"],
    )
    prepared["AF"] = pd.to_numeric(prepared["AF"], errors="coerce")
    prepared["AC"] = pd.to_numeric(prepared["AC"], errors="coerce")
    prepared["AN"] = pd.to_numeric(prepared["AN"], errors="coerce")

    prepared = prepared.dropna(subset=["_match_key"])
    prepared = prepared.drop_duplicates(subset=["_match_key"], keep="first")
    return prepared


def merge_clinvar_gnomad(clinvar_df: pd.DataFrame, gnomad_df: pd.DataFrame) -> pd.DataFrame:
    merged = clinvar_df.merge(
        gnomad_df[["_match_key", "AF", "AC", "AN"]],
        on="_match_key",
        how="left",
        indicator=True,
    )
    merged["gnomad_match"] = merged["_merge"].eq("both")
    return merged.drop(columns=["_merge", "_match_key"])


def print_series(title: str, series: pd.Series) -> None:
    print(f"\n{title}")
    if series.empty:
        print("  None")
        return
    for label, value in series.items():
        print(f"  {label}: {value}")


def print_af_summary(merged_df: pd.DataFrame) -> None:
    matched_af = merged_df.loc[merged_df["gnomad_match"], "AF"].dropna()
    print("\nAF distribution for matched variants")
    if matched_af.empty:
        print("  No matched variants with non-missing AF values.")
        return

    describe = matched_af.describe(percentiles=[0.25, 0.5, 0.75])
    for stat_name, value in describe.items():
        print(f"  {stat_name}: {value}")


def save_merged_dataset(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)


def main() -> None:
    args = parse_args()

    print(f"Loading ClinVar variant-level data from: {CLINVAR_PATH}")
    clinvar_df = load_table(CLINVAR_PATH)
    print(f"Loading gnomAD data from: {args.gnomad}")
    gnomad_df = load_table(args.gnomad)

    clinvar_prepared = prepare_clinvar(clinvar_df)
    gnomad_prepared = prepare_gnomad(gnomad_df)
    merged = merge_clinvar_gnomad(clinvar_prepared, gnomad_prepared)

    save_merged_dataset(merged, OUTPUT_PATH)

    matched_count = int(merged["gnomad_match"].sum())
    unmatched_count = int((~merged["gnomad_match"]).sum())

    print(f"\nSaved merged dataset to: {OUTPUT_PATH}")
    print(f"Total ClinVar variants: {len(merged)}")
    print(f"Variants with gnomAD match: {matched_count}")
    print(f"Variants without gnomAD match: {unmatched_count}")
    print_af_summary(merged)


if __name__ == "__main__":
    main()
