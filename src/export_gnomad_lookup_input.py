from __future__ import annotations

from pathlib import Path

try:
    import pandas as pd
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires pandas. Install it with: python -m pip install pandas"
    ) from exc


BASE_DIR = Path(__file__).resolve().parents[1]
INPUT_PATH = BASE_DIR / "data" / "processed" / "clinvar_variant_level.csv"
OUTPUT_PATH = BASE_DIR / "data" / "processed" / "gnomad_lookup_input.csv"

REQUIRED_COLUMNS = (
    "variant_key",
    "Chromosome",
    "Start",
    "ReferenceAllele",
    "AlternateAllele",
    "GeneSymbol",
    "best_review_status",
    "is_high_confidence",
)

GENOMIC_COLUMNS = ("Chromosome", "Start", "ReferenceAllele", "AlternateAllele")
MISSING_TOKENS = {"", "na", "nan", "none", "null", "n/a"}


def load_variants(input_path: Path) -> pd.DataFrame:
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    return pd.read_csv(input_path, dtype=str, keep_default_na=False)


def validate_required_columns(df: pd.DataFrame) -> None:
    missing = [column for column in REQUIRED_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")


def normalize_text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def is_missing_value(series: pd.Series) -> pd.Series:
    return normalize_text(series).str.lower().isin(MISSING_TOKENS)


def keep_one_row_per_variant(df: pd.DataFrame) -> pd.DataFrame:
    deduplicated = df.drop_duplicates(subset=["variant_key"], keep="first").copy()
    return deduplicated.reset_index(drop=True)


def filter_lookup_ready_rows(df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    missing_mask = pd.Series(False, index=df.index)
    for column in GENOMIC_COLUMNS:
        missing_mask = missing_mask | is_missing_value(df[column])

    retained = df.loc[~missing_mask, REQUIRED_COLUMNS].copy().reset_index(drop=True)
    excluded_count = int(missing_mask.sum())
    return retained, excluded_count


def save_lookup_table(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)


def main() -> None:
    variants = load_variants(INPUT_PATH)
    validate_required_columns(variants)

    total_variants = len(variants)
    variants = keep_one_row_per_variant(variants)
    lookup_ready, excluded_count = filter_lookup_ready_rows(variants)

    save_lookup_table(lookup_ready, OUTPUT_PATH)

    print(f"Total variants in input: {total_variants}")
    print(f"Number retained for lookup: {len(lookup_ready)}")
    print(f"Number excluded for missing genomic fields: {excluded_count}")


if __name__ == "__main__":
    main()
