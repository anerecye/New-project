from __future__ import annotations

from pathlib import Path

try:
    import pandas as pd
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires pandas. Install it with: python -m pip install pandas"
    ) from exc


BASE_DIR = Path(__file__).resolve().parents[1]
INPUT_PATH = BASE_DIR / "data" / "processed" / "gnomad_lookup_input.csv"
OUTPUT_PATH = BASE_DIR / "data" / "processed" / "gnomad_snv_sanity_sample.csv"

REQUIRED_COLUMNS = (
    "variant_key",
    "Chromosome",
    "Start",
    "ReferenceAllele",
    "AlternateAllele",
    "GeneSymbol",
)

VALID_BASES = {"A", "C", "G", "T"}
SAMPLE_SIZE = 10
RANDOM_SEED = 42


def load_lookup_input(input_path: Path) -> pd.DataFrame:
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    return pd.read_csv(input_path, dtype=str, keep_default_na=False)


def validate_required_columns(df: pd.DataFrame) -> None:
    missing = [column for column in REQUIRED_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")


def filter_snv_rows(df: pd.DataFrame) -> pd.DataFrame:
    snv_df = df.copy()
    snv_df["ReferenceAllele"] = snv_df["ReferenceAllele"].astype(str).str.strip().str.upper()
    snv_df["AlternateAllele"] = snv_df["AlternateAllele"].astype(str).str.strip().str.upper()

    ref_is_snv = snv_df["ReferenceAllele"].str.len().eq(1) & snv_df["ReferenceAllele"].isin(
        VALID_BASES
    )
    alt_is_snv = snv_df["AlternateAllele"].str.len().eq(1) & snv_df["AlternateAllele"].isin(
        VALID_BASES
    )
    return snv_df.loc[ref_is_snv & alt_is_snv].copy().reset_index(drop=True)


def sample_rows(df: pd.DataFrame) -> pd.DataFrame:
    if len(df) < SAMPLE_SIZE:
        raise ValueError(
            f"Not enough SNV rows to sample {SAMPLE_SIZE} variants. Found {len(df)}."
        )
    return df.sample(n=SAMPLE_SIZE, random_state=RANDOM_SEED).reset_index(drop=True)


def save_sample(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)


def print_sample(df: pd.DataFrame) -> None:
    columns = [
        "variant_key",
        "Chromosome",
        "Start",
        "ReferenceAllele",
        "AlternateAllele",
        "GeneSymbol",
    ]
    print(df.loc[:, columns].to_string(index=False))


def main() -> None:
    lookup_df = load_lookup_input(INPUT_PATH)
    validate_required_columns(lookup_df)

    snv_df = filter_snv_rows(lookup_df)
    sample_df = sample_rows(snv_df)
    save_sample(sample_df, OUTPUT_PATH)
    print_sample(sample_df)


if __name__ == "__main__":
    main()
