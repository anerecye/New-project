from __future__ import annotations

from pathlib import Path

try:
    import pandas as pd
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires pandas. Install it with: python -m pip install pandas"
    ) from exc


BASE_DIR = Path(__file__).resolve().parents[1]
INPUT_PATH = BASE_DIR / "data" / "processed" / "clinvar_arrhythmia_clean.csv"
WITH_KEYS_PATH = BASE_DIR / "data" / "processed" / "clinvar_arrhythmia_with_keys.csv"
HIGH_CONFIDENCE_PATH = (
    BASE_DIR / "data" / "processed" / "clinvar_arrhythmia_high_confidence.csv"
)
DUPLICATE_ROWS_PATH = BASE_DIR / "data" / "processed" / "duplicate_variant_rows.csv"
DUPLICATE_GROUPS_PATH = (
    BASE_DIR / "data" / "processed" / "duplicate_variant_group_counts.csv"
)

REQUIRED_COLUMNS = (
    "GeneSymbol",
    "ClinicalSignificance",
    "ReviewStatus",
    "Chromosome",
    "Start",
    "ReferenceAllele",
    "AlternateAllele",
    "Assembly",
)

HIGH_CONFIDENCE_REVIEW_STATUSES = {
    "criteria provided, multiple submitters, no conflicts",
    "reviewed by expert panel",
}


def load_variants(input_path: Path) -> pd.DataFrame:
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    return pd.read_csv(input_path, dtype=str, keep_default_na=False)


def validate_required_columns(df: pd.DataFrame) -> None:
    missing = [column for column in REQUIRED_COLUMNS if column not in df.columns]
    if missing:
        missing_str = ", ".join(missing)
        raise ValueError(f"Missing required columns: {missing_str}")


def standardize_chromosome_names(df: pd.DataFrame) -> pd.DataFrame:
    standardized = df.copy()
    standardized["Chromosome"] = (
        standardized["Chromosome"]
        .astype(str)
        .str.strip()
        .str.replace(r"^chr", "", regex=True, case=False)
    )
    return standardized


def add_variant_key(df: pd.DataFrame) -> pd.DataFrame:
    keyed = df.copy()
    key_columns = ["Chromosome", "Start", "ReferenceAllele", "AlternateAllele"]
    keyed[key_columns] = keyed[key_columns].apply(lambda column: column.astype(str).str.strip())
    keyed["variant_key"] = keyed[key_columns].agg(":".join, axis=1)
    return keyed


def get_duplicate_rows(df: pd.DataFrame) -> pd.DataFrame:
    duplicate_mask = df["variant_key"].duplicated(keep=False)
    duplicate_rows = df.loc[duplicate_mask].copy()
    return duplicate_rows.sort_values(["variant_key", "GeneSymbol", "Start"]).reset_index(
        drop=True
    )


def get_duplicate_group_counts(df: pd.DataFrame) -> pd.DataFrame:
    group_counts = (
        df.groupby("variant_key", dropna=False)
        .size()
        .reset_index(name="row_count")
        .sort_values(["row_count", "variant_key"], ascending=[False, True])
    )
    duplicate_groups = group_counts.loc[group_counts["row_count"] > 1].copy()
    return duplicate_groups.reset_index(drop=True)


def filter_high_confidence(df: pd.DataFrame) -> pd.DataFrame:
    review_status = df["ReviewStatus"].str.strip().str.lower()
    mask = review_status.isin(HIGH_CONFIDENCE_REVIEW_STATUSES)
    return df.loc[mask].copy().reset_index(drop=True)


def summarize_variants(df: pd.DataFrame, label: str) -> None:
    duplicate_rows = get_duplicate_rows(df)
    duplicate_groups = get_duplicate_group_counts(df)

    print(f"\n{label}")
    print(f"  Total rows: {len(df)}")
    print(f"  Unique variant_key values: {df['variant_key'].nunique()}")
    print(f"  Duplicated variant_key rows: {len(duplicate_rows)}")
    print(f"  Duplicated variant groups: {len(duplicate_groups)}")


def save_dataframe(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)


def main() -> None:
    variants = load_variants(INPUT_PATH)
    validate_required_columns(variants)

    variants = standardize_chromosome_names(variants)
    variants = add_variant_key(variants)

    duplicate_rows = get_duplicate_rows(variants)
    duplicate_group_counts = get_duplicate_group_counts(variants)
    high_confidence = filter_high_confidence(variants)

    save_dataframe(variants, WITH_KEYS_PATH)
    save_dataframe(duplicate_rows, DUPLICATE_ROWS_PATH)
    save_dataframe(duplicate_group_counts, DUPLICATE_GROUPS_PATH)
    save_dataframe(high_confidence, HIGH_CONFIDENCE_PATH)

    summarize_variants(variants, "All variants")
    summarize_variants(high_confidence, "High-confidence subset")


if __name__ == "__main__":
    main()
