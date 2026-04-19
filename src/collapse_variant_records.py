from __future__ import annotations

from pathlib import Path

try:
    import pandas as pd
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires pandas. Install it with: python -m pip install pandas"
    ) from exc


BASE_DIR = Path(__file__).resolve().parents[1]
INPUT_PATH = BASE_DIR / "data" / "processed" / "clinvar_arrhythmia_with_keys.csv"
OUTPUT_PATH = BASE_DIR / "data" / "processed" / "clinvar_variant_level.csv"

REQUIRED_COLUMNS = (
    "variant_key",
    "GeneSymbol",
    "ClinicalSignificance",
    "ReviewStatus",
    "Chromosome",
    "Start",
    "ReferenceAllele",
    "AlternateAllele",
)

REVIEW_STATUS_PRIORITY = (
    "reviewed by expert panel",
    "criteria provided, multiple submitters, no conflicts",
    "criteria provided, single submitter",
    "no assertion criteria provided",
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
        raise ValueError(f"Missing required columns: {', '.join(missing)}")


def ordered_unique(values: pd.Series) -> list[str]:
    seen: set[str] = set()
    unique_values: list[str] = []
    for value in values.fillna("").astype(str).str.strip():
        if not value or value in seen:
            continue
        seen.add(value)
        unique_values.append(value)
    return unique_values


def join_unique_values(values: pd.Series) -> str:
    unique_values = ordered_unique(values)
    return ";".join(unique_values)


def choose_best_review_status(values: pd.Series) -> str:
    unique_values = ordered_unique(values)
    unique_value_set = set(unique_values)

    for status in REVIEW_STATUS_PRIORITY:
        if status in unique_value_set:
            return status

    if unique_values:
        return unique_values[0]
    return ""


def is_high_confidence(values: pd.Series) -> bool:
    review_statuses = {value.lower() for value in ordered_unique(values)}
    return any(status in review_statuses for status in HIGH_CONFIDENCE_REVIEW_STATUSES)


def collapse_variant_group(group: pd.DataFrame) -> pd.Series:
    return pd.Series(
        {
            "variant_key": group.name,
            "GeneSymbol": join_unique_values(group["GeneSymbol"]),
            "Chromosome": group["Chromosome"].iat[0],
            "Start": group["Start"].iat[0],
            "ReferenceAllele": group["ReferenceAllele"].iat[0],
            "AlternateAllele": group["AlternateAllele"].iat[0],
            "n_clinvar_rows": len(group),
            "ClinicalSignificance_values": join_unique_values(
                group["ClinicalSignificance"]
            ),
            "ReviewStatus_values": join_unique_values(group["ReviewStatus"]),
            "best_review_status": choose_best_review_status(group["ReviewStatus"]),
            "is_high_confidence": is_high_confidence(group["ReviewStatus"]),
        }
    )


def collapse_variants(df: pd.DataFrame) -> pd.DataFrame:
    collapsed = (
        df.groupby("variant_key", sort=True, dropna=False)
        .apply(collapse_variant_group, include_groups=False)
        .reset_index(drop=True)
    )
    collapsed["n_clinvar_rows"] = collapsed["n_clinvar_rows"].astype(int)
    collapsed["is_high_confidence"] = collapsed["is_high_confidence"].astype(bool)
    return collapsed


def save_dataframe(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)


def print_counts(title: str, counts: pd.Series) -> None:
    print(f"\n{title}")
    if counts.empty:
        print("  None")
        return
    for label, count in counts.items():
        print(f"  {label}: {count}")


def print_summary(df: pd.DataFrame) -> None:
    print(f"Total unique variants: {len(df)}")
    print_counts(
        "Distribution of n_clinvar_rows",
        df["n_clinvar_rows"].value_counts().sort_index(),
    )
    print_counts(
        "Counts by best_review_status",
        df["best_review_status"].replace("", "Missing").value_counts(),
    )
    print(f"\nHigh-confidence variants: {int(df['is_high_confidence'].sum())}")


def main() -> None:
    variants = load_variants(INPUT_PATH)
    validate_required_columns(variants)

    collapsed = collapse_variants(variants)
    save_dataframe(collapsed, OUTPUT_PATH)
    print_summary(collapsed)


if __name__ == "__main__":
    main()
