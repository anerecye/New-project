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
OUTPUT_PATH = BASE_DIR / "data" / "processed" / "gnomad_lookup.bed"

REQUIRED_COLUMNS = ("Chromosome", "Start")
MISSING_TOKENS = {"", "na", "nan", "none", "null", "n/a"}


def load_lookup_input(input_path: Path) -> pd.DataFrame:
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    return pd.read_csv(input_path, dtype=str, keep_default_na=False)


def validate_required_columns(df: pd.DataFrame) -> None:
    missing = [column for column in REQUIRED_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")


def normalize_text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def validate_required_values(df: pd.DataFrame) -> None:
    chromosome = normalize_text(df["Chromosome"])
    start = normalize_text(df["Start"])

    missing_mask = chromosome.str.lower().isin(MISSING_TOKENS) | start.str.lower().isin(
        MISSING_TOKENS
    )
    if missing_mask.any():
        raise ValueError("Found rows with missing Chromosome or Start values.")


def standardize_chromosome(series: pd.Series) -> pd.Series:
    normalized = normalize_text(series).str.replace(r"^chr", "", regex=True, case=False)
    return "chr" + normalized


def make_bed_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    validate_required_values(df)

    bed = pd.DataFrame()
    bed["Chromosome"] = standardize_chromosome(df["Chromosome"])
    bed["bed_end"] = pd.to_numeric(df["Start"], errors="raise").astype(int)
    bed["bed_start"] = bed["bed_end"] - 1

    if (bed["bed_start"] < 0).any():
        raise ValueError("Encountered Start positions below 1, which are invalid for BED.")

    bed = bed.loc[:, ["Chromosome", "bed_start", "bed_end"]]
    return bed.drop_duplicates().sort_values(["Chromosome", "bed_start", "bed_end"]).reset_index(
        drop=True
    )


def save_bed(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", header=False, index=False)


def main() -> None:
    lookup_df = load_lookup_input(INPUT_PATH)
    validate_required_columns(lookup_df)

    bed_df = make_bed_dataframe(lookup_df)
    save_bed(bed_df, OUTPUT_PATH)

    print(f"Total input rows: {len(lookup_df)}")
    print(f"Number of BED rows: {len(bed_df)}")


if __name__ == "__main__":
    main()
