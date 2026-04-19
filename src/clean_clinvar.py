from __future__ import annotations

import argparse
import re
from pathlib import Path

try:
    import pandas as pd
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires pandas. Install it with: python -m pip install pandas"
    ) from exc


TARGET_GENES = (
    "RYR2",
    "KCNQ1",
    "KCNH2",
    "SCN5A",
    "CASQ2",
    "CALM1",
    "CALM2",
    "CALM3",
    "TRDN",
)

BASE_DIR = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = BASE_DIR / "data" / "raw" / "variant_summary.txt.gz"
DEFAULT_OUTPUT = BASE_DIR / "data" / "processed" / "clinvar_arrhythmia_clean.csv"
CHUNK_SIZE = 100_000

REFERENCE_ALLELE_CANDIDATES = ("ReferenceAlleleVCF", "ReferenceAllele")
ALTERNATE_ALLELE_CANDIDATES = ("AlternateAlleleVCF", "AlternateAllele")
REQUIRED_COLUMNS = (
    "GeneSymbol",
    "ClinicalSignificance",
    "ReviewStatus",
    "Chromosome",
    "Start",
    "Assembly",
)
MISSING_VALUE_TOKENS = {"", "na"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Clean ClinVar variant_summary data for inherited arrhythmia genes."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Path to ClinVar variant_summary.txt.gz (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Path for cleaned CSV output (default: {DEFAULT_OUTPUT})",
    )
    return parser.parse_args()


def find_available_columns(
    columns: pd.Index, candidates: tuple[str, ...], label: str
) -> list[str]:
    available = [candidate for candidate in candidates if candidate in columns]
    if available:
        return available

    candidate_list = ", ".join(candidates)
    raise KeyError(f"Could not find a {label} column. Expected one of: {candidate_list}")


def detect_input_columns(input_path: Path) -> tuple[list[str], list[str], list[str]]:
    if not input_path.exists():
        raise FileNotFoundError(
            f"Input file not found: {input_path}\n"
            "Expected ClinVar data at data/raw/variant_summary.txt.gz."
        )

    header = pd.read_csv(
        input_path,
        sep="\t",
        compression="gzip",
        nrows=0,
    )
    header.columns = header.columns.str.strip()

    reference_cols = find_available_columns(
        header.columns, REFERENCE_ALLELE_CANDIDATES, "reference allele"
    )
    alternate_cols = find_available_columns(
        header.columns, ALTERNATE_ALLELE_CANDIDATES, "alternate allele"
    )
    validate_required_columns(header.columns)

    usecols = list(dict.fromkeys([*REQUIRED_COLUMNS, *reference_cols, *alternate_cols]))
    return usecols, reference_cols, alternate_cols


def validate_required_columns(columns: pd.Index) -> None:
    missing = [column for column in REQUIRED_COLUMNS if column not in columns]

    if missing:
        missing_str = ", ".join(sorted(set(missing)))
        raise KeyError(f"Missing required columns in input file: {missing_str}")


def load_variant_summary(
    input_path: Path,
    usecols: list[str],
    reference_cols: list[str],
    alternate_cols: list[str],
    chunksize: int = CHUNK_SIZE,
) -> tuple[pd.DataFrame, int, int]:
    usecols_set = set(usecols)
    filtered_chunks: list[pd.DataFrame] = []
    gene_filtered_total = 0
    retained_total = 0

    chunk_reader = pd.read_csv(
        input_path,
        sep="\t",
        compression="gzip",
        dtype=str,
        keep_default_na=False,
        low_memory=False,
        usecols=lambda column: column.strip() in usecols_set,
        chunksize=chunksize,
    )

    for chunk_number, chunk in enumerate(chunk_reader, start=1):
        chunk.columns = chunk.columns.str.strip()
        chunk = build_final_allele_columns(chunk, reference_cols, alternate_cols)

        gene_filtered_chunk = filter_target_genes(chunk)
        gene_filtered_total += len(gene_filtered_chunk)

        retained_chunk = filter_clinical_significance(gene_filtered_chunk)
        retained_total += len(retained_chunk)

        if not retained_chunk.empty:
            filtered_chunks.append(retained_chunk)

        print(
            f"Processed chunk {chunk_number}: retained {retained_total} rows cumulatively"
        )

    if filtered_chunks:
        filtered_df = pd.concat(filtered_chunks, ignore_index=True)
    else:
        filtered_df = pd.DataFrame(
            columns=[
                *usecols,
                "ReferenceAllele_final",
                "AlternateAllele_final",
                "_MatchedGenes",
            ]
        )

    return filtered_df, gene_filtered_total, retained_total


def normalize_missing_values(series: pd.Series) -> pd.Series:
    normalized = series.fillna("").astype(str).str.strip()
    return normalized.mask(normalized.str.lower().isin(MISSING_VALUE_TOKENS))


def coalesce_columns(df: pd.DataFrame, columns: list[str]) -> pd.Series:
    final_values = pd.Series(pd.NA, index=df.index, dtype="object")
    for column in columns:
        if column in df.columns:
            final_values = final_values.combine_first(normalize_missing_values(df[column]))
    return final_values.fillna("na")


def build_final_allele_columns(
    df: pd.DataFrame, reference_cols: list[str], alternate_cols: list[str]
) -> pd.DataFrame:
    final_df = df.copy()
    final_df["ReferenceAllele_final"] = coalesce_columns(final_df, reference_cols)
    final_df["AlternateAllele_final"] = coalesce_columns(final_df, alternate_cols)
    return final_df


def split_gene_symbols(gene_symbol: str) -> list[str]:
    if not gene_symbol:
        return []
    parts = [part.strip() for part in re.split(r"[;,/|]", gene_symbol) if part.strip()]
    return parts or [gene_symbol.strip()]


def extract_target_genes(gene_symbol: str) -> list[str]:
    matches: list[str] = []
    for gene in split_gene_symbols(gene_symbol):
        if gene in TARGET_GENES and gene not in matches:
            matches.append(gene)
    return matches


def filter_target_genes(df: pd.DataFrame) -> pd.DataFrame:
    matched_genes = df["GeneSymbol"].map(extract_target_genes)
    filtered = df.loc[matched_genes.map(bool)].copy()
    filtered["_MatchedGenes"] = matched_genes[matched_genes.map(bool)]
    return filtered


def filter_clinical_significance(df: pd.DataFrame) -> pd.DataFrame:
    significance = df["ClinicalSignificance"].fillna("").str.strip().str.lower()
    has_pathogenic_label = significance.str.contains("pathogenic", na=False)
    has_conflict = significance.str.contains("conflicting", na=False)
    return df.loc[has_pathogenic_label & ~has_conflict].copy()


def validate_no_conflicting_significance(df: pd.DataFrame) -> None:
    has_conflict = (
        df["ClinicalSignificance"]
        .fillna("")
        .str.strip()
        .str.contains("conflicting", case=False, na=False)
    )
    if has_conflict.any():
        raise ValueError(
            "Clinical significance filtering failed: conflicting records remain after filtering."
        )


def apply_assembly_filter(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series, bool]:
    assembly = df["Assembly"].str.strip()
    assembly_counts = assembly.replace("", "Missing").value_counts()
    grch38_mask = assembly.str.upper().eq("GRCH38")

    if grch38_mask.any():
        return df.loc[grch38_mask].copy(), assembly_counts, True
    return df.copy(), assembly_counts, False


def select_output_columns(df: pd.DataFrame) -> pd.DataFrame:
    output = df.loc[
        :,
        [
            "GeneSymbol",
            "ClinicalSignificance",
            "ReviewStatus",
            "Chromosome",
            "Start",
            "ReferenceAllele_final",
            "AlternateAllele_final",
            "Assembly",
        ],
    ].copy()
    return output.rename(
        columns={
            "ReferenceAllele_final": "ReferenceAllele",
            "AlternateAllele_final": "AlternateAllele",
        }
    )


def save_cleaned_dataset(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)


def gene_counts(df: pd.DataFrame) -> pd.Series:
    if df.empty:
        return pd.Series(0, index=list(TARGET_GENES), dtype="int64")

    counts = df["_MatchedGenes"].explode().value_counts()
    return counts.reindex(TARGET_GENES, fill_value=0)


def printable_value_counts(series: pd.Series) -> pd.Series:
    if series.empty:
        return series
    cleaned = series.replace("", "Missing")
    return cleaned.value_counts()


def allele_summary(cleaned_df: pd.DataFrame) -> tuple[int, int, int]:
    reference = normalize_missing_values(cleaned_df["ReferenceAllele"])
    alternate = normalize_missing_values(cleaned_df["AlternateAllele"])
    rows_with_missing_alleles = (reference.isna() | alternate.isna()).sum()
    return int(reference.notna().sum()), int(alternate.notna().sum()), int(
        rows_with_missing_alleles
    )


def print_counts(title: str, counts: pd.Series) -> None:
    print(f"\n{title}")
    if counts.empty:
        print("  None")
        return

    for label, count in counts.items():
        print(f"  {label}: {count}")


def print_summary(
    cleaned_df: pd.DataFrame,
    filtered_df: pd.DataFrame,
    assembly_counts: pd.Series,
    kept_only_grch38: bool,
    output_path: Path,
) -> None:
    if kept_only_grch38:
        assembly_message = "GRCh38 rows were available, so only GRCh38 variants were retained."
    else:
        assembly_message = "No GRCh38 rows were available after the other filters, so all assemblies were retained."

    print(f"Saved cleaned data to: {output_path}")
    print(f"Assembly handling: {assembly_message}")
    print(f"\nTotal variants: {len(cleaned_df)}")

    non_missing_reference, non_missing_alternate, rows_with_missing_alleles = (
        allele_summary(cleaned_df)
    )
    print(f"Rows with non-missing final ReferenceAllele: {non_missing_reference}")
    print(f"Rows with non-missing final AlternateAllele: {non_missing_alternate}")
    print(f"Rows with missing/na final allele values: {rows_with_missing_alleles}")

    print_counts("Assembly counts before final assembly selection", assembly_counts)
    print_counts("Variants per gene", gene_counts(filtered_df))
    print_counts(
        "Variants per ClinicalSignificance",
        printable_value_counts(cleaned_df["ClinicalSignificance"]),
    )
    print_counts(
        "Variants per ReviewStatus",
        printable_value_counts(cleaned_df["ReviewStatus"]),
    )


def main() -> None:
    args = parse_args()

    print(f"Loading ClinVar data from: {args.input}")
    usecols, reference_cols, alternate_cols = detect_input_columns(args.input)
    filtered, gene_filtered_rows, retained_rows = load_variant_summary(
        args.input, usecols, reference_cols, alternate_cols
    )
    print(f"Rows after gene filtering: {gene_filtered_rows}")
    print(f"Rows after pathogenic/non-conflicting filtering: {retained_rows}")
    validate_no_conflicting_significance(filtered)
    filtered, assembly_counts, kept_only_grch38 = apply_assembly_filter(filtered)
    print(f"Rows after assembly filtering: {len(filtered)}")

    cleaned = select_output_columns(filtered)
    cleaned = cleaned.reset_index(drop=True)
    save_cleaned_dataset(cleaned, args.output)
    print_summary(cleaned, filtered, assembly_counts, kept_only_grch38, args.output)


if __name__ == "__main__":
    main()
