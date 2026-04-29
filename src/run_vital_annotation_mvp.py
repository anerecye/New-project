from __future__ import annotations

import argparse
import gzip
from pathlib import Path

try:
    import pandas as pd
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires pandas. Install dependencies with: "
        "python -m pip install -r requirements.txt"
    ) from exc

from vital_standard import infer_alert_state, infer_certification, infer_public_use, infer_sv_required


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"

DEFAULT_SCORES = DATA_DIR / "arrhythmia_vital_scores.csv"
DEFAULT_RECONCILIATION = DATA_DIR / "vital_tiered_match_reconciliation_detail.csv"
DEFAULT_LOOKUP = DATA_DIR / "arrhythmia_vital_mvp_lookup.tsv"
DEFAULT_ANNOVAR = DATA_DIR / "arrhythmia_vital_mvp_annovar.tsv"

RECESSIVE_CONTEXT_GENES = {"CASQ2", "TRDN"}
TIER2_RECONCILIATION_TIERS = {
    "tier3_locus_observed_no_exact_af",
    "tier3_regional_context_only",
}
REVIEW_TRIGGER_LABEL = "POPULATION_REVIEW_TRIGGER_MAX_AF_GT_1E-5"
DOMINANT_MODEL_TRIGGER_LABEL = "UNQUALIFIED_DOMINANT_HIGH_PENETRANCE_GLOBAL_AF_GT_1E-4_AC_GE_20"
TIER2_TRIGGER_LABEL = "USABLE_ALLELE_RESOLVED_AF_REQUIRED"
UNEVALUABLE_TRIGGER_LABEL = "GRCH38_NORMALIZED_ALLELE_REQUIRED"


def truthy(value: object) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def numeric(value: object) -> float | None:
    if value is None or pd.isna(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def detect_sep(path: Path) -> str:
    return "\t" if path.suffix.lower() in {".tsv", ".tab", ".txt"} else ","


def is_vcf_path(path: Path) -> bool:
    suffixes = [suffix.lower() for suffix in path.suffixes]
    return suffixes[-1:] == [".vcf"] or suffixes[-2:] == [".vcf", ".gz"]


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix.lower() == ".gz" else path.open("r", encoding="utf-8")


def read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing input table: {path}")
    return pd.read_csv(path, sep=detect_sep(path), low_memory=False)


def read_vcf(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing input VCF: {path}")
    header: list[str] | None = None
    records: list[list[str]] = []
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.split("\t")
                continue
            if header is None:
                raise SystemExit("VCF header not found. Expected a #CHROM line.")
            values = line.split("\t")
            if len(values) < len(header):
                values.extend([""] * (len(header) - len(values)))
            records.append(values[: len(header)])
    if header is None:
        raise SystemExit("VCF header not found. Expected a #CHROM line.")
    return pd.DataFrame(records, columns=header, dtype="string")


def read_variant_input(path: Path) -> pd.DataFrame:
    return read_vcf(path) if is_vcf_path(path) else read_table(path)


def save_table(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep=detect_sep(path), index=False, lineterminator="\n")
    resolved = path.resolve()
    try:
        display_path = resolved.relative_to(BASE_DIR)
    except ValueError:
        display_path = resolved
    print(f"Saved {display_path} ({len(df)} rows)")


def normalize_chrom(value: object) -> str:
    chrom = str(value).strip()
    if not chrom:
        return chrom
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    chrom = chrom.upper()
    if chrom == "M":
        return "MT"
    return chrom


def trim_variant(pos: int, ref: object, alt: object) -> tuple[int, str, str]:
    ref_text = str(ref).strip().upper()
    alt_text = str(alt).strip().upper()
    current_pos = int(pos)

    while len(ref_text) > 1 and len(alt_text) > 1 and ref_text[-1] == alt_text[-1]:
        ref_text = ref_text[:-1]
        alt_text = alt_text[:-1]

    while len(ref_text) > 1 and len(alt_text) > 1 and ref_text[0] == alt_text[0]:
        ref_text = ref_text[1:]
        alt_text = alt_text[1:]
        current_pos += 1

    return current_pos, ref_text, alt_text


def build_vcf_lookup_key(chrom: object, pos: object, ref: object, alt: object) -> str:
    alt_text = str(alt).strip().upper()
    if "," in alt_text:
        raise SystemExit(
            "Multiallelic ALT values are not supported for VITAL lookup. "
            "Normalize and split with: bcftools norm -f GRCh38.fa -m -both ..."
        )
    normalized_pos, normalized_ref, normalized_alt = trim_variant(int(pos), ref, alt)
    return (
        f"{normalize_chrom(chrom)}:{normalized_pos}:{normalized_ref}:{normalized_alt}"
    )


def vcf_to_annovar_fields(pos: int, ref: str, alt: str) -> tuple[int, int, str, str]:
    if len(ref) == 1 and len(alt) == 1:
        return pos, pos, ref, alt

    if len(alt) > len(ref) and alt.startswith(ref):
        inserted = alt[len(ref):]
        return pos + len(ref) - 1, pos + len(ref) - 1, "-", inserted

    if len(ref) > len(alt) and ref.startswith(alt):
        deleted = ref[len(alt):]
        start = pos + len(alt)
        end = pos + len(ref) - 1
        return start, end, deleted, "-"

    return pos, pos + len(ref) - 1, ref, alt


def build_annovar_lookup_key(
    chrom: object,
    start: object,
    end: object,
    ref: object,
    alt: object,
) -> str:
    return f"{normalize_chrom(chrom)}:{int(start)}:{int(end)}:{str(ref).strip().upper()}:{str(alt).strip().upper()}"


def infer_evaluability(row: pd.Series) -> str:
    frequency_status = str(row.get("frequency_evidence_status", "")).strip()
    tier = str(row.get("reconciliation_tier", "")).strip()
    if frequency_status == "frequency_observed":
        return "TIER_1"
    if tier in {"tier4_still_unevaluable", "tier4_query_error"}:
        return "UNEVALUABLE"
    if tier in TIER2_RECONCILIATION_TIERS or frequency_status in {
        "exact_match_without_AF",
        "allele_discordance_no_exact_AF",
    }:
        return "TIER_2"
    return "UNEVALUABLE"


def infer_context(gene: object) -> str:
    return "RECESSIVE_CONTEXT" if str(gene).strip().upper() in RECESSIVE_CONTEXT_GENES else "STANDARD"


def infer_regime(row: pd.Series) -> str:
    if row["VITAL_evaluability"] != "TIER_1":
        return "."
    if not truthy(row.get("standard_acmg_frequency_flag")):
        return "OK"
    if row["VITAL_context"] == "RECESSIVE_CONTEXT":
        return "BOUNDARY"
    global_af = numeric(row.get("global_af")) or 0.0
    if truthy(row.get("frequency_signal_ac_ge_20")) and global_af > 1e-4:
        return "DOMINANT_INCOMPATIBLE"
    return "BOUNDARY"


def infer_flag(row: pd.Series) -> str:
    if row["VITAL_evaluability"] != "TIER_1":
        return "."
    if row["VITAL_regime"] == "OK":
        return "OK"
    if row["VITAL_regime"] == "DOMINANT_INCOMPATIBLE":
        return "MODEL_CONFLICT"
    return "CHECK_POPMAX"


def infer_action(row: pd.Series) -> str:
    if row["VITAL_evaluability"] == "TIER_2":
        return "REVIEW_NO_EXACT_MATCH"
    if row["VITAL_evaluability"] == "UNEVALUABLE":
        return "NORMALIZE_AND_REQUERY"
    if row["VITAL_flag"] == "OK":
        return "PASS"
    if row["VITAL_flag"] == "MODEL_CONFLICT":
        return "REQUIRE_MODEL_SPECIFICATION"
    if row["VITAL_context"] == "RECESSIVE_CONTEXT":
        return "REQUIRE_MODEL_SPECIFICATION"
    return "REVIEW_POPMAX"


def infer_reason(row: pd.Series) -> str:
    if row["VITAL_evaluability"] == "TIER_2":
        return "EVAL_LOCI_ONLY"
    if row["VITAL_evaluability"] == "UNEVALUABLE":
        return "EVAL_NO_ALLELE_MATCH"
    if row["VITAL_flag"] == "OK":
        return "EXPERT_REVIEW_RETAINED"
    if row["VITAL_flag"] == "MODEL_CONFLICT":
        return "MODEL_CONFLICT_DOMINANT_HP"
    if row["VITAL_context"] == "RECESSIVE_CONTEXT":
        return "RECESSIVE_COMPATIBLE_NOT_DOMINANT"
    if not truthy(row.get("frequency_signal_ac_ge_20")):
        return "CHECK_LOW_COUNT"
    return "POPMAX_EXCEEDS_MCAF"


def infer_threshold(row: pd.Series) -> str:
    if row["VITAL_evaluability"] == "TIER_2":
        return TIER2_TRIGGER_LABEL
    if row["VITAL_evaluability"] == "UNEVALUABLE":
        return UNEVALUABLE_TRIGGER_LABEL
    if row["VITAL_flag"] == "MODEL_CONFLICT":
        return DOMINANT_MODEL_TRIGGER_LABEL
    return REVIEW_TRIGGER_LABEL


def build_lookup_table(scores: pd.DataFrame, reconciliation: pd.DataFrame) -> pd.DataFrame:
    keep_scores = [
        "variant_key",
        "gene",
        "clinvar_id",
        "variation_id",
        "title",
        "review_status",
        "match_category",
        "frequency_evidence_status",
        "global_af",
        "global_ac",
        "popmax_af",
        "popmax_ac",
        "popmax_population",
        "qualifying_frequency_ac",
        "standard_acmg_frequency_flag",
        "popmax_only_frequency_flag",
        "frequency_signal_ac_ge_20",
        "vital_red_flag",
        "vital_score",
        "vital_band",
    ]
    merged = reconciliation.merge(
        scores.loc[:, [column for column in keep_scores if column in scores.columns]],
        on="variant_key",
        how="left",
        suffixes=("_recon", ""),
    ).copy()

    for column in ["chrom", "pos", "ref", "alt"]:
        if column not in merged.columns:
            raise ValueError(f"Reconciliation table is missing required column: {column}")

    merged["chrom"] = merged["chrom"].map(normalize_chrom)
    normalized_fields = merged.apply(
        lambda row: trim_variant(int(row["pos"]), row["ref"], row["alt"]),
        axis=1,
        result_type="expand",
    )
    normalized_fields.columns = ["normalized_pos", "normalized_ref", "normalized_alt"]
    merged[["normalized_pos", "normalized_ref", "normalized_alt"]] = normalized_fields
    merged["lookup_key_vcf"] = merged.apply(
        lambda row: build_vcf_lookup_key(
            row["chrom"],
            row["pos"],
            row["ref"],
            row["alt"],
        ),
        axis=1,
    )
    annovar_fields = merged.apply(
        lambda row: vcf_to_annovar_fields(
            int(row["normalized_pos"]),
            str(row["normalized_ref"]),
            str(row["normalized_alt"]),
        ),
        axis=1,
        result_type="expand",
    )
    annovar_fields.columns = ["annovar_start", "annovar_end", "annovar_ref", "annovar_alt"]
    merged[["annovar_start", "annovar_end", "annovar_ref", "annovar_alt"]] = annovar_fields
    merged["lookup_key_annovar"] = merged.apply(
        lambda row: build_annovar_lookup_key(
            row["chrom"],
            row["annovar_start"],
            row["annovar_end"],
            row["annovar_ref"],
            row["annovar_alt"],
        ),
        axis=1,
    )
    merged["VITAL_evaluability"] = merged.apply(infer_evaluability, axis=1)
    merged["VITAL_context"] = merged["gene"].map(infer_context)
    merged["VITAL_regime"] = merged.apply(infer_regime, axis=1)
    merged["VITAL_flag"] = merged.apply(infer_flag, axis=1)
    merged["VITAL_action"] = merged.apply(infer_action, axis=1)
    merged["VITAL_reason"] = merged.apply(infer_reason, axis=1)
    merged["VITAL_threshold"] = merged.apply(infer_threshold, axis=1)
    merged["VITAL_sv_required"] = merged.apply(infer_sv_required, axis=1)
    merged["VITAL_certification"] = merged.apply(infer_certification, axis=1)
    merged["VITAL_alert"] = merged.apply(infer_alert_state, axis=1)
    merged["VITAL_public_use"] = merged.apply(infer_public_use, axis=1)
    merged["VITAL_source"] = "ARRHYTHMIA_MVP_V1"

    output_columns = [
        "chrom",
        "normalized_pos",
        "normalized_ref",
        "normalized_alt",
        "annovar_start",
        "annovar_end",
        "annovar_ref",
        "annovar_alt",
        "lookup_key_vcf",
        "lookup_key_annovar",
        "variant_key",
        "clinvar_id",
        "variation_id",
        "gene",
        "title",
        "reconciliation_tier",
        "VITAL_evaluability",
        "VITAL_context",
        "VITAL_flag",
        "VITAL_regime",
        "VITAL_certification",
        "VITAL_alert",
        "VITAL_public_use",
        "VITAL_sv_required",
        "VITAL_action",
        "VITAL_reason",
        "VITAL_threshold",
        "VITAL_source",
        "global_af",
        "global_ac",
        "popmax_af",
        "popmax_ac",
        "popmax_population",
        "qualifying_frequency_ac",
        "frequency_evidence_status",
        "standard_acmg_frequency_flag",
        "popmax_only_frequency_flag",
        "frequency_signal_ac_ge_20",
        "vital_red_flag",
        "vital_score",
        "vital_band",
    ]
    lookup = merged.loc[:, output_columns].rename(
        columns={
            "chrom": "chr",
            "normalized_pos": "pos",
            "normalized_ref": "ref",
            "normalized_alt": "alt",
            "annovar_start": "Start",
            "annovar_end": "End",
            "annovar_ref": "Ref",
            "annovar_alt": "Alt",
            "global_af": "VITAL_global_af",
            "global_ac": "VITAL_global_ac",
            "popmax_af": "VITAL_popmax_af",
            "popmax_ac": "VITAL_popmax_ac",
            "popmax_population": "VITAL_popmax_population",
            "qualifying_frequency_ac": "VITAL_qualifying_frequency_ac",
            "frequency_evidence_status": "VITAL_frequency_evidence_status",
            "standard_acmg_frequency_flag": "VITAL_frequency_triggered",
            "popmax_only_frequency_flag": "VITAL_popmax_only_trigger",
            "frequency_signal_ac_ge_20": "VITAL_ac_ge_20",
            "vital_red_flag": "VITAL_red_flag",
            "vital_score": "VITAL_score",
            "vital_band": "VITAL_band",
        }
    )
    lookup = lookup.sort_values(["chr", "pos", "ref", "alt", "clinvar_id"]).reset_index(drop=True)
    if lookup["lookup_key_vcf"].duplicated().any():
        duplicates = lookup.loc[lookup["lookup_key_vcf"].duplicated(), "lookup_key_vcf"].tolist()
        raise ValueError(f"Lookup table contains duplicate normalized keys: {duplicates[:5]}")
    return lookup


def write_lookup_outputs(lookup: pd.DataFrame, lookup_path: Path, annovar_path: Path) -> None:
    save_table(lookup, lookup_path)
    annovar_columns = [
        "chr",
        "Start",
        "End",
        "Ref",
        "Alt",
        "VITAL_evaluability",
        "VITAL_flag",
        "VITAL_regime",
        "VITAL_certification",
        "VITAL_alert",
        "VITAL_public_use",
        "VITAL_sv_required",
        "VITAL_popmax_af",
        "VITAL_global_af",
        "VITAL_threshold",
        "VITAL_reason",
        "VITAL_action",
        "VITAL_context",
        "VITAL_qualifying_frequency_ac",
        "clinvar_id",
        "gene",
    ]
    annovar_view = lookup.loc[:, annovar_columns].rename(columns={"chr": "Chr"})
    save_table(annovar_view, annovar_path)


def infer_input_key_mode(df: pd.DataFrame) -> str:
    columns = set(df.columns)
    if {"Chr", "Start", "End", "Ref", "Alt"}.issubset(columns):
        return "annovar"
    if {"#CHROM", "POS", "REF", "ALT"}.issubset(columns):
        return "vcf"
    vcf_like_groups = [
        {"chr", "pos", "ref", "alt"},
        {"chrom", "pos", "ref", "alt"},
        {"Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF"},
    ]
    if any(group.issubset(columns) for group in vcf_like_groups):
        return "vcf"
    raise SystemExit(
        "Could not infer variant columns. Use either VCF-style columns "
        "(chr/pos/ref/alt, chrom/pos/ref/alt, #CHROM/POS/REF/ALT, or "
        "Chromosome/PositionVCF/ReferenceAlleleVCF/AlternateAlleleVCF) "
        "or ANNOVAR-style columns (Chr/Start/End/Ref/Alt)."
    )


def build_input_lookup_keys(df: pd.DataFrame, mode: str) -> pd.Series:
    if mode == "annovar":
        return df.apply(
            lambda row: build_annovar_lookup_key(
                row["Chr"],
                row["Start"],
                row["End"],
                row["Ref"],
                row["Alt"],
            ),
            axis=1,
        )

    if {"#CHROM", "POS", "REF", "ALT"}.issubset(df.columns):
        chrom_col, pos_col, ref_col, alt_col = "#CHROM", "POS", "REF", "ALT"
    elif {"chrom", "pos", "ref", "alt"}.issubset(df.columns):
        chrom_col, pos_col, ref_col, alt_col = "chrom", "pos", "ref", "alt"
    elif {"chr", "pos", "ref", "alt"}.issubset(df.columns):
        chrom_col, pos_col, ref_col, alt_col = "chr", "pos", "ref", "alt"
    else:
        chrom_col, pos_col, ref_col, alt_col = (
            "Chromosome",
            "PositionVCF",
            "ReferenceAlleleVCF",
            "AlternateAlleleVCF",
        )
    return df.apply(
        lambda row: build_vcf_lookup_key(
            row[chrom_col],
            row[pos_col],
            row[ref_col],
            row[alt_col],
        ),
        axis=1,
    )


def annovar_export_view(input_df: pd.DataFrame, mode: str) -> pd.DataFrame:
    table = input_df.copy()
    if mode == "annovar":
        table = table.rename(columns={"Chr": "chr", "Ref": "ref", "Alt": "alt"})
        table["pos"] = table["Start"]
        table = table.rename(columns={"Start": "Start", "End": "End"})
        return table

    if {"#CHROM", "POS", "REF", "ALT"}.issubset(table.columns):
        chrom_col, pos_col, ref_col, alt_col = "#CHROM", "POS", "REF", "ALT"
    elif {"chrom", "pos", "ref", "alt"}.issubset(table.columns):
        chrom_col, pos_col, ref_col, alt_col = "chrom", "pos", "ref", "alt"
    elif {"chr", "pos", "ref", "alt"}.issubset(table.columns):
        chrom_col, pos_col, ref_col, alt_col = "chr", "pos", "ref", "alt"
    else:
        chrom_col, pos_col, ref_col, alt_col = (
            "Chromosome",
            "PositionVCF",
            "ReferenceAlleleVCF",
            "AlternateAlleleVCF",
        )
    export = pd.DataFrame()
    export["chr"] = table[chrom_col].map(normalize_chrom)
    export["pos"] = pd.to_numeric(table[pos_col], errors="raise").astype(int)
    export["ref"] = table[ref_col].astype("string").str.upper()
    export["alt"] = table[alt_col].astype("string").str.upper()
    annovar_fields = export.apply(
        lambda row: vcf_to_annovar_fields(int(row["pos"]), str(row["ref"]), str(row["alt"])),
        axis=1,
        result_type="expand",
    )
    annovar_fields.columns = ["Start", "End", "Ref", "Alt"]
    export[["Start", "End", "Ref", "Alt"]] = annovar_fields
    return export


def maybe_fill_defaults(df: pd.DataFrame, defaults: dict[str, object]) -> None:
    for column, default in defaults.items():
        if column in df.columns:
            df[column] = df[column].fillna(default)


def annotate_table(
    input_path: Path,
    output_path: Path,
    lookup: pd.DataFrame,
    annovar_export_output: Path | None = None,
) -> None:
    input_df = read_variant_input(input_path)
    mode = infer_input_key_mode(input_df)
    input_df["VITAL_lookup_key"] = build_input_lookup_keys(input_df, mode)
    join_key = "lookup_key_annovar" if mode == "annovar" else "lookup_key_vcf"
    annotation_columns = [
        join_key,
        "VITAL_flag",
        "VITAL_evaluability",
        "VITAL_regime",
        "VITAL_action",
        "VITAL_context",
        "VITAL_reason",
        "VITAL_threshold",
        "clinvar_id",
        "gene",
        "VITAL_popmax_af",
        "VITAL_global_af",
        "VITAL_qualifying_frequency_ac",
    ]
    annotated = input_df.merge(
        lookup.loc[:, annotation_columns],
        left_on="VITAL_lookup_key",
        right_on=join_key,
        how="left",
    )
    annotated = annotated.drop(columns=[join_key])
    maybe_fill_defaults(
        annotated,
        {
            "VITAL_flag": ".",
            "VITAL_evaluability": "NOT_IN_VITAL_DB",
            "VITAL_regime": ".",
            "VITAL_action": ".",
            "VITAL_context": ".",
            "VITAL_reason": "not_in_vital_db",
            "VITAL_threshold": ".",
        },
    )
    save_table(annotated, output_path)

    if annovar_export_output is not None:
        annovar_view = annovar_export_view(input_df, mode)
        annovar_view["VITAL_lookup_key"] = annotated["VITAL_lookup_key"].values
        annovar_annotated = annovar_view.merge(
            annotated[
                [
                    "VITAL_lookup_key",
                    "VITAL_evaluability",
                    "VITAL_flag",
                    "VITAL_regime",
                    "VITAL_popmax_af",
                    "VITAL_global_af",
                    "VITAL_threshold",
                    "VITAL_reason",
                ]
            ],
            on="VITAL_lookup_key",
            how="left",
        )
        maybe_fill_defaults(
            annovar_annotated,
            {
                "VITAL_flag": ".",
                "VITAL_evaluability": "NOT_IN_VITAL_DB",
                "VITAL_regime": ".",
                "VITAL_reason": "not_in_vital_db",
                "VITAL_threshold": ".",
            },
        )
        annovar_annotated = annovar_annotated[
            [
                "chr",
                "Start",
                "End",
                "Ref",
                "Alt",
                "VITAL_evaluability",
                "VITAL_flag",
                "VITAL_regime",
                "VITAL_popmax_af",
                "VITAL_global_af",
                "VITAL_threshold",
                "VITAL_reason",
            ]
        ].rename(columns={"chr": "Chr"})
        save_table(annovar_annotated, annovar_export_output)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a compact VITAL MVP lookup table and optionally annotate an external "
            "CSV/TSV table with machine-readable VITAL columns."
        )
    )
    parser.add_argument("--scores", type=Path, default=DEFAULT_SCORES, help="VITAL score table")
    parser.add_argument(
        "--reconciliation",
        type=Path,
        default=DEFAULT_RECONCILIATION,
        help="Tiered reconciliation detail table",
    )
    parser.add_argument(
        "--lookup-out",
        type=Path,
        default=DEFAULT_LOOKUP,
        help="Output TSV for normalized chr/pos/ref/alt lookup",
    )
    parser.add_argument(
        "--annovar-out",
        type=Path,
        default=DEFAULT_ANNOVAR,
        help="Output TSV in ANNOVAR-style Chr/Start/End/Ref/Alt layout",
    )
    parser.add_argument(
        "--lookup",
        type=Path,
        help="Optional precomputed VITAL lookup TSV. If omitted, rebuild from --scores and --reconciliation.",
    )
    parser.add_argument("--input", type=Path, help="Optional VCF/CSV/TSV to annotate")
    parser.add_argument("--output", type=Path, help="Optional annotated CSV/TSV output path")
    parser.add_argument(
        "--annovar-export-output",
        type=Path,
        help="Optional annotated ANNOVAR-style export for the provided --input.",
    )
    parser.add_argument(
        "--skip-db-write",
        action="store_true",
        help="Build the lookup in memory without writing the lookup TSV outputs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.lookup is not None:
        lookup = read_table(args.lookup)
    else:
        scores = read_table(args.scores)
        reconciliation = read_table(args.reconciliation)
        lookup = build_lookup_table(scores, reconciliation)

    if args.lookup is None and not args.skip_db_write:
        write_lookup_outputs(lookup, args.lookup_out, args.annovar_out)

    if args.input is not None:
        if args.output is None:
            raise SystemExit("--output is required when --input is provided.")
        annotate_table(
            args.input,
            args.output,
            lookup,
            annovar_export_output=args.annovar_export_output,
        )


if __name__ == "__main__":
    main()
