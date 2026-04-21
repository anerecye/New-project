from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires pandas and numpy. Install them with: "
        "python -m pip install -r requirements.txt"
    ) from exc

sys.path.insert(0, str(Path(__file__).resolve().parent))

import advanced_variant_analyses as ava
from validate_vital_reclassification import ARRHYTHMIA_GENES, clinical_group


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"


def data_path(prefix: str, name: str) -> Path:
    return DATA_DIR / ava.prefixed_name(prefix, name)


def save_table(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = output_path.with_name(f".{output_path.name}.tmp")
    last_error: PermissionError | None = None
    for attempt in range(1, 7):
        try:
            df.to_csv(temp_path, index=False)
            temp_path.replace(output_path)
            print(f"Saved {output_path} ({len(df)} rows)")
            return
        except PermissionError as exc:
            last_error = exc
            delay = 1.5 * attempt
            print(f"Save retry {attempt}/6 for {output_path} in {delay:.1f}s: {exc}")
            time.sleep(delay)
    if last_error is not None:
        raise last_error


def make_vcv_accession(variation_id: object) -> str:
    try:
        return f"VCV{int(str(variation_id)):09d}"
    except ValueError:
        return f"VCV{variation_id}"


def collapse_first_nonempty(values: pd.Series) -> str:
    cleaned = [str(value) for value in values.dropna() if str(value) and str(value) != "-"]
    return cleaned[0] if cleaned else ""


def collapse_unique(values: pd.Series) -> str:
    cleaned = sorted({str(value) for value in values.dropna() if str(value) and str(value) != "-"})
    return "|".join(cleaned)


def load_plp_snapshot(
    variant_summary_path: Path,
    genes: list[str] | None,
    max_variants: int | None = None,
    sample_variants: int | None = None,
    random_seed: int = 13,
) -> pd.DataFrame:
    usecols = [
        "Type",
        "Name",
        "GeneSymbol",
        "ClinicalSignificance",
        "ReviewStatus",
        "NumberSubmitters",
        "Assembly",
        "Chromosome",
        "PositionVCF",
        "ReferenceAlleleVCF",
        "AlternateAlleleVCF",
        "VariationID",
    ]
    chunks: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        variant_summary_path,
        sep="\t",
        compression="infer",
        usecols=lambda column: column in usecols,
        dtype=str,
        chunksize=250_000,
        low_memory=False,
    ):
        if genes is not None:
            chunk = chunk[chunk["GeneSymbol"].isin(genes)].copy()
        else:
            chunk = chunk.copy()
        if "Assembly" in chunk.columns:
            chunk = chunk[chunk["Assembly"].eq("GRCh38")]
        if chunk.empty:
            continue
        chunk["clinical_group"] = chunk["ClinicalSignificance"].map(clinical_group)
        chunk = chunk[chunk["clinical_group"].eq("P_LP")].copy()
        if chunk.empty:
            continue
        chunks.append(chunk)
    if not chunks:
        return pd.DataFrame()

    df = pd.concat(chunks, ignore_index=True)
    for column in ["Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF", "VariationID"]:
        df[column] = df[column].fillna("").astype(str)
    df = df[
        df["Chromosome"].ne("")
        & df["PositionVCF"].ne("")
        & df["ReferenceAlleleVCF"].ne("")
        & df["AlternateAlleleVCF"].ne("")
        & ~df["ReferenceAlleleVCF"].isin({"na", "-"})
        & ~df["AlternateAlleleVCF"].isin({"na", "-"})
    ].copy()
    df["variant_key"] = (
        df["Chromosome"]
        + ":"
        + df["PositionVCF"]
        + ":"
        + df["ReferenceAlleleVCF"]
        + ":"
        + df["AlternateAlleleVCF"]
    )
    grouped = (
        df.groupby(["VariationID", "variant_key"], dropna=False)
        .agg(
            gene=("GeneSymbol", collapse_first_nonempty),
            chrom=("Chromosome", collapse_first_nonempty),
            pos=("PositionVCF", collapse_first_nonempty),
            ref=("ReferenceAlleleVCF", collapse_first_nonempty),
            alt=("AlternateAlleleVCF", collapse_first_nonempty),
            clinsig=("ClinicalSignificance", collapse_unique),
            review_status=("ReviewStatus", collapse_unique),
            title=("Name", collapse_first_nonempty),
            variant_type=("Type", collapse_first_nonempty),
            submitter_count=("NumberSubmitters", "max"),
        )
        .reset_index()
        .rename(columns={"VariationID": "variation_id"})
    )
    grouped["clinvar_id"] = grouped["variation_id"].map(make_vcv_accession)
    grouped["gnomad_af"] = np.nan
    grouped["gnomad_ac"] = np.nan
    grouped["gnomad_an"] = np.nan
    grouped["match_category"] = "exact_match"
    grouped["gnomad_status"] = ""
    grouped["gnomad_variants_at_pos_for_match"] = ""
    grouped["gnomad_region_query_status"] = ""
    grouped["freq_class"] = ""
    grouped["submitter_count"] = pd.to_numeric(grouped["submitter_count"], errors="coerce")
    grouped = grouped.sort_values(["gene", "variation_id", "variant_key"]).reset_index(drop=True)
    if sample_variants is not None and len(grouped) > sample_variants:
        grouped = (
            grouped.sample(n=sample_variants, random_state=random_seed)
            .sort_values(["gene", "variation_id", "variant_key"])
            .reset_index(drop=True)
        )
    elif max_variants is not None:
        grouped = grouped.head(max_variants).copy()
    return grouped


def attach_gnomad_exact_matches(
    annotated: pd.DataFrame,
    output_prefix: str,
    dataset: str,
    pause: float,
    fetch: bool,
    force: bool,
) -> pd.DataFrame:
    exome_genome = ava.load_or_fetch_exome_genome_af(
        annotated,
        data_path(output_prefix, "exome_genome_af_comparison.csv"),
        dataset=dataset,
        pause=pause,
        fetch=fetch,
        force=force,
    )
    if exome_genome.empty:
        result = annotated.copy()
        result["match_category"] = "no_gnomad_record"
        return result
    keep = [
        "variant_key",
        "exome_af",
        "exome_ac",
        "exome_an",
        "exome_genome_query_status",
    ]
    merged = annotated.merge(
        exome_genome.loc[:, [column for column in keep if column in exome_genome.columns]],
        on="variant_key",
        how="left",
    )
    merged["gnomad_af"] = pd.to_numeric(merged.get("exome_af"), errors="coerce")
    merged["gnomad_ac"] = pd.to_numeric(merged.get("exome_ac"), errors="coerce")
    merged["gnomad_an"] = pd.to_numeric(merged.get("exome_an"), errors="coerce")
    query_status = merged.get(
        "exome_genome_query_status",
        pd.Series("", index=merged.index, dtype="object"),
    ).fillna("").astype(str)
    variant_not_found = query_status.str.contains("Variant not found", na=False)
    merged["match_category"] = np.select(
        [
            merged["gnomad_af"].notna(),
            variant_not_found,
            query_status.str.startswith("error:"),
        ],
        [
            "exact_match",
            "no_gnomad_record",
            "gnomad_query_error",
        ],
        default="no_gnomad_record",
    )
    merged["gnomad_status"] = query_status
    return merged


def score_snapshot(
    variant_summary_path: Path,
    output_prefix: str,
    genes: list[str] | None,
    dataset: str,
    pause: float,
    fetch_gnomad: bool,
    force_gnomad: bool,
    max_variants: int | None,
    sample_variants: int | None = None,
    random_seed: int = 13,
) -> dict[str, Path]:
    snapshot = load_plp_snapshot(
        variant_summary_path,
        genes=genes,
        max_variants=max_variants,
        sample_variants=sample_variants,
        random_seed=random_seed,
    )
    if snapshot.empty:
        raise ValueError(f"No P/LP GRCh38 variants found in {variant_summary_path}.")
    save_table(snapshot, data_path(output_prefix, "clinvar_variants.csv"))
    annotated = ava.prepare_annotations(snapshot)
    annotated = attach_gnomad_exact_matches(
        annotated,
        output_prefix=output_prefix,
        dataset=dataset,
        pause=pause,
        fetch=fetch_gnomad,
        force=force_gnomad,
    )
    annotated = ava.prepare_annotations(annotated)
    save_table(annotated, data_path(output_prefix, "gnomad_matched.csv"))

    population_df = ava.load_or_fetch_population_af(
        annotated,
        data_path(output_prefix, "population_af.csv"),
        dataset=dataset,
        pause=pause,
        fetch=fetch_gnomad,
        force=force_gnomad,
    )
    vital_scores, vital_summary, vital_predictions, gene_constraint, type_detectability = ava.make_vital_score_tables(
        annotated,
        population_df,
    )
    component_breakdown, component_summary = ava.make_vital_component_breakdown(vital_scores)
    method_comparison, validation_curves, retention_summary = ava.make_vital_benchmark_tables(vital_scores)
    threshold_sweep = ava.make_vital_threshold_sweep(vital_scores)
    ac_threshold_sensitivity = ava.make_vital_ac_threshold_sensitivity(vital_scores)
    acmg_disagreement = ava.make_vital_acmg_disagreement_table(vital_scores)
    top_suspicious = ava.make_vital_top_suspicious_table(vital_scores)
    absence_detectability_bias = ava.make_absence_detectability_bias_table(vital_scores)
    review_fragility = ava.make_review_fragility_summary(vital_scores)

    outputs = {
        "vital_scores": data_path(output_prefix, "vital_scores.csv"),
        "vital_summary": data_path(output_prefix, "vital_summary.csv"),
        "vital_predictions": data_path(output_prefix, "vital_predictions.csv"),
        "gene_frequency_constraint": data_path(output_prefix, "gene_frequency_constraint.csv"),
        "variant_type_detectability": data_path(output_prefix, "variant_type_detectability.csv"),
        "vital_component_breakdown": data_path(output_prefix, "vital_component_breakdown.csv"),
        "vital_component_summary": data_path(output_prefix, "vital_component_summary.csv"),
        "vital_method_comparison": data_path(output_prefix, "vital_method_comparison.csv"),
        "vital_validation_curves": data_path(output_prefix, "vital_validation_curves.csv"),
        "vital_retention_summary": data_path(output_prefix, "vital_retention_summary.csv"),
        "vital_threshold_sweep": data_path(output_prefix, "vital_threshold_sweep.csv"),
        "vital_ac_threshold_sensitivity": data_path(output_prefix, "vital_ac_threshold_sensitivity.csv"),
        "vital_acmg_disagreement": data_path(output_prefix, "vital_acmg_disagreement.csv"),
        "vital_top_suspicious": data_path(output_prefix, "vital_top_suspicious.csv"),
        "vital_absence_detectability_bias": data_path(output_prefix, "vital_absence_detectability_bias.csv"),
        "review_fragility_summary": data_path(output_prefix, "review_fragility_summary.csv"),
    }
    save_table(vital_scores, outputs["vital_scores"])
    save_table(vital_summary, outputs["vital_summary"])
    save_table(vital_predictions, outputs["vital_predictions"])
    save_table(gene_constraint, outputs["gene_frequency_constraint"])
    save_table(type_detectability, outputs["variant_type_detectability"])
    save_table(component_breakdown, outputs["vital_component_breakdown"])
    save_table(component_summary, outputs["vital_component_summary"])
    save_table(method_comparison, outputs["vital_method_comparison"])
    save_table(validation_curves, outputs["vital_validation_curves"])
    save_table(retention_summary, outputs["vital_retention_summary"])
    save_table(threshold_sweep, outputs["vital_threshold_sweep"])
    save_table(ac_threshold_sensitivity, outputs["vital_ac_threshold_sensitivity"])
    save_table(acmg_disagreement, outputs["vital_acmg_disagreement"])
    save_table(top_suspicious, outputs["vital_top_suspicious"])
    save_table(absence_detectability_bias, outputs["vital_absence_detectability_bias"])
    save_table(review_fragility, outputs["review_fragility_summary"])
    return outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate VITAL scores from an archived ClinVar variant_summary snapshot."
    )
    parser.add_argument("--variant-summary", type=Path, required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--genes", nargs="+", default=ARRHYTHMIA_GENES)
    parser.add_argument(
        "--all-genes",
        action="store_true",
        help="Sample from all ClinVar genes instead of filtering to --genes.",
    )
    parser.add_argument("--dataset", default=ava.GNOMAD_DATASET)
    parser.add_argument("--gnomad-pause", type=float, default=0.25)
    parser.add_argument("--no-fetch-gnomad", action="store_true")
    parser.add_argument("--force-gnomad-fetch", action="store_true")
    parser.add_argument("--max-variants", type=int)
    parser.add_argument(
        "--sample-variants",
        type=int,
        help="Randomly sample this many collapsed P/LP variants after filtering.",
    )
    parser.add_argument("--random-seed", type=int, default=13)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = score_snapshot(
        variant_summary_path=args.variant_summary,
        output_prefix=args.output_prefix,
        genes=None if args.all_genes else args.genes,
        dataset=args.dataset,
        pause=args.gnomad_pause,
        fetch_gnomad=not args.no_fetch_gnomad,
        force_gnomad=args.force_gnomad_fetch,
        max_variants=args.max_variants,
        sample_variants=args.sample_variants,
        random_seed=args.random_seed,
    )
    print("VITAL snapshot scoring outputs")
    for name, path in outputs.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
