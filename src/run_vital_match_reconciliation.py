from __future__ import annotations

import math
import time
from pathlib import Path

import pandas as pd
import requests

from arrhythmia_variant_pipeline import GNOMAD_DATASET, fetch_gnomad_region


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
SUPPLEMENT_DIR = BASE_DIR / "supplementary_tables"

MATCHED_IN = DATA_DIR / "arrhythmia_gnomad_matched.csv"
EXOME_GENOME_IN = DATA_DIR / "arrhythmia_exome_genome_af_comparison.csv"

DETAIL_OUT = DATA_DIR / "vital_tiered_match_reconciliation_detail.csv"
SUMMARY_OUT = DATA_DIR / "vital_tiered_match_reconciliation_summary.csv"
EVIDENCE_LAYERS_OUT = DATA_DIR / "vital_tiered_match_reconciliation_layers.csv"

DETAIL_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S47_tiered_match_reconciliation_detail.tsv"
SUMMARY_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S48_tiered_match_reconciliation_summary.tsv"
LAYERS_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S49_tiered_match_reconciliation_layers.tsv"

WINDOW_PAD = 2
REQUEST_PAUSE = 0.10
REQUEST_PAUSE_SECOND_PASS = 0.35


def normalize_variant(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    pos = int(pos)
    ref = str(ref)
    alt = str(alt)
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1
    return pos, ref, alt


def event_signature(pos: int, ref: str, alt: str) -> tuple[str, str, int] | None:
    pos, ref, alt = normalize_variant(pos, ref, alt)
    if len(ref) == len(alt):
        return None
    if len(ref) < len(alt):
        return ("ins", alt[1:], pos)
    return ("del", ref[1:], pos)


def decompose_substitution(pos: int, ref: str, alt: str) -> set[tuple[int, str, str]] | None:
    if len(ref) != len(alt) or len(ref) <= 1:
        return None
    components = set()
    for offset, (r, a) in enumerate(zip(ref, alt)):
        if r != a:
            components.add((int(pos) + offset, r, a))
    return components or None


def overlap_window(row: pd.Series) -> tuple[int, int]:
    pos = int(row["pos"])
    span = max(len(str(row["ref"])), len(str(row["alt"])))
    return max(1, pos - WINDOW_PAD), pos + span + WINDOW_PAD


def candidate_source(candidate: dict[str, object]) -> str:
    exome = candidate.get("exome") or {}
    genome = candidate.get("genome") or {}
    exome_present = exome.get("af") is not None or exome.get("ac") is not None
    genome_present = genome.get("af") is not None or genome.get("ac") is not None
    if exome_present and genome_present:
        return "exome_and_genome"
    if exome_present:
        return "exome"
    if genome_present:
        return "genome"
    return "locus_only"


def evaluate_unmatched(row: pd.Series, region_variants: list[dict[str, object]]) -> dict[str, object]:
    pos = int(row["pos"])
    ref = str(row["ref"])
    alt = str(row["alt"])
    npos, nref, nalt = normalize_variant(pos, ref, alt)
    signature = event_signature(pos, ref, alt)
    components = decompose_substitution(pos, ref, alt)
    exact_pos_variants = [v for v in region_variants if int(v.get("pos", -1)) == pos]

    for candidate in region_variants:
        cpos = int(candidate["pos"])
        cref = str(candidate["ref"])
        calt = str(candidate["alt"])
        if (cpos, cref, calt) == (pos, ref, alt):
            return {
                "reconciliation_tier": "tier1_exact_found_by_region",
                "candidate_variant_id": candidate.get("variantId", ""),
                "candidate_source": candidate_source(candidate),
                "reason": "Exact allele recovered in local region query despite original variant-id miss.",
            }

    for candidate in region_variants:
        cpos, cref, calt = normalize_variant(int(candidate["pos"]), str(candidate["ref"]), str(candidate["alt"]))
        if (cpos, cref, calt) == (npos, nref, nalt):
            return {
                "reconciliation_tier": "tier2_trim_normalized_equivalent",
                "candidate_variant_id": candidate.get("variantId", ""),
                "candidate_source": candidate_source(candidate),
                "reason": "Equivalent after trimming shared prefix/suffix without full left-alignment.",
            }

    if signature is not None:
        kind, payload, spos = signature
        for candidate in region_variants:
            csig = event_signature(int(candidate["pos"]), str(candidate["ref"]), str(candidate["alt"]))
            if csig is None:
                continue
            ckind, cpayload, cspos = csig
            if kind == ckind and payload == cpayload and abs(spos - cspos) <= WINDOW_PAD:
                return {
                    "reconciliation_tier": "tier2_event_equivalent_window",
                    "candidate_variant_id": candidate.get("variantId", ""),
                    "candidate_source": candidate_source(candidate),
                    "reason": "Equivalent indel payload observed in local +/-2 bp window.",
                }

    if components is not None:
        observed_components = set()
        for candidate in region_variants:
            if len(str(candidate["ref"])) == 1 and len(str(candidate["alt"])) == 1:
                observed_components.add((int(candidate["pos"]), str(candidate["ref"]), str(candidate["alt"])))
        if components.issubset(observed_components):
            return {
                "reconciliation_tier": "tier2_decomposed_equivalent_unphased",
                "candidate_variant_id": ";".join(
                    [
                        str(v.get("variantId", ""))
                        for v in region_variants
                        if len(str(v["ref"])) == 1 and len(str(v["alt"])) == 1
                        and (int(v["pos"]), str(v["ref"]), str(v["alt"])) in components
                    ]
                ),
                "candidate_source": "component_snv_set",
                "reason": "All decomposed substitution components are present, but phase/haplotype is unconfirmed.",
            }

    if exact_pos_variants:
        return {
            "reconciliation_tier": "tier3_locus_observed_no_exact_af",
            "candidate_variant_id": ";".join(str(v.get("variantId", "")) for v in exact_pos_variants[:5]),
            "candidate_source": "locus_only",
            "reason": "Same genomic position observed in gnomAD, but not the queried allele.",
        }

    if region_variants:
        return {
            "reconciliation_tier": "tier3_regional_context_only",
            "candidate_variant_id": ";".join(str(v.get("variantId", "")) for v in region_variants[:5]),
            "candidate_source": "regional_only",
            "reason": "Nearby population variation observed within +/-2 bp window, but no exact or equivalent allele recovered.",
        }

    return {
        "reconciliation_tier": "tier4_still_unevaluable",
        "candidate_variant_id": "",
        "candidate_source": "none",
        "reason": "No exact, equivalent, or locus-level population evidence recovered.",
    }


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    matched = pd.read_csv(MATCHED_IN)
    exg = pd.read_csv(EXOME_GENOME_IN)
    for df in [matched, exg]:
        for col in ["pos", "gnomad_af", "gnomad_ac", "gnomad_an", "genome_af", "genome_ac", "genome_an"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
    return matched, exg


def base_detail(matched: pd.DataFrame, exg: pd.DataFrame) -> pd.DataFrame:
    df = matched.copy()
    df = df.merge(
        exg[["variant_key", "genome_af", "genome_ac", "genome_an"]],
        on="variant_key",
        how="left",
    )
    df["reconciliation_tier"] = ""
    df["candidate_variant_id"] = ""
    df["candidate_source"] = ""
    df["reason"] = ""
    exact_mask = df["match_category"].eq("exact_match")
    exact_exome_mask = exact_mask & df["gnomad_status"].eq("matched")
    exact_no_exome_mask = exact_mask & df["gnomad_status"].eq("no_exome_data")
    df.loc[exact_exome_mask, "reconciliation_tier"] = "tier1_exact_exome"
    df.loc[exact_exome_mask, "candidate_source"] = "exome"
    df.loc[exact_exome_mask, "reason"] = "Original chr-pos-ref-alt exact exome match."
    df.loc[exact_no_exome_mask, "reconciliation_tier"] = "tier1_exact_genome_or_no_exome"
    df.loc[exact_no_exome_mask, "candidate_source"] = "genome_or_missing_exome"
    df.loc[exact_no_exome_mask, "reason"] = (
        "Original chr-pos-ref-alt exact match without exome block; genome query retained separately."
    )
    if df["variant_key"].duplicated().any():
        dup_keys = df.loc[df["variant_key"].duplicated(), "variant_key"].tolist()
        raise ValueError(f"Matched analytical cohort contains duplicate variant_key values: {dup_keys[:5]}")
    return df


def reconcile_unmatched(detail: pd.DataFrame) -> pd.DataFrame:
    unresolved = detail["reconciliation_tier"].eq("")
    if not unresolved.any():
        return detail
    windows = {}
    for _, row in detail.loc[unresolved].iterrows():
        start, stop = overlap_window(row)
        windows[(str(row["chrom"]), start, stop)] = None

    with requests.Session() as session:
        session.headers.update({"User-Agent": "vital-tiered-reconciliation/1.0"})
        for key in windows:
            chrom, start, stop = key
            windows[key] = fetch_gnomad_region(session, chrom=chrom, start=start, stop=stop, dataset=GNOMAD_DATASET)
            time.sleep(REQUEST_PAUSE)
        failed_windows = [key for key, value in windows.items() if value is None]
        if failed_windows:
            for key in failed_windows:
                chrom, start, stop = key
                windows[key] = fetch_gnomad_region(
                    session,
                    chrom=chrom,
                    start=start,
                    stop=stop,
                    dataset=GNOMAD_DATASET,
                )
                time.sleep(REQUEST_PAUSE_SECOND_PASS)

    tiers = []
    candidate_ids = []
    candidate_sources = []
    reasons = []
    window_variant_counts = []
    for _, row in detail.loc[unresolved].iterrows():
        start, stop = overlap_window(row)
        region_variants = windows[(str(row["chrom"]), start, stop)]
        if region_variants is None:
            result = {
                "reconciliation_tier": "tier4_query_error",
                "candidate_variant_id": "",
                "candidate_source": "query_error",
                "reason": "Region query failed.",
            }
            window_variant_counts.append(math.nan)
        else:
            result = evaluate_unmatched(row, region_variants)
            window_variant_counts.append(len(region_variants))
        tiers.append(result["reconciliation_tier"])
        candidate_ids.append(result["candidate_variant_id"])
        candidate_sources.append(result["candidate_source"])
        reasons.append(result["reason"])

    detail.loc[unresolved, "reconciliation_tier"] = tiers
    detail.loc[unresolved, "candidate_variant_id"] = candidate_ids
    detail.loc[unresolved, "candidate_source"] = candidate_sources
    detail.loc[unresolved, "reason"] = reasons
    detail.loc[unresolved, "window_variants_plusminus2bp"] = window_variant_counts
    return detail


def build_summary(detail: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    total = len(detail)
    rows = []
    tier_order = [
        "tier1_exact_exome",
        "tier1_exact_genome_or_no_exome",
        "tier1_exact_found_by_region",
        "tier2_trim_normalized_equivalent",
        "tier2_event_equivalent_window",
        "tier2_decomposed_equivalent_unphased",
        "tier3_locus_observed_no_exact_af",
        "tier3_regional_context_only",
        "tier4_still_unevaluable",
        "tier4_query_error",
    ]
    for tier in tier_order:
        sub = detail[detail["reconciliation_tier"].eq(tier)]
        rows.append(
            {
                "reconciliation_tier": tier,
                "variants": len(sub),
                "percent_of_total": 100 * len(sub) / total if total else math.nan,
                "interpretation": tier_interpretation(tier),
            }
        )
    summary = pd.DataFrame(rows)

    layer_rows = [
        {
            "evidence_layer": "original_exact_exome",
            "variants": int(detail["reconciliation_tier"].eq("tier1_exact_exome").sum()),
            "percent_of_total": 100 * detail["reconciliation_tier"].eq("tier1_exact_exome").sum() / total,
            "meaning": "Original strict exome AF-observed space used in the manuscript core analysis.",
        },
        {
            "evidence_layer": "original_exact_any_dataset",
            "variants": int(
                detail["reconciliation_tier"]
                .isin(["tier1_exact_exome", "tier1_exact_genome_or_no_exome", "tier1_exact_found_by_region"])
                .sum()
            ),
            "percent_of_total": 100
            * detail["reconciliation_tier"]
            .isin(["tier1_exact_exome", "tier1_exact_genome_or_no_exome", "tier1_exact_found_by_region"])
            .sum()
            / total,
            "meaning": "Original strict exact-match space including genome-only/no-exome exact records.",
        },
        {
            "evidence_layer": "reconciled_equivalent",
            "variants": int(
                detail["reconciliation_tier"]
                .isin(
                    [
                        "tier1_exact_exome",
                        "tier1_exact_genome_or_no_exome",
                        "tier1_exact_found_by_region",
                        "tier2_trim_normalized_equivalent",
                        "tier2_event_equivalent_window",
                        "tier2_decomposed_equivalent_unphased",
                    ]
                )
                .sum()
            ),
            "percent_of_total": 100
            * detail["reconciliation_tier"]
            .isin(
                [
                    "tier1_exact_exome",
                    "tier1_exact_genome_or_no_exome",
                    "tier1_exact_found_by_region",
                    "tier2_trim_normalized_equivalent",
                    "tier2_event_equivalent_window",
                    "tier2_decomposed_equivalent_unphased",
                ]
            )
            .sum()
            / total,
            "meaning": "Strict exact plus equivalent-allele reconciliation tiers; should not be mixed with locus-only evidence.",
        },
        {
            "evidence_layer": "locus_or_regional_observed",
            "variants": int(
                detail["reconciliation_tier"]
                .isin(["tier3_locus_observed_no_exact_af", "tier3_regional_context_only"])
                .sum()
            ),
            "percent_of_total": 100
            * detail["reconciliation_tier"]
            .isin(["tier3_locus_observed_no_exact_af", "tier3_regional_context_only"])
            .sum()
            / total,
            "meaning": "Region observed, but no exact/equivalent allele recovered; not AF evidence.",
        },
        {
            "evidence_layer": "still_unevaluable_after_reconciliation",
            "variants": int(
                detail["reconciliation_tier"].isin(["tier4_still_unevaluable", "tier4_query_error"]).sum()
            ),
            "percent_of_total": 100
            * detail["reconciliation_tier"].isin(["tier4_still_unevaluable", "tier4_query_error"]).sum()
            / total,
            "meaning": "Substantial representation gap remains even after aggressive local-window reconciliation.",
        },
    ]
    return summary, pd.DataFrame(layer_rows)


def tier_interpretation(tier: str) -> str:
    mapping = {
        "tier1_exact_exome": "Gold-standard exact exome match.",
        "tier1_exact_genome_or_no_exome": "Gold-standard exact match without exome AF block; genome evidence separate.",
        "tier1_exact_found_by_region": "Exact allele recovered by local region query despite original miss.",
        "tier2_trim_normalized_equivalent": "Equivalent after trim-normalization only; not full left-alignment.",
        "tier2_event_equivalent_window": "Equivalent indel payload in local window; moderate-confidence equivalent evidence.",
        "tier2_decomposed_equivalent_unphased": "All MNV/substitution components observed, but phase is unconfirmed.",
        "tier3_locus_observed_no_exact_af": "Same position observed with different allele; locus-level evidence only.",
        "tier3_regional_context_only": "Nearby regional variation present, but no same-locus or equivalent allele.",
        "tier4_still_unevaluable": "No exact, equivalent, or locus-level evidence recovered.",
        "tier4_query_error": "Reconciliation query failed.",
    }
    return mapping[tier]


def save(df: pd.DataFrame, path: Path, sep: str = ",") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep=sep, index=False, lineterminator="\n")
    print(f"Saved {path.relative_to(BASE_DIR)} ({len(df)} rows)")


def main() -> None:
    matched, exg = load_inputs()
    detail = base_detail(matched, exg)
    detail = reconcile_unmatched(detail)
    if len(detail) != len(matched):
        raise ValueError(
            f"Reconciled detail row count {len(detail)} does not match analytical cohort size {len(matched)}."
        )
    summary, layers = build_summary(detail)

    save(detail, DETAIL_OUT)
    save(summary, SUMMARY_OUT)
    save(layers, EVIDENCE_LAYERS_OUT)

    save(detail, DETAIL_SUPP, sep="\t")
    save(summary, SUMMARY_SUPP, sep="\t")
    save(layers, LAYERS_SUPP, sep="\t")


if __name__ == "__main__":
    main()
