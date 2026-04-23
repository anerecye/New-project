from __future__ import annotations

from pathlib import Path

import pandas as pd


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
SUPPLEMENT_DIR = BASE_DIR / "supplementary_tables"

RECONCILIATION_IN = DATA_DIR / "vital_tiered_match_reconciliation_detail.csv"
SCORES_IN = DATA_DIR / "arrhythmia_vital_scores.csv"

TIER2_SUBCATEGORY_OUT = DATA_DIR / "vital_tier2_subcategory_summary.csv"
TIER1_TIER2_CLASS_OUT = DATA_DIR / "vital_tier2_vs_tier1_variant_class_summary.csv"

TIER2_SUBCATEGORY_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S54_tier2_subcategory_summary.tsv"
TIER1_TIER2_CLASS_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S55_tier2_vs_tier1_variant_class_summary.tsv"

TIER1_EQUIVALENT = {
    "tier1_exact_exome",
    "tier1_exact_genome_or_no_exome",
    "tier1_exact_found_by_region",
    "tier2_trim_normalized_equivalent",
    "tier2_event_equivalent_window",
    "tier2_decomposed_equivalent_unphased",
}
TIER2_CONTEXT = {
    "tier3_locus_observed_no_exact_af",
    "tier3_regional_context_only",
}
INDEL_TYPES = {"deletion", "duplication", "insertion"}


def save_table(df: pd.DataFrame, path: Path, sep: str = ",") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep=sep, index=False, lineterminator="\n")
    print(f"Saved {path.relative_to(BASE_DIR)} ({len(df)} rows)")


def load_detail() -> pd.DataFrame:
    detail = pd.read_csv(RECONCILIATION_IN, low_memory=False)
    scores = pd.read_csv(SCORES_IN, usecols=["clinvar_id", "variant_type"])
    detail = detail.merge(scores, on="clinvar_id", how="left", validate="many_to_one")
    if detail["variant_type"].isna().any():
        missing = int(detail["variant_type"].isna().sum())
        raise ValueError(f"Missing variant_type for {missing} reconciliation rows")
    return detail


def class_counts(sub: pd.DataFrame) -> dict[str, float | int]:
    n = len(sub)
    indel_n = int(sub["variant_type"].isin(INDEL_TYPES).sum())
    dup_n = int(sub["variant_type"].eq("duplication").sum())
    snv_n = int(sub["variant_type"].eq("SNV").sum())
    return {
        "n_variants": n,
        "indel_dup_ins_n": indel_n,
        "indel_dup_ins_percent": 100 * indel_n / n if n else 0.0,
        "duplication_n": dup_n,
        "duplication_percent": 100 * dup_n / n if n else 0.0,
        "snv_n": snv_n,
        "snv_percent": 100 * snv_n / n if n else 0.0,
    }


def build_tier2_subcategory_summary(detail: pd.DataFrame) -> pd.DataFrame:
    tier2 = detail.loc[detail["reconciliation_tier"].isin(TIER2_CONTEXT)].copy()
    total = len(tier2)
    rows: list[dict[str, object]] = []
    specs = [
        (
            "allele_discordance_at_locus",
            "tier3_locus_observed_no_exact_af",
            "Same locus observed in gnomAD with a different allele; locus-level context is present, but no allele-resolved AF exists for the queried variant.",
        ),
        (
            "no_same_locus_record_regional_only",
            "tier3_regional_context_only",
            "No same-locus or equivalent allele recovered; only nearby regional variation is visible in the local window.",
        ),
    ]
    for label, tier_name, interpretation in specs:
        sub = tier2.loc[tier2["reconciliation_tier"].eq(tier_name)].copy()
        row = {
            "tier2_subcategory": label,
            "reconciliation_tier": tier_name,
            "percent_of_tier2": 100 * len(sub) / total if total else 0.0,
            "interpretation": interpretation,
        }
        row.update(class_counts(sub))
        rows.append(row)
    return pd.DataFrame(rows)


def build_tier1_tier2_class_summary(detail: pd.DataFrame) -> pd.DataFrame:
    tier1 = detail.loc[detail["reconciliation_tier"].isin(TIER1_EQUIVALENT)].copy()
    tier2 = detail.loc[detail["reconciliation_tier"].isin(TIER2_CONTEXT)].copy()

    tier1_indel = int(tier1["variant_type"].isin(INDEL_TYPES).sum())
    tier1_non = len(tier1) - tier1_indel
    tier2_indel = int(tier2["variant_type"].isin(INDEL_TYPES).sum())
    tier2_non = len(tier2) - tier2_indel
    indel_or = (tier2_indel * tier1_non) / (tier2_non * tier1_indel)

    tier1_dup = int(tier1["variant_type"].eq("duplication").sum())
    tier1_nondup = len(tier1) - tier1_dup
    tier2_dup = int(tier2["variant_type"].eq("duplication").sum())
    tier2_nondup = len(tier2) - tier2_dup
    dup_or = (tier2_dup * tier1_nondup) / (tier2_nondup * tier1_dup)

    rows = []
    for label, sub, interpretation in [
        (
            "tier1_allele_resolved_exact_or_equivalent",
            tier1,
            "Allele-resolved exact or equivalent population context suitable for direct AF interpretation.",
        ),
        (
            "tier2_locus_or_regional_context_only",
            tier2,
            "Structured context without allele-resolved AF; includes same-locus discordance and regional-only observations.",
        ),
    ]:
        row = {
            "comparison_group": label,
            "interpretation": interpretation,
            "indel_dup_ins_or_vs_tier1": 1.0 if label.startswith("tier1") else indel_or,
            "duplication_or_vs_tier1": 1.0 if label.startswith("tier1") else dup_or,
        }
        row.update(class_counts(sub))
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    detail = load_detail()
    tier2_subcategories = build_tier2_subcategory_summary(detail)
    tier1_tier2_classes = build_tier1_tier2_class_summary(detail)

    save_table(tier2_subcategories, TIER2_SUBCATEGORY_OUT)
    save_table(tier1_tier2_classes, TIER1_TIER2_CLASS_OUT)
    save_table(tier2_subcategories, TIER2_SUBCATEGORY_SUPP, sep="\t")
    save_table(tier1_tier2_classes, TIER1_TIER2_CLASS_SUPP, sep="\t")

    print("\nTier 2 subcategories")
    print(tier2_subcategories.to_string(index=False))
    print("\nTier 1 vs Tier 2 class comparison")
    print(tier1_tier2_classes.to_string(index=False))


if __name__ == "__main__":
    main()
