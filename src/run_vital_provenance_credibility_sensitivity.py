from __future__ import annotations

from pathlib import Path

import pandas as pd


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
SUPPLEMENT_DIR = BASE_DIR / "supplementary_tables"

SCORES_IN = DATA_DIR / "arrhythmia_vital_scores.csv"
RECON_IN = DATA_DIR / "vital_tiered_match_reconciliation_detail.csv"

SUMMARY_OUT = DATA_DIR / "vital_provenance_credibility_filter_summary.csv"
EXCLUDED_OUT = DATA_DIR / "vital_provenance_credibility_filter_excluded_cases.csv"

SUMMARY_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S50_provenance_credibility_filter_summary.tsv"
EXCLUDED_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S51_provenance_credibility_filter_excluded_cases.tsv"


def save(df: pd.DataFrame, path: Path, sep: str = ",") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep=sep, index=False, lineterminator="\n")
    print(f"Saved {path.relative_to(BASE_DIR)} ({len(df)} rows)")


def main() -> None:
    scores = pd.read_csv(SCORES_IN)
    recon = pd.read_csv(RECON_IN)[["variant_key", "reconciliation_tier"]]
    df = scores.merge(recon, on="variant_key", how="left")

    low_credibility_mask = (
        df["review_strength"].isin(["single_submitter", "weak_or_no_assertion"])
        & (df["max_frequency_signal"].fillna(0) > 1e-4)
        & (df["qualifying_frequency_ac"].fillna(0) >= 20)
    )

    excluded = (
        df.loc[
            low_credibility_mask,
            [
                "variant_key",
                "gene",
                "clinvar_id",
                "title",
                "review_status",
                "review_strength",
                "submitter_count",
                "max_frequency_signal",
                "max_frequency_source",
                "qualifying_frequency_ac",
                "vital_score",
                "vital_band",
                "vital_red_flag",
            ],
        ]
        .sort_values(["max_frequency_signal", "qualifying_frequency_ac"], ascending=[False, False])
        .reset_index(drop=True)
    )

    filtered = df.loc[~low_credibility_mask].copy()

    exact_equivalent_tiers = {
        "tier1_exact_exome",
        "tier1_exact_genome_or_no_exome",
        "tier1_exact_found_by_region",
        "tier2_trim_normalized_equivalent",
        "tier2_event_equivalent_window",
        "tier2_decomposed_equivalent_unphased",
    }
    locus_regional_tiers = {"tier3_locus_observed_no_exact_af", "tier3_regional_context_only"}
    af_observed_filtered = filtered[filtered["frequency_evidence_status"].eq("frequency_observed")].copy()

    summary_rows = [
        {
            "metric": "excluded_low_credibility_extreme_frequency_records",
            "value": int(low_credibility_mask.sum()),
            "note": "Weak/no-assertion or single-submitters with max_frequency_signal >1e-4 and qualifying_frequency_ac >=20.",
        },
        {
            "metric": "filtered_total_records",
            "value": int(len(filtered)),
            "note": "Remaining arrhythmia P/LP assertions after the provenance credibility filter.",
        },
        {
            "metric": "filtered_exact_or_equivalent_records",
            "value": int(filtered["reconciliation_tier"].isin(exact_equivalent_tiers).sum()),
            "note": "Tier 1 exact/equivalent context after filtering.",
        },
        {
            "metric": "filtered_locus_or_regional_records",
            "value": int(filtered["reconciliation_tier"].isin(locus_regional_tiers).sum()),
            "note": "Tier 2 locus/regional context after filtering.",
        },
        {
            "metric": "filtered_still_unevaluable_records",
            "value": int(filtered["reconciliation_tier"].eq("tier4_still_unevaluable").sum()),
            "note": "Tier 3 genuinely unevaluable context after filtering.",
        },
        {
            "metric": "filtered_af_observed_records",
            "value": int(len(af_observed_filtered)),
            "note": "Canonical exact AF-observed subset after filtering.",
        },
        {
            "metric": "filtered_global_af_gt_1e_5",
            "value": int((af_observed_filtered["global_af"].fillna(0) > 1e-5).sum()),
            "note": "Global-only frequency alerts in the filtered exact AF subset.",
        },
        {
            "metric": "filtered_popmax_or_global_af_gt_1e_5",
            "value": int(af_observed_filtered["standard_acmg_frequency_flag"].sum()),
            "note": "Ancestry-aware frequency alerts in the filtered exact AF subset.",
        },
        {
            "metric": "filtered_popmax_only_alerts",
            "value": int(af_observed_filtered["popmax_only_frequency_flag"].sum()),
            "note": "Alerts missed by global AF but visible under popmax in the filtered exact AF subset.",
        },
        {
            "metric": "filtered_ac_supported_frequency_records",
            "value": int(af_observed_filtered["frequency_signal_ac_ge_20"].sum()),
            "note": "Remaining AC-supported frequency signals after filtering.",
        },
        {
            "metric": "filtered_vital_red_records",
            "value": int(af_observed_filtered["vital_red_flag"].sum()),
            "note": "Red-priority records remaining after filtering.",
        },
    ]

    summary = pd.DataFrame(summary_rows)
    save(summary, SUMMARY_OUT)
    save(excluded, EXCLUDED_OUT)
    save(summary, SUMMARY_SUPP, sep="\t")
    save(excluded, EXCLUDED_SUPP, sep="\t")


if __name__ == "__main__":
    main()
