from __future__ import annotations

from pathlib import Path

import pandas as pd


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
SUPPLEMENTARY_DIR = BASE_DIR / "supplementary_tables"

HISTORICAL_PREDICTIONS = (
    DATA_DIR / "vital_cross_disease_3000_2023_01_to_current_vital_historical_predictions.csv"
)


def bool_series(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def pct(numerator: int | float, denominator: int | float) -> float:
    return float(100 * numerator / denominator) if denominator else 0.0


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    SUPPLEMENTARY_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(HISTORICAL_PREDICTIONS)
    df["strict_event"] = bool_series(df["strict_revised_to_b_or_vus"])
    df["red_priority"] = bool_series(df["vital_red_flag"])
    df["expert_panel_current"] = df["review_strength"].eq("expert_panel")
    df["frequency_observed"] = df["frequency_evidence_status"].eq("frequency_observed")

    strict = df[df["strict_event"]].copy()
    red = df[df["red_priority"]].copy()
    red_strict = strict[strict["red_priority"]].copy()
    expert_strict = strict[strict["expert_panel_current"]].copy()
    frequency_strict = strict[strict["frequency_observed"]].copy()

    summary = pd.DataFrame(
        [
            {
                "audit_layer": "real_strict_reclassification_events",
                "denominator": len(df),
                "event_count": len(strict),
                "event_definition": "ClinVar P/LP in January 2023 to VUS/B/LB by April 2026",
                "vital_red_count": len(red),
                "vital_red_event_count": len(red_strict),
                "vital_red_precision_percent": pct(len(red_strict), len(red)),
                "vital_red_recall_percent": pct(len(red_strict), len(strict)),
                "clinical_interpretation": (
                    "This is the real classification-change endpoint. VITAL red is not sensitive for all future "
                    "downgrades, but it did capture one high-frequency weak-review P/LP assertion that later moved to VUS."
                ),
            },
            {
                "audit_layer": "frequency_observed_real_strict_events",
                "denominator": int(df["frequency_observed"].sum()),
                "event_count": len(frequency_strict),
                "event_definition": "Strict future downgrade among records with observed frequency evidence",
                "vital_red_count": len(red),
                "vital_red_event_count": len(red_strict),
                "vital_red_precision_percent": pct(len(red_strict), len(red)),
                "vital_red_recall_percent": pct(len(red_strict), len(frequency_strict)),
                "clinical_interpretation": (
                    "Restricting to variants where VITAL can read population frequency still leaves low recall; "
                    "the endpoint validates scope rather than broad prediction."
                ),
            },
            {
                "audit_layer": "expert_panel_curated_strict_downgrade",
                "denominator": len(strict),
                "event_count": len(expert_strict),
                "event_definition": "Strict future downgrade with current expert-panel review status",
                "vital_red_count": int(expert_strict["red_priority"].sum()),
                "vital_red_event_count": int(expert_strict["red_priority"].sum()),
                "vital_red_precision_percent": pct(int(expert_strict["red_priority"].sum()), max(len(expert_strict), 1)),
                "vital_red_recall_percent": pct(int(expert_strict["red_priority"].sum()), len(expert_strict)),
                "clinical_interpretation": (
                    "External curated downgrade stress-test: the expert-panel RYR1 downgrade had no frequency tension, "
                    "so VITAL appropriately did not flag it. VITAL is not a general reclassification detector."
                ),
            },
        ]
    )

    case_cols = [
        "gene_baseline",
        "clinvar_id",
        "name",
        "followup_clinical_group",
        "followup_clinical_significance",
        "review_status_score",
        "review_strength",
        "vital_score",
        "vital_band",
        "red_priority",
        "max_frequency_signal",
        "qualifying_frequency_ac",
        "vital_signal_reason",
    ]
    case_cols = [c for c in case_cols if c in df.columns]

    real_red_cases = red_strict.loc[:, case_cols].copy()
    real_red_cases["case_role"] = "real_strict_reclassification_captured_by_VITAL_red"
    real_red_cases["clinical_impact"] = (
        "Public ClinVar classification changed from P/LP to VUS/B/LB; patient-level downstream outcomes are not available."
    )

    curated_cases = expert_strict.loc[:, case_cols].copy()
    curated_cases["case_role"] = "expert_panel_curated_strict_downgrade_not_frequency_tension"
    curated_cases["clinical_impact"] = (
        "External curated downgrade event; absence of frequency tension marks the boundary of VITAL's intended scope."
    )

    all_cases = pd.concat([real_red_cases, curated_cases], ignore_index=True, sort=False)

    summary.to_csv(DATA_DIR / "vital_real_reclassification_clinical_impact_audit.csv", index=False)
    all_cases.to_csv(DATA_DIR / "vital_real_reclassification_case_examples.csv", index=False)

    supplement = pd.concat(
        [
            summary.assign(table_section="real_reclassification_summary"),
            all_cases.assign(table_section="real_and_curated_case_examples"),
        ],
        ignore_index=True,
        sort=False,
    )
    supplement.fillna("NA").to_csv(
        SUPPLEMENTARY_DIR / "Supplementary_Table_S23_real_reclassification_truth_audit.tsv",
        sep="\t",
        index=False,
    )
    print("Saved real reclassification and curated truth audit outputs")


if __name__ == "__main__":
    main()
