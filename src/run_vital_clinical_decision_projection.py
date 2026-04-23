from __future__ import annotations

from pathlib import Path

import pandas as pd


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
SUPPLEMENT_DIR = BASE_DIR / "supplementary_tables"

SCORES_IN = DATA_DIR / "arrhythmia_vital_scores.csv"
MANUAL_IN = DATA_DIR / "vital_red_manual_review_live_check.csv"
BOUNDARY_IN = DATA_DIR / "vital_public_evidence_boundary_audit.csv"

PROJECTION_OUT = DATA_DIR / "vital_red_clinical_decision_projection.csv"
GUIDELINE_OUT = DATA_DIR / "vital_guideline_tension_audit.csv"
PROJECTION_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S40_clinical_decision_projection.tsv"
GUIDELINE_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S41_guideline_tension_audit.tsv"


CLINVAR_CLASSIFICATION_DOC = "https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/"
CLINVAR_REVIEW_STATUS_DOC = "https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/"
ACMG_AMP_PMID = "https://pubmed.ncbi.nlm.nih.gov/25741868/"


CASE_PROJECTIONS = {
    "VCV000440850": {
        "state_aware_frame": "susceptibility_or_haplotype_context",
        "mendelian_diagnosis_pathway": "Can be read as support for Brugada syndrome diagnostic closure if accepted as generic high-penetrance P/LP.",
        "mendelian_cascade_pathway": "Can motivate broad familial variant testing as a monogenic disease allele.",
        "mendelian_management_pathway": "Can influence arrhythmia surveillance, drug-avoidance counseling, and phenotype-triggered ICD/risk-stratification discussions.",
        "state_aware_diagnosis_pathway": "Does not establish monogenic Brugada diagnosis by itself; phenotype, ECG/drug-trigger context, ancestry, and haplotype evidence remain decisive.",
        "state_aware_cascade_pathway": "Cascade testing should be limited or qualified as susceptibility/haplotype tracking rather than deterministic Mendelian segregation.",
        "state_aware_management_pathway": "Management should be phenotype- and context-driven, not label-driven.",
        "one_label_two_pathways": "Generic Pathogenic label can imply Mendelian disease, while the frequency/mechanism frame implies susceptibility or modifier biology.",
    },
    "VCV001325231": {
        "state_aware_frame": "carrier_or_biallelic_affected_state",
        "mendelian_diagnosis_pathway": "Can be misread as sufficient support for CPVT5 diagnosis under a dominant-like P/LP interpretation.",
        "mendelian_cascade_pathway": "Can motivate dominant-risk cascade counseling if phase and second-allele status are not checked.",
        "mendelian_management_pathway": "Can influence surveillance or restriction decisions before biallelic affected-state evidence is established.",
        "state_aware_diagnosis_pathway": "Heterozygous carrier observation is not an affected-state diagnosis without phase, second allele, and phenotype fit.",
        "state_aware_cascade_pathway": "Family testing should prioritize carrier/phase clarification and second-allele search.",
        "state_aware_management_pathway": "Management should follow phenotype and biallelic/recessive evidence, not a generic dominant P/LP read.",
        "one_label_two_pathways": "Generic Likely pathogenic label can imply disease causality, while the mechanism-consistent reading is carrier-compatible.",
    },
    "VCV004535537": {
        "state_aware_frame": "borderline_annotation_inflation_or_mechanism_qualified_review",
        "mendelian_diagnosis_pathway": "Can be read as support for Long QT syndrome diagnostic closure if accepted as unqualified LP.",
        "mendelian_cascade_pathway": "Can motivate familial testing and segregation assumptions before allele-specific splice evidence is resolved.",
        "mendelian_management_pathway": "Can influence QT surveillance or medication-avoidance decisions in a constrained dominant gene context.",
        "state_aware_diagnosis_pathway": "Does not establish high-penetrance LQTS causality from public evidence alone; splice effect, same-site ambiguity, and ancestry frequency need expert adjudication.",
        "state_aware_cascade_pathway": "Cascade testing should be considered provisional and linked to phenotype/segregation evidence.",
        "state_aware_management_pathway": "Management should remain phenotype-driven while variant-level evidence is adjudicated.",
        "one_label_two_pathways": "Generic Likely pathogenic label can imply a dominant LQTS allele, while the evidence supports a borderline mechanism-review state.",
    },
}


def read_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    for path in [SCORES_IN, MANUAL_IN, BOUNDARY_IN]:
        if not path.exists():
            raise FileNotFoundError(f"Missing input: {path}")
    scores = pd.read_csv(SCORES_IN)
    manual = pd.read_csv(MANUAL_IN)
    boundary = pd.read_csv(BOUNDARY_IN)
    for col in ["global_af", "global_ac", "popmax_af", "popmax_ac", "qualifying_frequency_ac", "vital_score"]:
        if col in scores.columns:
            scores[col] = pd.to_numeric(scores[col], errors="coerce")
    return scores, manual, boundary


def build_projection(scores: pd.DataFrame, manual: pd.DataFrame, boundary: pd.DataFrame) -> pd.DataFrame:
    red = scores[scores["vital_red_flag"].astype(str).str.lower().isin(["true", "1"])].copy()
    merged = red.merge(
        manual[
            [
                "gene",
                "clinvar_id",
                "live_clinvar_status",
                "condition",
                "review_support",
                "clinvar_url",
            ]
        ],
        on=["gene", "clinvar_id"],
        how="left",
    ).merge(
        boundary[
            [
                "gene",
                "clinvar_id",
                "structured_patient_level_phenotype_linkage_available",
                "public_penetrance_estimate_available",
                "public_segregation_or_phase_evidence_available",
                "documented_cascade_testing_or_therapy_outcome_available",
            ]
        ],
        on=["gene", "clinvar_id"],
        how="left",
    )

    rows = []
    for _, row in merged.sort_values("vital_score", ascending=False).iterrows():
        case = CASE_PROJECTIONS[str(row["clinvar_id"])]
        rows.append(
            {
                "gene": row["gene"],
                "clinvar_id": row["clinvar_id"],
                "condition": row.get("condition", ""),
                "current_public_label": row.get("live_clinvar_status", row.get("clinsig", "")),
                "review_support": row.get("review_support", row.get("review_status", "")),
                "frequency_signal": (
                    f"global_AF={row['global_af']:.3g}; global_AC={row['global_ac']:.0f}; "
                    f"popmax_AF={row['popmax_af']:.3g}; popmax_AC={row['popmax_ac']:.0f} "
                    f"({row.get('popmax_population', '')})"
                ),
                "vital_score": round(float(row["vital_score"]), 1),
                "state_aware_frame": case["state_aware_frame"],
                "if_generic_PLP_read_as_Mendelian_diagnosis": case["mendelian_diagnosis_pathway"],
                "if_generic_PLP_read_as_Mendelian_cascade_testing": case["mendelian_cascade_pathway"],
                "if_generic_PLP_read_as_Mendelian_management": case["mendelian_management_pathway"],
                "if_state_aware_read_diagnosis": case["state_aware_diagnosis_pathway"],
                "if_state_aware_read_cascade_testing": case["state_aware_cascade_pathway"],
                "if_state_aware_read_management": case["state_aware_management_pathway"],
                "one_label_two_incompatible_pathways": case["one_label_two_pathways"],
                "public_patient_level_evidence_available": public_evidence_summary(row),
                "clinical_projection_boundary": (
                    "Decision projection only; public data do not prove actual patient misdiagnosis, "
                    "cascade testing, ICD placement, or management change."
                ),
            }
        )
    return pd.DataFrame(rows)


def public_evidence_summary(row: pd.Series) -> str:
    unavailable = []
    for col, label in [
        ("structured_patient_level_phenotype_linkage_available", "phenotype linkage"),
        ("public_penetrance_estimate_available", "penetrance estimate"),
        ("public_segregation_or_phase_evidence_available", "segregation/phase"),
        ("documented_cascade_testing_or_therapy_outcome_available", "documented downstream outcome"),
    ]:
        if not bool(row.get(col, False)):
            unavailable.append(label)
    if unavailable:
        return "Not public: " + ", ".join(unavailable)
    return "Public evidence present for all tracked fields"


def build_guideline_tension(projection: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in projection.iterrows():
        if row["gene"] == "SCN5A":
            guideline_gap = (
                "ClinVar/ClinGen terminology allows risk-allele and low-penetrance labels, "
                "but this record is publicly aggregated as generic Pathogenic with no assertion criteria."
            )
        elif row["gene"] == "TRDN":
            guideline_gap = (
                "ACMG/ClinGen germline terms classify variant-condition assertions, but the generic label "
                "does not encode carrier versus biallelic affected-state use."
            )
        else:
            guideline_gap = (
                "Single-submitter LP status plus AC-supported frequency tension creates a guideline-facing "
                "need for expert adjudication, but the public label does not expose splice certainty or penetrance state."
            )
        rows.append(
            {
                "gene": row["gene"],
                "clinvar_id": row["clinvar_id"],
                "current_public_label": row["current_public_label"],
                "expert_panel_or_practice_guideline_review_present": "no",
                "clinvar_state_terms_relevant": "Pathogenic, low penetrance; Likely pathogenic, low penetrance; Established risk allele; Likely risk allele",
                "guideline_tension": guideline_gap,
                "why_this_is_not_a_ClinGen_disagreement_claim": (
                    "No red-priority case is asserted here to conflict with a ClinGen expert-panel classification; "
                    "the tension is between a generic public P/LP label and the state-aware interpretation required "
                    "for clinical use."
                ),
                "source_clinvar_classification_terms": CLINVAR_CLASSIFICATION_DOC,
                "source_clinvar_review_status": CLINVAR_REVIEW_STATUS_DOC,
                "source_acmg_amp_framework": ACMG_AMP_PMID,
            }
        )
    return pd.DataFrame(rows)


def save(df: pd.DataFrame, path: Path, sep: str = ",") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep=sep, index=False, lineterminator="\n")
    print(f"Saved {path.relative_to(BASE_DIR)} ({len(df)} rows)")


def main() -> None:
    scores, manual, boundary = read_inputs()
    projection = build_projection(scores, manual, boundary)
    guideline = build_guideline_tension(projection)

    save(projection, PROJECTION_OUT)
    save(guideline, GUIDELINE_OUT)
    save(projection, PROJECTION_SUPP, sep="\t")
    save(guideline, GUIDELINE_SUPP, sep="\t")


if __name__ == "__main__":
    main()
