from __future__ import annotations

from pathlib import Path

try:
    import pandas as pd
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires pandas. Install dependencies with: "
        "python -m pip install -r requirements.txt"
    ) from exc


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
SUPPLEMENT_DIR = BASE_DIR / "supplementary_tables"

MANUAL_REVIEW = DATA_DIR / "vital_red_manual_review_live_check.csv"
BOUNDARY_TABLE = DATA_DIR / "vital_public_evidence_boundary_audit.csv"
BOUNDARY_SUMMARY = DATA_DIR / "vital_public_evidence_boundary_summary.csv"
SUPPLEMENT_TABLE = SUPPLEMENT_DIR / "Supplementary_Table_S29_public_evidence_boundary_audit.tsv"


def has_public_citation(value: object) -> bool:
    text = "" if pd.isna(value) else str(value).strip()
    if not text:
        return False
    return not text.lower().startswith("no clinvar germline citations")


def weak_review(value: object) -> bool:
    text = "" if pd.isna(value) else str(value).lower()
    return "single" in text or "no assertion" in text or "weak" in text


def build_boundary_table(red_review: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, row in red_review.iterrows():
        rows.append(
            {
                "gene": row.get("gene"),
                "clinvar_id": row.get("clinvar_id"),
                "live_clinvar_status": row.get("live_clinvar_status"),
                "public_condition_label": row.get("condition"),
                "weak_or_single_submitter_public_review": weak_review(row.get("review_support")),
                "public_citation_or_mention_present": has_public_citation(
                    row.get("clinvar_publications_or_mentions")
                ),
                "structured_patient_level_phenotype_linkage_available": False,
                "public_penetrance_estimate_available": False,
                "public_segregation_or_phase_evidence_available": False,
                "documented_cascade_testing_or_therapy_outcome_available": False,
                "vital_supported_claim": (
                    "public frequency-assertion tension and decision-risk exposure"
                ),
                "vital_not_supported_claim": (
                    "proven patient misdiagnosis, penetrance estimate, phenotype linkage, "
                    "segregation, cascade-testing outcome, or therapy outcome"
                ),
                "evidence_boundary_interpretation": (
                    "public data support urgent expert review, not causal proof of downstream harm"
                ),
            }
        )
    return pd.DataFrame(rows)


def summarize_boundary(boundary: pd.DataFrame) -> pd.DataFrame:
    denominator = len(boundary)
    items = [
        (
            "P/LP public status and disease/condition label",
            boundary["live_clinvar_status"].notna() & boundary["public_condition_label"].notna(),
            "The public label is capable of entering diagnostic workflows.",
        ),
        (
            "Weak or single-submitter public review support",
            boundary["weak_or_single_submitter_public_review"],
            "The record is fragile enough to justify review prioritization.",
        ),
        (
            "Public citation or ClinVar literature mention",
            boundary["public_citation_or_mention_present"],
            "Literature exists for some cases, but this is not patient-level outcome evidence.",
        ),
        (
            "Structured patient-level phenotype linkage sufficient for causality",
            boundary["structured_patient_level_phenotype_linkage_available"],
            "Not available in the public scoring layer.",
        ),
        (
            "Public penetrance estimate",
            boundary["public_penetrance_estimate_available"],
            "Not available in the public scoring layer.",
        ),
        (
            "Public segregation or allelic phase evidence",
            boundary["public_segregation_or_phase_evidence_available"],
            "Not available in the public scoring layer.",
        ),
        (
            "Documented cascade-testing or therapy outcome",
            boundary["documented_cascade_testing_or_therapy_outcome_available"],
            "Not available in the public scoring layer.",
        ),
        (
            "Measurable public decision-risk exposure proxy",
            pd.Series([True] * denominator),
            "Available as record-level and submitter-exposure counts.",
        ),
    ]
    rows = []
    for evidence_question, mask, interpretation in items:
        count = int(mask.sum())
        rows.append(
            {
                "public_evidence_question": evidence_question,
                "red_priority_cases_with_evidence": count,
                "denominator": denominator,
                "percent": 100 * count / denominator if denominator else 0.0,
                "interpretation": interpretation,
            }
        )
    return pd.DataFrame(rows)


def save_table(df: pd.DataFrame, path: Path, sep: str = ",") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep=sep, index=False, lineterminator="\n")
    print(f"Saved {path} ({len(df)} rows)")


def main() -> None:
    if not MANUAL_REVIEW.exists():
        raise FileNotFoundError(f"Missing manual review table: {MANUAL_REVIEW}")
    red_review = pd.read_csv(MANUAL_REVIEW)
    boundary = build_boundary_table(red_review)
    summary = summarize_boundary(boundary)
    save_table(boundary, BOUNDARY_TABLE)
    save_table(summary, BOUNDARY_SUMMARY)
    save_table(boundary, SUPPLEMENT_TABLE, sep="\t")


if __name__ == "__main__":
    main()
