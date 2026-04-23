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

SCORES = DATA_DIR / "arrhythmia_vital_scores.csv"
MANUAL_REVIEW = DATA_DIR / "vital_red_manual_review_live_check.csv"
OUTPUT = DATA_DIR / "vital_logical_incompatibility_audit.csv"
SUMMARY_OUTPUT = DATA_DIR / "vital_logical_incompatibility_summary.csv"
SUPPLEMENT_OUTPUT = SUPPLEMENT_DIR / "Supplementary_Table_S31_logical_incompatibility_audit.tsv"


CASE_RULES = {
    "VCV000440850": {
        "mechanism_frame": "susceptibility_or_haplotype_state",
        "high_penetrance_mendelian_expectation": (
            "A high-penetrance dominant Brugada/arrhythmia allele should not reach a "
            "population-specific frequency near 0.6% with AC-supported observation."
        ),
        "logical_conclusion": (
            "The unqualified high-penetrance Mendelian reading is logically incompatible "
            "with the observed AFR popmax frequency unless the biological state is changed."
        ),
        "state_required_for_consistency": (
            "context-dependent susceptibility, haplotype, low-penetrance, or modifier state"
        ),
    },
    "VCV001325231": {
        "mechanism_frame": "carrier_or_biallelic_affected_state",
        "high_penetrance_mendelian_expectation": (
            "A dominant high-penetrance CPVT assertion would require stronger rarity than "
            "the observed AC-supported carrier-compatible frequency."
        ),
        "logical_conclusion": (
            "The dominant Mendelian reading is logically incompatible; the observation can "
            "remain coherent only under carrier, phase-aware, or biallelic affected-state framing."
        ),
        "state_required_for_consistency": (
            "carrier-compatible recessive or affected-state-specific assertion"
        ),
    },
    "VCV004535537": {
        "mechanism_frame": "annotation_inflation_boundary_state",
        "high_penetrance_mendelian_expectation": (
            "A constrained dominant LQTS splice assertion should combine strong rarity with "
            "robust allele-specific evidence."
        ),
        "logical_conclusion": (
            "The unqualified high-penetrance reading is not supportable from public data alone; "
            "frequency tension plus weak review makes this a boundary case requiring expert adjudication."
        ),
        "state_required_for_consistency": (
            "unresolved annotation-inflation or mechanism-qualified review state"
        ),
    },
}


def format_af(value: object) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "NA"
    return f"{numeric:.3g}"


def format_ac(value: object) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "NA"
    return str(int(numeric))


def build_audit(scores: pd.DataFrame, manual: pd.DataFrame) -> pd.DataFrame:
    red = scores.loc[scores["vital_red_flag"].astype(str).str.lower().isin({"true", "1"})].copy()
    manual_small = manual[
        [
            "clinvar_id",
            "live_clinvar_status",
            "condition",
            "review_support",
            "clinvar_publications_or_mentions",
            "why_suspicious_but_not_resolved",
        ]
    ].copy()
    merged = red.merge(manual_small, on="clinvar_id", how="left")
    rows: list[dict[str, object]] = []
    for _, row in merged.iterrows():
        rules = CASE_RULES.get(str(row["clinvar_id"]), {})
        pop = str(row.get("popmax_population", "NA")).upper()
        population_reality = (
            f"global AF={format_af(row.get('global_af'))}, global AC={format_ac(row.get('global_ac'))}; "
            f"{pop} popmax AF={format_af(row.get('popmax_af'))}, "
            f"popmax AC={format_ac(row.get('popmax_ac'))}; "
            f"qualifying AC={format_ac(row.get('qualifying_frequency_ac'))}"
        )
        clinvar_assertion = (
            f"{row.get('live_clinvar_status', row.get('clinsig'))} for "
            f"{row.get('condition', 'NA')}; review={row.get('review_support', row.get('review_strength'))}"
        )
        rows.append(
            {
                "gene": row["gene"],
                "clinvar_id": row["clinvar_id"],
                "variant": row["title"],
                "clinvar_assertion_public_claim": clinvar_assertion,
                "population_reality": population_reality,
                "mechanism_frame": rules.get("mechanism_frame", "NA"),
                "high_penetrance_mendelian_expectation": rules.get(
                    "high_penetrance_mendelian_expectation", "NA"
                ),
                "logical_conclusion": rules.get("logical_conclusion", "NA"),
                "state_required_for_consistency": rules.get(
                    "state_required_for_consistency", "NA"
                ),
                "not_claimed": (
                    "This audit does not prove benignity, misdiagnosis, penetrance, or patient outcome; "
                    "it proves that the unqualified P/LP label is insufficient or incompatible with a "
                    "high-penetrance Mendelian reading."
                ),
                "vital_score": row["vital_score"],
                "vital_signal_reason": row["vital_signal_reason"],
                "supporting_public_note": row.get("why_suspicious_but_not_resolved", "NA"),
            }
        )
    return pd.DataFrame(rows).sort_values("vital_score", ascending=False)


def build_summary(audit: pd.DataFrame) -> pd.DataFrame:
    rows = [
        {
            "summary_statement": "red_cases_with_logical_incompatibility_or_state_insufficiency",
            "count": len(audit),
            "denominator": len(audit),
            "interpretation": (
                "Every red case requires a state change or state qualification before the public P/LP label can be read coherently."
            ),
        },
        {
            "summary_statement": "susceptibility_or_haplotype_case",
            "count": int(audit["mechanism_frame"].eq("susceptibility_or_haplotype_state").sum()),
            "denominator": len(audit),
            "interpretation": "High population frequency breaks the high-penetrance Mendelian reading.",
        },
        {
            "summary_statement": "carrier_or_biallelic_case",
            "count": int(audit["mechanism_frame"].eq("carrier_or_biallelic_affected_state").sum()),
            "denominator": len(audit),
            "interpretation": "Carrier-compatible frequency can be coherent only with phase-aware disease framing.",
        },
        {
            "summary_statement": "annotation_inflation_boundary_case",
            "count": int(audit["mechanism_frame"].eq("annotation_inflation_boundary_state").sum()),
            "denominator": len(audit),
            "interpretation": "Severe annotation does not substitute for robust allele-specific evidence.",
        },
    ]
    return pd.DataFrame(rows)


def save_table(df: pd.DataFrame, path: Path, sep: str = ",") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep=sep, index=False, lineterminator="\n")
    print(f"Saved {path} ({len(df)} rows)")


def main() -> None:
    if not SCORES.exists():
        raise FileNotFoundError(f"Missing VITAL score table: {SCORES}")
    if not MANUAL_REVIEW.exists():
        raise FileNotFoundError(f"Missing manual review table: {MANUAL_REVIEW}")
    scores = pd.read_csv(SCORES, low_memory=False)
    manual = pd.read_csv(MANUAL_REVIEW)
    audit = build_audit(scores, manual)
    summary = build_summary(audit)
    save_table(audit, OUTPUT)
    save_table(summary, SUMMARY_OUTPUT)
    save_table(audit, SUPPLEMENT_OUTPUT, sep="\t")


if __name__ == "__main__":
    main()
