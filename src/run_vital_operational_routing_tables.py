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

LABEL_STATE_OUT = DATA_DIR / "vital_label_state_transition_operationalization.csv"
NON_EVALUABLE_OUT = DATA_DIR / "vital_non_evaluable_variant_routing.csv"
LABEL_STATE_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S32_label_state_transition_operationalization.tsv"
NON_EVALUABLE_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S33_non_evaluable_variant_routing.tsv"


LABEL_STATE_ROWS = [
    {
        "label_state_class": "High-penetrance Mendelian assertion",
        "trigger": "P/LP label remains frequency-consistent with mechanism-specific support.",
        "required_evidence_before_transition": (
            "Phenotype match, inheritance fit, segregation or de novo evidence when relevant, "
            "and no AC-supported population contradiction."
        ),
        "output_format": "P/LP - high-penetrance Mendelian assertion; disease and inheritance specified.",
    },
    {
        "label_state_class": "Susceptibility or context-dependent risk allele",
        "trigger": "AC-supported popmax/global AF incompatible with high-penetrance Mendelian disease, especially for haplotype, drug-triggered, or modifier contexts.",
        "required_evidence_before_transition": (
            "Literature or functional context supporting susceptibility/modifier action; ancestry-aware AF and AC; "
            "exclusion of straightforward Mendelian penetrance."
        ),
        "output_format": "LP/P risk allele with susceptibility/context qualifier; not monogenic diagnostic P/LP.",
    },
    {
        "label_state_class": "Carrier-compatible recessive assertion",
        "trigger": "Severe allele with frequency compatible with heterozygous carrier state in a recessive or biallelic architecture gene.",
        "required_evidence_before_transition": (
            "Inheritance model, phase, second allele status, and affected-state evidence."
        ),
        "output_format": "Carrier-state or affected-state-specific assertion; phase confirmation required.",
    },
    {
        "label_state_class": "Founder or low-penetrance allele",
        "trigger": "AC-supported high popmax in a plausible founder or reduced-penetrance context.",
        "required_evidence_before_transition": (
            "Ancestry-specific penetrance estimate, founder evidence, disease mechanism, and phenotype context."
        ),
        "output_format": "P/LP with penetrance and ancestry qualifier, or risk allele with penetrance qualifier.",
    },
    {
        "label_state_class": "Annotation-inflated assertion",
        "trigger": "Severe annotation plus weak review and AC-supported frequency contradiction without mechanism-specific support.",
        "required_evidence_before_transition": (
            "Allele-specific functional/splice data, phenotype linkage, segregation, transcript/NMD context, and expert review."
        ),
        "output_format": "Urgent expert review; candidate VUS/LB or mechanism-qualified assertion.",
    },
    {
        "label_state_class": "Representation-uncertain gray state",
        "trigger": "No exact usable AF because of no record, allele discordance, exact match without AF, or query failure.",
        "required_evidence_before_transition": (
            "Normalization review, genome-level query when available, orthogonal validation for complex alleles, and future freeze recheck."
        ),
        "output_format": "Gray no-frequency-evidence state; do not apply rarity inference by default.",
    },
]


NON_EVALUABLE_ROWS = [
    {
        "non_evaluable_variant_type": "Indel or duplication with no gnomAD record",
        "recommended_next_step": "Genome-level query plus orthogonal validation when clinically material.",
        "interpretation_guardrail": "Do not interpret absence as rarity; treat as representation-sensitive.",
    },
    {
        "non_evaluable_variant_type": "Allele-discordant SNV or MNV",
        "recommended_next_step": "Normalize, left-align, verify transcript/genomic allele, and requery exact allele.",
        "interpretation_guardrail": "Same-site observation is not exact-allele frequency evidence.",
    },
    {
        "non_evaluable_variant_type": "Exact match without usable AF block",
        "recommended_next_step": "Retain gray status and recheck in future gnomAD/ClinVar freeze.",
        "interpretation_guardrail": "Do not convert missing AF to AF=0.",
    },
    {
        "non_evaluable_variant_type": "No record SNV in constrained gene",
        "recommended_next_step": "Deferred review with phenotype, segregation, and mechanism context.",
        "interpretation_guardrail": "Absence supports neither frequency consistency nor contradiction by itself.",
    },
    {
        "non_evaluable_variant_type": "Complex allele or haplotype",
        "recommended_next_step": "Resolve allele representation, phase, and haplotype-specific literature before frequency interpretation.",
        "interpretation_guardrail": "Do not transfer single-site frequency assumptions to multi-allelic/haplotype assertions.",
    },
]


def save_table(df: pd.DataFrame, path: Path, sep: str = ",") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep=sep, index=False, lineterminator="\n")
    print(f"Saved {path} ({len(df)} rows)")


def main() -> None:
    label_state = pd.DataFrame(LABEL_STATE_ROWS)
    non_evaluable = pd.DataFrame(NON_EVALUABLE_ROWS)
    save_table(label_state, LABEL_STATE_OUT)
    save_table(non_evaluable, NON_EVALUABLE_OUT)
    save_table(label_state, LABEL_STATE_SUPP, sep="\t")
    save_table(non_evaluable, NON_EVALUABLE_SUPP, sep="\t")


if __name__ == "__main__":
    main()
