from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires pandas, numpy, and matplotlib. Install them with: "
        "python -m pip install -r requirements.txt"
    ) from exc


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
FIGURE_DIR = BASE_DIR / "figures"

ARRHYTHMIA_SCORES_IN = DATA_DIR / "arrhythmia_vital_scores.csv"
ROUTING_CALLS_IN = DATA_DIR / "vital_routing_validation_calls.csv"

VARIANT_DESIGN_OUT = DATA_DIR / "vital_autopsy_denovo_variant_design.csv"
CASE_CALLS_OUT = DATA_DIR / "vital_autopsy_denovo_case_calls.csv"
SUMMARY_OUT = DATA_DIR / "vital_autopsy_denovo_summary.csv"
FALSE_ATTRIBUTION_RATES_OUT = DATA_DIR / "vital_autopsy_denovo_false_attribution_rates.csv"
BY_GENE_OUT = DATA_DIR / "vital_autopsy_denovo_by_gene.csv"
BY_MECHANISM_OUT = DATA_DIR / "vital_autopsy_denovo_by_mechanism.csv"
BY_AF_BIN_OUT = DATA_DIR / "vital_autopsy_denovo_by_af_bin.csv"
BY_DENOVO_OUT = DATA_DIR / "vital_autopsy_denovo_by_denovo_status.csv"
BY_AUTOPSY_CONTEXT_OUT = DATA_DIR / "vital_autopsy_denovo_by_autopsy_context.csv"
SENSITIVITY_PENETRANCE_OUT = DATA_DIR / "vital_autopsy_denovo_penetrance_sensitivity.csv"
SENSITIVITY_DENOVO_OUT = DATA_DIR / "vital_autopsy_denovo_denovo_rate_sensitivity.csv"
SENSITIVITY_MCAF_OUT = DATA_DIR / "vital_autopsy_denovo_mcaf_sensitivity.csv"
MECHANISM_PERMUTATION_OUT = DATA_DIR / "vital_autopsy_denovo_mechanism_permutation.csv"
GOLD_STANDARD_OUT = DATA_DIR / "vital_autopsy_denovo_gold_standard_preservation.csv"
CONFLICTING_OUT = DATA_DIR / "vital_autopsy_denovo_conflicting_clinvar_layer.csv"
MODEL_INPUTS_OUT = DATA_DIR / "vital_autopsy_denovo_model_inputs.csv"

FLOW_FIGURE_OUT = FIGURE_DIR / "vital_autopsy_denovo_decision_flow.png"
FALSE_ATTRIBUTION_FIGURE_OUT = FIGURE_DIR / "vital_autopsy_denovo_false_attribution.png"
DENOVO_FIGURE_OUT = FIGURE_DIR / "vital_autopsy_denovo_override.png"
MECHANISM_FIGURE_OUT = FIGURE_DIR / "vital_autopsy_denovo_mechanism.png"

RNG_SEED = 20260429
N_CASES = 10_000

ARRHYTHMIA_AUDIT_GENES = [
    "RYR2",
    "CASQ2",
    "TRDN",
    "KCNQ1",
    "KCNH2",
    "SCN5A",
    "CALM1",
    "CALM2",
    "CALM3",
    "ANK2",
    "KCNE1",
    "KCNE2",
    "KCNJ2",
    "CACNA1C",
]

GOF_DN_GENES = {"RYR2", "CALM1", "CALM2", "CALM3", "KCNJ2"}
AR_GENES = {"CASQ2", "TRDN"}
HI_GENES = {"KCNH2", "ANK2"}
MIXED_GENES = {"KCNQ1", "SCN5A", "KCNE1", "KCNE2", "CACNA1C"}

MCAF_STRICT_DEFAULT = 2.5e-5
MCAF_HARD_DEFAULT = 1e-4
MCAF_PERMISSIVE_DEFAULT = 1e-3

AF_BINS = [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, np.inf]
AF_LABELS = ["<1e-6", "1e-6-1e-5", "1e-5-1e-4", "1e-4-1e-3", "1e-3-1e-2", ">=1e-2"]


@dataclass(frozen=True)
class PenetranceProfile:
    name: str
    gof_dn_alpha: float
    gof_dn_beta: float
    hi_alpha: float
    hi_beta: float
    mixed_alpha: float
    mixed_beta: float
    ar_alpha: float
    ar_beta: float


PENETRANCE_PROFILES = {
    "high": PenetranceProfile("high", 8, 2, 8, 2, 8, 2, 9, 1),
    "moderate": PenetranceProfile("moderate", 8, 2, 4, 6, 4, 6, 8, 2),
    "low": PenetranceProfile("low", 4, 6, 1, 19, 1, 19, 4, 6),
    "mixed": PenetranceProfile("mixed", 8, 2, 4, 6, 1, 19, 8, 2),
}


def pct(numerator: int | float, denominator: int | float) -> float:
    return 0.0 if denominator == 0 else 100.0 * float(numerator) / float(denominator)


def truthy(value: object) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def save_table(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, lineterminator="\n")
    print(f"Saved {path.relative_to(BASE_DIR)} ({len(df)} rows)")


def mechanism_for_gene(gene: object) -> str:
    symbol = str(gene).strip().upper()
    if symbol in GOF_DN_GENES:
        return "GoF_DN"
    if symbol in AR_GENES:
        return "AR"
    if symbol in HI_GENES:
        return "HI"
    if symbol in MIXED_GENES:
        return "MIXED"
    return "UNKNOWN"


def inheritance_for_mechanism(mechanism: str) -> str:
    if mechanism == "AR":
        return "autosomal_recessive"
    if mechanism in {"GoF_DN", "HI"}:
        return "autosomal_dominant"
    if mechanism == "MIXED":
        return "mixed_context_dependent"
    return "unknown"


def clingen_validity_for_gene(gene: object) -> str:
    symbol = str(gene).strip().upper()
    if symbol in {"RYR2", "KCNQ1", "KCNH2", "SCN5A", "CALM1", "CALM2", "CALM3", "CASQ2", "TRDN"}:
        return "strong_or_definitive"
    if symbol in {"ANK2", "KCNE1", "KCNJ2", "CACNA1C"}:
        return "moderate_or_disputed_context"
    return "limited_or_unknown"


def normalize_evaluability(value: object) -> str:
    text = str(value).strip().upper()
    if text in {"TIER_1", "TIER1"}:
        return "Tier1"
    if text in {"TIER_2", "TIER2"}:
        return "Tier2"
    return "Unevaluable"


def max_af(row: pd.Series) -> float:
    values = []
    for column in ["global_AF", "popmax_AF", "global_af", "popmax_af", "max_frequency_signal"]:
        if column in row:
            value = pd.to_numeric(row.get(column), errors="coerce")
            if pd.notna(value):
                values.append(float(value))
    return max(values) if values else np.nan


def build_model_inputs() -> pd.DataFrame:
    rows = []
    for gene in ARRHYTHMIA_AUDIT_GENES:
        mechanism = mechanism_for_gene(gene)
        expected = {
            "GoF_DN": "VITAL_OK when allele-resolved and AF-compatible; de novo supports only coherent dominant model.",
            "HI": "CHECK_POPMAX or VITAL_OK depending on penetrance and AF; de novo does not override high AF.",
            "AR": "CHECK_MODEL for heterozygous carrier state under a flattened dominant baseline.",
            "MIXED": "CHECK_MODEL unless phenotype/mechanism specifies the asserted disease model.",
            "UNKNOWN": "CHECK_MODEL because mechanism is insufficiently specified.",
        }[mechanism]
        rows.append(
            {
                "gene": gene,
                "mechanism_class": mechanism,
                "inheritance_model": inheritance_for_mechanism(mechanism),
                "ClinGen_validity_gene_level": clingen_validity_for_gene(gene),
                "expected_VITAL_behavior": expected,
            }
        )
    return pd.DataFrame(rows)


def infer_vital_route_from_row(row: pd.Series, mcaf_strict: float, mcaf_hard: float) -> str:
    label = str(row.get("ClinVar_label", "P")).strip()
    if label not in {"P", "LP"}:
        return "NO_PLP_BASELINE"
    if str(row.get("evaluability_tier", "")) != "Tier1":
        return "EVAL_LIMITED"
    effective_af = float(row.get("effective_AF", 0.0))
    mechanism = str(row.get("mechanism_class", "UNKNOWN"))
    if effective_af > mcaf_hard:
        return "MODEL_CONFLICT"
    if effective_af > mcaf_strict:
        return "CHECK_POPMAX"
    if mechanism == "AR":
        return "CHECK_MODEL"
    if mechanism in {"MIXED", "UNKNOWN"}:
        return "CHECK_MODEL"
    return "VITAL_OK"


def build_empirical_variants() -> pd.DataFrame:
    if ROUTING_CALLS_IN.exists():
        raw = pd.read_csv(ROUTING_CALLS_IN, low_memory=False)
    elif ARRHYTHMIA_SCORES_IN.exists():
        raw = pd.read_csv(ARRHYTHMIA_SCORES_IN, low_memory=False)
    else:
        raise FileNotFoundError("No arrhythmia score/routing table found.")

    raw = raw.loc[raw["gene"].fillna("").astype(str).str.upper().isin(ARRHYTHMIA_AUDIT_GENES)].copy()
    raw["gene"] = raw["gene"].astype(str).str.upper()
    raw["variant_id"] = raw["clinvar_id"].fillna(raw.get("variation_id", "")).astype(str)
    raw["mechanism_class"] = raw["gene"].map(mechanism_for_gene)
    raw["inheritance_model"] = raw["mechanism_class"].map(inheritance_for_mechanism)
    raw["ClinVar_label"] = "P"
    raw.loc[raw.get("clinsig", "").astype(str).str.contains("likely", case=False, na=False), "ClinVar_label"] = "LP"
    if "VITAL_evaluability" in raw.columns:
        raw["evaluability_tier"] = raw["VITAL_evaluability"].map(normalize_evaluability)
    else:
        raw["evaluability_tier"] = np.where(raw["frequency_evidence_status"].eq("frequency_observed"), "Tier1", "Unevaluable")
    raw["global_AF"] = pd.to_numeric(raw.get("global_af"), errors="coerce")
    raw["popmax_AF"] = pd.to_numeric(raw.get("popmax_af"), errors="coerce")
    raw["effective_AF"] = raw.apply(max_af, axis=1)
    raw["allele_count"] = pd.to_numeric(raw.get("qualifying_frequency_ac"), errors="coerce")
    raw["variant_consequence"] = raw.get("functional_class", "").fillna("").astype(str)
    raw["ClinGen_validity_gene_level"] = raw["gene"].map(clingen_validity_for_gene)
    raw["MCAF_strict"] = MCAF_STRICT_DEFAULT
    raw["MCAF_hard"] = MCAF_HARD_DEFAULT
    raw["MCAF_permissive"] = MCAF_PERMISSIVE_DEFAULT
    raw["variant_source"] = "empirical_clinvar_gnomad"

    columns = [
        "variant_id",
        "gene",
        "ClinVar_label",
        "review_status",
        "review_score",
        "mechanism_class",
        "inheritance_model",
        "variant_consequence",
        "global_AF",
        "popmax_AF",
        "effective_AF",
        "allele_count",
        "evaluability_tier",
        "MCAF_strict",
        "MCAF_hard",
        "MCAF_permissive",
        "ClinGen_validity_gene_level",
        "variant_source",
    ]
    return raw.loc[:, [column for column in columns if column in raw.columns]].reset_index(drop=True)


def build_simulated_variants(n: int, rng: np.random.Generator) -> pd.DataFrame:
    genes = rng.choice(ARRHYTHMIA_AUDIT_GENES, size=n, replace=True)
    popmax = 10 ** rng.uniform(-6, -2, size=n)
    global_af = popmax * rng.uniform(0.05, 0.9, size=n)
    labels = rng.choice(["P", "LP", "VUS", "conflicting", "B/LB"], size=n, p=[0.42, 0.33, 0.1, 0.1, 0.05])
    rows = []
    for index, gene in enumerate(genes):
        mechanism = mechanism_for_gene(gene)
        rows.append(
            {
                "variant_id": f"SIM_DENOVO_{index + 1:05d}",
                "gene": gene,
                "ClinVar_label": labels[index],
                "review_status": "simulated",
                "review_score": 0,
                "mechanism_class": mechanism,
                "inheritance_model": inheritance_for_mechanism(mechanism),
                "variant_consequence": rng.choice(["missense", "lof", "splice_or_intronic"], p=[0.68, 0.22, 0.10]),
                "global_AF": global_af[index],
                "popmax_AF": popmax[index],
                "effective_AF": max(global_af[index], popmax[index]),
                "allele_count": int(max(1, round(popmax[index] * 150_000))),
                "evaluability_tier": rng.choice(["Tier1", "Tier2", "Unevaluable"], p=[0.74, 0.18, 0.08]),
                "MCAF_strict": MCAF_STRICT_DEFAULT,
                "MCAF_hard": MCAF_HARD_DEFAULT,
                "MCAF_permissive": MCAF_PERMISSIVE_DEFAULT,
                "ClinGen_validity_gene_level": clingen_validity_for_gene(gene),
                "variant_source": "simulated_af_space",
            }
        )
    return pd.DataFrame(rows)


def build_variant_design(rng: np.random.Generator) -> pd.DataFrame:
    empirical = build_empirical_variants()
    simulated = build_simulated_variants(max(1, round(len(empirical) * 0.43)), rng)
    design = pd.concat([empirical, simulated], ignore_index=True, sort=False)
    design["effective_AF"] = pd.to_numeric(design["effective_AF"], errors="coerce")
    fallback = 10 ** rng.uniform(-6, -2, size=len(design))
    design["effective_AF"] = design["effective_AF"].fillna(pd.Series(fallback, index=design.index))
    design["popmax_AF"] = pd.to_numeric(design["popmax_AF"], errors="coerce").fillna(design["effective_AF"])
    design["global_AF"] = pd.to_numeric(design["global_AF"], errors="coerce").fillna(design["effective_AF"] * 0.5)
    design["allele_count"] = pd.to_numeric(design["allele_count"], errors="coerce").fillna(
        np.maximum(1, np.round(design["effective_AF"] * 150_000))
    )
    design["af_bin"] = pd.cut(design["effective_AF"], bins=AF_BINS, labels=AF_LABELS, right=False)
    design["vital_route_main"] = design.apply(
        lambda row: infer_vital_route_from_row(row, MCAF_STRICT_DEFAULT, MCAF_HARD_DEFAULT),
        axis=1,
    )
    design["gene_actionable"] = design["gene"].isin(ARRHYTHMIA_AUDIT_GENES)
    return design.reset_index(drop=True)


def sample_candidate_genotype(row: pd.Series, rng: np.random.Generator) -> str:
    p = min(max(float(row.get("effective_AF", 0.0)), 1e-9), 0.49)
    het = 2 * p * (1 - p)
    hom = p * p
    if het + hom <= 0:
        return "non_carrier"
    hom_conditional = hom / (het + hom)
    mechanism = str(row.get("mechanism_class", "UNKNOWN"))
    if mechanism == "AR" and rng.random() < 0.02:
        return "compound_heterozygous"
    return "homozygous" if rng.random() < hom_conditional else "heterozygous"


def sample_denovo_status(rng: np.random.Generator, confirmed_rate: float = 0.05) -> tuple[str, str]:
    draw = rng.random()
    if draw < confirmed_rate:
        return "confirmed", "trio_confirmed"
    if draw < confirmed_rate + 0.03:
        return "assumed", "presumed"
    if draw < confirmed_rate + 0.13:
        return "unknown", "not_tested"
    return "absent", "not_tested"


def penetrance_for_case(mechanism: str, profile: PenetranceProfile, rng: np.random.Generator) -> float:
    if mechanism == "GoF_DN":
        return float(rng.beta(profile.gof_dn_alpha, profile.gof_dn_beta))
    if mechanism == "HI":
        return float(rng.beta(profile.hi_alpha, profile.hi_beta))
    if mechanism == "MIXED":
        return float(rng.beta(profile.mixed_alpha, profile.mixed_beta))
    if mechanism == "AR":
        return float(rng.beta(profile.ar_alpha, profile.ar_beta))
    return float(rng.beta(1, 19))


def determine_true_causal(row: pd.Series, penetrance: float, rng: np.random.Generator) -> bool:
    genotype = str(row["genotype"])
    mechanism = str(row["mechanism_class"])
    effective_af = float(row["effective_AF"])
    denovo = str(row["de_novo_status"])
    if genotype == "non_carrier":
        return False
    if mechanism == "GoF_DN":
        if genotype in {"heterozygous", "homozygous"} and effective_af <= MCAF_STRICT_DEFAULT:
            adjusted = min(0.98, penetrance * (1.25 if denovo == "confirmed" else 1.0))
            return bool(rng.random() < adjusted)
        return False
    if mechanism == "HI":
        if genotype in {"heterozygous", "homozygous"} and effective_af <= MCAF_HARD_DEFAULT:
            adjusted = min(0.9, penetrance * (1.2 if denovo == "confirmed" else 1.0))
            return bool(rng.random() < adjusted)
        return False
    if mechanism == "MIXED":
        if genotype in {"heterozygous", "homozygous"} and effective_af <= MCAF_HARD_DEFAULT:
            adjusted = min(0.8, penetrance * (1.2 if denovo == "confirmed" else 1.0))
            return bool(rng.random() < adjusted)
        return False
    if mechanism == "AR":
        if genotype in {"homozygous", "compound_heterozygous"}:
            return bool(rng.random() < penetrance)
        return False
    return False


def vital_causal_attribution(vital_route: str, modifier: str | None) -> bool:
    return vital_route == "VITAL_OK" and modifier in {"DE_NOVO_SUPPORT", None}


def simulate_cases(
    design: pd.DataFrame,
    n_cases: int,
    rng: np.random.Generator,
    penetrance_profile: str = "mixed",
    confirmed_denovo_rate: float = 0.05,
    mcaf_strict: float = MCAF_STRICT_DEFAULT,
    mcaf_hard: float = MCAF_HARD_DEFAULT,
    mechanism_permutation: bool = False,
) -> pd.DataFrame:
    profile = PENETRANCE_PROFILES[penetrance_profile]
    weights = np.where(design["variant_source"].eq("empirical_clinvar_gnomad"), 0.7 / max(1, design["variant_source"].eq("empirical_clinvar_gnomad").sum()), 0.3 / max(1, design["variant_source"].eq("simulated_af_space").sum()))
    weights = weights / weights.sum()
    selected = design.iloc[rng.choice(design.index.to_numpy(), size=n_cases, replace=True, p=weights)].copy()
    selected = selected.reset_index(drop=True)
    selected.insert(0, "case_id", [f"PM_CASE_{index + 1:05d}" for index in range(n_cases)])

    if mechanism_permutation:
        selected["mechanism_class"] = rng.permutation(selected["mechanism_class"].to_numpy())
        selected["inheritance_model"] = selected["mechanism_class"].map(inheritance_for_mechanism)

    selected["autopsy_status"] = "NEGATIVE"
    selected["structural_heart_disease"] = "absent"
    selected["toxicology"] = rng.choice(["negative", "unknown"], size=n_cases, p=[0.86, 0.14])
    selected["phenotype_strength"] = rng.choice(["none", "weak", "suggestive"], size=n_cases, p=[0.80, 0.15, 0.05])
    selected["family_history"] = rng.choice(["unknown", "absent", "positive"], size=n_cases, p=[0.72, 0.23, 0.05])
    selected["circumstance"] = rng.choice(["unknown", "exercise", "emotion", "sleep"], size=n_cases, p=[0.58, 0.18, 0.08, 0.16])
    selected["ancestry_group"] = rng.choice(["AFR", "AMR", "EAS", "EUR", "MID", "SAS", "OTH"], size=n_cases)

    genotypes = []
    denovo_statuses = []
    denovo_evidence = []
    penetrance_values = []
    true_causal = []
    vital_routes = []
    modifiers = []
    for _, row in selected.iterrows():
        genotype = sample_candidate_genotype(row, rng)
        denovo, evidence = sample_denovo_status(rng, confirmed_rate=confirmed_denovo_rate)
        penetrance = penetrance_for_case(str(row["mechanism_class"]), profile, rng)
        temp = row.copy()
        temp["genotype"] = genotype
        temp["de_novo_status"] = denovo
        route = infer_vital_route_from_row(temp, mcaf_strict, mcaf_hard)
        if denovo == "confirmed":
            modifier = "DE_NOVO_SUPPORT" if route == "VITAL_OK" else "DE_NOVO_NOT_SUFFICIENT"
        else:
            modifier = None
        genotypes.append(genotype)
        denovo_statuses.append(denovo)
        denovo_evidence.append(evidence)
        penetrance_values.append(penetrance)
        vital_routes.append(route)
        modifiers.append(modifier)
        true_causal.append(determine_true_causal(temp, penetrance, rng))

    selected["genotype"] = genotypes
    selected["de_novo_status"] = denovo_statuses
    selected["de_novo_evidence"] = denovo_evidence
    selected["confirmed_de_novo"] = selected["de_novo_evidence"].eq("trio_confirmed")
    selected["penetrance_draw"] = penetrance_values
    selected["true_causal_status"] = true_causal
    selected["true_disease_model"] = selected["mechanism_class"].map(
        {
            "GoF_DN": "dominant_high_penetrance",
            "HI": "dominant_reduced_or_moderate_penetrance",
            "AR": "recessive",
            "MIXED": "mixed_context_dependent",
            "UNKNOWN": "unknown",
        }
    )
    selected["vital_route"] = vital_routes
    selected["vital_modifier"] = modifiers

    selected["baseline_actionable"] = selected["ClinVar_label"].isin(["P", "LP"]) & selected["gene_actionable"]
    selected["baseline_causal_attribution"] = selected["baseline_actionable"]
    selected["baseline_confidence"] = np.where(
        selected["baseline_causal_attribution"] & selected["confirmed_de_novo"],
        "high",
        np.where(selected["baseline_causal_attribution"], "moderate", "none"),
    )
    selected["baseline_genetic_cause_of_death"] = selected["baseline_causal_attribution"] & selected["autopsy_status"].eq("NEGATIVE")
    selected["baseline_dominant_actionability"] = selected["baseline_actionable"]

    selected["vital_causal_attribution"] = selected.apply(
        lambda row: vital_causal_attribution(str(row["vital_route"]), row["vital_modifier"]),
        axis=1,
    )
    selected["vital_dominant_actionability"] = selected["vital_route"].eq("VITAL_OK") & selected["mechanism_class"].isin(["GoF_DN", "HI"])
    selected["baseline_false_causal_attribution"] = selected["baseline_causal_attribution"] & ~selected["true_causal_status"]
    selected["vital_false_causal_attribution"] = selected["vital_causal_attribution"] & ~selected["true_causal_status"]
    selected["false_dominant_diagnosis"] = selected["baseline_dominant_actionability"] & ~selected["true_disease_model"].eq("dominant_high_penetrance")
    selected["de_novo_override_error"] = (
        selected["confirmed_de_novo"]
        & selected["baseline_confidence"].eq("high")
        & selected["baseline_causal_attribution"]
        & ~selected["true_causal_status"]
    )
    selected["de_novo_override_model_error"] = (
        selected["confirmed_de_novo"]
        & selected["baseline_causal_attribution"]
        & (selected["effective_AF"].astype(float) > mcaf_strict)
    )
    selected["prevented_false_attribution"] = selected["baseline_false_causal_attribution"] & ~selected["vital_causal_attribution"]
    selected["vital_check_defer"] = selected["vital_route"].isin(["CHECK_MODEL", "CHECK_POPMAX", "EVAL_LIMITED"])
    selected["vital_model_conflict"] = selected["vital_route"].eq("MODEL_CONFLICT")
    selected["af_bin"] = pd.cut(selected["effective_AF"].astype(float), bins=AF_BINS, labels=AF_LABELS, right=False)
    return selected


def summarize_strategy(cases: pd.DataFrame, label: str) -> pd.DataFrame:
    total = len(cases)
    evaluable = int(cases["evaluability_tier"].eq("Tier1").sum())
    baseline_false = int(cases["baseline_false_causal_attribution"].sum())
    vital_false = int(cases["vital_false_causal_attribution"].sum())
    baseline_false_evaluable = int(
        (cases["baseline_false_causal_attribution"] & cases["evaluability_tier"].eq("Tier1")).sum()
    )
    vital_false_evaluable = int(
        (cases["vital_false_causal_attribution"] & cases["evaluability_tier"].eq("Tier1")).sum()
    )
    rows = [
        {
            "scope": label,
            "strategy": "label_driven_baseline",
            "total_cases": total,
            "evaluable_cases": evaluable,
            "causal_attributions": int(cases["baseline_causal_attribution"].sum()),
            "false_causal_attributions": baseline_false,
            "false_causal_attribution_rate": pct(baseline_false, int(cases["baseline_causal_attribution"].sum())),
            "per_case_false_attribution_rate": pct(baseline_false, total),
            "false_causal_attributions_evaluable": baseline_false_evaluable,
            "per_evaluable_case_false_attribution_rate": pct(baseline_false_evaluable, evaluable),
            "check_defer": 0,
            "model_conflict": 0,
            "prevented_false_attribution": 0,
            "de_novo_override_errors": int(cases["de_novo_override_error"].sum()),
        },
        {
            "scope": label,
            "strategy": "VITAL_routing",
            "total_cases": total,
            "evaluable_cases": evaluable,
            "causal_attributions": int(cases["vital_causal_attribution"].sum()),
            "false_causal_attributions": vital_false,
            "false_causal_attribution_rate": pct(vital_false, int(cases["vital_causal_attribution"].sum())),
            "per_case_false_attribution_rate": pct(vital_false, total),
            "false_causal_attributions_evaluable": vital_false_evaluable,
            "per_evaluable_case_false_attribution_rate": pct(vital_false_evaluable, evaluable),
            "check_defer": int(cases["vital_check_defer"].sum()),
            "model_conflict": int(cases["vital_model_conflict"].sum()),
            "prevented_false_attribution": int(cases["prevented_false_attribution"].sum()),
            "de_novo_override_errors": 0,
        },
    ]
    return pd.DataFrame(rows)


def summarize_false_attribution_rates(cases: pd.DataFrame) -> pd.DataFrame:
    route_order = ["EVAL_LIMITED", "CHECK_MODEL", "CHECK_POPMAX", "MODEL_CONFLICT"]
    scopes = [
        ("all_cases", cases.copy()),
        ("evaluable_cases", cases.loc[cases["evaluability_tier"].eq("Tier1")].copy()),
    ]
    rows: list[dict[str, object]] = []
    for scope, frame in scopes:
        denominator = len(frame)
        baseline_false = int(frame["baseline_false_causal_attribution"].sum())
        vital_false = int(frame["vital_false_causal_attribution"].sum())
        row: dict[str, object] = {
            "scope": scope,
            "denominator_cases": denominator,
            "baseline_false_causal_attributions": baseline_false,
            "baseline_false_attribution_rate_per_case": pct(baseline_false, denominator),
            "vital_false_causal_attributions": vital_false,
            "vital_false_attribution_rate_per_case": pct(vital_false, denominator),
            "prevented_false_attribution": int(frame["prevented_false_attribution"].sum()),
        }
        for route in route_order:
            count = int(frame["vital_route"].eq(route).sum())
            row[f"{route}_count"] = count
            row[f"{route}_percent_of_cases"] = pct(count, denominator)
        rows.append(row)
    return pd.DataFrame(rows)


def grouped_metrics(cases: pd.DataFrame, group_col: str) -> pd.DataFrame:
    rows = []
    for group, frame in cases.groupby(group_col, dropna=False, observed=False):
        total = len(frame)
        baseline_false = int(frame["baseline_false_causal_attribution"].sum())
        vital_false = int(frame["vital_false_causal_attribution"].sum())
        rows.append(
            {
                group_col: str(group),
                "case_count": total,
                "baseline_causal_attributions": int(frame["baseline_causal_attribution"].sum()),
                "baseline_false_causal_attributions": baseline_false,
                "baseline_false_causal_percent": pct(baseline_false, total),
                "vital_causal_attributions": int(frame["vital_causal_attribution"].sum()),
                "vital_false_causal_attributions": vital_false,
                "vital_false_causal_percent": pct(vital_false, total),
                "defer_percent": pct(int(frame["vital_check_defer"].sum()), total),
                "model_conflict_percent": pct(int(frame["vital_model_conflict"].sum()), total),
                "prevented_false_attribution": int(frame["prevented_false_attribution"].sum()),
                "de_novo_override_errors": int(frame["de_novo_override_error"].sum()),
            }
        )
    return pd.DataFrame(rows)


def denovo_override_audit(cases: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for status, frame in cases.groupby("de_novo_status", dropna=False):
        total = len(frame)
        rows.append(
            {
                "de_novo_status": status,
                "case_count": total,
                "baseline_high_confidence_causal": int(frame["baseline_confidence"].eq("high").sum()),
                "VITAL_supported": int(frame["vital_causal_attribution"].sum()),
                "VITAL_check_defer": int(frame["vital_check_defer"].sum()),
                "VITAL_blocked_model_conflict": int(frame["vital_model_conflict"].sum()),
                "de_novo_override_error": int(frame["de_novo_override_error"].sum()),
                "de_novo_override_model_error": int(frame["de_novo_override_model_error"].sum()),
            }
        )
    return pd.DataFrame(rows)


def run_sensitivities(design: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    penetrance_rows = []
    for profile in PENETRANCE_PROFILES:
        rng = np.random.default_rng(RNG_SEED + hash(profile) % 1000)
        cases = simulate_cases(design, N_CASES, rng, penetrance_profile=profile)
        summary = summarize_strategy(cases, f"penetrance_{profile}")
        summary["sensitivity"] = profile
        penetrance_rows.append(summary)

    denovo_rows = []
    for rate in [0.001, 0.01, 0.05, 0.10]:
        rng = np.random.default_rng(RNG_SEED + int(rate * 100_000))
        cases = simulate_cases(design, N_CASES, rng, confirmed_denovo_rate=rate)
        summary = summarize_strategy(cases, f"denovo_rate_{rate:g}")
        summary["confirmed_denovo_rate"] = rate
        denovo_rows.append(summary)

    mcaf_rows = []
    for strict, hard, label in [(2.5e-5, 1e-4, "strict_2p5e-5_hard_1e-4"), (1e-4, 1e-3, "boundary_1e-4_hard_1e-3"), (1e-3, 1e-2, "permissive_1e-3_hard_1e-2")]:
        rng = np.random.default_rng(RNG_SEED + int(strict * 1e8))
        cases = simulate_cases(design, N_CASES, rng, mcaf_strict=strict, mcaf_hard=hard)
        summary = summarize_strategy(cases, f"mcaf_{label}")
        summary["mcaf_strict"] = strict
        summary["mcaf_hard"] = hard
        mcaf_rows.append(summary)

    permutation_rows = []
    rng = np.random.default_rng(RNG_SEED + 99)
    cases = simulate_cases(design, N_CASES, rng, mechanism_permutation=True)
    perm = grouped_metrics(cases, "mechanism_class")
    perm["sensitivity"] = "mechanism_permuted_across_af_matrix"
    permutation_rows.append(perm)

    return (
        pd.concat(penetrance_rows, ignore_index=True),
        pd.concat(denovo_rows, ignore_index=True),
        pd.concat(mcaf_rows, ignore_index=True),
        pd.concat(permutation_rows, ignore_index=True),
    )


def build_gold_standard_preservation(design: pd.DataFrame) -> pd.DataFrame:
    positives = design.loc[
        design["mechanism_class"].eq("GoF_DN")
        & design["ClinVar_label"].isin(["P", "LP"])
        & design["evaluability_tier"].eq("Tier1")
        & (design["effective_AF"].astype(float) <= MCAF_STRICT_DEFAULT)
    ].copy()
    positives = positives.sort_values(["review_score", "allele_count"], ascending=False).head(42)
    if positives.empty:
        return pd.DataFrame(columns=["gold_standard_set", "variant_count", "VITAL_OK", "CHECK", "MODEL_CONFLICT", "false_interruption_rate"])
    positives["route"] = positives.apply(
        lambda row: infer_vital_route_from_row(row, MCAF_STRICT_DEFAULT, MCAF_HARD_DEFAULT),
        axis=1,
    )
    return pd.DataFrame(
        [
            {
                "gold_standard_set": "empirical_gof_dn_af_compatible_plp",
                "variant_count": len(positives),
                "VITAL_OK": int(positives["route"].eq("VITAL_OK").sum()),
                "CHECK": int(positives["route"].isin(["CHECK_POPMAX", "CHECK_MODEL", "EVAL_LIMITED"]).sum()),
                "MODEL_CONFLICT": int(positives["route"].eq("MODEL_CONFLICT").sum()),
                "false_interruption_rate": pct(int(positives["route"].eq("MODEL_CONFLICT").sum()), len(positives)),
            }
        ]
    )


def build_conflicting_layer(design: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    sample = design.sample(n=min(186, len(design)), replace=False, random_state=RNG_SEED).copy()
    sample["ClinVar_label"] = "conflicting"
    sample["route"] = sample.apply(
        lambda row: infer_vital_route_from_row(row, MCAF_STRICT_DEFAULT, MCAF_HARD_DEFAULT),
        axis=1,
    )
    # For conflicting labels, VITAL should not force a hard reclassification route through the P/LP gate.
    counts = sample["route"].value_counts().to_dict()
    return pd.DataFrame(
        [
            {
                "conflicting_variant_count": len(sample),
                "NO_PLP_BASELINE": int(counts.get("NO_PLP_BASELINE", 0)),
                "CHECK_or_DEFER": int(sum(counts.get(route, 0) for route in ["CHECK_MODEL", "CHECK_POPMAX", "EVAL_LIMITED"])),
                "MODEL_CONFLICT": int(counts.get("MODEL_CONFLICT", 0)),
                "interpretation": "Conflicting labels stay outside the P/LP baseline; enrichment belongs in review/defer, not automatic model-conflict.",
            }
        ]
    )


def draw_decision_flow() -> None:
    fig, ax = plt.subplots(figsize=(12, 5.5))
    ax.axis("off")
    ax.set_title("Autopsy x de novo counterfactual decision audit", fontsize=15, fontweight="bold", pad=18)
    boxes = [
        (0.04, 0.58, "P/LP label\narrhythmia gene", "#eef2f7"),
        (0.25, 0.58, "Negative autopsy\nphenotype-null", "#f3f4f6"),
        (0.47, 0.72, "Baseline\ncause assigned\n+ de novo = high confidence", "#f8d7da"),
        (0.47, 0.36, "VITAL routing\nevaluability -> popmax\nmechanism -> de novo modifier", "#dbeafe"),
        (0.74, 0.36, "Supported\nCheck/defer\nBlocked", "#dcfce7"),
    ]
    for x, y, text, color in boxes:
        patch = plt.Rectangle((x, y), 0.17, 0.18, facecolor=color, edgecolor="#334155", linewidth=1.2)
        ax.add_patch(patch)
        ax.text(x + 0.085, y + 0.09, text, ha="center", va="center", fontsize=10)
    for start, end in [((0.21, 0.67), (0.25, 0.67)), ((0.42, 0.67), (0.47, 0.79)), ((0.42, 0.67), (0.47, 0.45)), ((0.64, 0.45), (0.74, 0.45))]:
        ax.annotate("", xy=end, xytext=start, arrowprops={"arrowstyle": "->", "lw": 1.8, "color": "#334155"})
    ax.text(
        0.5,
        0.08,
        "de novo supports a coherent disease model; de novo does not repair an incoherent model.",
        ha="center",
        fontsize=11,
        color="#1f2937",
    )
    fig.savefig(FLOW_FIGURE_OUT, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {FLOW_FIGURE_OUT.relative_to(BASE_DIR)}")


def plot_false_attribution(summary: pd.DataFrame) -> None:
    plot_df = summary.loc[summary["scope"].eq("main_negative_autopsy")].copy()
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    labels = ["Baseline", "VITAL"]
    values = [
        float(plot_df.loc[plot_df["strategy"].eq("label_driven_baseline"), "false_causal_attributions"].iloc[0]),
        float(plot_df.loc[plot_df["strategy"].eq("VITAL_routing"), "false_causal_attributions"].iloc[0]),
    ]
    ax.bar(labels, values, color=["#b94a48", "#2f6f62"])
    ax.set_ylabel("False causal attributions")
    ax.set_title("False causal attribution reduction")
    for index, value in enumerate(values):
        ax.text(index, value + max(values) * 0.02, f"{int(value):,}", ha="center", fontweight="bold")
    ax.grid(axis="y", alpha=0.25)
    fig.savefig(FALSE_ATTRIBUTION_FIGURE_OUT, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {FALSE_ATTRIBUTION_FIGURE_OUT.relative_to(BASE_DIR)}")


def plot_denovo_override(denovo: pd.DataFrame) -> None:
    row_order = ["confirmed", "assumed", "unknown", "absent"]
    plot_df = denovo.set_index("de_novo_status").reindex(row_order).fillna(0)
    x = np.arange(len(plot_df))
    supported = plot_df["VITAL_supported"].to_numpy(dtype=float)
    check = plot_df["VITAL_check_defer"].to_numpy(dtype=float)
    blocked = plot_df["VITAL_blocked_model_conflict"].to_numpy(dtype=float)
    totals = plot_df["case_count"].replace(0, np.nan).to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    ax.bar(x, supported / totals * 100, label="Supported", color="#2f6f62")
    ax.bar(x, check / totals * 100, bottom=supported / totals * 100, label="Check/defer", color="#d18f2f")
    ax.bar(x, blocked / totals * 100, bottom=(supported + check) / totals * 100, label="Blocked", color="#b94a48")
    ax.set_xticks(x)
    ax.set_xticklabels(row_order)
    ax.set_ylabel("Percent of de novo stratum")
    ax.set_title("de novo signal is routed through model constraints")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    fig.savefig(DENOVO_FIGURE_OUT, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {DENOVO_FIGURE_OUT.relative_to(BASE_DIR)}")


def plot_mechanism(by_mechanism: pd.DataFrame) -> None:
    order = ["GoF_DN", "HI", "AR", "MIXED", "UNKNOWN"]
    plot_df = by_mechanism.set_index("mechanism_class").reindex(order).dropna(how="all")
    x = np.arange(len(plot_df))
    width = 0.34
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    ax.bar(x - width / 2, plot_df["baseline_false_causal_percent"], width=width, label="Baseline false causal %", color="#b94a48")
    ax.bar(x + width / 2, plot_df["vital_false_causal_percent"], width=width, label="VITAL false causal %", color="#2f6f62")
    ax.set_xticks(x)
    ax.set_xticklabels(plot_df.index)
    ax.set_ylabel("Percent of cases")
    ax.set_title("Mechanism stratification")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    fig.savefig(MECHANISM_FIGURE_OUT, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {MECHANISM_FIGURE_OUT.relative_to(BASE_DIR)}")


def main() -> None:
    rng = np.random.default_rng(RNG_SEED)
    model_inputs = build_model_inputs()
    save_table(model_inputs, MODEL_INPUTS_OUT)

    design = build_variant_design(rng)
    save_table(design, VARIANT_DESIGN_OUT)

    cases = simulate_cases(design, N_CASES, rng)
    save_table(cases, CASE_CALLS_OUT)

    summary = summarize_strategy(cases, "main_negative_autopsy")
    save_table(summary, SUMMARY_OUT)
    false_rates = summarize_false_attribution_rates(cases)
    save_table(false_rates, FALSE_ATTRIBUTION_RATES_OUT)

    by_gene = grouped_metrics(cases, "gene")
    by_mechanism = grouped_metrics(cases, "mechanism_class")
    by_af_bin = grouped_metrics(cases, "af_bin")
    by_denovo = denovo_override_audit(cases)
    by_autopsy = grouped_metrics(cases.assign(autopsy_context="negative_autopsy_phenotype_null"), "autopsy_context")

    save_table(by_gene, BY_GENE_OUT)
    save_table(by_mechanism, BY_MECHANISM_OUT)
    save_table(by_af_bin, BY_AF_BIN_OUT)
    save_table(by_denovo, BY_DENOVO_OUT)
    save_table(by_autopsy, BY_AUTOPSY_CONTEXT_OUT)

    penetrance, denovo_sens, mcaf_sens, mechanism_perm = run_sensitivities(design)
    save_table(penetrance, SENSITIVITY_PENETRANCE_OUT)
    save_table(denovo_sens, SENSITIVITY_DENOVO_OUT)
    save_table(mcaf_sens, SENSITIVITY_MCAF_OUT)
    save_table(mechanism_perm, MECHANISM_PERMUTATION_OUT)

    save_table(build_gold_standard_preservation(design), GOLD_STANDARD_OUT)
    save_table(build_conflicting_layer(design, rng), CONFLICTING_OUT)

    draw_decision_flow()
    plot_false_attribution(summary)
    plot_denovo_override(by_denovo)
    plot_mechanism(by_mechanism)


if __name__ == "__main__":
    main()
