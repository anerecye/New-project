from __future__ import annotations

from math import sqrt
from pathlib import Path

import pandas as pd


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
SUPPLEMENT_DIR = BASE_DIR / "supplementary_tables"

SCORES_IN = DATA_DIR / "arrhythmia_vital_scores.csv"

MAX_CREDIBLE_AF_OUT = DATA_DIR / "vital_max_credible_af_analysis.csv"
PREVALENCE_OUT = DATA_DIR / "vital_prevalence_constraint_analysis.csv"
BAYES_OUT = DATA_DIR / "vital_bayesian_acmg_frequency_tension.csv"
CARRIER_OUT = DATA_DIR / "vital_recessive_carrier_logic.csv"
AC_RELIABILITY_OUT = DATA_DIR / "vital_allele_count_reliability.csv"

MAX_CREDIBLE_AF_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S42_max_credible_af.tsv"
PREVALENCE_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S43_prevalence_constraint.tsv"
BAYES_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S44_bayesian_acmg_frequency_tension.tsv"
CARRIER_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S45_recessive_carrier_logic.tsv"
AC_RELIABILITY_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S46_allele_count_reliability.tsv"

PRIOR_PROBABILITY = 0.10
ACMG_LR = {
    "no_pathogenic_evidence": 1.0,
    "pathogenic_supporting": 2.08,
    "pathogenic_strong": 18.7,
    "pathogenic_very_strong": 350.0,
}


# Allelic contribution is the fraction of disease attributable to this allele,
# not gnomAD allele count. We spell this out in outputs to avoid AC ambiguity.
SCENARIOS = [
    {
        "scenario": "strict_high_penetrance_dominant",
        "applies_to": "all_red_cases",
        "disease_prevalence": 1 / 2000,
        "allelic_contribution_fraction": 0.10,
        "penetrance": 0.50,
        "model_factor": 2.0,
        "interpretation": "Stringent dominant model: rare disease, limited allelic contribution, moderate/high penetrance.",
    },
    {
        "scenario": "generous_dominant_low_penetrance",
        "applies_to": "all_red_cases",
        "disease_prevalence": 1 / 1000,
        "allelic_contribution_fraction": 0.20,
        "penetrance": 0.10,
        "model_factor": 2.0,
        "interpretation": "Generous dominant model allowing higher prevalence, higher allelic contribution, and reduced penetrance.",
    },
    {
        "scenario": "SCN5A_Brugada_high_penetrance_read",
        "applies_to": "VCV000440850",
        "disease_prevalence": 1 / 2000,
        "allelic_contribution_fraction": 0.10,
        "penetrance": 0.20,
        "model_factor": 2.0,
        "interpretation": "Illustrative Brugada-like Mendelian read for the SCN5A haplotype label.",
    },
    {
        "scenario": "KCNH2_LQTS_high_penetrance_read",
        "applies_to": "VCV004535537",
        "disease_prevalence": 1 / 2000,
        "allelic_contribution_fraction": 0.10,
        "penetrance": 0.50,
        "model_factor": 2.0,
        "interpretation": "Illustrative LQTS-like dominant high-penetrance read for KCNH2.",
    },
    {
        "scenario": "TRDN_CPVT_dominant_misread",
        "applies_to": "VCV001325231",
        "disease_prevalence": 1 / 10000,
        "allelic_contribution_fraction": 0.10,
        "penetrance": 0.50,
        "model_factor": 2.0,
        "interpretation": "Intentional dominant misread for a TRDN allele to show why carrier architecture matters.",
    },
]


def load_red_cases() -> pd.DataFrame:
    if not SCORES_IN.exists():
        raise FileNotFoundError(f"Missing input: {SCORES_IN}")
    df = pd.read_csv(SCORES_IN)
    for col in [
        "global_af",
        "global_ac",
        "global_an",
        "popmax_af",
        "popmax_ac",
        "qualifying_frequency_ac",
        "vital_score",
    ]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    red = df[df["vital_red_flag"].astype(str).str.lower().isin(["true", "1"])].copy()
    if red.empty:
        raise ValueError("No VITAL red-priority cases found.")
    return red.sort_values("vital_score", ascending=False)


def max_credible_af(prevalence: float, contribution: float, penetrance: float, model_factor: float) -> float:
    return (prevalence * contribution) / (penetrance * model_factor)


def required_penetrance(prevalence: float, contribution: float, observed_af: float, model_factor: float) -> float:
    return (prevalence * contribution) / (observed_af * model_factor)


def build_max_credible_af(red: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, variant in red.iterrows():
        for scenario in SCENARIOS:
            if scenario["applies_to"] not in {"all_red_cases", variant["clinvar_id"]}:
                continue
            threshold = max_credible_af(
                scenario["disease_prevalence"],
                scenario["allelic_contribution_fraction"],
                scenario["penetrance"],
                scenario["model_factor"],
            )
            observed = float(variant["popmax_af"])
            rows.append(
                {
                    "gene": variant["gene"],
                    "clinvar_id": variant["clinvar_id"],
                    "variant": variant["title"],
                    "scenario": scenario["scenario"],
                    "disease_prevalence": scenario["disease_prevalence"],
                    "allelic_contribution_fraction": scenario["allelic_contribution_fraction"],
                    "penetrance": scenario["penetrance"],
                    "model_factor": scenario["model_factor"],
                    "max_credible_af": threshold,
                    "observed_popmax_af": observed,
                    "popmax_population": variant["popmax_population"],
                    "observed_to_max_credible_af_ratio": observed / threshold if threshold else float("nan"),
                    "gnomad_popmax_allele_count": variant["popmax_ac"],
                    "gnomad_qualifying_allele_count": variant["qualifying_frequency_ac"],
                    "interpretation": scenario["interpretation"],
                    "constraint_call": constraint_call(observed / threshold if threshold else float("nan")),
                }
            )
    return pd.DataFrame(rows)


def constraint_call(ratio: float) -> str:
    if ratio >= 10:
        return "very_strong_frequency_incompatibility_under_this_model"
    if ratio > 1:
        return "BS1_like_frequency_incompatibility_under_this_model"
    return "frequency_not_above_model_ceiling_under_this_model"


def build_prevalence_constraint(red: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, variant in red.iterrows():
        observed = float(variant["popmax_af"])
        for prevalence, contribution, label in [
            (1 / 2000, 0.10, "strict_1_in_2000_10pct_contribution"),
            (1 / 1000, 0.20, "generous_1_in_1000_20pct_contribution"),
        ]:
            req_pen = required_penetrance(prevalence, contribution, observed, model_factor=2.0)
            rows.append(
                {
                    "gene": variant["gene"],
                    "clinvar_id": variant["clinvar_id"],
                    "scenario": label,
                    "observed_popmax_af": observed,
                    "disease_prevalence": prevalence,
                    "allelic_contribution_fraction": contribution,
                    "required_penetrance_for_observed_AF_to_fit_dominant_model": req_pen,
                    "expected_affected_prevalence_if_penetrance_10pct": 2 * observed * 0.10,
                    "expected_affected_prevalence_if_penetrance_50pct": 2 * observed * 0.50,
                    "logic": (
                        "For a rare dominant allele, expected affected prevalence is approximately "
                        "2 * AF * penetrance. If this exceeds disease prevalence times allelic "
                        "contribution, the high-penetrance Mendelian read is not credible."
                    ),
                }
            )
    return pd.DataFrame(rows)


def build_bayesian_table(max_af: pd.DataFrame) -> pd.DataFrame:
    disease_specific = max_af[
        max_af["scenario"].isin(
            [
                "SCN5A_Brugada_high_penetrance_read",
                "KCNH2_LQTS_high_penetrance_read",
                "TRDN_CPVT_dominant_misread",
            ]
        )
    ].copy()
    rows = []
    prior_odds = PRIOR_PROBABILITY / (1 - PRIOR_PROBABILITY)
    for _, row in disease_specific.iterrows():
        ratio = float(row["observed_to_max_credible_af_ratio"])
        benign_lr = 350.0 if ratio >= 10 else 18.7 if ratio > 1 else 1.0
        benign_weight = (
            "very_strong_frequency_conflict"
            if ratio >= 10
            else "strong_frequency_conflict"
            if ratio > 1
            else "no_frequency_conflict"
        )
        for evidence_label, pathogenic_lr in ACMG_LR.items():
            posterior_odds = prior_odds * pathogenic_lr / benign_lr
            posterior_probability = posterior_odds / (1 + posterior_odds)
            rows.append(
                {
                    "gene": row["gene"],
                    "clinvar_id": row["clinvar_id"],
                    "scenario": row["scenario"],
                    "frequency_conflict_weight": benign_weight,
                    "benign_frequency_likelihood_ratio_divisor": benign_lr,
                    "pathogenic_evidence_scenario": evidence_label,
                    "pathogenic_likelihood_ratio_multiplier": pathogenic_lr,
                    "prior_probability": PRIOR_PROBABILITY,
                    "posterior_probability_after_frequency_and_pathogenic_evidence": posterior_probability,
                    "interpretation": (
                        "Illustrative Tavtigian-style odds calculation; not a final ACMG classification."
                    ),
                }
            )
    return pd.DataFrame(rows)


def build_carrier_logic(red: pd.DataFrame) -> pd.DataFrame:
    rows = []
    trdn = red[red["gene"].eq("TRDN")]
    if trdn.empty:
        return pd.DataFrame()
    row = trdn.iloc[0]
    q = float(row["popmax_af"])
    global_q = float(row["global_af"])
    for q_label, q_value in [("popmax", q), ("global", global_q)]:
        rows.append(
            {
                "gene": row["gene"],
                "clinvar_id": row["clinvar_id"],
                "frequency_source": q_label,
                "allele_frequency_q": q_value,
                "heterozygote_carrier_frequency_approx_2q": 2 * q_value,
                "homozygote_frequency_q_squared": q_value**2,
                "compound_heterozygote_not_modeled": True,
                "dominant_expected_prevalence_if_penetrance_10pct": 2 * q_value * 0.10,
                "dominant_expected_prevalence_if_penetrance_50pct": 2 * q_value * 0.50,
                "interpretation": (
                    "Observed heterozygous frequency can be compatible with carrier-state biology "
                    "while remaining incompatible with a dominant high-penetrance read."
                ),
            }
        )
    return pd.DataFrame(rows)


def build_ac_reliability(red: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in red.iterrows():
        ac = float(row["qualifying_frequency_ac"])
        global_af = float(row["global_af"])
        global_an = float(row["global_an"])
        relative_se = 1 / sqrt(ac) if ac > 0 else float("nan")
        rows.append(
            {
                "gene": row["gene"],
                "clinvar_id": row["clinvar_id"],
                "gnomad_qualifying_allele_count": ac,
                "poisson_relative_standard_error_approx": relative_se,
                "approx_95pct_relative_margin": 1.96 * relative_se if ac > 0 else float("nan"),
                "global_af": global_af,
                "global_an": global_an,
                "binomial_global_af_standard_error": sqrt(global_af * (1 - global_af) / global_an)
                if global_an > 0
                else float("nan"),
                "interpretation": (
                    "Allele-count reliability is reported as a graded measurement; AC>=20 is an "
                    "operational floor, not a biological discontinuity."
                ),
            }
        )
    return pd.DataFrame(rows)


def save(df: pd.DataFrame, path: Path, sep: str = ",") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep=sep, index=False, lineterminator="\n")
    print(f"Saved {path.relative_to(BASE_DIR)} ({len(df)} rows)")


def main() -> None:
    red = load_red_cases()
    max_af = build_max_credible_af(red)
    prevalence = build_prevalence_constraint(red)
    bayes = build_bayesian_table(max_af)
    carrier = build_carrier_logic(red)
    ac_reliability = build_ac_reliability(red)

    save(max_af, MAX_CREDIBLE_AF_OUT)
    save(prevalence, PREVALENCE_OUT)
    save(bayes, BAYES_OUT)
    save(carrier, CARRIER_OUT)
    save(ac_reliability, AC_RELIABILITY_OUT)

    save(max_af, MAX_CREDIBLE_AF_SUPP, sep="\t")
    save(prevalence, PREVALENCE_SUPP, sep="\t")
    save(bayes, BAYES_SUPP, sep="\t")
    save(carrier, CARRIER_SUPP, sep="\t")
    save(ac_reliability, AC_RELIABILITY_SUPP, sep="\t")


if __name__ == "__main__":
    main()
