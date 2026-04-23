from __future__ import annotations

import math
from pathlib import Path

import pandas as pd


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
SUPPLEMENT_DIR = BASE_DIR / "supplementary_tables"

SCORES_IN = DATA_DIR / "arrhythmia_vital_scores.csv"

ACTION_CONTEXT_OUT = DATA_DIR / "vital_clinical_action_context_summary.csv"
REGIME_OUT = DATA_DIR / "vital_frequency_tension_regime_distribution.csv"

ACTION_CONTEXT_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S52_clinical_action_context_summary.tsv"
REGIME_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S53_frequency_tension_regime_distribution.tsv"

RECESSIVE_CONTEXT_GENES = {"CASQ2", "TRDN"}

ACTION_CONTEXTS: list[tuple[str, list[str], str]] = [
    (
        "guideline_driven_core",
        ["SCN5A", "KCNH2", "KCNQ1", "RYR2", "CACNA1C", "KCNE1", "KCNJ2", "CASQ2", "TRDN"],
        "Investigator-defined core arrhythmia genes where public P/LP labels can directly influence syndrome diagnosis, cascade testing, and management.",
    ),
    (
        "drug_restriction_context",
        ["SCN5A", "KCNH2", "KCNQ1", "KCNE1", "CACNA1C", "KCNJ2"],
        "Genes in which a pathogenic arrhythmia label can feed into medication avoidance or drug-restriction counseling.",
    ),
    (
        "device_or_intensive_surveillance_context",
        ["SCN5A", "KCNH2", "KCNQ1", "RYR2", "CACNA1C"],
        "Genes in which a pathogenic arrhythmia label can influence intensive surveillance or device-risk discussions.",
    ),
    (
        "cascade_testing_context",
        ["SCN5A", "KCNH2", "KCNQ1", "RYR2", "CACNA1C", "KCNE1", "KCNJ2", "CASQ2", "TRDN"],
        "Genes in which a public P/LP label can plausibly trigger family cascade testing or targeted familial follow-up.",
    ),
]


def truthy(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def wilson_ci(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total == 0:
        return math.nan, math.nan
    phat = successes / total
    denom = 1 + z * z / total
    center = (phat + z * z / (2 * total)) / denom
    margin = (
        z
        * math.sqrt((phat * (1 - phat) / total) + (z * z / (4 * total * total)))
        / denom
    )
    return 100 * max(0.0, center - margin), 100 * min(1.0, center + margin)


def save_table(df: pd.DataFrame, path: Path, sep: str = ",") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep=sep, index=False, lineterminator="\n")
    print(f"Saved {path.relative_to(BASE_DIR)} ({len(df)} rows)")


def load_scores() -> pd.DataFrame:
    df = pd.read_csv(SCORES_IN, low_memory=False)
    for column in [
        "standard_acmg_frequency_flag",
        "frequency_signal_ac_ge_20",
        "weak_review_signal",
        "vital_red_flag",
        "popmax_only_frequency_flag",
    ]:
        df[column] = truthy(df[column])
    for column in ["global_af", "max_frequency_signal", "vital_score"]:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    return df


def build_action_context_summary(freq: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    total = len(freq)
    for context_name, genes, interpretation in ACTION_CONTEXTS:
        mask = freq["gene"].isin(genes)
        n = int(mask.sum())
        ci_low, ci_high = wilson_ci(n, total)
        rows.append(
            {
                "clinical_context": context_name,
                "genes": "|".join(genes),
                "n_frequency_tension": n,
                "percent_of_tension": 100 * n / total if total else math.nan,
                "ci95_low_percent": ci_low,
                "ci95_high_percent": ci_high,
                "ac_supported_n": int((mask & freq["frequency_signal_ac_ge_20"]).sum()),
                "weak_review_n": int((mask & freq["weak_review_signal"]).sum()),
                "red_n": int((mask & freq["vital_red_flag"]).sum()),
                "popmax_only_n": int((mask & freq["popmax_only_frequency_flag"]).sum()),
                "clinical_interpretation": interpretation,
            }
        )
    return pd.DataFrame(rows)


def build_regime_distribution(freq: pd.DataFrame) -> pd.DataFrame:
    standard_context = ~freq["gene"].isin(RECESSIVE_CONTEXT_GENES)
    recessive_context = freq["gene"].isin(RECESSIVE_CONTEXT_GENES)

    hard = (
        standard_context
        & freq["frequency_signal_ac_ge_20"]
        & freq["global_af"].gt(1e-4)
    )
    boundary = standard_context & ~hard
    recessive = recessive_context

    regime_specs = [
        (
            "hard_dominant_incompatibility_candidate",
            hard,
            "Standard-context frequency-tension variants with AC support and global AF >1e-4; operational proxy for strongest unqualified dominant-model incompatibility.",
        ),
        (
            "boundary_or_monitoring_tension",
            boundary,
            "Standard-context frequency-tension variants without the hard-incompatibility signature; includes popmax-driven boundary cases and monitoring-priority signals.",
        ),
        (
            "recessive_carrier_compatible_tension",
            recessive,
            "Frequency-tension variants in CASQ2/TRDN where dominant over-reading competes with recessive/carrier-compatible architecture.",
        ),
    ]

    rows: list[dict[str, object]] = []
    total = len(freq)
    for regime_name, mask, definition in regime_specs:
        sub = freq.loc[mask].copy()
        n = len(sub)
        ci_low, ci_high = wilson_ci(n, total)
        rows.append(
            {
                "regime": regime_name,
                "n": n,
                "percent_of_tension": 100 * n / total if total else math.nan,
                "ci95_low_percent": ci_low,
                "ci95_high_percent": ci_high,
                "ac_ge_20_n": int(sub["frequency_signal_ac_ge_20"].sum()),
                "weak_review_n": int(sub["weak_review_signal"].sum()),
                "red_n": int(sub["vital_red_flag"].sum()),
                "popmax_only_n": int(sub["popmax_only_frequency_flag"].sum()),
                "af_gt_1e_4_n": int(sub["max_frequency_signal"].gt(1e-4).sum()),
                "genes": "|".join(sorted(sub["gene"].dropna().unique())),
                "operational_definition": definition,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    scores = load_scores()
    freq = scores.loc[scores["standard_acmg_frequency_flag"]].copy()
    if len(freq) == 0:
        raise SystemExit("No frequency-tension variants found in arrhythmia_vital_scores.csv")

    action_summary = build_action_context_summary(freq)
    regime_summary = build_regime_distribution(freq)

    save_table(action_summary, ACTION_CONTEXT_OUT)
    save_table(regime_summary, REGIME_OUT)
    save_table(action_summary, ACTION_CONTEXT_SUPP, sep="\t")
    save_table(regime_summary, REGIME_SUPP, sep="\t")

    print("\nClinical action contexts")
    print(action_summary.to_string(index=False))
    print("\nOperational regime distribution")
    print(regime_summary.to_string(index=False))


if __name__ == "__main__":
    main()
