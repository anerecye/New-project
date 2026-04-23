from __future__ import annotations

import argparse
import math
import time
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires pandas and numpy. Install them with: "
        "python -m pip install -r requirements.txt"
    ) from exc


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
SUPPLEMENT_DIR = BASE_DIR / "supplementary_tables"

ARRHYTHMIA_CURRENT = DATA_DIR / "arrhythmia_vital_scores.csv"
CROSS_DISEASE_CURRENT = DATA_DIR / "vital_cross_disease_3000_vital_scores.csv"
TEMPORAL_ROBUSTNESS = DATA_DIR / "arrhythmia_temporal_robustness_vital_scored_variants.tsv"

RECESSIVE_CONTEXT_GENES = {"CASQ2", "TRDN"}


def save_table(df: pd.DataFrame, output_path: Path, sep: str = ",") -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table = df.copy()
    if sep == "\t":
        object_columns = table.select_dtypes(include=["object", "string"]).columns
        table.loc[:, object_columns] = table.loc[:, object_columns].replace("", "NA")
    temp_path = output_path.with_name(f".{output_path.name}.tmp")
    last_error: PermissionError | None = None
    for attempt in range(1, 7):
        try:
            table.to_csv(temp_path, sep=sep, index=False, na_rep="NA", lineterminator="\n")
            temp_path.replace(output_path)
            print(f"Saved {output_path} ({len(df)} rows)")
            return
        except PermissionError as exc:
            last_error = exc
            delay = 1.5 * attempt
            print(f"Save retry {attempt}/6 for {output_path} in {delay:.1f}s: {exc}")
            time.sleep(delay)
    if last_error is not None:
        raise last_error


def truthy(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def load_scores(path: Path, sep: str = ",") -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing score table: {path}")
    df = pd.read_csv(path, sep=sep, low_memory=False)
    if "inheritance_flag" not in df.columns:
        df["inheritance_flag"] = np.where(
            df["gene"].isin(RECESSIVE_CONTEXT_GENES),
            "recessive_context_required",
            "standard",
        )
    for column in [
        "standard_acmg_frequency_flag",
        "frequency_signal_ac_ge_20",
        "weak_review_signal",
        "vital_red_flag",
    ]:
        if column not in df.columns:
            df[column] = False
        df[column] = truthy(df[column])
    for column in [
        "submitter_count",
        "global_af",
        "popmax_af",
        "max_frequency_signal",
        "qualifying_frequency_ac",
        "vital_score",
    ]:
        if column not in df.columns:
            df[column] = np.nan
        df[column] = pd.to_numeric(df[column], errors="coerce")
    df["submitter_exposure_units"] = df["submitter_count"].fillna(1).clip(lower=1)
    return df


def summarize_layer(
    df: pd.DataFrame,
    cohort: str,
    layer: str,
    mask: pd.Series,
    clinical_risk: str,
) -> dict[str, object]:
    sub = df.loc[mask].copy()
    denominator = len(df)
    observed = int(df["frequency_evidence_status"].fillna("").astype(str).eq("frequency_observed").sum())
    return {
        "cohort": cohort,
        "decision_risk_layer": layer,
        "records": len(sub),
        "percent_of_all_records": 100 * len(sub) / denominator if denominator else math.nan,
        "percent_of_frequency_observed": 100 * len(sub) / observed if observed else math.nan,
        "minimum_submitter_exposure_units": int(sub["submitter_exposure_units"].sum()),
        "standard_context_records": int(sub["inheritance_flag"].eq("standard").sum()),
        "recessive_context_required_records": int(
            sub["inheritance_flag"].eq("recessive_context_required").sum()
        ),
        "weak_review_records": int(sub["weak_review_signal"].sum()),
        "median_vital_score": float(sub["vital_score"].median()) if len(sub) else math.nan,
        "max_vital_score": float(sub["vital_score"].max()) if len(sub) else math.nan,
        "clinical_risk_proxy": clinical_risk,
    }


def cohort_summary(df: pd.DataFrame, cohort: str) -> pd.DataFrame:
    frequency_flag = df["standard_acmg_frequency_flag"]
    ac_supported = df["frequency_signal_ac_ge_20"]
    weak_frequency = frequency_flag & df["weak_review_signal"]
    red = df["vital_red_flag"]
    standard_red = red & df["inheritance_flag"].eq("standard")
    recessive_red = red & df["inheritance_flag"].eq("recessive_context_required")
    rows = [
        summarize_layer(
            df,
            cohort,
            "upper_bound_frequency_discordant_PLP_records",
            frequency_flag,
            "Could support a P/LP label despite ancestry-aware frequency contradiction if accepted uncritically.",
        ),
        summarize_layer(
            df,
            cohort,
            "AC_supported_frequency_discordant_records",
            ac_supported,
            "Frequency contradiction has allele-count support; suitable for active review rather than passive watchlist.",
        ),
        summarize_layer(
            df,
            cohort,
            "weak_review_frequency_discordant_records",
            weak_frequency,
            "Weak public assertions with frequency contradiction create higher diagnostic-closure and cascade-testing pressure.",
        ),
        summarize_layer(
            df,
            cohort,
            "VITAL_red_urgent_review_records",
            red,
            "Smallest urgent review queue: high score, weak review, and AC-supported frequency contradiction.",
        ),
        summarize_layer(
            df,
            cohort,
            "standard_context_red_records",
            standard_red,
            "Potential monogenic-disease label risk: diagnostic closure, cascade testing, or surveillance before adjudication.",
        ),
        summarize_layer(
            df,
            cohort,
            "recessive_context_red_records",
            recessive_red,
            "Carrier-architecture routing risk: inappropriate dominant-disease counseling if phase/inheritance is ignored.",
        ),
    ]
    return pd.DataFrame(rows)


def top_decision_risk_records(df: pd.DataFrame, cohort: str) -> pd.DataFrame:
    mask = df["vital_red_flag"] | (
        df["frequency_signal_ac_ge_20"] & df["weak_review_signal"] & df["vital_score"].ge(55)
    )
    columns = [
        "gene",
        "clinvar_id",
        "variation_id",
        "title",
        "review_strength",
        "inheritance_flag",
        "global_af",
        "global_ac",
        "popmax_af",
        "popmax_ac",
        "popmax_population",
        "qualifying_frequency_ac",
        "submitter_exposure_units",
        "vital_score",
        "vital_band",
        "vital_red_flag",
        "vital_signal_reason",
    ]
    result = df.loc[mask, [column for column in columns if column in df.columns]].copy()
    result.insert(0, "cohort", cohort)
    result["decision_risk_mode"] = np.select(
        [
            result["vital_red_flag"] & result["inheritance_flag"].eq("standard"),
            result["vital_red_flag"] & result["inheritance_flag"].eq("recessive_context_required"),
            result["inheritance_flag"].eq("recessive_context_required"),
        ],
        [
            "standard_context_urgent_review",
            "recessive_carrier_context_urgent_review",
            "recessive_carrier_context_watchlist",
        ],
        default="AC_supported_weak_review_watchlist",
    )
    return result.sort_values(["cohort", "vital_red_flag", "vital_score"], ascending=[True, False, False])


def temporal_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for snapshot, sub in df.groupby("snapshot", sort=False):
        summary = cohort_summary(sub, f"arrhythmia_{snapshot}")
        summary.insert(1, "snapshot", snapshot)
        rows.append(summary)
    return pd.concat(rows, ignore_index=True)


def run() -> dict[str, Path]:
    arrhythmia = load_scores(ARRHYTHMIA_CURRENT)
    cross_disease = load_scores(CROSS_DISEASE_CURRENT)
    temporal = load_scores(TEMPORAL_ROBUSTNESS, sep="\t")

    summary = pd.concat(
        [
            cohort_summary(arrhythmia, "arrhythmia_current"),
            cohort_summary(cross_disease, "cross_disease_3000_current"),
        ],
        ignore_index=True,
    )
    temporal_layers = temporal_summary(temporal)
    top_cases = pd.concat(
        [
            top_decision_risk_records(arrhythmia, "arrhythmia_current"),
            top_decision_risk_records(cross_disease, "cross_disease_3000_current"),
        ],
        ignore_index=True,
    )

    outputs = {
        "summary": DATA_DIR / "vital_clinical_decision_risk_proxy_summary.tsv",
        "temporal": DATA_DIR / "vital_clinical_decision_risk_proxy_temporal.tsv",
        "top_cases": DATA_DIR / "vital_clinical_decision_risk_proxy_top_cases.tsv",
        "supplement": SUPPLEMENT_DIR / "Supplementary_Table_S28_clinical_decision_risk_proxy.tsv",
    }
    save_table(summary, outputs["summary"], sep="\t")
    save_table(temporal_layers, outputs["temporal"], sep="\t")
    save_table(top_cases, outputs["top_cases"], sep="\t")
    save_table(summary, outputs["supplement"], sep="\t")

    print("\nClinical decision-risk proxy summary")
    print(summary.to_string(index=False))
    return outputs


def parse_args() -> argparse.Namespace:
    return argparse.ArgumentParser(
    description="Quantify public ClinVar P/LP decision-risk exposure using cached score layers."
    ).parse_args()


def main() -> None:
    parse_args()
    outputs = run()
    print("\nOutputs")
    for name, path in outputs.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
