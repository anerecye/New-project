from __future__ import annotations

import math
from pathlib import Path

try:
    import pandas as pd
    from scipy.stats import chi2_contingency
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires pandas and scipy. Install dependencies with: "
        "python -m pip install -r requirements.txt"
    ) from exc


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
SUPPLEMENT_DIR = BASE_DIR / "supplementary_tables"

SCORES = DATA_DIR / "arrhythmia_vital_scores.csv"
SUMMARY = DATA_DIR / "vital_frequency_evaluability_audit_summary.csv"
DIMENSION_AUDIT = DATA_DIR / "vital_frequency_evaluability_dimension_audit.csv"
SUPPLEMENT_TABLE = SUPPLEMENT_DIR / "Supplementary_Table_S30_frequency_evaluability_gap.tsv"

RECESSIVE_CONTEXT_GENES = {"CASQ2", "TRDN"}


def truthy(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def pct(numerator: float, denominator: float) -> float:
    if denominator == 0:
        return math.nan
    return 100 * numerator / denominator


def summarize_signal_subset(df: pd.DataFrame, label: str, interpretation: str) -> dict[str, object]:
    denominator = len(df)
    frequency_flag = truthy(df["standard_acmg_frequency_flag"])
    ac_supported = truthy(df["frequency_signal_ac_ge_20"])
    red = truthy(df["vital_red_flag"])
    global_flag = pd.to_numeric(df["global_af"], errors="coerce").gt(1e-5)
    popmax_only = truthy(df["popmax_only_frequency_flag"])
    return {
        "audit_layer": label,
        "variants": denominator,
        "percent_of_total_arrhythmia_space": math.nan,
        "af_observed_variants": denominator,
        "non_evaluable_variants": 0,
        "frequency_discordant_variants": int(frequency_flag.sum()),
        "frequency_discordant_percent": pct(frequency_flag.sum(), denominator),
        "global_af_gt_1e_5_variants": int(global_flag.sum()),
        "popmax_only_alerts": int(popmax_only.sum()),
        "ac_supported_frequency_variants": int(ac_supported.sum()),
        "vital_red_variants": int(red.sum()),
        "interpretation": interpretation,
    }


def build_summary(df: pd.DataFrame) -> pd.DataFrame:
    total = len(df)
    af_observed = df["frequency_evidence_status"].eq("frequency_observed")
    non_evaluable = ~af_observed
    rows: list[dict[str, object]] = [
        {
            "audit_layer": "full_arrhythmia_clinvar_plp_space",
            "variants": total,
            "percent_of_total_arrhythmia_space": 100.0,
            "af_observed_variants": int(af_observed.sum()),
            "non_evaluable_variants": int(non_evaluable.sum()),
            "frequency_discordant_variants": int(truthy(df["standard_acmg_frequency_flag"]).sum()),
            "frequency_discordant_percent": pct(
                truthy(df["standard_acmg_frequency_flag"]).sum(), int(af_observed.sum())
            ),
            "global_af_gt_1e_5_variants": int(pd.to_numeric(df["global_af"], errors="coerce").gt(1e-5).sum()),
            "popmax_only_alerts": int(truthy(df["popmax_only_frequency_flag"]).sum()),
            "ac_supported_frequency_variants": int(truthy(df["frequency_signal_ac_ge_20"]).sum()),
            "vital_red_variants": int(truthy(df["vital_red_flag"]).sum()),
            "interpretation": (
                "Most public P/LP assertions are not evaluable by exact AF matching; "
                "frequency conclusions are restricted to the AF-observed subset."
            ),
        },
        {
            "audit_layer": "af_observed_evaluable_space",
            "variants": int(af_observed.sum()),
            "percent_of_total_arrhythmia_space": pct(af_observed.sum(), total),
            "af_observed_variants": int(af_observed.sum()),
            "non_evaluable_variants": 0,
            "frequency_discordant_variants": int(
                truthy(df.loc[af_observed, "standard_acmg_frequency_flag"]).sum()
            ),
            "frequency_discordant_percent": pct(
                truthy(df.loc[af_observed, "standard_acmg_frequency_flag"]).sum(), af_observed.sum()
            ),
            "global_af_gt_1e_5_variants": int(
                pd.to_numeric(df.loc[af_observed, "global_af"], errors="coerce").gt(1e-5).sum()
            ),
            "popmax_only_alerts": int(truthy(df.loc[af_observed, "popmax_only_frequency_flag"]).sum()),
            "ac_supported_frequency_variants": int(
                truthy(df.loc[af_observed, "frequency_signal_ac_ge_20"]).sum()
            ),
            "vital_red_variants": int(truthy(df.loc[af_observed, "vital_red_flag"]).sum()),
            "interpretation": (
                "This is the only space used for frequency-discordance conclusions."
            ),
        },
        {
            "audit_layer": "non_evaluable_interpretability_gap",
            "variants": int(non_evaluable.sum()),
            "percent_of_total_arrhythmia_space": pct(non_evaluable.sum(), total),
            "af_observed_variants": 0,
            "non_evaluable_variants": int(non_evaluable.sum()),
            "frequency_discordant_variants": 0,
            "frequency_discordant_percent": math.nan,
            "global_af_gt_1e_5_variants": 0,
            "popmax_only_alerts": 0,
            "ac_supported_frequency_variants": 0,
            "vital_red_variants": 0,
            "interpretation": (
                "This majority is a separate representation and interpretability gap, "
                "not evidence of frequency consistency."
            ),
        },
    ]

    snv_observed = df.loc[af_observed & df["variant_type"].eq("SNV")]
    rows.append(
        summarize_signal_subset(
            snv_observed,
            "SNV_only_AF_observed_sensitivity",
            "Popmax frequency tension persists after restricting to exact-observed SNVs.",
        )
    )
    non_recessive_snv_observed = snv_observed.loc[~snv_observed["gene"].isin(RECESSIVE_CONTEXT_GENES)]
    rows.append(
        summarize_signal_subset(
            non_recessive_snv_observed,
            "non_recessive_SNV_AF_observed_sensitivity",
            "Popmax tension persists after removing CASQ2/TRDN and restricting to SNVs.",
        )
    )
    return pd.DataFrame(rows)


def dimension_audit(df: pd.DataFrame) -> pd.DataFrame:
    af_observed = df["frequency_evidence_status"].eq("frequency_observed")
    records: list[pd.DataFrame] = []
    for dimension in ["gene", "variant_type", "review_strength"]:
        table = pd.crosstab(df[dimension], af_observed)
        table = table.rename(columns={False: "non_evaluable", True: "af_observed"})
        for column in ["non_evaluable", "af_observed"]:
            if column not in table.columns:
                table[column] = 0
        table["total"] = table["non_evaluable"] + table["af_observed"]
        chi2_p = math.nan
        if table.shape[0] > 1 and table[["non_evaluable", "af_observed"]].to_numpy().sum() > 0:
            _, chi2_p, _, _ = chi2_contingency(table[["non_evaluable", "af_observed"]])
        out = table.reset_index().rename(columns={dimension: "category"})
        out.insert(0, "dimension", dimension)
        out["af_observed_percent_within_category"] = out.apply(
            lambda row: pct(row["af_observed"], row["total"]), axis=1
        )
        out["share_of_af_observed_space_percent"] = out["af_observed"].apply(
            lambda value: pct(value, af_observed.sum())
        )
        out["share_of_non_evaluable_space_percent"] = out["non_evaluable"].apply(
            lambda value: pct(value, (~af_observed).sum())
        )
        out["chi_square_p_for_dimension"] = chi2_p
        out["audit_interpretation"] = (
            "AF-observed and non-evaluable spaces differ systematically; "
            "this is a measured representation bias, not random missingness."
        )
        records.append(out)
    return pd.concat(records, ignore_index=True).sort_values(
        ["dimension", "total"], ascending=[True, False]
    )


def save_table(df: pd.DataFrame, path: Path, sep: str = ",") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep=sep, index=False, lineterminator="\n")
    print(f"Saved {path} ({len(df)} rows)")


def main() -> None:
    if not SCORES.exists():
        raise FileNotFoundError(f"Missing score table: {SCORES}")
    df = pd.read_csv(SCORES, low_memory=False)
    summary = build_summary(df)
    dimensions = dimension_audit(df)
    save_table(summary, SUMMARY)
    save_table(dimensions, DIMENSION_AUDIT)
    save_table(dimensions, SUPPLEMENT_TABLE, sep="\t")


if __name__ == "__main__":
    main()
