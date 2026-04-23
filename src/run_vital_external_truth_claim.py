from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
FIGURE_DIR = BASE_DIR / "figures"
SUPPLEMENTARY_DIR = BASE_DIR / "supplementary_tables"

ARRHYTHMIA_SCORES = DATA_DIR / "arrhythmia_vital_scores.csv"
ARRHYTHMIA_LOF_SUMMARY = DATA_DIR / "arrhythmia_vital_lof_subtype_discordance_summary.csv"
ARRHYTHMIA_CONTEXT = DATA_DIR / "arrhythmia_gene_biological_context.csv"
CROSS_DISEASE_SCORES = DATA_DIR / "vital_cross_disease_3000_vital_scores.csv"

CANONICAL_RECESSIVE_OR_BIALLELIC_ARRHYTHMIA_GENES = {"CASQ2", "TRDN"}


def bool_series(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def pct(numerator: int | float, denominator: int | float) -> float:
    return float(100 * numerator / denominator) if denominator else np.nan


def summarize_expert_panel_truth(scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    table = scores.copy()
    for column in ["vital_score", "max_frequency_signal", "qualifying_frequency_ac"]:
        table[column] = pd.to_numeric(table.get(column), errors="coerce")
    table["naive_af_flag"] = table["max_frequency_signal"].fillna(0).gt(1e-5) & table[
        "frequency_evidence_status"
    ].eq("frequency_observed")
    table["expert_panel_truth"] = table["review_strength"].eq("expert_panel")
    table["red_priority"] = bool_series(table["vital_red_flag"])
    table["score_ge70"] = table["vital_score"].ge(70)
    table["ac_supported_frequency"] = bool_series(table["frequency_signal_ac_ge_20"])
    table["severe_annotation"] = table["functional_class"].isin(["LOF", "splice_or_intronic"])

    rows = []
    for label, sub in [
        ("expert_panel_reviewed_PLP", table[table["expert_panel_truth"]]),
        ("non_expert_or_non_panel_PLP", table[~table["expert_panel_truth"]]),
        ("all_cross_disease_PLP", table),
    ]:
        rows.append(
            {
                "truth_comparator_group": label,
                "n": len(sub),
                "severe_annotation_count": int(sub["severe_annotation"].sum()),
                "severe_annotation_percent": pct(sub["severe_annotation"].sum(), len(sub)),
                "naive_af_flags": int(sub["naive_af_flag"].sum()),
                "naive_af_flag_percent": pct(sub["naive_af_flag"].sum(), len(sub)),
                "ac_supported_frequency_flags": int(
                    (sub["naive_af_flag"] & sub["ac_supported_frequency"]).sum()
                ),
                "ac_supported_frequency_percent": pct(
                    (sub["naive_af_flag"] & sub["ac_supported_frequency"]).sum(),
                    len(sub),
                ),
                "score_ge70_count": int(sub["score_ge70"].sum()),
                "score_ge70_percent": pct(sub["score_ge70"].sum(), len(sub)),
                "red_priority_count": int(sub["red_priority"].sum()),
                "red_priority_percent": pct(sub["red_priority"].sum(), len(sub)),
                "severe_naive_af_flags": int((sub["severe_annotation"] & sub["naive_af_flag"]).sum()),
                "severe_naive_af_flag_percent_of_severe": pct(
                    (sub["severe_annotation"] & sub["naive_af_flag"]).sum(),
                    sub["severe_annotation"].sum(),
                ),
                "severe_ac_supported_frequency_flags": int(
                    (sub["severe_annotation"] & sub["naive_af_flag"] & sub["ac_supported_frequency"]).sum()
                ),
                "severe_red_priority_count": int((sub["severe_annotation"] & sub["red_priority"]).sum()),
                "orange_high_tension_count": int(sub["vital_band"].eq("orange_high_tension").sum()),
                "max_vital_score": float(sub["vital_score"].max()) if len(sub) else np.nan,
                "interpretation": (
                    "Curated expert-panel P/LP assertions are treated as an external specificity comparator, "
                    "not as ClinVar churn. The review-routing layer should not convert high-frequency curated exceptions into red calls."
                    if label == "expert_panel_reviewed_PLP"
                    else "Comparator group."
                ),
            }
        )

    expert = table[table["expert_panel_truth"]].copy()
    top_expert = expert.sort_values(["vital_score", "max_frequency_signal"], ascending=[False, False]).head(12)
    top_cols = [
        "gene",
        "clinvar_id",
        "variation_id",
        "title",
        "functional_class",
        "variant_type",
        "review_status",
        "vital_score",
        "vital_band",
        "naive_af_flag",
        "ac_supported_frequency",
        "global_af",
        "popmax_af",
        "qualifying_frequency_ac",
        "known_founder_or_carrier_context_gene",
    ]
    return pd.DataFrame(rows), top_expert.loc[:, [c for c in top_cols if c in top_expert.columns]]


def prepare_observed_scores(scores: pd.DataFrame, context: pd.DataFrame) -> pd.DataFrame:
    observed = scores[scores["frequency_evidence_status"].eq("frequency_observed")].copy()
    context_cols = ["gene", "gnomad_lof_oe_ci_upper"]
    observed = observed.merge(context.loc[:, context_cols], on="gene", how="left")
    observed["naive_af_flag"] = pd.to_numeric(observed["max_frequency_signal"], errors="coerce").fillna(0).gt(1e-5)
    observed["ac_supported_frequency"] = bool_series(observed["frequency_signal_ac_ge_20"])
    observed["red_priority"] = bool_series(observed["vital_red_flag"])
    observed["severe_annotation"] = observed["functional_class"].isin(["LOF", "splice_or_intronic"])
    observed["missense_annotation"] = observed["functional_class"].eq("missense")
    observed["canonical_recessive_or_biallelic_gene"] = observed["gene"].isin(
        CANONICAL_RECESSIVE_OR_BIALLELIC_ARRHYTHMIA_GENES
    )
    observed["constrained_loeuf_lt_0_5"] = pd.to_numeric(
        observed["gnomad_lof_oe_ci_upper"], errors="coerce"
    ).lt(0.5)
    observed["constrained_loeuf_lt_0_6"] = pd.to_numeric(
        observed["gnomad_lof_oe_ci_upper"], errors="coerce"
    ).lt(0.6)
    return observed


def severe_annotation_summary(scores: pd.DataFrame, lof_summary: pd.DataFrame, context: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    observed = prepare_observed_scores(scores, context)

    rows = []
    for label, sub in [
        ("all_frequency_observed_arrhythmia", observed),
        ("severe_LOF_or_splice_annotation", observed[observed["severe_annotation"]]),
        (
            "severe_excluding_CASQ2_TRDN",
            observed[observed["severe_annotation"] & ~observed["canonical_recessive_or_biallelic_gene"]],
        ),
        (
            "severe_in_constrained_genes_LOEUF_lt_0_5",
            observed[observed["severe_annotation"] & observed["constrained_loeuf_lt_0_5"]],
        ),
        (
            "severe_in_constrained_genes_LOEUF_lt_0_6",
            observed[observed["severe_annotation"] & observed["constrained_loeuf_lt_0_6"]],
        ),
        ("missense_annotation", observed[observed["missense_annotation"]]),
        ("non_severe_or_missense_other", observed[~observed["severe_annotation"]]),
    ]:
        discordant = sub[sub["naive_af_flag"]]
        rows.append(
            {
                "group": label,
                "n": len(sub),
                "naive_af_flags": int(sub["naive_af_flag"].sum()),
                "naive_af_flag_percent": pct(sub["naive_af_flag"].sum(), len(sub)),
                "ac_supported_frequency_flags": int((sub["naive_af_flag"] & sub["ac_supported_frequency"]).sum()),
                "ac_supported_frequency_percent": pct(
                    (sub["naive_af_flag"] & sub["ac_supported_frequency"]).sum(), len(sub)
                ),
                "red_priority_count": int(sub["red_priority"].sum()),
                "red_priority_percent": pct(sub["red_priority"].sum(), len(sub)),
                "discordant_gene_count": int(discordant["gene"].nunique()),
                "discordant_genes": "|".join(sorted(discordant["gene"].dropna().astype(str).unique())),
                "claim": (
                    "Severe HGVS consequence is not sufficient evidence of high-penetrance pathogenicity when "
                    "population frequency is discordant; this is the severe-annotation frequency-discordance pattern."
                    if label == "severe_LOF_or_splice_annotation"
                    else ""
                ),
            }
        )

    subtype = lof_summary.copy()
    subtype["group"] = "lof_subtype:" + subtype["lof_subtype"].astype(str)
    subtype = subtype.rename(
        columns={
            "lof_variant_count": "n",
            "naive_af_flag_count": "naive_af_flags",
            "ac_supported_frequency_count": "ac_supported_frequency_flags",
            "vital_red_count": "red_priority_count",
        }
    )
    subtype["claim"] = "LOF subtype-specific severe-annotation frequency discordance audit"
    subtype_cols = [
        "group",
        "n",
        "naive_af_flags",
        "naive_af_flag_percent",
        "ac_supported_frequency_flags",
        "ac_supported_frequency_percent",
        "red_priority_count",
        "claim",
    ]
    subtype = subtype.loc[:, [c for c in subtype_cols if c in subtype.columns]].copy()
    if "vital_red_count" in subtype.columns:
        subtype = subtype.drop(columns=["vital_red_count"])
    subtype["red_priority_percent"] = subtype.apply(lambda row: pct(row["red_priority_count"], row["n"]), axis=1)
    combined = pd.concat([pd.DataFrame(rows), subtype], ignore_index=True, sort=False)
    gene_spread = (
        observed[observed["severe_annotation"] & observed["naive_af_flag"]]
        .groupby("gene", dropna=False)
        .agg(
            severe_discordant_count=("variation_id", "count"),
            ac_supported_severe_discordant_count=("ac_supported_frequency", "sum"),
            red_priority_count=("red_priority", "sum"),
            median_vital_score=("vital_score", "median"),
            max_vital_score=("vital_score", "max"),
            loeuf_ci_upper=("gnomad_lof_oe_ci_upper", "first"),
            canonical_recessive_or_biallelic_gene=("canonical_recessive_or_biallelic_gene", "first"),
        )
        .sort_values(["severe_discordant_count", "max_vital_score"], ascending=[False, False])
        .reset_index()
    )
    return combined, gene_spread


def mechanism_triage_summary(scores: pd.DataFrame, context: pd.DataFrame) -> pd.DataFrame:
    observed = prepare_observed_scores(scores, context)
    severe_discordant = observed[observed["severe_annotation"] & observed["naive_af_flag"]].copy()
    non_recessive = ~severe_discordant["canonical_recessive_or_biallelic_gene"]
    ac_supported = severe_discordant["ac_supported_frequency"]
    constrained = severe_discordant["constrained_loeuf_lt_0_5"]

    categories = [
        (
            "carrier-compatible recessive/biallelic architecture",
            severe_discordant["canonical_recessive_or_biallelic_gene"],
            "Frequency can be compatible with heterozygous carrier state; not evidence of misclassification by itself.",
        ),
        (
            "non-recessive AC-supported high-penetrance tension",
            non_recessive & ac_supported,
            "Highest-priority mechanism adjudication: population frequency has allele-count support outside canonical recessive genes.",
        ),
        (
            "non-recessive constrained-gene low-AC surveillance",
            non_recessive & ~ac_supported & constrained,
            "Biologically uncomfortable direction in constrained genes, but AC support is insufficient for urgent action.",
        ),
        (
            "non-recessive other low-AC or architecture-unresolved surveillance",
            non_recessive & ~ac_supported & ~constrained,
            "Frequency tension remains visible, but current public data support monitoring or routine mechanism review rather than urgent downgrade.",
        ),
    ]

    rows = []
    total = len(severe_discordant)
    for label, mask, interpretation in categories:
        sub = severe_discordant[mask]
        rows.append(
            {
                "mechanism_triage_class": label,
                "n": len(sub),
                "percent_of_severe_discordant": pct(len(sub), total),
                "genes": "|".join(sorted(sub["gene"].dropna().astype(str).unique())),
                "ac_supported_count": int(sub["ac_supported_frequency"].sum()),
                "red_priority_count": int(sub["red_priority"].sum()),
                "median_vital_score": float(sub["vital_score"].median()) if len(sub) else np.nan,
                "max_vital_score": float(sub["vital_score"].max()) if len(sub) else np.nan,
                "interpretation": interpretation,
            }
        )

    rows.append(
        {
            "mechanism_triage_class": "all severe-annotation frequency-discordant assertions",
            "n": total,
            "percent_of_severe_discordant": 100.0 if total else np.nan,
            "genes": "|".join(sorted(severe_discordant["gene"].dropna().astype(str).unique())),
            "ac_supported_count": int(severe_discordant["ac_supported_frequency"].sum()),
            "red_priority_count": int(severe_discordant["red_priority"].sum()),
            "median_vital_score": float(severe_discordant["vital_score"].median()) if total else np.nan,
            "max_vital_score": float(severe_discordant["vital_score"].max()) if total else np.nan,
            "interpretation": "This is the tension universe, not an error-rate estimate.",
        }
    )
    return pd.DataFrame(rows)


def fisher_claim_tests(severe_summary: pd.DataFrame, truth_summary: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    observed = scores[scores["frequency_evidence_status"].eq("frequency_observed")].copy()
    observed = observed.merge(
        pd.read_csv(ARRHYTHMIA_CONTEXT).loc[:, ["gene", "gnomad_lof_oe_ci_upper"]],
        on="gene",
        how="left",
    )
    observed["naive_af_flag"] = pd.to_numeric(observed["max_frequency_signal"], errors="coerce").fillna(0).gt(1e-5)
    observed["severe_annotation"] = observed["functional_class"].isin(["LOF", "splice_or_intronic"])
    observed["missense_annotation"] = observed["functional_class"].eq("missense")
    observed["canonical_recessive_or_biallelic_gene"] = observed["gene"].isin(
        CANONICAL_RECESSIVE_OR_BIALLELIC_ARRHYTHMIA_GENES
    )
    observed["constrained_loeuf_lt_0_5"] = pd.to_numeric(
        observed["gnomad_lof_oe_ci_upper"], errors="coerce"
    ).lt(0.5)

    severe = observed["severe_annotation"]
    missense = observed["missense_annotation"]
    naive = observed["naive_af_flag"]
    table = [
        [int((severe & naive).sum()), int((severe & ~naive).sum())],
        [int((~severe & naive).sum()), int((~severe & ~naive).sum())],
    ]
    severe_or, severe_p = fisher_exact(table)

    table_missense = [
        [int((severe & naive).sum()), int((severe & ~naive).sum())],
        [int((missense & naive).sum()), int((missense & ~naive).sum())],
    ]
    missense_or, missense_p = fisher_exact(table_missense)

    severe_non_recessive = severe & ~observed["canonical_recessive_or_biallelic_gene"]
    other_non_recessive = ~severe & ~observed["canonical_recessive_or_biallelic_gene"]
    table_non_recessive = [
        [int((severe_non_recessive & naive).sum()), int((severe_non_recessive & ~naive).sum())],
        [int((other_non_recessive & naive).sum()), int((other_non_recessive & ~naive).sum())],
    ]
    non_recessive_or, non_recessive_p = fisher_exact(table_non_recessive)

    severe_constrained = severe & observed["constrained_loeuf_lt_0_5"]
    other_constrained = ~severe & observed["constrained_loeuf_lt_0_5"]
    table_constrained = [
        [int((severe_constrained & naive).sum()), int((severe_constrained & ~naive).sum())],
        [int((other_constrained & naive).sum()), int((other_constrained & ~naive).sum())],
    ]
    constrained_or, constrained_p = fisher_exact(table_constrained)

    cross = pd.read_csv(CROSS_DISEASE_SCORES)
    cross["expert_panel_truth"] = cross["review_strength"].eq("expert_panel")
    cross["red_priority"] = bool_series(cross["vital_red_flag"])
    table_truth = [
        [
            int((cross["expert_panel_truth"] & cross["red_priority"]).sum()),
            int((cross["expert_panel_truth"] & ~cross["red_priority"]).sum()),
        ],
        [
            int((~cross["expert_panel_truth"] & cross["red_priority"]).sum()),
            int((~cross["expert_panel_truth"] & ~cross["red_priority"]).sum()),
        ],
    ]
    truth_or, truth_p = fisher_exact(table_truth)

    return pd.DataFrame(
        [
            {
                "test": "severe_annotation_vs_naive_frequency_tension",
                "odds_ratio": severe_or,
                "p_value": severe_p,
                "table": str(table),
                "interpretation": "Tests whether severe LOF/splice annotations are overrepresented among AF-positive tension records.",
            },
            {
                "test": "severe_annotation_vs_missense_frequency_tension",
                "odds_ratio": missense_or,
                "p_value": missense_p,
                "table": str(table_missense),
                "interpretation": "Tests the key claim that severe annotations are not cleaner than missense with respect to AF tension.",
            },
            {
                "test": "severe_annotation_excluding_CASQ2_TRDN",
                "odds_ratio": non_recessive_or,
                "p_value": non_recessive_p,
                "table": str(table_non_recessive),
                "interpretation": "Controls the obvious recessive/biallelic explanation by excluding CASQ2 and TRDN.",
            },
            {
                "test": "severe_annotation_in_LOEUF_lt_0_5_genes",
                "odds_ratio": constrained_or,
                "p_value": constrained_p,
                "table": str(table_constrained),
                "interpretation": "Controls the constraint objection by focusing on genes with LOEUF upper CI <0.5.",
            },
            {
                "test": "expert_panel_truth_vs_red_priority",
                "odds_ratio": truth_or,
                "p_value": truth_p,
                "table": str(table_truth),
                "interpretation": "Tests whether expert-panel curated P/LP assertions are protected from red-priority calls.",
            },
        ]
    )


def plot_claim(severe: pd.DataFrame, truth: pd.DataFrame) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(13.8, 5.0))

    subtype = severe[severe["group"].str.startswith("lof_subtype:", na=False)].copy()
    subtype["label"] = subtype["group"].str.replace("lof_subtype:", "", regex=False)
    x = np.arange(len(subtype))
    axes[0].bar(x - 0.2, subtype["naive_af_flag_percent"], width=0.2, label="Naive AF >1e-5", color="#5b8db8")
    axes[0].bar(
        x,
        subtype["ac_supported_frequency_percent"],
        width=0.2,
        label="AC-supported",
        color="#f4a261",
    )
    axes[0].bar(x + 0.2, subtype["red_priority_percent"], width=0.2, label="Urgent review", color="#b23a48")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(subtype["label"], rotation=20, ha="right")
    axes[0].set_ylabel("Percent of exact AF-observed subtype")
    axes[0].set_title("Severe-annotation frequency discordance", weight="bold")
    axes[0].grid(axis="y", alpha=0.25)
    axes[0].legend(frameon=False, fontsize=8)

    truth_plot = truth[truth["truth_comparator_group"].isin(["expert_panel_reviewed_PLP", "all_cross_disease_PLP"])].copy()
    labels = ["Expert-panel P/LP", "All cross-disease P/LP"]
    x = np.arange(len(truth_plot))
    axes[1].bar(x - 0.24, truth_plot["naive_af_flag_percent"], width=0.24, label="Naive AF", color="#5b8db8")
    axes[1].bar(
        x,
        truth_plot["ac_supported_frequency_percent"],
        width=0.24,
        label="AC-supported",
        color="#f4a261",
    )
    axes[1].bar(x + 0.24, truth_plot["red_priority_percent"], width=0.24, label="Urgent review", color="#b23a48")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=0)
    axes[1].set_ylabel("Percent")
    axes[1].set_title("Expert-panel truth comparator", weight="bold")
    axes[1].grid(axis="y", alpha=0.25)
    axes[1].legend(frameon=False, fontsize=8)

    fig.suptitle("External truth and biological pattern audit", weight="bold")
    fig.tight_layout()
    path = FIGURE_DIR / "vital_external_truth_biological_claim.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    print(f"Saved {path}")


def main() -> None:
    arrhythmia = pd.read_csv(ARRHYTHMIA_SCORES)
    lof_summary = pd.read_csv(ARRHYTHMIA_LOF_SUMMARY)
    context = pd.read_csv(ARRHYTHMIA_CONTEXT)
    cross = pd.read_csv(CROSS_DISEASE_SCORES)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    SUPPLEMENTARY_DIR.mkdir(parents=True, exist_ok=True)

    severe, gene_spread = severe_annotation_summary(arrhythmia, lof_summary, context)
    mechanism_triage = mechanism_triage_summary(arrhythmia, context)
    truth, top_expert = summarize_expert_panel_truth(cross)
    tests = fisher_claim_tests(severe, truth, arrhythmia)

    severe.to_csv(DATA_DIR / "vital_severe_annotation_frequency_discordance_summary.csv", index=False)
    gene_spread.to_csv(DATA_DIR / "vital_severe_annotation_gene_spread.csv", index=False)
    mechanism_triage.to_csv(DATA_DIR / "vital_severe_annotation_mechanism_triage.csv", index=False)
    truth.to_csv(DATA_DIR / "vital_expert_panel_truth_validation.csv", index=False)
    top_expert.to_csv(DATA_DIR / "vital_expert_panel_high_tension_examples.csv", index=False)
    tests.to_csv(DATA_DIR / "vital_external_truth_claim_tests.csv", index=False)

    supplement = pd.concat(
        [
            severe.assign(table_section="severe_annotation_frequency_discordance"),
            gene_spread.assign(table_section="severe_annotation_gene_spread"),
            mechanism_triage.assign(table_section="severe_annotation_mechanism_triage"),
            truth.assign(table_section="expert_panel_truth_validation"),
            top_expert.assign(table_section="expert_panel_high_tension_examples"),
            tests.assign(table_section="statistical_tests"),
        ],
        ignore_index=True,
        sort=False,
    )
    supplement.fillna("NA").to_csv(
        SUPPLEMENTARY_DIR / "Supplementary_Table_S22_external_truth_biological_claim.tsv",
        sep="\t",
        index=False,
    )
    plot_claim(severe, truth)
    print("Saved external truth and biological claim outputs")


if __name__ == "__main__":
    main()
