from __future__ import annotations

from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, PathPatch, Rectangle
    from matplotlib.path import Path as MplPath
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
ARRHYTHMIA_RECON_IN = DATA_DIR / "vital_tiered_match_reconciliation_detail.csv"
EXPERT_ARRHYTHMIA_IN = DATA_DIR / "arrhythmia_expert_consensus_validation_regime_calls.csv"
EXTERNAL_EXPERT_PANEL_IN = DATA_DIR / "vital_cross_disease_3000_expert_panel_regime_calls.csv"

ROUTING_CALLS_OUT = DATA_DIR / "vital_routing_validation_calls.csv"
ROUTING_SUMMARY_OUT = DATA_DIR / "vital_routing_validation_summary.csv"
ROUTING_CONTEXT_OUT = DATA_DIR / "vital_routing_validation_context_summary.csv"
ROUTING_TRANSITIONS_OUT = DATA_DIR / "vital_routing_validation_route_transitions.csv"
ROUTING_HIGH_REVIEW_OUT = DATA_DIR / "vital_routing_validation_high_review_summary.csv"
ROUTING_EXPERT_OUT = DATA_DIR / "vital_routing_validation_expert_concordance.csv"
ROUTING_VIGNETTES_OUT = DATA_DIR / "vital_routing_validation_case_vignettes.csv"
BASELINE_DECISION_MODEL_OUT = DATA_DIR / "vital_formal_baseline_decision_model.csv"
COUNTERFACTUAL_AUDIT_CALLS_OUT = DATA_DIR / "vital_counterfactual_decision_audit_calls.csv"
COUNTERFACTUAL_AUDIT_SUMMARY_OUT = DATA_DIR / "vital_counterfactual_decision_audit_summary.csv"
ADS_OUT = DATA_DIR / "vital_actionability_discordance_audit.csv"
SIMULATED_CDS_ALERT_CALLS_OUT = DATA_DIR / "vital_simulated_cds_alert_calls.csv"
SIMULATED_CDS_ALERT_SUMMARY_OUT = DATA_DIR / "vital_simulated_cds_alert_summary.csv"
REPAIR_LOGIC_OUT = DATA_DIR / "vital_minimal_repair_logic.csv"
REASON_CODE_DEFINITIONS_OUT = DATA_DIR / "vital_reason_code_definitions.csv"
ROUTING_FIGURE_OUT = FIGURE_DIR / "vital_routing_validation.png"
DECISION_DISRUPTION_FIGURE_OUT = FIGURE_DIR / "vital_decision_disruption.png"

RECESSIVE_CONTEXT_GENES = {"CASQ2", "TRDN"}

CLINICAL_CONTEXTS: list[tuple[str, set[str], str]] = [
    (
        "cascade_testing",
        {"SCN5A", "KCNH2", "KCNQ1", "RYR2", "CACNA1C", "KCNE1", "KCNJ2", "CASQ2", "TRDN"},
        "Cascade testing / family follow-up exposure.",
    ),
    (
        "drug_guidance",
        {"SCN5A", "KCNH2", "KCNQ1", "KCNE1", "CACNA1C", "KCNJ2"},
        "Medication restriction or drug-avoidance exposure.",
    ),
    (
        "surveillance_device",
        {"SCN5A", "KCNH2", "KCNQ1", "RYR2", "CACNA1C"},
        "Surveillance-intensity or device-related exposure.",
    ),
    (
        "carrier_recessive",
        RECESSIVE_CONTEXT_GENES,
        "Carrier-state or recessive/affected-state exposure.",
    ),
]

ROUTE_ORDER = [
    "VITAL_OK",
    "CHECK_POPMAX",
    "LOW_COUNT_UNCERTAIN",
    "MODEL_CONFLICT_REVIEW",
    "RECESSIVE_MODEL_REQUIRED",
    "RECESSIVE_COMPATIBLE",
    "NOT_EVALUABLE_NO_FREQUENCY_CLAIM",
]

ROUTE_LABELS = {
    "VITAL_OK": "VITAL_OK",
    "CHECK_POPMAX": "CHECK_POPMAX",
    "LOW_COUNT_UNCERTAIN": "LOW_COUNT_UNCERTAIN",
    "MODEL_CONFLICT_REVIEW": "MODEL_CONFLICT_REVIEW",
    "RECESSIVE_MODEL_REQUIRED": "RECESSIVE_MODEL_REQUIRED",
    "RECESSIVE_COMPATIBLE": "RECESSIVE_COMPATIBLE",
    "NOT_EVALUABLE_NO_FREQUENCY_CLAIM": "NOT_EVALUABLE",
    "CONTROL_REVIEW_FLAG": "CONTROL_REVIEW",
    "MODEL_CONFLICT_CONTROL_FLAG": "CONTROL_MODEL_CONFLICT",
}

CONSTRAINT_ROUTE_LABELS = {
    "VITAL_OK": "VITAL_OK",
    "CHECK_POPMAX": "CHECK_POPMAX",
    "CHECK_MODEL": "CHECK_MODEL",
    "MODEL_CONFLICT": "MODEL_CONFLICT",
    "EVAL_LIMITED": "EVAL_LIMITED",
}

SIMULATED_ALERT_LABELS = {
    "NONE": "No alert",
    "YELLOW": "Yellow alert: actionability requires population/mechanism review",
    "ORANGE": "Orange alert: allele-resolved population evidence unavailable",
    "RED": "Red alert: dominant high-penetrance model incompatible",
}

ROUTE_COLORS = {
    "VITAL_OK": "#2f6f62",
    "CHECK_POPMAX": "#d18f2f",
    "LOW_COUNT_UNCERTAIN": "#e2b565",
    "MODEL_CONFLICT_REVIEW": "#b94a48",
    "RECESSIVE_MODEL_REQUIRED": "#4b6cb7",
    "RECESSIVE_COMPATIBLE": "#809ad5",
    "NOT_EVALUABLE_NO_FREQUENCY_CLAIM": "#8b98a5",
    "CONTROL_REVIEW_FLAG": "#d18f2f",
    "MODEL_CONFLICT_CONTROL_FLAG": "#b94a48",
}

CHANGE_CLASS_LABELS = {
    "unchanged": "Unchanged",
    "risk_escalating": "Risk-escalating",
    "model_correcting": "Model-correcting",
    "evaluation_limiting": "Evaluation-limiting",
}

CHANGE_CLASS_COLORS = {
    "unchanged": "#2f6f62",
    "risk_escalating": "#d18f2f",
    "model_correcting": "#4b6cb7",
    "evaluation_limiting": "#8b98a5",
}

CASE_SPECS = [
    {
        "case_id": "Case 1",
        "source_dataset": "arrhythmia_main",
        "clinvar_id": "VCV000440850",
        "clinical_context": "cascade_testing|drug_guidance|surveillance_device",
        "baseline_route": "ROUTE_PLP_ACTIONABLE",
        "interpretation": "Requires explicit dominant-model specification before cascade testing or management reuse.",
    },
    {
        "case_id": "Case 2",
        "source_dataset": "arrhythmia_main",
        "clinvar_id": "VCV004535537",
        "clinical_context": "cascade_testing|drug_guidance|surveillance_device",
        "baseline_route": "ROUTE_PLP_ACTIONABLE",
        "interpretation": "Ancestry-aware review is required before medication or surveillance reuse.",
    },
    {
        "case_id": "Case 3",
        "source_dataset": "arrhythmia_main",
        "clinvar_id": "VCV001325231",
        "clinical_context": "carrier_recessive|cascade_testing",
        "baseline_route": "ROUTE_PLP_ACTIONABLE",
        "interpretation": "Dominant-actionable reuse should be replaced by recessive/carrier-state logic.",
    },
    {
        "case_id": "Case 4",
        "source_dataset": "arrhythmia_main",
        "clinvar_id": "VCV000405355",
        "clinical_context": "surveillance_device|cascade_testing",
        "baseline_route": "ROUTE_PLP_ACTIONABLE",
        "interpretation": "No allele-resolved frequency claim should be made before representation-sensitive review.",
    },
    {
        "case_id": "Case 5",
        "source_dataset": "expert_arrhythmia_plp",
        "clinvar_id": "VCV000003118",
        "clinical_context": "cascade_testing|drug_guidance|surveillance_device",
        "baseline_route": "ROUTE_PLP_ACTIONABLE",
        "interpretation": "Positive-control expert-panel assertion remains stable under the VITAL layer.",
    },
]


def truthy(value: object) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def pct(numerator: int, denominator: int) -> float | None:
    if denominator == 0:
        return None
    return 100.0 * numerator / denominator


def save_table(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, lineterminator="\n")
    print(f"Saved {path.relative_to(BASE_DIR)} ({len(df)} rows)")


def infer_evaluability(row: pd.Series) -> str:
    frequency_status = str(row.get("frequency_evidence_status", "")).strip()
    tier = str(row.get("reconciliation_tier", "")).strip()
    if frequency_status == "frequency_observed":
        return "TIER_1"
    if tier in {"tier4_still_unevaluable", "tier4_query_error"}:
        return "UNEVALUABLE"
    if tier in {"tier3_locus_observed_no_exact_af", "tier3_regional_context_only"} or frequency_status in {
        "exact_match_without_AF",
        "allele_discordance_no_exact_AF",
    }:
        return "TIER_2"
    return "UNEVALUABLE"


def infer_flag(row: pd.Series) -> str:
    if row["VITAL_evaluability"] != "TIER_1":
        return "."
    if not truthy(row.get("standard_acmg_frequency_flag")):
        return "OK"
    global_af = pd.to_numeric(row.get("global_af"), errors="coerce")
    if (
        str(row.get("gene", "")).upper() not in RECESSIVE_CONTEXT_GENES
        and truthy(row.get("frequency_signal_ac_ge_20"))
        and pd.notna(global_af)
        and float(global_af) > 1e-4
    ):
        return "MODEL_CONFLICT"
    return "CHECK_POPMAX"


def infer_routing_route(row: pd.Series) -> str:
    if row["VITAL_evaluability"] != "TIER_1":
        return "NOT_EVALUABLE_NO_FREQUENCY_CLAIM"
    gene = str(row.get("gene", "")).upper()
    flag = str(row.get("VITAL_flag", "")).strip()
    ac_supported = truthy(row.get("frequency_signal_ac_ge_20"))
    if gene in RECESSIVE_CONTEXT_GENES and flag != "OK":
        return "RECESSIVE_MODEL_REQUIRED" if ac_supported else "RECESSIVE_COMPATIBLE"
    if flag == "MODEL_CONFLICT":
        return "MODEL_CONFLICT_REVIEW"
    if flag == "CHECK_POPMAX":
        return "CHECK_POPMAX" if ac_supported else "LOW_COUNT_UNCERTAIN"
    return "VITAL_OK"


def infer_change_class(route: str) -> str:
    if route == "VITAL_OK":
        return "unchanged"
    if route in {"CHECK_POPMAX", "LOW_COUNT_UNCERTAIN", "MODEL_CONFLICT_REVIEW"}:
        return "risk_escalating"
    if route in {"RECESSIVE_MODEL_REQUIRED", "RECESSIVE_COMPATIBLE"}:
        return "model_correcting"
    return "evaluation_limiting"


def infer_constraint_route(row: pd.Series) -> str:
    route = str(row.get("vital_route", "")).strip()
    if route == "VITAL_OK":
        return "VITAL_OK"
    if route in {"CHECK_POPMAX", "LOW_COUNT_UNCERTAIN"}:
        return "CHECK_POPMAX"
    if route == "MODEL_CONFLICT_REVIEW":
        return "MODEL_CONFLICT"
    if route in {"RECESSIVE_MODEL_REQUIRED", "RECESSIVE_COMPATIBLE"}:
        return "CHECK_MODEL"
    return "EVAL_LIMITED"


def infer_constraint_violated(row: pd.Series) -> str:
    route = str(row.get("VITAL_constraint_route", "")).strip()
    if route == "VITAL_OK":
        return "none_detected_under_tested_model"
    if route == "CHECK_POPMAX":
        return "ancestry-aware frequency or allele-count review required"
    if route == "CHECK_MODEL":
        return "inheritance or mechanism model cannot be flattened into direct actionability"
    if route == "MODEL_CONFLICT":
        return "dominant high-penetrance model incompatible under current frequency constraints"
    return "allele-resolved population evaluability unavailable"


def infer_reason_code(row: pd.Series) -> str:
    route = str(row.get("vital_route", "")).strip()
    evaluability = str(row.get("VITAL_evaluability", "")).strip()
    if route == "VITAL_OK":
        return "EXPERT_REVIEW_RETAINED"
    if route == "LOW_COUNT_UNCERTAIN":
        return "CHECK_LOW_COUNT"
    if route == "CHECK_POPMAX":
        return "POPMAX_EXCEEDS_MCAF"
    if route == "MODEL_CONFLICT_REVIEW":
        return "MODEL_CONFLICT_DOMINANT_HP"
    if route in {"RECESSIVE_MODEL_REQUIRED", "RECESSIVE_COMPATIBLE"}:
        return "RECESSIVE_COMPATIBLE_NOT_DOMINANT"
    if evaluability == "TIER_2":
        return "EVAL_LOCI_ONLY"
    return "EVAL_NO_ALLELE_MATCH"


def infer_simulated_alert_level(row: pd.Series) -> str:
    route = str(row.get("VITAL_constraint_route", "")).strip()
    if route == "VITAL_OK":
        return "NONE"
    if route == "MODEL_CONFLICT":
        return "RED"
    if route == "EVAL_LIMITED":
        return "ORANGE"
    return "YELLOW"


def infer_repair_next_step(row: pd.Series) -> str:
    route = str(row.get("VITAL_constraint_route", "")).strip()
    reason = str(row.get("reason_code", "")).strip()
    if route == "VITAL_OK":
        return "Proceed with standard expert interpretation."
    if route == "CHECK_POPMAX":
        if reason == "CHECK_LOW_COUNT":
            return "Repeat ancestry-aware review with allele-count support and penetrance/model assumptions."
        return "Perform ancestry-aware frequency review and penetrance/model reassessment."
    if route == "MODEL_CONFLICT":
        return "Do not use flattened P/LP label as a dominant-actionability claim."
    if route == "CHECK_MODEL":
        return "Specify inheritance, phase, and mechanism before clinical actionability is routed."
    return "Defer direct actionability; require allele representation, callability, or orthogonal evidence."


def infer_contexts(gene: object) -> str:
    gene_symbol = str(gene).strip().upper()
    labels = [name for name, genes, _ in CLINICAL_CONTEXTS if gene_symbol in genes]
    return "|".join(labels) if labels else "syndrome_diagnosis"


def infer_primary_context(gene: object) -> str:
    gene_symbol = str(gene).strip().upper()
    if gene_symbol in RECESSIVE_CONTEXT_GENES:
        return "carrier_recessive"
    if gene_symbol in {"SCN5A", "KCNH2", "KCNQ1", "RYR2", "CACNA1C"}:
        return "surveillance_device"
    if gene_symbol in {"KCNE1", "KCNJ2"}:
        return "drug_guidance"
    if gene_symbol in {"SCN5A", "KCNH2", "KCNQ1", "RYR2", "CACNA1C", "KCNE1", "KCNJ2"}:
        return "cascade_testing"
    return "syndrome_diagnosis"


def baseline_decision_without_vital(row: pd.Series) -> str:
    clinical_group = str(row.get("clinical_group", "P_LP")).strip()
    if clinical_group == "B_LB":
        return "ROUTE_NO_ACTION"
    if clinical_group != "P_LP":
        return "ROUTE_REVIEW"
    return "ROUTE_PLP_ACTIONABLE"


def decision_with_vital(route: str) -> str:
    if route == "VITAL_OK":
        return "ACTIONABLE_MAINTAIN"
    if route in {"CHECK_POPMAX", "LOW_COUNT_UNCERTAIN", "MODEL_CONFLICT_REVIEW"}:
        return "ACTION_REVIEW_REQUIRED"
    if route in {"RECESSIVE_MODEL_REQUIRED", "RECESSIVE_COMPATIBLE"}:
        return "ACTION_MODEL_SPECIFIC_REROUTE"
    return "ACTION_DEFER_NO_FREQUENCY_CLAIM"


def counterfactual_audit_category(row: pd.Series) -> str:
    baseline = str(row.get("decision_without_vital", "")).strip()
    with_vital = str(row.get("decision_with_vital", "")).strip()
    if baseline == "ROUTE_NO_ACTION":
        return "control_no_action"
    if baseline != "ROUTE_PLP_ACTIONABLE":
        return "baseline_review"
    if with_vital == "ACTIONABLE_MAINTAIN":
        return "sustained_action"
    if with_vital == "ACTION_DEFER_NO_FREQUENCY_CLAIM":
        return "unjustified_action_evaluation_limiting"
    if with_vital == "ACTION_REVIEW_REQUIRED":
        return "unjustified_action_risk_escalating"
    if with_vital == "ACTION_MODEL_SPECIFIC_REROUTE":
        return "unjustified_action_model_correcting"
    return "baseline_review"


def build_formal_baseline_decision_model() -> pd.DataFrame:
    rows = [
        {
            "rule_id": "B1_LABEL_DRIVEN_BASELINE",
            "input_condition": "A variant is exported as public ClinVar P/LP.",
            "decision_without_vital": "ROUTE_PLP_ACTIONABLE",
            "what_counts_as_actionable": "Downstream systems may treat the label as eligible for direct actionability unless additional context is explicitly encoded.",
            "why_this_is_a_realistic_baseline": "Formalizes label-driven reuse as the null model being audited, not as a recommended clinical practice.",
        },
        {
            "rule_id": "B2_ACTIONABILITY_CONTEXT_BASELINE",
            "input_condition": "Public P/LP appears in a gene linked to intervention, surveillance, cascade testing, device consideration, or therapy selection.",
            "decision_without_vital": "ROUTE_PLP_ACTIONABLE",
            "what_counts_as_actionable": "The baseline route is direct-actionable in the relevant workflow context.",
            "why_this_is_a_realistic_baseline": "Captures the practical actionability proxy created when condition, inheritance, and evaluability fields are flattened.",
        },
        {
            "rule_id": "B3_NO_IMPLIED_BENIGNITY",
            "input_condition": "A public P/LP variant does not remain VITAL_OK after constraint routing.",
            "decision_without_vital": "NO_BENIGNITY_INFERENCE",
            "what_counts_as_actionable": "Absence of VITAL_OK does not imply benignity, incorrectness, or clinical irrelevance.",
            "why_this_is_a_realistic_baseline": "Separates pathogenicity classification from actionability routing.",
        },
        {
            "rule_id": "B4_CONTROL_NO_ACTION",
            "input_condition": "ClinVar classification in {Benign, Likely benign}",
            "decision_without_vital": "ROUTE_NO_ACTION",
            "what_counts_as_actionable": "No positive disease-action escalation.",
            "why_this_is_a_realistic_baseline": "Used for the B/LB control layer.",
        },
        {
            "rule_id": "B5_CONTEXT_INSUFFICIENT_REVIEW",
            "input_condition": "Any non-P/LP or context-insufficient record outside the defined reuse environments",
            "decision_without_vital": "ROUTE_REVIEW",
            "what_counts_as_actionable": "No direct action; manual review required.",
            "why_this_is_a_realistic_baseline": "Completes the formal decision model even though the main arrhythmia cohort is P/LP-enriched.",
        },
    ]
    return pd.DataFrame(rows)


def load_arrhythmia_routing_calls() -> pd.DataFrame:
    if not ARRHYTHMIA_SCORES_IN.exists():
        raise FileNotFoundError(f"Missing input: {ARRHYTHMIA_SCORES_IN}")
    if not ARRHYTHMIA_RECON_IN.exists():
        raise FileNotFoundError(f"Missing input: {ARRHYTHMIA_RECON_IN}")

    scores = pd.read_csv(ARRHYTHMIA_SCORES_IN, low_memory=False)
    reconciliation = pd.read_csv(
        ARRHYTHMIA_RECON_IN,
        usecols=["variant_key", "reconciliation_tier", "reason"],
        low_memory=False,
    )
    calls = scores.merge(reconciliation, on="variant_key", how="left")
    calls["VITAL_evaluability"] = calls.apply(infer_evaluability, axis=1)
    calls["VITAL_flag"] = calls.apply(infer_flag, axis=1)
    calls["primary_workflow_context"] = calls["gene"].map(infer_primary_context)
    calls["baseline_route"] = calls.apply(baseline_decision_without_vital, axis=1)
    calls["vital_route"] = calls.apply(infer_routing_route, axis=1)
    calls["VITAL_constraint_route"] = calls.apply(infer_constraint_route, axis=1)
    calls["constraint_violated"] = calls.apply(infer_constraint_violated, axis=1)
    calls["reason_code"] = calls.apply(infer_reason_code, axis=1)
    calls["simulated_cds_alert"] = calls.apply(infer_simulated_alert_level, axis=1)
    calls["simulated_cds_alert_label"] = calls["simulated_cds_alert"].map(SIMULATED_ALERT_LABELS)
    calls["recommended_next_step"] = calls.apply(infer_repair_next_step, axis=1)
    calls["routing_changed"] = calls["vital_route"].ne("VITAL_OK")
    calls["actionability_at_risk"] = calls["routing_changed"]
    calls["routing_change_class"] = calls["vital_route"].map(infer_change_class)
    calls["clinical_contexts"] = calls["gene"].map(infer_contexts)
    calls["decision_without_vital"] = calls["baseline_route"]
    calls["decision_with_vital"] = calls["vital_route"].map(decision_with_vital)
    calls["counterfactual_audit_category"] = calls.apply(counterfactual_audit_category, axis=1)
    calls["unjustified_action_without_vital"] = calls["counterfactual_audit_category"].isin(
        {
            "unjustified_action_evaluation_limiting",
            "unjustified_action_risk_escalating",
            "unjustified_action_model_correcting",
        }
    )
    calls["corrected_routing_by_vital"] = calls["unjustified_action_without_vital"]
    calls["review_score"] = pd.to_numeric(calls.get("review_score"), errors="coerce")
    return calls


def summarize_main_cohort(calls: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    route_counts = (
        calls.groupby(["baseline_route", "vital_route", "routing_change_class"], dropna=False)
        .size()
        .rename("variant_count")
        .reset_index()
    )
    route_counts["percent_of_arrhythmia_plp"] = route_counts["variant_count"] / len(calls) * 100

    summary_rows = [
        {
            "cohort": "arrhythmia_plp_all",
            "metric": "routing_changed",
            "numerator": int(calls["routing_changed"].sum()),
            "denominator": len(calls),
            "percent": pct(int(calls["routing_changed"].sum()), len(calls)),
            "note": "Baseline ROUTE_PLP_ACTIONABLE rerouted by the VITAL layer.",
        },
        {
            "cohort": "arrhythmia_plp_all",
            "metric": "actionability_at_risk",
            "numerator": int(calls["actionability_at_risk"].sum()),
            "denominator": len(calls),
            "percent": pct(int(calls["actionability_at_risk"].sum()), len(calls)),
            "note": "Public P/LP variants that would be directly reused as actionable but are flagged or limited by VITAL.",
        },
        {
            "cohort": "arrhythmia_plp_all",
            "metric": "vital_ok",
            "numerator": int(calls["vital_route"].eq("VITAL_OK").sum()),
            "denominator": len(calls),
            "percent": pct(int(calls["vital_route"].eq("VITAL_OK").sum()), len(calls)),
            "note": "No routing change under the routing-validation layer.",
        },
    ]
    for change_class in ["evaluation_limiting", "risk_escalating", "model_correcting"]:
        count = int(calls["routing_change_class"].eq(change_class).sum())
        summary_rows.append(
            {
                "cohort": "arrhythmia_plp_all",
                "metric": change_class,
                "numerator": count,
                "denominator": len(calls),
                "percent": pct(count, len(calls)),
                "note": CHANGE_CLASS_LABELS[change_class],
            }
        )
    summary = pd.DataFrame(summary_rows)

    context_rows: list[dict[str, object]] = []
    for context_name, genes, interpretation in CLINICAL_CONTEXTS:
        sub = calls.loc[calls["gene"].fillna("").astype(str).isin(genes)].copy()
        if sub.empty:
            continue
        total = len(sub)
        context_rows.append(
            {
                "clinical_context": context_name,
                "n_total": total,
                "n_actionability_at_risk": int(sub["actionability_at_risk"].sum()),
                "actionability_at_risk_percent": pct(int(sub["actionability_at_risk"].sum()), total),
                "n_vital_ok": int(sub["vital_route"].eq("VITAL_OK").sum()),
                "n_evaluation_limiting": int(sub["routing_change_class"].eq("evaluation_limiting").sum()),
                "n_risk_escalating": int(sub["routing_change_class"].eq("risk_escalating").sum()),
                "n_model_correcting": int(sub["routing_change_class"].eq("model_correcting").sum()),
                "interpretation": interpretation,
            }
        )
    context_summary = pd.DataFrame(context_rows)
    return route_counts, summary, context_summary


def summarize_high_review_subset(calls: pd.DataFrame) -> pd.DataFrame:
    high_review = calls.loc[calls["review_score"].ge(2)].copy()
    if high_review.empty:
        return pd.DataFrame()

    rows = [
        {
            "cohort": "arrhythmia_high_review_ge2",
            "metric": "total_variants",
            "numerator": len(high_review),
            "denominator": len(high_review),
            "percent": 100.0,
            "note": "ClinVar review score >= 2 subset.",
        },
        {
            "cohort": "arrhythmia_high_review_ge2",
            "metric": "tier1_evaluable",
            "numerator": int(high_review["VITAL_evaluability"].eq("TIER_1").sum()),
            "denominator": len(high_review),
            "percent": pct(int(high_review["VITAL_evaluability"].eq("TIER_1").sum()), len(high_review)),
            "note": "Allele-resolved AF evaluable within the high-review subset.",
        },
        {
            "cohort": "arrhythmia_high_review_ge2",
            "metric": "actionability_at_risk",
            "numerator": int(high_review["actionability_at_risk"].sum()),
            "denominator": len(high_review),
            "percent": pct(int(high_review["actionability_at_risk"].sum()), len(high_review)),
            "note": "Baseline actionable variants rerouted or limited by VITAL.",
        },
    ]
    for route in ROUTE_ORDER:
        route_count = int(high_review["vital_route"].eq(route).sum())
        if route_count == 0:
            continue
        rows.append(
            {
                "cohort": "arrhythmia_high_review_ge2",
                "metric": route,
                "numerator": route_count,
                "denominator": len(high_review),
                "percent": pct(route_count, len(high_review)),
                "note": ROUTE_LABELS[route],
            }
        )
    return pd.DataFrame(rows)


def summarize_counterfactual_decision_audit(calls: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    def add_scope(scope_name: str, frame: pd.DataFrame, note: str) -> None:
        total = len(frame)
        if total == 0:
            return
        unjustified = int(frame["unjustified_action_without_vital"].sum())
        sustained = int(frame["counterfactual_audit_category"].eq("sustained_action").sum())
        rows.append(
            {
                "scope": scope_name,
                "metric": "incorrectly_escalated_without_vital",
                "numerator": unjustified,
                "denominator": total,
                "percent": pct(unjustified, total),
                "note": note,
            }
        )
        rows.append(
            {
                "scope": scope_name,
                "metric": "sustained_action_after_vital",
                "numerator": sustained,
                "denominator": total,
                "percent": pct(sustained, total),
                "note": "Baseline-actionable decisions preserved after VITAL.",
            }
        )
        for category in [
            "unjustified_action_evaluation_limiting",
            "unjustified_action_risk_escalating",
            "unjustified_action_model_correcting",
        ]:
            count = int(frame["counterfactual_audit_category"].eq(category).sum())
            rows.append(
                {
                    "scope": scope_name,
                    "metric": category,
                    "numerator": count,
                    "denominator": total,
                    "percent": pct(count, total),
                    "note": "Counterfactual audit subtype.",
                }
            )

    baseline_actionable = calls.loc[calls["decision_without_vital"].eq("ROUTE_PLP_ACTIONABLE")].copy()
    add_scope(
        "all_baseline_actionable",
        baseline_actionable,
        "Formal counterfactual decision audit across the full ClinVar-driven actionable baseline.",
    )
    add_scope(
        "tier1_evaluable_baseline_actionable",
        baseline_actionable.loc[baseline_actionable["VITAL_evaluability"].eq("TIER_1")].copy(),
        "Restricted to variants with allele-resolved AF available.",
    )
    add_scope(
        "high_review_baseline_actionable",
        baseline_actionable.loc[baseline_actionable["review_score"].ge(2)].copy(),
        "Restricted to ClinVar review score >= 2.",
    )
    add_scope(
        "high_review_tier1_evaluable",
        baseline_actionable.loc[
            baseline_actionable["review_score"].ge(2) & baseline_actionable["VITAL_evaluability"].eq("TIER_1")
        ].copy(),
        "High-review subset with allele-resolved AF available.",
    )
    return pd.DataFrame(rows)


def expert_route_group(row: pd.Series, control_mode: bool = False) -> str:
    if control_mode:
        if str(row.get("VITAL_flag", "")).strip() == "MODEL_CONFLICT":
            return "MODEL_CONFLICT_CONTROL_FLAG"
        if str(row.get("VITAL_flag", "")).strip() == "CHECK_POPMAX":
            return "CONTROL_REVIEW_FLAG"
        return "CONTROL_NO_ACTION"
    if str(row.get("VITAL_evaluability", "")).strip() != "TIER_1":
        return "NOT_EVALUABLE_NO_FREQUENCY_CLAIM"
    gene = str(row.get("gene", "")).strip().upper()
    flag = str(row.get("VITAL_flag", "")).strip()
    if gene in RECESSIVE_CONTEXT_GENES and flag != "OK":
        return "CONTROL_REVIEW_FLAG"
    if flag == "MODEL_CONFLICT":
        return "MODEL_CONFLICT_REVIEW"
    if flag == "CHECK_POPMAX":
        return "CONTROL_REVIEW_FLAG"
    return "VITAL_OK"


def summarize_expert_concordance() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if not EXPERT_ARRHYTHMIA_IN.exists():
        raise FileNotFoundError(f"Missing input: {EXPERT_ARRHYTHMIA_IN}")
    if not EXTERNAL_EXPERT_PANEL_IN.exists():
        raise FileNotFoundError(f"Missing input: {EXTERNAL_EXPERT_PANEL_IN}")

    expert_arrhythmia = pd.read_csv(EXPERT_ARRHYTHMIA_IN, low_memory=False)
    external_expert = pd.read_csv(EXTERNAL_EXPERT_PANEL_IN, low_memory=False)

    arr_plp = expert_arrhythmia.loc[expert_arrhythmia["clinical_group"].eq("P_LP")].copy()
    arr_blb = expert_arrhythmia.loc[expert_arrhythmia["clinical_group"].eq("B_LB")].copy()
    arr_plp["expert_route_group"] = arr_plp.apply(expert_route_group, axis=1)
    arr_blb["expert_route_group"] = arr_blb.apply(lambda row: expert_route_group(row, control_mode=True), axis=1)
    external_expert["expert_route_group"] = external_expert.apply(expert_route_group, axis=1)

    dataset_specs = [
        (
            "arrhythmia_expert_plp",
            arr_plp,
            "Expert-curated arrhythmia positive-control set.",
        ),
        (
            "external_expert_plp",
            external_expert,
            "Independent non-arrhythmia expert-panel stress set.",
        ),
        (
            "arrhythmia_blb_controls",
            arr_blb,
            "Benign/likely benign control layer for dominant-model compatibility.",
        ),
        (
            "combined_expert_plp",
            pd.concat([arr_plp, external_expert], ignore_index=True),
            "Combined evaluability and concordance view across expert-curated P/LP positives.",
        ),
    ]

    summary_rows: list[dict[str, object]] = []
    for dataset_name, frame, note in dataset_specs:
        total = len(frame)
        tier1 = int(frame["VITAL_evaluability"].eq("TIER_1").sum()) if "VITAL_evaluability" in frame.columns else 0
        for route_group in [
            "VITAL_OK",
            "CONTROL_REVIEW_FLAG",
            "MODEL_CONFLICT_REVIEW",
            "NOT_EVALUABLE_NO_FREQUENCY_CLAIM",
            "MODEL_CONFLICT_CONTROL_FLAG",
        ]:
            if route_group not in set(frame["expert_route_group"]):
                continue
            count = int(frame["expert_route_group"].eq(route_group).sum())
            summary_rows.append(
                {
                    "dataset": dataset_name,
                    "route_group": route_group,
                    "count": count,
                    "denominator": total,
                    "percent_of_dataset": pct(count, total),
                    "tier1_evaluable_n": tier1,
                    "dataset_note": note,
                }
            )

    summary = pd.DataFrame(summary_rows)
    return summary, arr_plp, arr_blb


def build_case_vignettes(
    arrhythmia_calls: pd.DataFrame,
    expert_arrhythmia_plp: pd.DataFrame,
) -> pd.DataFrame:
    lookup_frames = [
        arrhythmia_calls.assign(source_dataset="arrhythmia_main"),
        expert_arrhythmia_plp.assign(source_dataset="expert_arrhythmia_plp", vital_route=expert_arrhythmia_plp["expert_route_group"]),
    ]
    lookup = pd.concat(lookup_frames, ignore_index=True, sort=False)
    lookup["clinvar_id"] = lookup["clinvar_id"].fillna("").astype(str)
    rows: list[dict[str, object]] = []
    for spec in CASE_SPECS:
        match = lookup.loc[
            lookup["source_dataset"].eq(spec["source_dataset"])
            & lookup["clinvar_id"].eq(spec["clinvar_id"])
        ]
        if match.empty:
            continue
        row = match.iloc[0]
        vital_route = row.get("vital_route", "")
        constraint_route = row.get("VITAL_constraint_route", "")
        if pd.isna(constraint_route) or str(constraint_route).strip() == "":
            constraint_route = infer_constraint_route(pd.Series({"vital_route": vital_route}))
        reason_code = row.get("reason_code", "")
        if pd.isna(reason_code) or str(reason_code).strip() == "":
            reason_code = infer_reason_code(
                pd.Series(
                    {
                        "vital_route": vital_route,
                        "VITAL_evaluability": row.get("VITAL_evaluability", ""),
                    }
                )
            )
        constraint_violated = row.get("constraint_violated", "")
        if pd.isna(constraint_violated) or str(constraint_violated).strip() == "":
            constraint_violated = infer_constraint_violated(
                pd.Series({"VITAL_constraint_route": constraint_route})
            )
        recommended_next_step = row.get("recommended_next_step", "")
        if pd.isna(recommended_next_step) or str(recommended_next_step).strip() == "":
            recommended_next_step = infer_repair_next_step(
                pd.Series({"VITAL_constraint_route": constraint_route, "reason_code": reason_code})
            )
        rows.append(
            {
                "case_id": spec["case_id"],
                "source_dataset": spec["source_dataset"],
                "clinical_context": spec["clinical_context"],
                "baseline_route": spec["baseline_route"],
                "vital_route": vital_route,
                "VITAL_constraint_route": constraint_route,
                "constraint_violated": constraint_violated,
                "reason_code": reason_code,
                "recommended_next_step": recommended_next_step,
                "gene": row.get("gene", ""),
                "clinvar_id": row.get("clinvar_id", ""),
                "title": row.get("title", ""),
                "review_status": row.get("review_status", ""),
                "frequency_evidence_status": row.get("frequency_evidence_status", ""),
                "reconciliation_tier": row.get("reconciliation_tier", ""),
                "global_af": row.get("global_af", np.nan),
                "popmax_af": row.get("popmax_af", np.nan),
                "popmax_population": row.get("popmax_population", ""),
                "qualifying_frequency_ac": row.get("qualifying_frequency_ac", np.nan),
                "interpretation": spec["interpretation"],
            }
        )
    return pd.DataFrame(rows)


def ads_row_from_series(row: pd.Series, source: str, interpretation: str | None = None) -> dict[str, object]:
    title = str(row.get("title", "")).strip()
    variant_id = str(row.get("clinvar_id", "")).strip()
    gene = str(row.get("gene", "")).strip()
    variant = f"{gene} {variant_id}".strip()
    if title and title.lower() != "nan":
        variant = f"{variant} ({title})"
    return {
        "Variant": variant,
        "Gene": gene,
        "Disease domain": "arrhythmia",
        "Public label": "P/LP",
        "Action context": row.get("clinical_context", row.get("clinical_contexts", "")),
        "Baseline route": row.get("baseline_route", "ROUTE_PLP_ACTIONABLE"),
        "VITAL route": row.get("VITAL_constraint_route", ""),
        "Constraint violated": row.get("constraint_violated", ""),
        "Reason code": row.get("reason_code", ""),
        "Interpretation": interpretation
        or "Direct-actionability no longer supported by the flattened label alone.",
        "Source": source,
    }


def build_actionability_discordance_audit(
    arrhythmia_calls: pd.DataFrame,
    case_vignettes: pd.DataFrame,
    target_n: int = 24,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    seen: set[str] = set()
    for _, row in case_vignettes.iterrows():
        constraint_route = str(row.get("VITAL_constraint_route", row.get("vital_route", ""))).strip()
        if constraint_route == "VITAL_OK":
            continue
        out = ads_row_from_series(
            row,
            source="repository_case_vignette",
            interpretation=(
                "Direct-actionability no longer supported by the flattened label alone; "
                f"{row.get('interpretation', '')}"
            ).strip(),
        )
        key = str(row.get("clinvar_id", row.get("Variant", ""))).strip()
        rows.append(out)
        seen.add(key)

    candidates = arrhythmia_calls.loc[
        arrhythmia_calls["baseline_route"].eq("ROUTE_PLP_ACTIONABLE")
        & arrhythmia_calls["VITAL_constraint_route"].ne("VITAL_OK")
    ].copy()
    candidates["popmax_af_numeric"] = pd.to_numeric(candidates.get("popmax_af"), errors="coerce").fillna(-1)
    candidates["qualifying_ac_numeric"] = pd.to_numeric(
        candidates.get("qualifying_frequency_ac"), errors="coerce"
    ).fillna(-1)
    candidates["review_score_numeric"] = pd.to_numeric(candidates.get("review_score"), errors="coerce").fillna(-1)

    route_targets = [
        ("MODEL_CONFLICT", 3),
        ("CHECK_POPMAX", 7),
        ("CHECK_MODEL", 7),
        ("EVAL_LIMITED", 7),
    ]
    for route, limit in route_targets:
        subset = candidates.loc[candidates["VITAL_constraint_route"].eq(route)].copy()
        if subset.empty:
            continue
        subset = subset.sort_values(
            ["popmax_af_numeric", "qualifying_ac_numeric", "review_score_numeric", "clinvar_id"],
            ascending=[False, False, False, True],
        )
        added = 0
        for _, row in subset.iterrows():
            key = str(row.get("clinvar_id", row.get("variant_key", ""))).strip()
            if key in seen:
                continue
            rows.append(ads_row_from_series(row, source="vital_routing_validation_calls"))
            seen.add(key)
            added += 1
            if added >= limit or len(rows) >= target_n:
                break
        if len(rows) >= target_n:
            break

    ads = pd.DataFrame(rows).head(target_n).reset_index(drop=True)
    if not ads.empty:
        ads.insert(0, "ADS_id", [f"ADS-{index + 1:02d}" for index in range(len(ads))])
    return ads


def build_minimal_repair_logic() -> pd.DataFrame:
    rows = [
        {
            "VITAL_route": "VITAL_OK",
            "meaning": "Label remains compatible with the tested actionability model.",
            "recommended_next_step": "Proceed with standard expert interpretation.",
        },
        {
            "VITAL_route": "CHECK_POPMAX",
            "meaning": "Label may remain valid, but direct actionability requires ancestry-aware frequency review.",
            "recommended_next_step": "Perform ancestry-aware review and penetrance/model reassessment.",
        },
        {
            "VITAL_route": "CHECK_MODEL",
            "meaning": "Label may remain valid, but inheritance, phase, or mechanism is not portable from the flattened label.",
            "recommended_next_step": "Specify inheritance/mechanism before action; use model-specific routing.",
        },
        {
            "VITAL_route": "MODEL_CONFLICT",
            "meaning": "Flattened P/LP label is incompatible with the tested dominant high-penetrance disease model.",
            "recommended_next_step": "Do not use the flattened P/LP label as a dominant-actionability claim.",
        },
        {
            "VITAL_route": "EVAL_LIMITED",
            "meaning": "No allele-resolved evaluability; direct actionability cannot be justified from population evidence.",
            "recommended_next_step": "Defer direct actionability; require representation, callability, or orthogonal evidence.",
        },
        {
            "VITAL_route": "RECESSIVE_COMPATIBLE",
            "meaning": "Dominant model unsupported but carrier or recessive model may remain coherent.",
            "recommended_next_step": "Specify inheritance and phase before routing clinical actionability.",
        },
    ]
    return pd.DataFrame(rows)


def build_reason_code_definitions() -> pd.DataFrame:
    rows = [
        ("EVAL_NO_ALLELE_MATCH", "No allele-resolved population match after normalization/reconciliation."),
        ("EVAL_LOCI_ONLY", "Only locus or regional context is available; allele-level AF is not usable."),
        ("POPMAX_EXCEEDS_MCAF", "Ancestry-aware population maximum exceeds the operational model/review ceiling."),
        ("RECESSIVE_COMPATIBLE_NOT_DOMINANT", "Carrier/recessive logic may remain coherent, but direct dominant actionability is unsupported."),
        ("EXPERT_REVIEW_RETAINED", "No VITAL constraint displaced the label under the tested model."),
        ("MODEL_CONFLICT_DOMINANT_HP", "Observed frequency is incompatible with an unqualified dominant high-penetrance model."),
        ("CHECK_LOW_COUNT", "A frequency signal exists but allele-count support is limited; review is required."),
        ("CHECK_MECHANISM_AMBIGUOUS", "Mechanism or inheritance cannot be inferred from the flattened label."),
    ]
    return pd.DataFrame(rows, columns=["reason_code", "definition"])


def summarize_simulated_cds_alert_layer(
    calls: pd.DataFrame,
    expert_summary: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    alert_calls = calls.loc[calls["baseline_route"].eq("ROUTE_PLP_ACTIONABLE")].copy()

    rows: list[dict[str, object]] = []

    def add_scope(scope: str, frame: pd.DataFrame, note: str) -> None:
        total = len(frame)
        if total == 0:
            return
        alerted = frame.loc[frame["simulated_cds_alert"].ne("NONE")].copy()
        rows.append(
            {
                "scope": scope,
                "metric": "alert_rate",
                "numerator": len(alerted),
                "denominator": total,
                "percent": pct(len(alerted), total),
                "note": note,
            }
        )
        for alert_level in ["YELLOW", "ORANGE", "RED"]:
            count = int(frame["simulated_cds_alert"].eq(alert_level).sum())
            rows.append(
                {
                    "scope": scope,
                    "metric": f"{alert_level.lower()}_alert_rate",
                    "numerator": count,
                    "denominator": total,
                    "percent": pct(count, total),
                    "note": SIMULATED_ALERT_LABELS[alert_level],
                }
            )
        for route in ["CHECK_POPMAX", "CHECK_MODEL", "MODEL_CONFLICT", "EVAL_LIMITED"]:
            count = int(frame["VITAL_constraint_route"].eq(route).sum())
            rows.append(
                {
                    "scope": scope,
                    "metric": f"{route.lower()}_rate",
                    "numerator": count,
                    "denominator": total,
                    "percent": pct(count, total),
                    "note": CONSTRAINT_ROUTE_LABELS[route],
                }
            )

    add_scope(
        "all_arrhythmia_p_lp",
        alert_calls,
        "Simulated clinical decision-support layer across all baseline direct-actionable arrhythmia P/LP labels.",
    )
    add_scope(
        "high_review_arrhythmia_p_lp",
        alert_calls.loc[alert_calls["review_score"].ge(2)].copy(),
        "Restricted to ClinVar review score >= 2.",
    )

    for context_name, genes, _ in CLINICAL_CONTEXTS:
        add_scope(
            f"context_{context_name}",
            alert_calls.loc[alert_calls["gene"].fillna("").astype(str).isin(genes)].copy(),
            f"Domain-specific alert burden for {context_name}.",
        )

    combined = expert_summary.loc[expert_summary["dataset"].eq("combined_expert_plp")].copy()
    if not combined.empty:
        tier1_values = pd.to_numeric(combined["tier1_evaluable_n"], errors="coerce").dropna()
        tier1 = int(tier1_values.iloc[0]) if not tier1_values.empty else 0
        hard_conflict = int(
            pd.to_numeric(
                combined.loc[combined["route_group"].eq("MODEL_CONFLICT_REVIEW"), "count"],
                errors="coerce",
            )
            .fillna(0)
            .sum()
        )
        retained_or_review = int(
            pd.to_numeric(
                combined.loc[combined["route_group"].isin(["VITAL_OK", "CONTROL_REVIEW_FLAG"]), "count"],
                errors="coerce",
            )
            .fillna(0)
            .sum()
        )
        rows.append(
            {
                "scope": "combined_expert_curated_evaluable_plp",
                "metric": "hard_conflict_rate",
                "numerator": hard_conflict,
                "denominator": tier1,
                "percent": pct(hard_conflict, tier1),
                "note": "Hard model-conflict rate among evaluable expert-curated P/LP positives.",
            }
        )
        rows.append(
            {
                "scope": "combined_expert_curated_evaluable_plp",
                "metric": "review_level_retention_rate",
                "numerator": retained_or_review,
                "denominator": tier1,
                "percent": pct(retained_or_review, tier1),
                "note": "Expert-curated positives retained as VITAL_OK or review-level rather than model-conflict.",
            }
        )

    return alert_calls, pd.DataFrame(rows)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.08,
        1.05,
        label,
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


def draw_box(ax: plt.Axes, xy: tuple[float, float], width: float, height: float, text: str, color: str) -> None:
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.018,rounding_size=0.03",
        linewidth=1.0,
        edgecolor="#2d3748",
        facecolor=color,
        alpha=0.95,
    )
    ax.add_patch(patch)
    ax.text(
        xy[0] + width / 2,
        xy[1] + height / 2,
        text,
        ha="center",
        va="center",
        fontsize=8.5,
        color="#101820",
        wrap=True,
    )


def plot_decision_disruption_figure(calls: pd.DataFrame) -> None:
    plot_counts = (
        calls.groupby("VITAL_constraint_route")
        .size()
        .rename("variant_count")
        .reindex(["VITAL_OK", "CHECK_POPMAX", "CHECK_MODEL", "MODEL_CONFLICT", "EVAL_LIMITED"])
        .fillna(0)
        .astype(int)
    )
    total = int(plot_counts.sum())
    nonpass = total - int(plot_counts.get("VITAL_OK", 0))
    nonpass_pct = pct(nonpass, total) or 0.0

    fig, ax = plt.subplots(figsize=(13.2, 6.6))
    ax.axis("off")
    ax.set_title(
        "Label-to-actionability routing audit",
        fontsize=16,
        fontweight="bold",
        pad=16,
    )

    draw_box(ax, (0.03, 0.58), 0.18, 0.18, "ClinVar\nP/LP label", "#eef2f7")
    draw_box(ax, (0.29, 0.58), 0.20, 0.18, "Baseline\nDIRECT-ACTIONABLE", "#c6f6d5")
    draw_box(
        ax,
        (0.57, 0.58),
        0.22,
        0.18,
        "VITAL constraint layer\npopulation + mechanism\n+ evaluability",
        "#e6f0ff",
    )

    ax.annotate("", xy=(0.29, 0.67), xytext=(0.21, 0.67), arrowprops={"arrowstyle": "->", "lw": 2.0, "color": "#2d3748"})
    ax.annotate("", xy=(0.57, 0.67), xytext=(0.49, 0.67), arrowprops={"arrowstyle": "->", "lw": 2.0, "color": "#2d3748"})

    route_specs = [
        ("VITAL_OK", "Proceed", "#2f6f62"),
        ("CHECK_POPMAX", "Review", "#d18f2f"),
        ("CHECK_MODEL", "Model repair", "#4b6cb7"),
        ("MODEL_CONFLICT", "Do not direct-route", "#b94a48"),
        ("EVAL_LIMITED", "Defer direct actionability", "#8b98a5"),
    ]
    y_positions = [0.78, 0.61, 0.44, 0.27, 0.10]
    for (route, action, color), y in zip(route_specs, y_positions):
        count = int(plot_counts.get(route, 0))
        label = f"{route}\n{count}/{total} ({count / total * 100:.1f}%)\n{action}"
        draw_box(ax, (0.84, y), 0.13, 0.12, label, color)
        ax.annotate("", xy=(0.84, y + 0.06), xytext=(0.79, 0.67), arrowprops={"arrowstyle": "->", "lw": 1.5, "color": "#2d3748"})

    ax.text(
        0.5,
        0.02,
        f"Across the arrhythmia cohort, VITAL rerouted {nonpass}/{total} ({nonpass_pct:.1f}%) "
        "label-driven interpretations away from direct actionability; the six-domain mean nonpass rate was 87.7%.",
        ha="center",
        va="bottom",
        fontsize=11,
        color="#1f2933",
    )

    DECISION_DISRUPTION_FIGURE_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(DECISION_DISRUPTION_FIGURE_OUT, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {DECISION_DISRUPTION_FIGURE_OUT.relative_to(BASE_DIR)}")


def draw_workflow_panel(ax: plt.Axes) -> None:
    ax.axis("off")
    add_panel_label(ax, "A")
    ax.set_title("Routing-validation design", fontsize=12, pad=10)
    ax.text(0.16, 0.94, "ClinVar-only baseline", fontsize=10, fontweight="bold", ha="center")
    ax.text(0.73, 0.94, "VITAL-augmented routing", fontsize=10, fontweight="bold", ha="center")

    draw_box(ax, (0.02, 0.63), 0.28, 0.18, "Input:\nvariant + ClinVar P/LP\n+ gene / condition if present", "#edf2f7")
    draw_box(ax, (0.02, 0.28), 0.28, 0.18, "Output:\nROUTE_PLP_ACTIONABLE", "#c6f6d5")
    ax.annotate("", xy=(0.16, 0.47), xytext=(0.16, 0.63), arrowprops={"arrowstyle": "->", "lw": 1.6, "color": "#2d3748"})

    draw_box(
        ax,
        (0.45, 0.63),
        0.48,
        0.18,
        "Input:\nnormalized allele + evaluability\n+ global AF + popmax AF\n+ inheritance/model-aware routing layer",
        "#e6f0ff",
    )
    draw_box(ax, (0.42, 0.31), 0.17, 0.14, "VITAL_OK", "#c6f6d5")
    draw_box(ax, (0.61, 0.31), 0.17, 0.14, "CHECK_POPMAX\nLOW_COUNT", "#f6e0a8")
    draw_box(ax, (0.80, 0.31), 0.17, 0.14, "RECESSIVE /\nMODEL-SPEC", "#cbd5f6")
    draw_box(ax, (0.61, 0.08), 0.17, 0.14, "NOT_EVALUABLE", "#d6dde5")
    for x in [0.505, 0.695, 0.885]:
        ax.annotate("", xy=(x, 0.45), xytext=(0.69, 0.63), arrowprops={"arrowstyle": "->", "lw": 1.4, "color": "#2d3748"})
    ax.annotate("", xy=(0.695, 0.22), xytext=(0.69, 0.63), arrowprops={"arrowstyle": "->", "lw": 1.4, "color": "#2d3748"})


def draw_one_to_many_sankey(ax: plt.Axes, route_counts: pd.DataFrame) -> None:
    ax.axis("off")
    add_panel_label(ax, "B")
    ax.set_title("Counterfactual decision audit", fontsize=12, pad=10)
    counts = (
        route_counts.groupby("vital_route")["variant_count"]
        .sum()
        .reindex(ROUTE_ORDER)
        .fillna(0)
        .astype(int)
    )
    counts = counts[counts > 0]
    total = int(counts.sum())
    gap = 0.008
    usable_height = 0.84 - gap * (len(counts) - 1)
    left_x0 = 0.12
    left_x1 = 0.24
    right_x0 = 0.76
    right_x1 = 0.90
    top = 0.92

    ax.add_patch(Rectangle((left_x0, top - 0.84), left_x1 - left_x0, 0.84, facecolor="#d8f3dc", edgecolor="#2d3748"))
    ax.text(0.02, 0.955, "Decision WITHOUT VITAL", ha="left", va="top", fontsize=9.5, fontweight="bold")
    ax.text(0.02, 0.92, f"ROUTE_PLP_ACTIONABLE\nn={total}", ha="left", va="top", fontsize=10, fontweight="bold")
    ax.text(0.62, 0.955, "Decision WITH VITAL", ha="left", va="top", fontsize=9.5, fontweight="bold")

    left_segments: dict[str, tuple[float, float]] = {}
    right_segments: dict[str, tuple[float, float]] = {}
    current_top = top
    for route, count in counts.items():
        height = usable_height * count / total
        y1 = current_top
        y0 = y1 - height
        left_segments[route] = (y0, y1)
        current_top = y0 - gap

    current_top = top
    label_targets: list[tuple[str, float, float, int]] = []
    for route, count in counts.items():
        height = usable_height * count / total
        y1 = current_top
        y0 = y1 - height
        right_segments[route] = (y0, y1)
        current_top = y0 - gap
        ax.add_patch(
            Rectangle(
                (right_x0, y0),
                right_x1 - right_x0,
                height,
                facecolor=ROUTE_COLORS[route],
                edgecolor="#2d3748",
                linewidth=0.8,
            )
        )
        label_targets.append((route, (y0 + y1) / 2, height, count))

    def ribbon(left_bounds: tuple[float, float], right_bounds: tuple[float, float], color: str) -> None:
        l0, l1 = left_bounds
        r0, r1 = right_bounds
        verts = [
            (left_x1, l1),
            (0.45, l1),
            (0.55, r1),
            (right_x0, r1),
            (right_x0, r0),
            (0.55, r0),
            (0.45, l0),
            (left_x1, l0),
            (left_x1, l1),
        ]
        codes = [
            MplPath.MOVETO,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.LINETO,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CLOSEPOLY,
        ]
        patch = PathPatch(MplPath(verts, codes), facecolor=color, edgecolor="none", alpha=0.7)
        ax.add_patch(patch)

    for route in counts.index:
        ribbon(left_segments[route], right_segments[route], ROUTE_COLORS[route])

    large_routes = {"VITAL_OK", "NOT_EVALUABLE_NO_FREQUENCY_CLAIM"}
    minor_lines: list[str] = []
    for route, desired, _, count in label_targets:
        if route in large_routes:
            ax.plot([right_x1, right_x1 + 0.015], [desired, desired], color="#4a5568", linewidth=0.7)
            ax.text(
                right_x1 + 0.022,
                desired,
                f"{ROUTE_LABELS[route]}\n{count} ({count / total * 100:.1f}%)",
                fontsize=8.8,
                ha="left",
                va="center",
            )
        else:
            minor_lines.append(f"{ROUTE_LABELS[route]}  {count} ({count / total * 100:.1f}%)")
    if minor_lines:
        ax.text(
            0.59,
            0.86,
            "Minor routing branches\n" + "\n".join(minor_lines),
            fontsize=8.5,
            ha="left",
            va="top",
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "#ffffff", "edgecolor": "#cbd5e0", "alpha": 0.95},
        )


def draw_context_panel(ax: plt.Axes, context_summary: pd.DataFrame) -> None:
    add_panel_label(ax, "C")
    ax.set_title("Incorrectly escalated decisions by workflow context", fontsize=12, pad=10)
    plot_df = context_summary.loc[
        context_summary["clinical_context"].isin(["cascade_testing", "drug_guidance", "surveillance_device"])
    ].copy()
    plot_df = plot_df.set_index("clinical_context").loc[
        ["cascade_testing", "drug_guidance", "surveillance_device"]
    ]
    x = np.arange(len(plot_df))
    total = plot_df["n_total"].to_numpy()
    unchanged = plot_df["n_vital_ok"].to_numpy() / total * 100
    evaluation_limiting = plot_df["n_evaluation_limiting"].to_numpy() / total * 100
    risk_escalating = plot_df["n_risk_escalating"].to_numpy() / total * 100
    model_correcting = plot_df["n_model_correcting"].to_numpy() / total * 100
    ax.bar(x, unchanged, color=CHANGE_CLASS_COLORS["unchanged"], label="Unchanged")
    ax.bar(x, evaluation_limiting, bottom=unchanged, color=CHANGE_CLASS_COLORS["evaluation_limiting"], label="Evaluation-limiting")
    ax.bar(
        x,
        risk_escalating,
        bottom=unchanged + evaluation_limiting,
        color=CHANGE_CLASS_COLORS["risk_escalating"],
        label="Risk-escalating",
    )
    ax.bar(
        x,
        model_correcting,
        bottom=unchanged + evaluation_limiting + risk_escalating,
        color=CHANGE_CLASS_COLORS["model_correcting"],
        label="Model-correcting",
    )
    for index, risk_percent in enumerate(plot_df["actionability_at_risk_percent"].to_numpy()):
        ax.text(index, 101.5, f"{risk_percent:.1f}%", ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(["Cascade", "Drug", "Surveillance/device"])
    ax.set_ylabel("Percent of baseline actionable decisions")
    ax.set_ylim(0, 108)
    ax.grid(axis="y", color="#d8dde6", linewidth=0.8, alpha=0.7)
    ax.legend(frameon=False, loc="upper right", fontsize=8)


def draw_expert_panel(ax: plt.Axes, expert_summary: pd.DataFrame) -> None:
    add_panel_label(ax, "D")
    ax.set_title("Expert-panel concordance and control behavior", fontsize=12, pad=10)
    dataset_order = ["arrhythmia_expert_plp", "external_expert_plp", "arrhythmia_blb_controls"]
    labels = ["Arrhythmia expert P/LP", "External expert P/LP", "B/LB controls"]
    route_groups = [
        ("VITAL_OK", "VITAL_OK", "#2f6f62"),
        ("CONTROL_REVIEW_FLAG", "Review", "#d18f2f"),
        ("MODEL_CONFLICT_REVIEW", "Model conflict", "#b94a48"),
        ("NOT_EVALUABLE_NO_FREQUENCY_CLAIM", "Not evaluable", "#8b98a5"),
        ("MODEL_CONFLICT_CONTROL_FLAG", "Model conflict", "#b94a48"),
    ]
    x = np.arange(len(dataset_order))
    bottoms = np.zeros(len(dataset_order), dtype=float)
    for route_group, label, color in route_groups:
        heights = []
        for dataset in dataset_order:
            row = expert_summary.loc[
                expert_summary["dataset"].eq(dataset) & expert_summary["route_group"].eq(route_group),
                "percent_of_dataset",
            ]
            heights.append(float(row.iloc[0]) if not row.empty else 0.0)
        if max(heights) == 0:
            continue
        ax.bar(x, heights, bottom=bottoms, color=color, label=label)
        bottoms += np.array(heights)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Percent within validation set")
    ax.set_ylim(0, 105)
    ax.grid(axis="y", color="#d8dde6", linewidth=0.8, alpha=0.7)
    handles, legend_labels = ax.get_legend_handles_labels()
    unique_handles: dict[str, object] = {}
    for handle, legend_label in zip(handles, legend_labels):
        unique_handles.setdefault(legend_label, handle)
    ax.legend(unique_handles.values(), unique_handles.keys(), frameon=False, loc="upper right", fontsize=8)


def plot_routing_validation_figure(
    route_counts: pd.DataFrame,
    context_summary: pd.DataFrame,
    expert_summary: pd.DataFrame,
) -> None:
    fig = plt.figure(figsize=(14.2, 10.6))
    gs = fig.add_gridspec(2, 2, hspace=0.28, wspace=0.2)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    draw_workflow_panel(ax_a)
    draw_one_to_many_sankey(ax_b, route_counts)
    draw_context_panel(ax_c, context_summary)
    draw_expert_panel(ax_d, expert_summary)

    fig.suptitle(
        "Counterfactual decision audit: VITAL changes downstream routing of flattened P/LP labels",
        fontsize=15,
        fontweight="bold",
        y=0.98,
    )
    ROUTING_FIGURE_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(ROUTING_FIGURE_OUT, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {ROUTING_FIGURE_OUT.relative_to(BASE_DIR)}")


def main() -> None:
    baseline_model = build_formal_baseline_decision_model()
    save_table(baseline_model, BASELINE_DECISION_MODEL_OUT)

    arrhythmia_calls = load_arrhythmia_routing_calls()
    save_table(arrhythmia_calls, ROUTING_CALLS_OUT)
    save_table(arrhythmia_calls, COUNTERFACTUAL_AUDIT_CALLS_OUT)

    route_counts, routing_summary, context_summary = summarize_main_cohort(arrhythmia_calls)
    save_table(route_counts, ROUTING_TRANSITIONS_OUT)
    save_table(routing_summary, ROUTING_SUMMARY_OUT)
    save_table(context_summary, ROUTING_CONTEXT_OUT)

    counterfactual_summary = summarize_counterfactual_decision_audit(arrhythmia_calls)
    save_table(counterfactual_summary, COUNTERFACTUAL_AUDIT_SUMMARY_OUT)

    high_review_summary = summarize_high_review_subset(arrhythmia_calls)
    if not high_review_summary.empty:
        save_table(high_review_summary, ROUTING_HIGH_REVIEW_OUT)

    expert_summary, expert_arrhythmia_plp, _ = summarize_expert_concordance()
    save_table(expert_summary, ROUTING_EXPERT_OUT)

    cds_alert_calls, cds_alert_summary = summarize_simulated_cds_alert_layer(arrhythmia_calls, expert_summary)
    save_table(cds_alert_calls, SIMULATED_CDS_ALERT_CALLS_OUT)
    save_table(cds_alert_summary, SIMULATED_CDS_ALERT_SUMMARY_OUT)

    case_vignettes = build_case_vignettes(arrhythmia_calls, expert_arrhythmia_plp)
    save_table(case_vignettes, ROUTING_VIGNETTES_OUT)
    ads = build_actionability_discordance_audit(arrhythmia_calls, case_vignettes)
    save_table(ads, ADS_OUT)

    save_table(build_minimal_repair_logic(), REPAIR_LOGIC_OUT)
    save_table(build_reason_code_definitions(), REASON_CODE_DEFINITIONS_OUT)

    plot_decision_disruption_figure(arrhythmia_calls)
    plot_routing_validation_figure(route_counts, context_summary, expert_summary)


if __name__ == "__main__":
    main()
