from __future__ import annotations

from pathlib import Path

import pandas as pd


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
SUPPLEMENT_DIR = BASE_DIR / "supplementary_tables"

SCORES_IN = DATA_DIR / "arrhythmia_vital_scores.csv"
TEMPORAL_IN = DATA_DIR / "arrhythmia_temporal_robustness_vital_scored_variants.tsv"

PILOT_KEY_OUT = DATA_DIR / "vital_blinded_expert_pilot_key.csv"
PILOT_FORM_OUT = DATA_DIR / "vital_blinded_expert_review_form.csv"
PILOT_PROTOCOL_OUT = DATA_DIR / "vital_independent_validation_protocol.csv"
AC_STRATA_OUT = DATA_DIR / "vital_ac_reliability_strata.csv"
AC_SENSITIVITY_OUT = DATA_DIR / "vital_ac_threshold_sensitivity_extended.csv"
DETECTABILITY_OUT = DATA_DIR / "vital_detectability_required_fields.csv"
TEMPORAL_TRACKING_OUT = DATA_DIR / "vital_temporal_red_tracking_summary.csv"
CONTEXT_OVERLAYS_OUT = DATA_DIR / "vital_contextual_overlay_modules.csv"

PILOT_KEY_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S34_blinded_expert_pilot_key.tsv"
PILOT_FORM_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S35_blinded_expert_review_form.tsv"
AC_STRATA_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S36_ac_reliability_strata.tsv"
DETECTABILITY_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S37_detectability_required_fields.tsv"
TEMPORAL_TRACKING_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S38_temporal_red_tracking.tsv"
CONTEXT_OVERLAYS_SUPP = SUPPLEMENT_DIR / "Supplementary_Table_S39_contextual_overlay_modules.tsv"


def read_scores() -> pd.DataFrame:
    if not SCORES_IN.exists():
        raise FileNotFoundError(f"Missing input: {SCORES_IN}")
    df = pd.read_csv(SCORES_IN)
    for col in [
        "vital_score",
        "global_af",
        "global_ac",
        "popmax_af",
        "popmax_ac",
        "qualifying_frequency_ac",
        "submitter_count",
    ]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in [
        "vital_red_flag",
        "standard_acmg_frequency_flag",
        "frequency_signal_ac_ge_20",
        "weak_review_signal",
    ]:
        if col in df.columns:
            df[col] = df[col].astype(str).str.lower().isin(["true", "1", "yes"])
    return df


def ac_stratum(ac: float | int | None) -> str:
    if pd.isna(ac) or float(ac) < 10:
        return "AC <10 or unavailable"
    ac = float(ac)
    if ac < 20:
        return "AC 10-19"
    if ac < 50:
        return "AC 20-49"
    if ac < 200:
        return "AC 50-199"
    return "AC >=200"


def public_variant_columns(df: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "gene",
        "clinvar_id",
        "variation_id",
        "title",
        "clinsig",
        "review_status",
        "review_strength",
        "match_category",
        "frequency_evidence_status",
        "variant_type",
        "functional_class",
        "global_af",
        "global_ac",
        "popmax_af",
        "popmax_ac",
        "popmax_population",
        "qualifying_frequency_ac",
    ]
    return df[[c for c in cols if c in df.columns]].copy()


def select_blinded_pilot(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    red = df[df["vital_red_flag"]].sort_values(["vital_score"], ascending=False).copy()
    red["pilot_stratum"] = "red_priority"

    gray_pool = df[df["vital_band"].eq("gray_no_frequency_evidence")].copy()
    gray_priority = gray_pool.assign(
        complex_rank=gray_pool["variant_type"].isin(["duplication", "deletion", "insertion"]).astype(int),
        gene_rank=gray_pool["gene"].isin(["KCNH2", "SCN5A", "KCNQ1", "RYR2"]).astype(int),
    ).sort_values(
        ["complex_rank", "gene_rank", "gene", "variant_type", "clinvar_id"],
        ascending=[False, False, True, True, True],
    )
    gray = select_diverse(gray_priority, n=10, max_per_gene=2).copy()
    gray["pilot_stratum"] = "gray_no_frequency_evidence"

    controls_pool = df[
        df["frequency_evidence_status"].eq("frequency_observed")
        & ~df["vital_red_flag"]
        & ~df["standard_acmg_frequency_flag"]
    ].copy()
    controls_priority = controls_pool.sort_values(
        ["review_strength", "gene", "variant_type", "clinvar_id"],
        ascending=[True, True, True, True],
    )
    controls = select_diverse(controls_priority, n=10, max_per_gene=2).copy()
    controls["pilot_stratum"] = "non_red_frequency_consistent_control"

    pilot = pd.concat([red, gray, controls], ignore_index=True)
    pilot.insert(0, "pilot_variant_id", [f"PILOT-{i:03d}" for i in range(1, len(pilot) + 1)])
    pilot["pre_specified_endpoint"] = "expert_requires_re_review"
    pilot["expert_consensus_rule"] = "two_of_three_reviewers_or_adjudication_panel"
    pilot["planned_metrics"] = (
        "sensitivity_specificity_ppv_for_red_and_gray_routing; "
        "cohen_or_fleiss_kappa; discordance_review"
    )
    pilot["vital_route_for_analysis"] = pilot["pilot_stratum"].map(
        {
            "red_priority": "requires_immediate_state_aware_re_review",
            "gray_no_frequency_evidence": "requires_detectability_or_deferred_review",
            "non_red_frequency_consistent_control": "no_frequency_triggered_re_review_expected",
        }
    )

    key_cols = [
        "pilot_variant_id",
        "pilot_stratum",
        "vital_route_for_analysis",
        "pre_specified_endpoint",
        "expert_consensus_rule",
        "planned_metrics",
        "vital_score",
        "vital_band",
        "vital_red_flag",
        "vital_signal_reason",
    ]
    key = pd.concat([pilot[key_cols], public_variant_columns(pilot)], axis=1)

    form = public_variant_columns(pilot)
    form.insert(0, "pilot_variant_id", pilot["pilot_variant_id"])
    form["expert_requires_re_review"] = ""
    form["expert_route"] = ""
    form["expert_confidence_1_to_5"] = ""
    form["expert_comments"] = ""
    form["discordance_adjudication_notes"] = ""
    return key, form


def select_diverse(df: pd.DataFrame, n: int, max_per_gene: int) -> pd.DataFrame:
    """Deterministically prefer cross-gene diversity for small expert-review pilots."""
    selected_indices: list[int] = []
    gene_counts: dict[str, int] = {}
    type_counts: dict[str, int] = {}

    for idx, row in df.iterrows():
        gene = str(row.get("gene", ""))
        variant_type = str(row.get("variant_type", ""))
        if gene_counts.get(gene, 0) >= max_per_gene:
            continue
        if type_counts.get(variant_type, 0) >= max(2, n // 2):
            continue
        selected_indices.append(idx)
        gene_counts[gene] = gene_counts.get(gene, 0) + 1
        type_counts[variant_type] = type_counts.get(variant_type, 0) + 1
        if len(selected_indices) == n:
            return df.loc[selected_indices]

    for idx, row in df.iterrows():
        if idx in selected_indices:
            continue
        gene = str(row.get("gene", ""))
        if gene_counts.get(gene, 0) >= max_per_gene + 1:
            continue
        selected_indices.append(idx)
        gene_counts[gene] = gene_counts.get(gene, 0) + 1
        if len(selected_indices) == n:
            break

    return df.loc[selected_indices]


def build_protocol() -> pd.DataFrame:
    rows = [
        {
            "protocol_element": "external_endpoint",
            "pre_specified_value": "requires re-review / does not require re-review",
            "rationale": "Separates external adjudication from internal score construction.",
        },
        {
            "protocol_element": "pilot_set_size",
            "pre_specified_value": "23 variants: 3 red, 10 gray, 10 non-red controls",
            "rationale": "Small blinded pilot designed to test routing plausibility before diagnostic claims.",
        },
        {
            "protocol_element": "reviewer_blinding",
            "pre_specified_value": "review form excludes score, band, red flag, and pilot stratum",
            "rationale": "Prevents circular confirmation of the route by reviewers.",
        },
        {
            "protocol_element": "consensus_rule",
            "pre_specified_value": "two of three reviewers, with adjudication for discordance",
            "rationale": "Produces a reproducible external route label.",
        },
        {
            "protocol_element": "metrics",
            "pre_specified_value": "sensitivity, specificity, PPV for red/gray routing; Cohen/Fleiss kappa; discordance review",
            "rationale": "Treats the score as review-routing support, not as a diagnostic classifier.",
        },
        {
            "protocol_element": "external_retrospective_sources",
            "pre_specified_value": "ClinGen/VCEP reclassification histories, LOVD/registry cases, national consortium candidates when available",
            "rationale": "Defines the non-internal validation path without changing the current core score.",
        },
    ]
    return pd.DataFrame(rows)


def build_ac_tables(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = df.copy()
    df["ac_reliability_stratum"] = df["qualifying_frequency_ac"].apply(ac_stratum)
    frequency_flagged = df[df["standard_acmg_frequency_flag"]].copy()
    order = ["AC <10 or unavailable", "AC 10-19", "AC 20-49", "AC 50-199", "AC >=200"]
    rows = []
    for stratum in order:
        sub = frequency_flagged[frequency_flagged["ac_reliability_stratum"].eq(stratum)]
        red = sub[sub["vital_red_flag"]]
        rows.append(
            {
                "ac_reliability_stratum": stratum,
                "frequency_flagged_variants": len(sub),
                "vital_red_variants": len(red),
                "red_variant_ids": "; ".join(red["gene"] + ":" + red["clinvar_id"]),
                "interpretation": interpretation_for_ac_stratum(stratum),
            }
        )
    strata = pd.DataFrame(rows)

    sensitivity_rows = []
    for threshold in [10, 20, 50, 100, 200]:
        red = df[
            df["standard_acmg_frequency_flag"]
            & df["weak_review_signal"]
            & (df["vital_score"] >= 70)
            & (df["qualifying_frequency_ac"] >= threshold)
        ].copy()
        sensitivity_rows.append(
            {
                "ac_threshold": f"AC >= {threshold}",
                "red_queue_count": len(red),
                "red_queue_variants": "; ".join(
                    red.sort_values("vital_score", ascending=False)["gene"] + ":" + red["clinvar_id"]
                ),
                "composition_note": composition_note(threshold, red),
            }
        )
    sensitivity = pd.DataFrame(sensitivity_rows)
    return strata, sensitivity


def interpretation_for_ac_stratum(stratum: str) -> str:
    if stratum == "AC <10 or unavailable":
        return "Low reliability; do not use as urgent actionability evidence without additional support."
    if stratum == "AC 10-19":
        return "Intermediate-low reliability; useful for watchlist but not default red gate."
    if stratum == "AC 20-49":
        return "Default actionability floor; enough to avoid one-allele noise but still near-boundary."
    if stratum == "AC 50-199":
        return "Stronger frequency reliability; supports higher confidence in frequency tension."
    return "Very strong frequency reliability; frequency signal should dominate review prioritization."


def composition_note(threshold: int, red: pd.DataFrame) -> str:
    ids = set(red["clinvar_id"].astype(str))
    if threshold < 20:
        return "Lower than default; admits lower-AC borderline calls such as CACNB2 when present."
    if threshold == 20:
        return "Default compromise gate retaining SCN5A, TRDN, and borderline KCNH2."
    if threshold >= 50 and ids == {"VCV000440850"}:
        return "High-reliability gate retains only SCN5A and drops TRDN/KCNH2 near-boundary AC cases."
    if len(red) == 0:
        return "Very stringent gate removes all current red-priority cases."
    return "Stringent gate narrows the queue to high-AC cases."


def build_detectability_table() -> pd.DataFrame:
    rows = [
        {
            "field_name": "gnomAD_exome_present",
            "required_for": "all gray indel/duplication or exact-match-without-AF cases",
            "allowed_values": "yes/no/unknown",
            "action_if_unavailable": "retain gray state; do not impute AF=0",
        },
        {
            "field_name": "gnomAD_genome_present",
            "required_for": "all gray indel/duplication or structural/complex alleles",
            "allowed_values": "yes/no/unknown",
            "action_if_unavailable": "request genome-level requery or defer to future freeze",
        },
        {
            "field_name": "coverage_sufficient",
            "required_for": "red or gray calls where absence/low AF is clinically material",
            "allowed_values": "yes/no/unknown/not_applicable",
            "action_if_unavailable": "route to detectability review before using absence as rarity evidence",
        },
        {
            "field_name": "mapping_quality_concern",
            "required_for": "indels, duplications, repeats, homopolymers, and complex alleles",
            "allowed_values": "yes/no/unknown",
            "action_if_unavailable": "treat as representation-sensitive; consider orthogonal validation",
        },
        {
            "field_name": "variant_representation_ambiguity",
            "required_for": "allele-discordant, same-site, multiallelic, and haplotype assertions",
            "allowed_values": "yes/no/unknown",
            "action_if_unavailable": "normalize/left-align/phase before frequency interpretation",
        },
        {
            "field_name": "detectability_assessment_complete",
            "required_for": "gray-to-actionable routing",
            "allowed_values": "yes/no",
            "action_if_unavailable": "variant cannot leave gray state on frequency evidence alone",
        },
    ]
    return pd.DataFrame(rows)


def build_temporal_tracking() -> pd.DataFrame:
    if not TEMPORAL_IN.exists():
        return pd.DataFrame(
            [
                {
                    "tracking_status": "template_only",
                    "reason": f"Missing {TEMPORAL_IN.name}",
                }
            ]
        )
    temporal = pd.read_csv(TEMPORAL_IN, sep="\t")
    red_or_history = temporal[
        temporal.get("is_red", False).astype(str).str.lower().isin(["true", "1"])
        | temporal.get("clinvar_id", "").astype(str).isin(["VCV000440850", "VCV001325231", "VCV004535537", "VCV001202620"])
    ].copy()
    cols = [
        "snapshot",
        "gene",
        "clinvar_id",
        "title",
        "clinsig",
        "vital_score",
        "vital_band",
        "is_red",
        "global_ac",
        "popmax_ac",
        "qualifying_frequency_ac",
        "review_status",
        "inheritance_flag",
    ]
    red_or_history = red_or_history[[c for c in cols if c in red_or_history.columns]].copy()
    red_or_history["exit_reason_if_removed_from_red_queue"] = red_or_history.apply(exit_reason, axis=1)
    red_or_history["tracking_requirement"] = (
        "quarterly ClinVar/gnomAD freeze check; record status, AC, review status, representation, and transcript changes"
    )
    return red_or_history.sort_values(["gene", "clinvar_id", "snapshot"])


def exit_reason(row: pd.Series) -> str:
    if str(row.get("is_red", "")).lower() in ["true", "1"]:
        return "still_in_red_queue_at_snapshot"
    return "not_red_at_snapshot; classify cause as ClinVar_status_change, gnomAD_AC_change, representation_change, transcript_change, or review_evidence_change"


def build_context_overlays() -> pd.DataFrame:
    rows = [
        {
            "overlay": "functional_flag",
            "examples": "MAVE/MaveDB, patch clamp, splicing assay, minigene",
            "role": "interpretive overlay for red cases",
            "score_policy": "not part of the core score",
        },
        {
            "overlay": "segregation_flag",
            "examples": "de novo, cosegregation, phase, second allele",
            "role": "case-level adjudication after routing",
            "score_policy": "not part of the core score",
        },
        {
            "overlay": "transcript_flag",
            "examples": "MANE Select, clinically relevant transcript, alternative transcript only, NMD rescue",
            "role": "mechanism qualification and annotation-inflation review",
            "score_policy": "not part of the core score",
        },
    ]
    return pd.DataFrame(rows)


def save(df: pd.DataFrame, path: Path, sep: str = ",") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep=sep, index=False, lineterminator="\n")
    print(f"Saved {path.relative_to(BASE_DIR)} ({len(df)} rows)")


def main() -> None:
    df = read_scores()
    key, form = select_blinded_pilot(df)
    protocol = build_protocol()
    strata, sensitivity = build_ac_tables(df)
    detectability = build_detectability_table()
    temporal = build_temporal_tracking()
    overlays = build_context_overlays()

    save(key, PILOT_KEY_OUT)
    save(form, PILOT_FORM_OUT)
    save(protocol, PILOT_PROTOCOL_OUT)
    save(strata, AC_STRATA_OUT)
    save(sensitivity, AC_SENSITIVITY_OUT)
    save(detectability, DETECTABILITY_OUT)
    save(temporal, TEMPORAL_TRACKING_OUT)
    save(overlays, CONTEXT_OVERLAYS_OUT)

    save(key, PILOT_KEY_SUPP, sep="\t")
    save(form, PILOT_FORM_SUPP, sep="\t")
    save(strata, AC_STRATA_SUPP, sep="\t")
    save(detectability, DETECTABILITY_SUPP, sep="\t")
    save(temporal, TEMPORAL_TRACKING_SUPP, sep="\t")
    save(overlays, CONTEXT_OVERLAYS_SUPP, sep="\t")


if __name__ == "__main__":
    main()
