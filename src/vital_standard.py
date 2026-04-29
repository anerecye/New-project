from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"

TIER2_RECONCILIATION_TIERS = {
    "tier3_locus_observed_no_exact_af",
    "tier3_regional_context_only",
}

HIGH_CONFIDENCE_RECESSIVE_GENES = {
    "CASQ2",
    "TRDN",
    "PAH",
    "GAA",
    "GALT",
    "SLC22A5",
    "ACADVL",
    "BTD",
    "ASS1",
    "MMACHC",
    "ASL",
    "HMGCL",
    "LPL",
    "SLC26A4",
    "OTOF",
    "CDH23",
    "USH2A",
    "MYO15A",
}

SV_PRIORITY_GENES = {
    "SMN1": "copy-number and paralog resolution required",
    "STRC": "copy-number and gene-conversion resolution required",
    "GJB6": "deletion/haplotype resolution required",
}

SV_KEYWORDS = (
    "copy number",
    "gene conversion",
    "whole gene",
    "whole-gene",
    "exon deletion",
    "exon duplication",
    "deletion encompassing",
)


@dataclass(frozen=True)
class DatasetSpec:
    dataset: str
    domain: str
    path: Path
    include_in_meta: bool = True


DEFAULT_DATASET_SPECS = (
    DatasetSpec(
        dataset="arrhythmia",
        domain="arrhythmia",
        path=DATA_DIR / "arrhythmia_vital_scores.csv",
    ),
    DatasetSpec(
        dataset="cardiomyopathy",
        domain="cardiomyopathy",
        path=DATA_DIR / "vital_domain_cardiomyopathy_vital_scores.csv",
    ),
    DatasetSpec(
        dataset="epilepsy",
        domain="epilepsy",
        path=DATA_DIR / "vital_domain_epilepsy_vital_scores.csv",
    ),
    DatasetSpec(
        dataset="metabolic_ar",
        domain="metabolic_ar",
        path=DATA_DIR / "vital_domain_metabolic_ar_vital_scores.csv",
    ),
    DatasetSpec(
        dataset="primary_immunodeficiency",
        domain="primary_immunodeficiency",
        path=DATA_DIR / "vital_domain_primary_immunodeficiency_vital_scores.csv",
    ),
    DatasetSpec(
        dataset="hearing_loss",
        domain="hearing_loss",
        path=DATA_DIR / "vital_domain_hearing_loss_vital_scores.csv",
    ),
)


def truthy(value: object) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def parse_gene_argument(raw: object) -> list[str]:
    if raw is None:
        return []
    if isinstance(raw, (list, tuple, set)):
        tokens = [str(item) for item in raw]
    else:
        tokens = str(raw).replace(";", ",").split(",")
    genes = []
    for token in tokens:
        normalized = token.strip().upper()
        if normalized:
            genes.append(normalized)
    return sorted(set(genes))


def infer_evaluability(row: pd.Series) -> str:
    frequency_status = str(row.get("frequency_evidence_status", "")).strip()
    tier = str(row.get("reconciliation_tier", "")).strip()
    match_category = str(row.get("match_category", "")).strip()

    if frequency_status == "frequency_observed":
        return "TIER_1"
    if tier in {"tier4_still_unevaluable", "tier4_query_error"}:
        return "UNEVALUABLE"
    if tier in TIER2_RECONCILIATION_TIERS:
        return "TIER_2"
    if frequency_status in {"exact_match_without_AF", "allele_discordance_no_exact_AF"}:
        return "TIER_2"
    if match_category == "allele_discordance":
        return "TIER_2"
    if frequency_status in {"not_observed_in_gnomAD", "gnomad_query_error_no_frequency_evidence"}:
        return "UNEVALUABLE"
    if match_category in {"no_gnomad_record", "gnomad_query_error"}:
        return "UNEVALUABLE"
    return "UNEVALUABLE"


def infer_flag(row: pd.Series) -> str:
    evaluability = row.get("VITAL_evaluability", infer_evaluability(row))
    if evaluability != "TIER_1":
        return "."
    if not truthy(row.get("standard_acmg_frequency_flag")):
        return "OK"
    gene = str(row.get("gene", "")).strip().upper()
    global_af = pd.to_numeric(row.get("global_af"), errors="coerce")
    if (
        gene not in HIGH_CONFIDENCE_RECESSIVE_GENES
        and truthy(row.get("frequency_signal_ac_ge_20"))
        and pd.notna(global_af)
        and float(global_af) > 1e-4
    ):
        return "MODEL_CONFLICT"
    return "CHECK_POPMAX"


def infer_sv_required(row: pd.Series) -> bool:
    gene = str(row.get("gene", "")).strip().upper()
    if gene in SV_PRIORITY_GENES:
        return True
    title = str(row.get("title", "")).strip().lower()
    return any(keyword in title for keyword in SV_KEYWORDS)


def infer_certification(row: pd.Series) -> str:
    evaluability = row.get("VITAL_evaluability", infer_evaluability(row))
    flag = row.get("VITAL_flag", infer_flag(row))
    gene = str(row.get("gene", "")).strip().upper()

    if evaluability == "TIER_2":
        return "VITAL-2"
    if evaluability != "TIER_1":
        return "VITAL-0"
    if gene in HIGH_CONFIDENCE_RECESSIVE_GENES and flag != "OK":
        return "VITAL-X"
    if flag == "OK":
        return "VITAL-1"
    return "VITAL-ALERT"


def infer_alert_state(row: pd.Series) -> str:
    certification = row.get("VITAL_certification", infer_certification(row))
    flag = row.get("VITAL_flag", infer_flag(row))

    if truthy(row.get("VITAL_sv_required", infer_sv_required(row))):
        return "SV_LAYER_REQUIRED"
    if certification == "VITAL-1":
        return "NONE"
    if certification == "VITAL-2":
        return "LOCUS_LEVEL_ONLY"
    if certification == "VITAL-0":
        return "ALLELE_UNRESOLVED"
    if certification == "VITAL-X":
        return "RECESSIVE_ROUTE_REQUIRED"
    if flag == "MODEL_CONFLICT":
        return "MODEL_CONFLICT"
    if truthy(row.get("frequency_signal_ac_ge_20")):
        return "POPMAX_REVIEW"
    return "POPULATION_SIGNAL_CHECK"


def infer_public_use(row: pd.Series) -> str:
    alert_state = row.get("VITAL_alert", infer_alert_state(row))
    certification = row.get("VITAL_certification", infer_certification(row))
    if alert_state == "SV_LAYER_REQUIRED":
        return "specialized_sv_review_required"
    if certification == "VITAL-1":
        return "clinical_portable"
    if certification == "VITAL-2":
        return "research_only_locus_level"
    if certification == "VITAL-X":
        return "recessive_or_phase_defined_route"
    return "normalize_or_review_before_public_reuse"


def build_standard_view(frame: pd.DataFrame, dataset: str, domain: str) -> pd.DataFrame:
    result = frame.copy()
    result["dataset"] = dataset
    result["domain"] = domain
    if "gene" in result.columns:
        gene_series = result["gene"]
    else:
        gene_series = pd.Series("", index=result.index, dtype="object")
    result["gene"] = gene_series.fillna("").astype(str).str.upper()
    result["VITAL_evaluability"] = result.apply(infer_evaluability, axis=1)
    result["VITAL_flag"] = result.apply(infer_flag, axis=1)
    result["VITAL_sv_required"] = result.apply(infer_sv_required, axis=1)
    result["VITAL_certification"] = result.apply(infer_certification, axis=1)
    result["VITAL_alert"] = result.apply(infer_alert_state, axis=1)
    result["VITAL_public_use"] = result.apply(infer_public_use, axis=1)
    result["VITAL_pass"] = result["VITAL_certification"].eq("VITAL-1") & ~result["VITAL_sv_required"]
    result["VITAL_alerted"] = result["VITAL_alert"].ne("NONE")
    result["VITAL_problem_evaluation"] = result["VITAL_certification"].isin(["VITAL-2", "VITAL-0"])
    result["VITAL_problem_popmax"] = result["VITAL_alert"].isin(
        ["POPMAX_REVIEW", "POPULATION_SIGNAL_CHECK", "MODEL_CONFLICT"]
    )
    result["VITAL_problem_model"] = result["VITAL_alert"].eq("MODEL_CONFLICT")
    result["VITAL_problem_recessive"] = result["VITAL_certification"].eq("VITAL-X")
    result["VITAL_problem_sv"] = result["VITAL_alert"].eq("SV_LAYER_REQUIRED")
    return result


def load_processed_prediction_tables(include_auxiliary: bool = False) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    specs = list(DEFAULT_DATASET_SPECS)
    if include_auxiliary:
        specs.extend(
            [
                DatasetSpec(
                    dataset="random_clinvar_plp",
                    domain="random_clinvar_plp",
                    path=DATA_DIR / "vital_random_clinvar_plp_vital_scores.csv",
                    include_in_meta=False,
                ),
                DatasetSpec(
                    dataset="expert_arrhythmia",
                    domain="expert_arrhythmia",
                    path=DATA_DIR / "arrhythmia_expert_consensus_validation_vital_scores.csv",
                    include_in_meta=False,
                ),
            ]
        )

    for spec in specs:
        if not spec.path.exists():
            continue
        frame = pd.read_csv(spec.path, low_memory=False)
        standard = build_standard_view(frame, dataset=spec.dataset, domain=spec.domain)
        standard["include_in_meta"] = spec.include_in_meta
        frames.append(standard)

    if not frames:
        raise FileNotFoundError("No processed VITAL prediction tables were found.")
    return pd.concat(frames, ignore_index=True)


def filter_standard_view(frame: pd.DataFrame, genes: object) -> pd.DataFrame:
    selected_genes = parse_gene_argument(genes)
    if not selected_genes:
        return frame.copy()
    return frame.loc[frame["gene"].isin(selected_genes)].copy()


def summarize_standard_view(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if frame.empty:
        return pd.DataFrame(
            columns=[
                "dataset",
                "domain",
                "variant_count",
                "gene_count",
                "vital_1_count",
                "vital_alert_count",
                "vital_x_count",
                "vital_2_count",
                "vital_0_count",
                "alert_count",
                "pass_rate",
                "nonpass_rate",
            ]
        )

    for (dataset, domain), sub in frame.groupby(["dataset", "domain"], dropna=False):
        variant_count = len(sub)
        pass_count = int(sub["VITAL_pass"].sum())
        alert_count = int(sub["VITAL_alerted"].sum())
        rows.append(
            {
                "dataset": dataset,
                "domain": domain,
                "variant_count": variant_count,
                "gene_count": int(sub["gene"].nunique()),
                "vital_1_count": int(sub["VITAL_certification"].eq("VITAL-1").sum()),
                "vital_alert_count": int(sub["VITAL_certification"].eq("VITAL-ALERT").sum()),
                "vital_x_count": int(sub["VITAL_certification"].eq("VITAL-X").sum()),
                "vital_2_count": int(sub["VITAL_certification"].eq("VITAL-2").sum()),
                "vital_0_count": int(sub["VITAL_certification"].eq("VITAL-0").sum()),
                "sv_required_count": int(sub["VITAL_problem_sv"].sum()),
                "alert_count": alert_count,
                "pass_count": pass_count,
                "nonpass_count": variant_count - pass_count,
                "pass_rate": pass_count / variant_count if variant_count else 0.0,
                "nonpass_rate": (variant_count - pass_count) / variant_count if variant_count else 0.0,
                "evaluation_limited_count": int(sub["VITAL_problem_evaluation"].sum()),
                "popmax_review_count": int(sub["VITAL_problem_popmax"].sum()),
                "model_conflict_count": int(sub["VITAL_problem_model"].sum()),
                "recessive_route_count": int(sub["VITAL_problem_recessive"].sum()),
            }
        )
    return pd.DataFrame(rows).sort_values(["dataset", "domain"]).reset_index(drop=True)


def build_gene_problem_matrix(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(
            columns=[
                "gene",
                "datasets_present",
                "variant_count",
                "evaluation_limited",
                "popmax_review",
                "model_conflict",
                "recessive_route",
                "sv_required",
            ]
        )

    grouped = frame.groupby("gene", dropna=False)
    matrix = grouped.agg(
        datasets_present=("dataset", lambda values: "|".join(sorted(set(map(str, values))))),
        variant_count=("gene", "size"),
        evaluation_limited=("VITAL_problem_evaluation", "any"),
        popmax_review=("VITAL_problem_popmax", "any"),
        model_conflict=("VITAL_problem_model", "any"),
        recessive_route=("VITAL_problem_recessive", "any"),
        sv_required=("VITAL_problem_sv", "any"),
    )
    return matrix.reset_index().sort_values(["variant_count", "gene"], ascending=[False, True]).reset_index(drop=True)


def build_problem_intersections(gene_matrix: pd.DataFrame) -> pd.DataFrame:
    if gene_matrix.empty:
        return pd.DataFrame(
            columns=[
                "evaluation_limited",
                "popmax_review",
                "model_conflict",
                "recessive_route",
                "sv_required",
                "gene_count",
                "genes",
            ]
        )

    keys = ["evaluation_limited", "popmax_review", "model_conflict", "recessive_route", "sv_required"]
    grouped = (
        gene_matrix.groupby(keys, dropna=False)
        .agg(
            gene_count=("gene", "size"),
            genes=("gene", lambda values: "|".join(sorted(map(str, values)))),
        )
        .reset_index()
        .sort_values("gene_count", ascending=False)
        .reset_index(drop=True)
    )
    return grouped
