from __future__ import annotations

import argparse
import math
import time
from pathlib import Path

try:
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires pandas, numpy, and matplotlib. Install them with: "
        "python -m pip install -r requirements.txt"
    ) from exc

from validate_vital_reclassification import clinical_group


BASE_DIR = Path(__file__).resolve().parents[1]
RAW_DIR = BASE_DIR / "data" / "raw"
DATA_DIR = BASE_DIR / "data" / "processed"
FIGURE_DIR = BASE_DIR / "figures"
SUPPLEMENT_DIR = BASE_DIR / "supplementary_tables"

DEFAULT_HISTORICAL_PREDICTIONS = (
    DATA_DIR / "vital_cross_disease_3000_2023_01_to_current_vital_historical_predictions.csv"
)
DEFAULT_OUTPUT_PREFIX = "vital_cross_disease_3000"

SNAPSHOTS = [
    ("2023-01", RAW_DIR / "variant_summary_2023-01.txt.gz"),
    ("2024-01", RAW_DIR / "variant_summary_2024-01.txt.gz"),
    ("2025-01", RAW_DIR / "variant_summary_2025-01.txt.gz"),
    ("2026-04", RAW_DIR / "variant_summary_2026-04.txt.gz"),
]

STRICT_DOWNGRADE_GROUPS = {"B_LB", "VUS"}
BROAD_INSTABILITY_GROUPS = {"B_LB", "VUS", "conflicting", "other", "missing"}
GROUP_ORDER = ["P_LP", "VUS", "B_LB", "conflicting", "other", "missing"]
GROUP_PRIORITY = {group: index for index, group in enumerate(GROUP_ORDER)}


def save_csv(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = output_path.with_name(f".{output_path.name}.tmp")
    last_error: PermissionError | None = None
    for attempt in range(1, 7):
        try:
            df.to_csv(temp_path, index=False)
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


def save_tsv(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = output_path.with_name(f".{output_path.name}.tmp")
    last_error: PermissionError | None = None
    for attempt in range(1, 7):
        try:
            df.to_csv(temp_path, sep="\t", index=False)
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


def normalize_id(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.endswith(".0"):
        text = text[:-2]
    return text


def to_vcv_id(variation_id: object) -> str:
    text = normalize_id(variation_id)
    if not text:
        return ""
    if text.upper().startswith("VCV"):
        return text.upper()
    try:
        return f"VCV{int(text):09d}"
    except ValueError:
        return f"VCV{text}"


def parse_bool_series(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def safe_divide(numerator: int | float, denominator: int | float) -> float:
    return float(numerator / denominator) if denominator else np.nan


def wilson_ci(successes: int | float, total: int | float, z: float = 1.96) -> tuple[float, float]:
    successes = int(successes)
    total = int(total)
    if total <= 0:
        return np.nan, np.nan
    proportion = successes / total
    denominator = 1 + (z**2 / total)
    center = (proportion + (z**2 / (2 * total))) / denominator
    margin = (
        z
        * math.sqrt((proportion * (1 - proportion) / total) + (z**2 / (4 * total**2)))
        / denominator
    )
    low = max(0.0, center - margin)
    high = min(1.0, center + margin)
    if successes == 0:
        low = 0.0
    if successes == total:
        high = 1.0
    return low, high


def collapse_clinical_groups(groups: pd.Series) -> str:
    unique = {group for group in groups.dropna().astype(str) if group and group != "nan"}
    if not unique:
        return "missing"
    if "conflicting" in unique:
        return "conflicting"
    if "P_LP" in unique and (unique & {"B_LB", "VUS", "other"}):
        return "conflicting"
    return sorted(unique, key=lambda group: GROUP_PRIORITY.get(group, 99))[0]


def unique_join(values: pd.Series) -> str:
    unique = sorted({str(value) for value in values.dropna() if str(value) and str(value) != "nan"})
    return "|".join(unique)


def max_numeric(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce")
    return float(numeric.max()) if numeric.notna().any() else np.nan


def read_snapshot_for_ids(path: Path, target_ids: set[str]) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing ClinVar snapshot: {path}")

    usecols = [
        "#AlleleID",
        "Type",
        "Name",
        "GeneSymbol",
        "ClinicalSignificance",
        "ReviewStatus",
        "NumberSubmitters",
        "Assembly",
        "Chromosome",
        "PositionVCF",
        "ReferenceAlleleVCF",
        "AlternateAlleleVCF",
        "VariationID",
    ]
    chunks: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        usecols=lambda column: column in usecols,
        dtype=str,
        chunksize=250_000,
        low_memory=False,
    ):
        if "VariationID" not in chunk.columns:
            continue
        chunk["variation_id"] = chunk["VariationID"].map(normalize_id)
        chunk = chunk[chunk["variation_id"].isin(target_ids)].copy()
        if chunk.empty:
            continue
        chunk["clinical_group"] = chunk["ClinicalSignificance"].map(clinical_group)
        for column in ["Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF"]:
            if column not in chunk.columns:
                chunk[column] = ""
        chunk["variant_key"] = (
            chunk["Chromosome"].fillna("").astype(str)
            + ":"
            + chunk["PositionVCF"].fillna("").astype(str)
            + ":"
            + chunk["ReferenceAlleleVCF"].fillna("").astype(str)
            + ":"
            + chunk["AlternateAlleleVCF"].fillna("").astype(str)
        )
        chunks.append(chunk)

    if not chunks:
        return pd.DataFrame(
            columns=[
                "variation_id",
                "snapshot_gene",
                "snapshot_name",
                "snapshot_variant_type",
                "clinical_group",
                "clinical_significance",
                "review_status",
                "number_submitters",
                "assemblies",
                "grch38_present",
                "variant_key",
                "snapshot_row_count",
            ]
        )

    df = pd.concat(chunks, ignore_index=True)
    df["grch38_flag"] = df.get("Assembly", "").fillna("").astype(str).eq("GRCh38")
    grouped = (
        df.groupby("variation_id", dropna=False)
        .agg(
            snapshot_gene=("GeneSymbol", unique_join),
            snapshot_name=("Name", "first"),
            snapshot_variant_type=("Type", "first"),
            clinical_group=("clinical_group", collapse_clinical_groups),
            clinical_significance=("ClinicalSignificance", unique_join),
            review_status=("ReviewStatus", unique_join),
            number_submitters=("NumberSubmitters", max_numeric),
            assemblies=("Assembly", unique_join),
            grch38_present=("grch38_flag", "max"),
            variant_key=("variant_key", "first"),
            snapshot_row_count=("variation_id", "size"),
        )
        .reset_index()
    )
    return grouped


def priority_group(row: pd.Series) -> str:
    band = str(row.get("vital_band", "")).lower()
    red = bool(row.get("vital_red_flag", False))
    if red or "red" in band:
        return "VITAL_red"
    if "orange" in band:
        return "orange_high_tension"
    if "yellow" in band:
        return "yellow_watch"
    if "gray" in band or "not_found" in band or "no_frequency" in band:
        return "gray_no_frequency"
    return "green_low_tension"


def normalize_review_status(value: object) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).lower().replace("|", " ").split())


def first_snapshot(row: pd.Series, snapshots: list[str], predicate) -> str:
    for snapshot in snapshots[1:]:
        value = row.get(f"clinical_group_{snapshot}", "missing")
        if predicate(value):
            return snapshot
    return ""


def trajectory_string(row: pd.Series, snapshots: list[str]) -> str:
    return "->".join(str(row.get(f"clinical_group_{snapshot}", "missing")) for snapshot in snapshots)


def has_stepwise_label_drift(row: pd.Series, snapshots: list[str]) -> bool:
    groups = [str(row.get(f"clinical_group_{snapshot}", "missing")) for snapshot in snapshots]
    collapsed: list[str] = []
    for group in groups:
        if not collapsed or collapsed[-1] != group:
            collapsed.append(group)
    if len(collapsed) >= 3 and collapsed[0] == "P_LP":
        return True
    if "P_LP" in collapsed and "VUS" in collapsed and "B_LB" in collapsed:
        return collapsed.index("P_LP") < collapsed.index("VUS") < collapsed.index("B_LB")
    return False


def add_rate_ci(row: dict[str, object], prefix: str, successes: int, total: int) -> None:
    low, high = wilson_ci(successes, total)
    row[f"{prefix}_count"] = successes
    row[f"{prefix}_rate"] = safe_divide(successes, total)
    row[f"{prefix}_ci_low"] = low
    row[f"{prefix}_ci_high"] = high


def build_long_table(baseline: pd.DataFrame, snapshots: list[tuple[str, Path]]) -> pd.DataFrame:
    target_ids = set(baseline["variation_id"].dropna().map(normalize_id))
    baseline = baseline.copy()
    if "gene_baseline" not in baseline.columns and "gene" in baseline.columns:
        baseline["gene_baseline"] = baseline["gene"]
    if "name" not in baseline.columns and "title" in baseline.columns:
        baseline["name"] = baseline["title"]
    meta_columns = [
        "variation_id",
        "clinvar_id",
        "gene_baseline",
        "name",
        "clinical_group",
        "followup_clinical_group",
        "vital_score",
        "vital_band",
        "vital_red_flag",
        "max_frequency_signal",
        "qualifying_frequency_ac",
        "vital_signal_reason",
    ]
    meta_columns = [column for column in meta_columns if column in baseline.columns]
    meta = baseline[meta_columns].copy()
    meta["variation_id"] = meta["variation_id"].map(normalize_id)
    if "clinvar_id" not in meta.columns:
        meta["clinvar_id"] = meta["variation_id"].map(to_vcv_id)
    meta["vcv_id"] = meta["variation_id"].map(to_vcv_id)
    if "vital_red_flag" in meta.columns:
        meta["vital_red_flag"] = parse_bool_series(meta["vital_red_flag"])
    else:
        meta["vital_red_flag"] = False
    meta["priority_group"] = meta.apply(priority_group, axis=1)

    long_parts: list[pd.DataFrame] = []
    for order, (snapshot, path) in enumerate(snapshots):
        print(f"Reading {snapshot}: {path}")
        snapshot_df = read_snapshot_for_ids(path, target_ids)
        merged = meta.merge(snapshot_df, on="variation_id", how="left")
        snapshot_group_column = "clinical_group_y" if "clinical_group_y" in merged.columns else "clinical_group"
        merged["snapshot"] = snapshot
        merged["snapshot_order"] = order
        merged["status_present"] = merged[snapshot_group_column].notna()
        merged["snapshot_clinical_group"] = merged[snapshot_group_column].fillna("missing")
        merged["snapshot_clinical_significance"] = merged["clinical_significance"].fillna("")
        merged["snapshot_review_status"] = merged["review_status"].fillna("")
        merged["snapshot_number_submitters"] = merged["number_submitters"]
        merged["snapshot_gene"] = merged["snapshot_gene"].fillna("")
        merged["snapshot_name"] = merged["snapshot_name"].fillna("")
        merged["snapshot_variant_type"] = merged["snapshot_variant_type"].fillna("")
        merged["assemblies"] = merged["assemblies"].fillna("")
        merged["grch38_present"] = merged["grch38_present"].fillna(False).astype(bool)
        merged["variant_key_snapshot"] = merged["variant_key"].fillna("")
        merged["snapshot_row_count"] = merged["snapshot_row_count"].fillna(0).astype(int)
        merged["missing_reason"] = np.where(
            merged["status_present"], "present", "not_found_in_variant_summary"
        )
        long_parts.append(
            merged[
                [
                    "variation_id",
                    "vcv_id",
                    "clinvar_id",
                    "gene_baseline",
                    "name",
                    "vital_score",
                    "vital_band",
                    "vital_red_flag",
                    "priority_group",
                    "max_frequency_signal",
                    "qualifying_frequency_ac",
                    "vital_signal_reason",
                    "snapshot",
                    "snapshot_order",
                    "status_present",
                    "missing_reason",
                    "snapshot_gene",
                    "snapshot_name",
                    "snapshot_variant_type",
                    "snapshot_clinical_group",
                    "snapshot_clinical_significance",
                    "snapshot_review_status",
                    "snapshot_number_submitters",
                    "assemblies",
                    "grch38_present",
                    "variant_key_snapshot",
                    "snapshot_row_count",
                ]
            ]
        )
    return pd.concat(long_parts, ignore_index=True)


def build_wide_table(long_df: pd.DataFrame, snapshots: list[str]) -> pd.DataFrame:
    meta_columns = [
        "variation_id",
        "vcv_id",
        "clinvar_id",
        "gene_baseline",
        "name",
        "vital_score",
        "vital_band",
        "vital_red_flag",
        "priority_group",
        "max_frequency_signal",
        "qualifying_frequency_ac",
        "vital_signal_reason",
    ]
    meta = long_df[meta_columns].drop_duplicates("variation_id").copy()
    wide = meta.copy()

    for value_column, prefix in [
        ("snapshot_clinical_group", "clinical_group"),
        ("snapshot_clinical_significance", "clinical_significance"),
        ("snapshot_review_status", "review_status"),
        ("status_present", "status_present"),
    ]:
        pivot = (
            long_df.pivot(index="variation_id", columns="snapshot", values=value_column)
            .rename(columns={snapshot: f"{prefix}_{snapshot}" for snapshot in snapshots})
            .reset_index()
        )
        wide = wide.merge(pivot, on="variation_id", how="left")

    for snapshot in snapshots:
        wide[f"clinical_group_{snapshot}"] = wide[f"clinical_group_{snapshot}"].fillna("missing")
        wide[f"review_status_{snapshot}"] = wide[f"review_status_{snapshot}"].fillna("")
        wide[f"status_present_{snapshot}"] = wide[f"status_present_{snapshot}"].fillna(False).astype(bool)

    post_snapshots = snapshots[1:]
    wide["clinical_group_trajectory"] = wide.apply(lambda row: trajectory_string(row, snapshots), axis=1)
    wide["ever_strict_downgrade"] = wide.apply(
        lambda row: any(str(row[f"clinical_group_{snapshot}"]) in STRICT_DOWNGRADE_GROUPS for snapshot in post_snapshots),
        axis=1,
    )
    wide["ever_broad_instability"] = wide.apply(
        lambda row: any(str(row[f"clinical_group_{snapshot}"]) in BROAD_INSTABILITY_GROUPS for snapshot in post_snapshots),
        axis=1,
    )
    wide["first_strict_downgrade_snapshot"] = wide.apply(
        lambda row: first_snapshot(row, snapshots, lambda value: str(value) in STRICT_DOWNGRADE_GROUPS),
        axis=1,
    )
    wide["first_broad_instability_snapshot"] = wide.apply(
        lambda row: first_snapshot(row, snapshots, lambda value: str(value) in BROAD_INSTABILITY_GROUPS),
        axis=1,
    )

    baseline_review_column = f"review_status_{snapshots[0]}"
    for snapshot in snapshots:
        wide[f"normalized_review_status_{snapshot}"] = wide[f"review_status_{snapshot}"].map(
            normalize_review_status
        )
    wide["ever_review_status_changed"] = wide.apply(
        lambda row: any(
            row[f"normalized_review_status_{snapshot}"]
            and row[f"normalized_review_status_{snapshot}"] != row[baseline_review_column.replace("review_status", "normalized_review_status")]
            for snapshot in post_snapshots
        ),
        axis=1,
    )
    wide["first_review_status_change_snapshot"] = wide.apply(
        lambda row: next(
            (
                snapshot
                for snapshot in post_snapshots
                if row[f"normalized_review_status_{snapshot}"]
                and row[f"normalized_review_status_{snapshot}"]
                != row[baseline_review_column.replace("review_status", "normalized_review_status")]
            ),
            "",
        ),
        axis=1,
    )
    wide["ever_expanded_instability"] = (
        wide["ever_broad_instability"].astype(bool) | wide["ever_review_status_changed"].astype(bool)
    )
    wide["stepwise_label_drift"] = wide.apply(lambda row: has_stepwise_label_drift(row, snapshots), axis=1)
    intermediate_snapshots = snapshots[1:-1]
    wide["flipper_plp_vus_back_to_plp"] = wide.apply(
        lambda row: str(row[f"clinical_group_{snapshots[0]}"]) == "P_LP"
        and any(str(row[f"clinical_group_{snapshot}"]) == "VUS" for snapshot in intermediate_snapshots)
        and str(row[f"clinical_group_{snapshots[-1]}"]) == "P_LP",
        axis=1,
    )
    wide["flipper_plp_nonplp_back_to_plp"] = wide.apply(
        lambda row: str(row[f"clinical_group_{snapshots[0]}"]) == "P_LP"
        and any(
            str(row[f"clinical_group_{snapshot}"]) in BROAD_INSTABILITY_GROUPS
            for snapshot in intermediate_snapshots
        )
        and str(row[f"clinical_group_{snapshots[-1]}"]) == "P_LP",
        axis=1,
    )
    wide["intermediate_only_strict_downgrade"] = wide.apply(
        lambda row: str(row[f"clinical_group_{snapshots[0]}"]) == "P_LP"
        and any(
            str(row[f"clinical_group_{snapshot}"]) in STRICT_DOWNGRADE_GROUPS
            for snapshot in intermediate_snapshots
        )
        and str(row[f"clinical_group_{snapshots[-1]}"]) == "P_LP",
        axis=1,
    )
    wide["intermediate_only_broad_instability"] = wide.apply(
        lambda row: str(row[f"clinical_group_{snapshots[0]}"]) == "P_LP"
        and any(
            str(row[f"clinical_group_{snapshot}"]) in BROAD_INSTABILITY_GROUPS
            for snapshot in intermediate_snapshots
        )
        and str(row[f"clinical_group_{snapshots[-1]}"]) == "P_LP",
        axis=1,
    )
    wide["transient_missing"] = wide.apply(
        lambda row: any(not bool(row[f"status_present_{snapshot}"]) for snapshot in snapshots[1:-1])
        and bool(row[f"status_present_{snapshots[-1]}"]),
        axis=1,
    )

    def current_missing_category(row: pd.Series) -> str:
        if not bool(row[f"status_present_{snapshots[0]}"]):
            if bool(row[f"status_present_{snapshots[-1]}"]):
                return "not_present_at_2023_baseline"
            if any(bool(row[f"status_present_{snapshot}"]) for snapshot in snapshots[1:-1]):
                return "intermediate_only_presence"
            return "mapping_failure_or_id_mismatch"
        if bool(row[f"status_present_{snapshots[-1]}"]):
            return "current_present"
        if any(bool(row[f"status_present_{snapshot}"]) for snapshot in snapshots[1:-1]):
            return "disappeared_after_prior_presence"
        return "disappeared_after_2023"

    wide["current_missing_category"] = wide.apply(current_missing_category, axis=1)
    wide["event_timing_strict"] = wide["first_strict_downgrade_snapshot"].replace("", "none")
    wide["event_timing_broad"] = wide["first_broad_instability_snapshot"].replace("", "none")
    return wide


def make_summary_tables(wide: pd.DataFrame, long_df: pd.DataFrame, snapshots: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    cohorts = [("all", pd.Series(True, index=wide.index))]
    for group in ["VITAL_red", "orange_high_tension", "yellow_watch", "gray_no_frequency", "green_low_tension"]:
        cohorts.append((group, wide["priority_group"].eq(group)))

    for cohort_name, mask in cohorts:
        sub = wide[mask].copy()
        n = len(sub)
        if n == 0:
            continue
        row: dict[str, object] = {
            "priority_group": cohort_name,
            "n_variants": n,
            "median_vital_score": float(pd.to_numeric(sub["vital_score"], errors="coerce").median()),
        }
        for prefix, column in [
            ("ever_strict_downgrade", "ever_strict_downgrade"),
            ("ever_broad_instability", "ever_broad_instability"),
            ("ever_expanded_instability", "ever_expanded_instability"),
            ("ever_review_status_changed", "ever_review_status_changed"),
            ("stepwise_label_drift", "stepwise_label_drift"),
            ("flipper_plp_vus_back_to_plp", "flipper_plp_vus_back_to_plp"),
            ("flipper_plp_nonplp_back_to_plp", "flipper_plp_nonplp_back_to_plp"),
            ("intermediate_only_strict_downgrade", "intermediate_only_strict_downgrade"),
            ("intermediate_only_broad_instability", "intermediate_only_broad_instability"),
            ("transient_missing", "transient_missing"),
        ]:
            add_rate_ci(row, prefix, int(sub[column].sum()), n)
        current_missing = int((~sub[f"status_present_{snapshots[-1]}"].astype(bool)).sum())
        add_rate_ci(row, "current_missing", current_missing, n)
        rows.append(row)

    snapshot_rows: list[dict[str, object]] = []
    for (cohort_name, mask) in cohorts:
        variant_ids = set(wide.loc[mask, "variation_id"])
        if not variant_ids:
            continue
        sub_long = long_df[long_df["variation_id"].isin(variant_ids)]
        for snapshot in snapshots:
            sub = sub_long[sub_long["snapshot"].eq(snapshot)]
            n = len(sub)
            row = {"priority_group": cohort_name, "snapshot": snapshot, "n_variants": n}
            for prefix, condition in [
                ("p_lp", sub["snapshot_clinical_group"].eq("P_LP")),
                ("strict_downgrade_status", sub["snapshot_clinical_group"].isin(STRICT_DOWNGRADE_GROUPS)),
                ("broad_non_plp_or_missing", sub["snapshot_clinical_group"].isin(BROAD_INSTABILITY_GROUPS)),
                ("missing", ~sub["status_present"].astype(bool)),
            ]:
                add_rate_ci(row, prefix, int(condition.sum()), n)
            snapshot_rows.append(row)

    missing_breakdown = (
        wide.groupby(["priority_group", "current_missing_category"], dropna=False)
        .size()
        .reset_index(name="n_variants")
        .sort_values(["priority_group", "current_missing_category"])
    )
    return pd.DataFrame(rows), pd.DataFrame(snapshot_rows), missing_breakdown


def make_retention_plot(snapshot_summary: pd.DataFrame, output_path: Path) -> None:
    plot_groups = ["VITAL_red", "orange_high_tension", "yellow_watch", "gray_no_frequency", "green_low_tension"]
    labels = {
        "VITAL_red": "Urgent review",
        "orange_high_tension": "Orange",
        "yellow_watch": "Yellow",
        "gray_no_frequency": "Gray",
        "green_low_tension": "Green",
    }
    colors = {
        "VITAL_red": "#b22a2a",
        "orange_high_tension": "#d87325",
        "yellow_watch": "#c6a700",
        "gray_no_frequency": "#6f7782",
        "green_low_tension": "#2d7a46",
    }
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    for group in plot_groups:
        sub = snapshot_summary[snapshot_summary["priority_group"].eq(group)].copy()
        if sub.empty:
            continue
        x = np.arange(len(sub))
        y = sub["p_lp_rate"].astype(float).to_numpy()
        low = sub["p_lp_ci_low"].astype(float).to_numpy()
        high = sub["p_lp_ci_high"].astype(float).to_numpy()
        ax.plot(x, y, marker="o", linewidth=2, label=labels[group], color=colors[group])
        ax.fill_between(x, low, high, color=colors[group], alpha=0.12, linewidth=0)
    ax.set_xticks(np.arange(snapshot_summary["snapshot"].nunique()))
    ax.set_xticklabels(sorted(snapshot_summary["snapshot"].unique()))
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("Fraction retained as P/LP")
    ax.set_xlabel("ClinVar snapshot")
    ax.set_title("ClinVar P/LP status retention across archived snapshots")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, ncol=2, fontsize=9)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    print(f"Saved {output_path}")


def run(historical_predictions: Path, output_prefix: str, snapshots: list[tuple[str, Path]]) -> None:
    baseline = pd.read_csv(historical_predictions, dtype={"variation_id": str})
    baseline["variation_id"] = baseline["variation_id"].map(normalize_id)
    if "vital_red_flag" in baseline.columns:
        baseline["vital_red_flag"] = parse_bool_series(baseline["vital_red_flag"])
    print(f"Loaded {historical_predictions} ({len(baseline)} variants)")

    duplicate_count = int(baseline["variation_id"].duplicated().sum())
    if duplicate_count:
        print(f"WARNING: {duplicate_count} duplicate variation_id rows in historical predictions; keeping first.")
        baseline = baseline.drop_duplicates("variation_id", keep="first")

    snapshot_labels = [label for label, _ in snapshots]
    long_df = build_long_table(baseline, snapshots)
    wide = build_wide_table(long_df, snapshot_labels)
    summary, snapshot_summary, missing_breakdown = make_summary_tables(wide, long_df, snapshot_labels)

    red_queue = wide[wide["priority_group"].eq("VITAL_red")].copy()
    red_queue = red_queue.sort_values(["vital_score", "variation_id"], ascending=[False, True])
    flippers = wide[
        wide["flipper_plp_vus_back_to_plp"].astype(bool)
        | wide["flipper_plp_nonplp_back_to_plp"].astype(bool)
        | wide["intermediate_only_strict_downgrade"].astype(bool)
        | wide["intermediate_only_broad_instability"].astype(bool)
    ].copy()
    flippers = flippers.sort_values(["priority_group", "vital_score", "variation_id"], ascending=[True, False, True])

    save_csv(long_df, DATA_DIR / f"{output_prefix}_clinvar_time_series_long.csv")
    save_csv(wide, DATA_DIR / f"{output_prefix}_clinvar_time_series_wide.csv")
    save_csv(summary, DATA_DIR / f"{output_prefix}_clinvar_time_series_summary.csv")
    save_csv(snapshot_summary, DATA_DIR / f"{output_prefix}_clinvar_time_series_snapshot_summary.csv")
    save_csv(missing_breakdown, DATA_DIR / f"{output_prefix}_clinvar_time_series_missing_breakdown.csv")
    save_csv(red_queue, DATA_DIR / f"{output_prefix}_clinvar_time_series_red_queue.csv")
    save_csv(flippers, DATA_DIR / f"{output_prefix}_clinvar_time_series_flippers.csv")

    save_tsv(
        wide,
        SUPPLEMENT_DIR / f"Supplementary_Table_S18_{output_prefix}_ClinVar_time_series_trajectories.tsv",
    )
    save_tsv(
        red_queue,
        SUPPLEMENT_DIR / f"Supplementary_Table_S19_{output_prefix}_ClinVar_time_series_red_queue.tsv",
    )
    save_tsv(
        flippers,
        SUPPLEMENT_DIR / f"Supplementary_Table_S20_{output_prefix}_ClinVar_time_series_flippers.tsv",
    )
    make_retention_plot(
        snapshot_summary,
        FIGURE_DIR / f"{output_prefix}_clinvar_time_series_retention.png",
    )

    print("\nTrajectory summary:")
    print(summary.to_string(index=False))
    print("\nMissing/current-presence breakdown:")
    print(missing_breakdown.to_string(index=False))
    if not red_queue.empty:
        cols = [
            "variation_id",
            "vcv_id",
            "gene_baseline",
            "vital_score",
            "clinical_group_trajectory",
            "ever_strict_downgrade",
            "ever_expanded_instability",
            "first_strict_downgrade_snapshot",
            "current_missing_category",
        ]
        print("\nUrgent-review trajectories:")
        print(red_queue[cols].to_string(index=False))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build VCV-level ClinVar classification trajectories across archived variant_summary snapshots."
    )
    parser.add_argument(
        "--historical-predictions",
        type=Path,
        default=DEFAULT_HISTORICAL_PREDICTIONS,
        help="Historical predictions table to annotate with archived ClinVar trajectories.",
    )
    parser.add_argument(
        "--output-prefix",
        default=DEFAULT_OUTPUT_PREFIX,
        help="Prefix for generated processed-data and figure files.",
    )
    parser.add_argument(
        "--snapshot",
        action="append",
        default=None,
        help="Optional snapshot override as LABEL=PATH. Repeat for multiple snapshots.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    snapshots = SNAPSHOTS
    if args.snapshot:
        snapshots = []
        for item in args.snapshot:
            if "=" not in item:
                raise SystemExit(f"Invalid --snapshot value {item!r}; expected LABEL=PATH")
            label, path_text = item.split("=", 1)
            snapshots.append((label, Path(path_text)))
    run(args.historical_predictions, args.output_prefix, snapshots)


if __name__ == "__main__":
    main()
