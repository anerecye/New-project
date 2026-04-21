from __future__ import annotations

import argparse
import math
import time
from pathlib import Path

try:
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from scipy.stats import fisher_exact
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires pandas, numpy, scipy, and matplotlib. Install them with: "
        "python -m pip install -r requirements.txt"
    ) from exc


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
FIGURE_DIR = BASE_DIR / "figures"
AF_ULTRA_RARE = 1e-5
AF_RARE = 1e-4
MIN_RECLASSIFICATION_AC = 20
VITAL_ACTION_THRESHOLD = 70

ARRHYTHMIA_GENES = [
    "KCNQ1",
    "KCNH2",
    "SCN5A",
    "KCNE1",
    "KCNE2",
    "RYR2",
    "CASQ2",
    "TRDN",
    "CALM1",
    "CALM2",
    "CALM3",
    "ANK2",
    "SCN4B",
    "KCNJ2",
    "HCN4",
    "CACNA1C",
    "CACNB2",
    "CACNA2D1",
    "AKAP9",
    "SNTA1",
]


def prefixed_name(prefix: str, name: str) -> str:
    return f"{prefix}_{name}" if prefix else name


def data_path(prefix: str, name: str) -> Path:
    return DATA_DIR / prefixed_name(prefix, name)


def figure_path(prefix: str, name: str) -> Path:
    return FIGURE_DIR / prefixed_name(prefix, name)


def save_table(df: pd.DataFrame, output_path: Path) -> None:
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


def clinical_group(value: object) -> str:
    text = str(value).lower()
    if not text or text == "nan" or text == "-":
        return "missing"
    if "conflicting" in text:
        return "conflicting"
    has_pathogenic = "pathogenic" in text
    has_benign = "benign" in text
    if has_pathogenic and not has_benign:
        return "P_LP"
    if has_benign and not has_pathogenic:
        return "B_LB"
    if "uncertain" in text or "vus" in text:
        return "VUS"
    return "other"


def collapse_groups(groups: pd.Series) -> str:
    unique = set(groups.dropna().astype(str))
    for group in ["P_LP", "B_LB", "VUS", "conflicting", "other", "missing"]:
        if group in unique:
            return group
    return "missing"


def read_variant_summary(path: Path, genes: list[str] | None) -> pd.DataFrame:
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
        if genes is not None:
            chunk = chunk[chunk["GeneSymbol"].isin(genes)].copy()
        else:
            chunk = chunk.copy()
        if "Assembly" in chunk.columns:
            chunk = chunk[chunk["Assembly"].eq("GRCh38")]
        if chunk.empty:
            continue
        chunk["clinical_group"] = chunk["ClinicalSignificance"].map(clinical_group)
        chunks.append(chunk)
    if not chunks:
        return pd.DataFrame(columns=["variation_id", "clinical_group"])

    df = pd.concat(chunks, ignore_index=True)
    df["variation_id"] = df["VariationID"].astype(str)
    for column in ["Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF"]:
        if column not in df.columns:
            df[column] = ""
    df["variant_key"] = (
        df["Chromosome"].fillna("").astype(str)
        + ":"
        + df["PositionVCF"].fillna("").astype(str)
        + ":"
        + df["ReferenceAlleleVCF"].fillna("").astype(str)
        + ":"
        + df["AlternateAlleleVCF"].fillna("").astype(str)
    )
    grouped = (
        df.groupby("variation_id", dropna=False)
        .agg(
            gene=("GeneSymbol", "first"),
            variant_key=("variant_key", "first"),
            name=("Name", "first"),
            variant_type=("Type", "first"),
            clinical_group=("clinical_group", collapse_groups),
            clinical_significance=("ClinicalSignificance", lambda values: "|".join(sorted(set(values.dropna().astype(str))))),
            review_status=("ReviewStatus", lambda values: "|".join(sorted(set(values.dropna().astype(str))))),
            number_submitters=("NumberSubmitters", "max"),
        )
        .reset_index()
    )
    return grouped


def normalize_change_text(series: pd.Series) -> pd.Series:
    return (
        series.fillna("")
        .astype(str)
        .str.lower()
        .str.replace(r"\s+", " ", regex=True)
        .str.strip()
    )


def historical_targets(evaluated: pd.DataFrame) -> dict[str, pd.Series]:
    targets = {
        "strict_PLP_to_BLB_or_VUS": evaluated["strict_revised_to_b_or_vus"].astype(bool),
        "broad_PLP_to_non_PLP_or_missing": evaluated["broad_revised_or_destabilized"].astype(bool),
    }
    if "expanded_revised_or_review_changed" in evaluated.columns:
        targets["expanded_PLP_to_non_PLP_or_review_change"] = evaluated[
            "expanded_revised_or_review_changed"
        ].astype(bool)
    return targets


def roc_auc(labels: pd.Series, scores: pd.Series) -> float:
    labels = labels.astype(bool)
    n_pos = int(labels.sum())
    n_neg = int((~labels).sum())
    if n_pos == 0 or n_neg == 0:
        return np.nan
    ranks = scores.rank(method="average", ascending=True)
    rank_sum_pos = ranks[labels].sum()
    return float((rank_sum_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def average_precision(labels: pd.Series, scores: pd.Series) -> float:
    ordered = pd.DataFrame({"label": labels.astype(bool), "score": scores}).sort_values(
        "score", ascending=False
    )
    n_pos = int(ordered["label"].sum())
    if n_pos == 0:
        return np.nan
    ordered["cum_pos"] = ordered["label"].cumsum()
    ordered["rank"] = np.arange(1, len(ordered) + 1)
    return float((ordered.loc[ordered["label"], "cum_pos"] / ordered.loc[ordered["label"], "rank"]).sum() / n_pos)


def precision_at_fraction(labels: pd.Series, scores: pd.Series, fraction: float) -> float:
    n = max(1, int(np.ceil(len(labels) * fraction)))
    top = pd.DataFrame({"label": labels.astype(bool), "score": scores}).sort_values(
        "score", ascending=False
    ).head(n)
    return float(top["label"].mean()) if len(top) else np.nan


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
    return max(0.0, center - margin), min(1.0, center + margin)


def log_risk_ratio_ci(
    events_a: int,
    total_a: int,
    events_b: int,
    total_b: int,
    z: float = 1.96,
) -> tuple[float, float]:
    if total_a <= 0 or total_b <= 0:
        return np.nan, np.nan
    # Haldane-Anscombe correction keeps intervals finite for sparse validation cells.
    a = events_a + 0.5
    b = total_a - events_a + 0.5
    c = events_b + 0.5
    d = total_b - events_b + 0.5
    risk_a = a / (a + b)
    risk_b = c / (c + d)
    if risk_b <= 0:
        return np.nan, np.nan
    log_rr = math.log(risk_a / risk_b)
    se = math.sqrt((1 / a) - (1 / (a + b)) + (1 / c) - (1 / (c + d)))
    return math.exp(log_rr - z * se), math.exp(log_rr + z * se)


def add_binary_metric_cis(row: dict[str, object], tp: int, fp: int, fn: int, tn: int) -> dict[str, object]:
    for metric, successes, total in [
        ("precision", tp, tp + fp),
        ("recall", tp, tp + fn),
        ("specificity", tn, tn + fp),
        ("false_positive_rate", fp, fp + tn),
        ("false_negative_rate", fn, fn + tp),
    ]:
        low, high = wilson_ci(successes, total)
        row[f"{metric}_ci_low"] = low
        row[f"{metric}_ci_high"] = high
    return row


def make_curve_points(labels: pd.Series, scores: pd.Series, target: str) -> pd.DataFrame:
    curve = pd.DataFrame(
        {
            "label": labels.astype(bool).to_numpy(),
            "score": pd.to_numeric(scores, errors="coerce").fillna(0).to_numpy(),
        }
    ).sort_values("score", ascending=False)
    n_positive = int(curve["label"].sum())
    n_negative = int((~curve["label"]).sum())
    rows: list[dict[str, object]] = []
    thresholds = [np.inf, *curve["score"].drop_duplicates().tolist(), -np.inf]
    for threshold in thresholds:
        predicted = curve["score"] >= threshold
        tp = int((predicted & curve["label"]).sum())
        fp = int((predicted & ~curve["label"]).sum())
        fn = int((~predicted & curve["label"]).sum())
        tn = int((~predicted & ~curve["label"]).sum())
        rows.append(
            {
                "target": target,
                "threshold": threshold,
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "tn": tn,
                "tpr_recall": safe_divide(tp, n_positive),
                "fpr": safe_divide(fp, n_negative),
                "precision": safe_divide(tp, tp + fp),
                "specificity": safe_divide(tn, tn + fp),
            }
        )
    return pd.DataFrame(rows)


def method_flags(df: pd.DataFrame) -> dict[str, pd.Series]:
    numeric = df.copy()
    for column in ["global_af", "popmax_af", "max_frequency_signal"]:
        if column not in numeric.columns:
            numeric[column] = np.nan
        numeric[column] = pd.to_numeric(numeric[column], errors="coerce")
    numeric["qualifying_frequency_ac"] = pd.to_numeric(
        numeric.get("qualifying_frequency_ac", np.nan), errors="coerce"
    )
    if "frequency_signal_ac_ge_20" in numeric.columns:
        ac_supported = parse_bool_series(numeric["frequency_signal_ac_ge_20"])
    elif "qualifying_frequency_ac" in numeric.columns:
        ac_supported = pd.to_numeric(numeric["qualifying_frequency_ac"], errors="coerce").fillna(0) >= MIN_RECLASSIFICATION_AC
    else:
        ac_supported = pd.Series(False, index=numeric.index)
    vital_red = (
        parse_bool_series(numeric["vital_red_flag"])
        if "vital_red_flag" in numeric.columns
        else pd.Series(False, index=numeric.index)
    )
    return {
        "global_AF_gt_1e-5": numeric["global_af"] > AF_ULTRA_RARE,
        "global_AF_gt_1e-4": numeric["global_af"] > AF_RARE,
        "popmax_or_global_AF_gt_1e-5": numeric["max_frequency_signal"] > AF_ULTRA_RARE,
        "popmax_or_global_AF_gt_1e-4": numeric["max_frequency_signal"] > AF_RARE,
        "AF_gt_1e-5_and_AC_ge_5": (numeric["max_frequency_signal"] > AF_ULTRA_RARE)
        & (numeric["qualifying_frequency_ac"].fillna(0) >= 5),
        "AF_gt_1e-5_and_AC_ge_10": (numeric["max_frequency_signal"] > AF_ULTRA_RARE)
        & (numeric["qualifying_frequency_ac"].fillna(0) >= 10),
        "AF_gt_1e-5_and_AC_ge_20": (numeric["max_frequency_signal"] > AF_ULTRA_RARE) & ac_supported,
        "AF_gt_1e-5_and_AC_ge_50": (numeric["max_frequency_signal"] > AF_ULTRA_RARE)
        & (numeric["qualifying_frequency_ac"].fillna(0) >= 50),
        "VITAL_red": vital_red,
    }


def make_historical_method_comparison(evaluated: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    targets = historical_targets(evaluated)
    flags = method_flags(evaluated)
    rows: list[dict[str, object]] = []
    curves: list[pd.DataFrame] = []
    for target_name, labels in targets.items():
        positives = labels
        negatives = ~labels
        for method, predicted in flags.items():
            predicted = predicted.fillna(False).astype(bool)
            tp = int((predicted & positives).sum())
            fp = int((predicted & negatives).sum())
            fn = int((~predicted & positives).sum())
            tn = int((~predicted & negatives).sum())
            rows.append(
                add_binary_metric_cis(
                    {
                        "target": target_name,
                        "method": method,
                        "true_positives": tp,
                        "false_positives": fp,
                        "false_negatives": fn,
                        "true_negatives": tn,
                        "precision": safe_divide(tp, tp + fp),
                        "recall": safe_divide(tp, tp + fn),
                        "specificity": safe_divide(tn, tn + fp),
                        "false_positive_rate": safe_divide(fp, fp + tn),
                        "false_negative_rate": safe_divide(fn, fn + tp),
                        "confidence_interval_method": "Wilson 95%",
                    },
                    tp,
                    fp,
                    fn,
                    tn,
                )
            )
        curves.append(make_curve_points(labels, evaluated["vital_score"], target_name))
    return pd.DataFrame(rows), pd.concat(curves, ignore_index=True) if curves else pd.DataFrame()


def vital_score_threshold_sweep(evaluated: pd.DataFrame) -> pd.DataFrame:
    targets = historical_targets(evaluated)
    score = pd.to_numeric(evaluated["vital_score"], errors="coerce").fillna(-np.inf)
    if "frequency_evidence_status" in evaluated.columns:
        observed_scope = evaluated["frequency_evidence_status"].fillna("").astype(str).eq("frequency_observed")
    else:
        observed_scope = pd.Series(True, index=evaluated.index)
    scopes = {
        "all_baseline_plp": pd.Series(True, index=evaluated.index),
        "frequency_observed": observed_scope,
    }
    thresholds = sorted({*range(40, 96, 5), VITAL_ACTION_THRESHOLD})
    rows: list[dict[str, object]] = []
    for scope_name, scope_mask in scopes.items():
        scoped_score = score.loc[scope_mask]
        if scoped_score.empty:
            continue
        for target_name, labels in targets.items():
            scoped_labels = labels.loc[scope_mask]
            positives = scoped_labels.astype(bool)
            for threshold in thresholds:
                predicted = scoped_score >= threshold
                tp = int((predicted & positives).sum())
                fp = int((predicted & ~positives).sum())
                fn = int((~predicted & positives).sum())
                tn = int((~predicted & ~positives).sum())
                rows.append(
                    {
                        "scope": scope_name,
                        "target": target_name,
                        "score_threshold": threshold,
                        "flagged_count": int(predicted.sum()),
                        "true_positives": tp,
                        "false_positives": fp,
                        "false_negatives": fn,
                        "true_negatives": tn,
                        "precision": safe_divide(tp, tp + fp),
                        "recall": safe_divide(tp, tp + fn),
                        "specificity": safe_divide(tn, tn + fp),
                        "false_positive_rate": safe_divide(fp, fp + tn),
                    }
                )
    return pd.DataFrame(rows)


def make_reclassification_enrichment(evaluated: pd.DataFrame) -> pd.DataFrame:
    targets = historical_targets(evaluated)
    red = parse_bool_series(evaluated["vital_red_flag"]) if "vital_red_flag" in evaluated.columns else pd.Series(False, index=evaluated.index)
    if "frequency_evidence_status" in evaluated.columns:
        observed_scope = evaluated["frequency_evidence_status"].fillna("").astype(str).eq("frequency_observed")
    else:
        observed_scope = pd.Series(True, index=evaluated.index)
    scopes = {
        "all_baseline_plp": pd.Series(True, index=evaluated.index),
        "frequency_observed": observed_scope,
    }
    rows: list[dict[str, object]] = []
    for scope_name, scope_mask in scopes.items():
        scoped_red = red.loc[scope_mask]
        for target_name, labels in targets.items():
            scoped_labels = labels.loc[scope_mask]
            total_n = int(scope_mask.sum())
            event_count = int(scoped_labels.sum())
            red_n = int(scoped_red.sum())
            red_events = int((scoped_red & scoped_labels).sum())
            nonred_n = total_n - red_n
            nonred_events = event_count - red_events
            overall_rate = safe_divide(event_count, total_n)
            red_rate = safe_divide(red_events, red_n)
            nonred_rate = safe_divide(nonred_events, nonred_n)
            odds_ratio, fisher_p = (
                fisher_exact([[red_events, red_n - red_events], [nonred_events, nonred_n - nonred_events]])
                if red_n and nonred_n
                else (np.nan, np.nan)
            )
            overall_ci_low, overall_ci_high = wilson_ci(event_count, total_n)
            red_ci_low, red_ci_high = wilson_ci(red_events, red_n)
            nonred_ci_low, nonred_ci_high = wilson_ci(nonred_events, nonred_n)
            rr_nonred_ci_low, rr_nonred_ci_high = log_risk_ratio_ci(
                red_events,
                red_n,
                nonred_events,
                nonred_n,
            )
            rows.append(
                {
                    "scope": scope_name,
                    "target": target_name,
                    "baseline_variant_count": total_n,
                    "event_count": event_count,
                    "overall_event_rate": overall_rate,
                    "overall_event_rate_ci_low": overall_ci_low,
                    "overall_event_rate_ci_high": overall_ci_high,
                    "vital_red_count": red_n,
                    "vital_red_event_count": red_events,
                    "vital_red_event_rate": red_rate,
                    "vital_red_event_rate_ci_low": red_ci_low,
                    "vital_red_event_rate_ci_high": red_ci_high,
                    "nonred_event_rate": nonred_rate,
                    "nonred_event_rate_ci_low": nonred_ci_low,
                    "nonred_event_rate_ci_high": nonred_ci_high,
                    "red_enrichment_vs_overall": safe_divide(red_rate, overall_rate),
                    "red_enrichment_vs_nonred": safe_divide(red_rate, nonred_rate),
                    "red_enrichment_vs_nonred_ci_low": rr_nonred_ci_low,
                    "red_enrichment_vs_nonred_ci_high": rr_nonred_ci_high,
                    "fisher_exact_odds_ratio": odds_ratio,
                    "fisher_exact_p": fisher_p,
                    "confidence_interval_method": "Wilson 95% for rates; Haldane-corrected log risk ratio for enrichment",
                }
            )
    return pd.DataFrame(rows)


def make_historical_ac_threshold_sensitivity(evaluated: pd.DataFrame) -> pd.DataFrame:
    targets = historical_targets(evaluated)
    table = evaluated.copy()
    if "variant_key" not in table.columns:
        table["variant_key"] = ""
        for source_column in ["variant_key_score", "variant_key_baseline", "variation_id"]:
            if source_column in table.columns:
                table["variant_key"] = table["variant_key"].where(
                    table["variant_key"].fillna("").astype(str).ne(""),
                    table[source_column].fillna("").astype(str),
                )
    if "gene" not in table.columns:
        table["gene"] = ""
        for source_column in ["gene_score", "gene_baseline"]:
            if source_column in table.columns:
                table["gene"] = table["gene"].where(
                    table["gene"].fillna("").astype(str).ne(""),
                    table[source_column].fillna("").astype(str),
                )
    for column in ["weak_review_signal"]:
        if column not in table.columns:
            table[column] = False
        table[column] = parse_bool_series(table[column])
    for column in ["max_frequency_signal", "qualifying_frequency_ac", "vital_score"]:
        table[column] = pd.to_numeric(table.get(column), errors="coerce")
    if "frequency_evidence_status" in table.columns:
        observed_scope = table["frequency_evidence_status"].fillna("").astype(str).eq("frequency_observed")
    else:
        observed_scope = table["max_frequency_signal"].notna()
    scopes = {
        "all_baseline_plp": pd.Series(True, index=table.index),
        "frequency_observed": observed_scope,
    }

    def variant_summary(mask: pd.Series) -> str:
        subset = table.loc[mask].copy()
        if subset.empty:
            return ""
        labels = []
        for _, row in subset.sort_values(["gene", "variation_id", "variant_key"]).iterrows():
            gene = str(row.get("gene", ""))
            variation = str(row.get("clinvar_id", row.get("variation_id", "")))
            score = row.get("vital_score", np.nan)
            score_label = f"{score:.1f}" if pd.notna(score) else "NA"
            labels.append(f"{gene}:{variation}:VITAL={score_label}")
        return "; ".join(labels)

    rows: list[dict[str, object]] = []
    for scope_name, scope_mask in scopes.items():
        frequency_signal = scope_mask & (table["max_frequency_signal"] > AF_ULTRA_RARE)
        reference_ac_supported = frequency_signal & (table["qualifying_frequency_ac"].fillna(0) >= MIN_RECLASSIFICATION_AC)
        reference_red_mask = (
            (table["vital_score"] >= VITAL_ACTION_THRESHOLD)
            & table["weak_review_signal"]
            & reference_ac_supported
        )
        reference_red_keys = set(table.loc[reference_red_mask, "variant_key"].fillna("").astype(str))
        for ac_threshold in [5, 10, 20, 50]:
            ac_supported = frequency_signal & (table["qualifying_frequency_ac"].fillna(0) >= ac_threshold)
            vital_red_mask = (
                (table["vital_score"] >= VITAL_ACTION_THRESHOLD)
                & table["weak_review_signal"]
                & ac_supported
            )
            threshold_red_keys = set(table.loc[vital_red_mask, "variant_key"].fillna("").astype(str))
            shared_keys = threshold_red_keys & reference_red_keys
            gained_keys = threshold_red_keys - reference_red_keys
            lost_keys = reference_red_keys - threshold_red_keys
            red_composition = {
                "vital_red_count_at_ac_threshold": int(vital_red_mask.sum()),
                "vital_red_variants": variant_summary(vital_red_mask),
                "vital_red_shared_with_ac_ge_20_count": len(shared_keys),
                "vital_red_shared_with_ac_ge_20_variants": variant_summary(
                    table["variant_key"].astype(str).isin(shared_keys)
                ),
                "vital_red_gained_vs_ac_ge_20_count": len(gained_keys),
                "vital_red_gained_vs_ac_ge_20_variants": variant_summary(
                    table["variant_key"].astype(str).isin(gained_keys)
                ),
                "vital_red_lost_vs_ac_ge_20_count": len(lost_keys),
                "vital_red_lost_vs_ac_ge_20_variants": variant_summary(
                    table["variant_key"].astype(str).isin(lost_keys)
                ),
            }
            predictions = {
                "naive_AF_gt_1e-5": frequency_signal,
                f"AF_gt_1e-5_and_AC_ge_{ac_threshold}": ac_supported,
                f"VITAL_red_AC_ge_{ac_threshold}": vital_red_mask,
            }
            for target_name, labels in targets.items():
                scoped_labels = labels.loc[scope_mask]
                target_events = int(scoped_labels.sum())
                for method, predicted_all in predictions.items():
                    predicted = predicted_all.loc[scope_mask].fillna(False).astype(bool)
                    tp = int((predicted & scoped_labels).sum())
                    fp = int((predicted & ~scoped_labels).sum())
                    fn = int((~predicted & scoped_labels).sum())
                    tn = int((~predicted & ~scoped_labels).sum())
                    row = {
                        "scope": scope_name,
                        "target": target_name,
                        "ac_threshold": ac_threshold,
                        "method": method,
                        "ac_threshold_role": "actionability_filter_not_scoring_function",
                        "variant_count_scope": int(scope_mask.sum()),
                        "target_event_count": target_events,
                        "naive_af_flag_count_scope": int(frequency_signal.sum()),
                        "ac_supported_frequency_flag_count_scope": int(ac_supported.sum()),
                        **red_composition,
                        "true_positives": tp,
                        "false_positives": fp,
                        "false_negatives": fn,
                        "true_negatives": tn,
                        "precision": safe_divide(tp, tp + fp),
                        "recall": safe_divide(tp, tp + fn),
                        "specificity": safe_divide(tn, tn + fp),
                        "false_positive_rate": safe_divide(fp, fp + tn),
                        "false_negative_rate": safe_divide(fn, fn + tp),
                        "confidence_interval_method": "Wilson 95%",
                    }
                    rows.append(add_binary_metric_cis(row, tp, fp, fn, tn))
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    naive_fp = (
        result.loc[result["method"].eq("naive_AF_gt_1e-5"), ["scope", "target", "ac_threshold", "false_positives"]]
        .rename(columns={"false_positives": "naive_false_positives"})
    )
    result = result.merge(naive_fp, on=["scope", "target", "ac_threshold"], how="left")
    result["false_positive_reduction_vs_naive"] = result["naive_false_positives"] - result["false_positives"]
    result["false_positive_reduction_vs_naive_percent"] = result.apply(
        lambda row: safe_divide(row["false_positive_reduction_vs_naive"], row["naive_false_positives"]) * 100
        if row["naive_false_positives"]
        else np.nan,
        axis=1,
    )
    return result


def plot_historical_curves(curves: pd.DataFrame, output_path: Path) -> None:
    if curves.empty:
        return
    output_path.parent.mkdir(parents=True, exist_ok=True)
    targets = list(curves["target"].drop_duplicates())
    fig, axes = plt.subplots(len(targets), 2, figsize=(10, 4 * len(targets)))
    if len(targets) == 1:
        axes = np.array([axes])
    for row_index, target in enumerate(targets):
        target_df = curves[curves["target"].eq(target)].copy()
        roc_df = target_df.sort_values("fpr")
        roc_auc = np.trapezoid(roc_df["tpr_recall"].fillna(0), roc_df["fpr"].fillna(0))
        pr_df = target_df.sort_values("tpr_recall")
        pr_auc = np.trapezoid(pr_df["precision"].fillna(0), pr_df["tpr_recall"].fillna(0))
        axes[row_index, 0].plot(roc_df["fpr"], roc_df["tpr_recall"], color="#2563EB", linewidth=2)
        axes[row_index, 0].plot([0, 1], [0, 1], color="#9CA3AF", linestyle="--", linewidth=1)
        axes[row_index, 0].set_title(f"{target} ROC (AUC={roc_auc:.2f})")
        axes[row_index, 0].set_xlabel("False positive rate")
        axes[row_index, 0].set_ylabel("Recall")
        axes[row_index, 1].plot(pr_df["tpr_recall"], pr_df["precision"].fillna(1), color="#DC2626", linewidth=2)
        axes[row_index, 1].set_title(f"{target} PR (AUC={pr_auc:.2f})")
        axes[row_index, 1].set_xlabel("Recall")
        axes[row_index, 1].set_ylabel("Precision")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    print(f"Saved {output_path}")


def evaluate_reclassification(
    score_table: Path,
    baseline_summary: Path,
    followup_summary: Path,
    output_prefix: str,
    genes: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    scores = pd.read_csv(score_table)
    scores["variation_id"] = scores["variation_id"].astype(str)
    scores["vital_score"] = pd.to_numeric(scores["vital_score"], errors="coerce")
    if "vital_red_flag" in scores.columns:
        scores["vital_red_flag"] = scores["vital_red_flag"].astype(str).str.lower().eq("true")
    else:
        scores["vital_red_flag"] = False

    baseline = read_variant_summary(baseline_summary, genes)
    followup = read_variant_summary(followup_summary, genes)
    baseline_plp = baseline[baseline["clinical_group"].eq("P_LP")].copy()
    followup_groups = followup.loc[:, ["variation_id", "clinical_group", "clinical_significance", "review_status"]].rename(
        columns={
            "clinical_group": "followup_clinical_group",
            "clinical_significance": "followup_clinical_significance",
            "review_status": "followup_review_status",
        }
    )
    evaluated = baseline_plp.merge(followup_groups, on="variation_id", how="left")
    evaluated["followup_clinical_group"] = evaluated["followup_clinical_group"].fillna("missing")
    evaluated["strict_revised_to_b_or_vus"] = evaluated["followup_clinical_group"].isin(["B_LB", "VUS"])
    evaluated["broad_revised_or_destabilized"] = evaluated["followup_clinical_group"].isin(
        ["B_LB", "VUS", "conflicting", "other", "missing"]
    )
    followup_present = ~evaluated["followup_clinical_group"].eq("missing")
    evaluated["clinical_significance_changed"] = followup_present & (
        normalize_change_text(evaluated["clinical_significance"])
        != normalize_change_text(evaluated["followup_clinical_significance"])
    )
    evaluated["review_status_changed"] = followup_present & (
        normalize_change_text(evaluated["review_status"])
        != normalize_change_text(evaluated["followup_review_status"])
    )
    evaluated["expanded_revised_or_review_changed"] = (
        evaluated["broad_revised_or_destabilized"]
        | evaluated["clinical_significance_changed"]
        | evaluated["review_status_changed"]
    )
    evaluated = evaluated.merge(
        scores,
        on="variation_id",
        how="inner",
        suffixes=("_baseline", "_score"),
    )

    targets = historical_targets(evaluated)
    score = evaluated["vital_score"].fillna(0)
    red = evaluated["vital_red_flag"].astype(bool)

    def metric_rows(target: str, labels: pd.Series) -> list[dict[str, object]]:
        positives = int(labels.sum())
        red_total = int(red.sum())
        red_tp = int((red & labels).sum())
        precision_low, precision_high = wilson_ci(red_tp, red_total)
        recall_low, recall_high = wilson_ci(red_tp, positives)
        return [
            {"target": target, "metric": "scored_baseline_plp_count", "value": len(evaluated)},
            {"target": target, "metric": "revision_event_count", "value": positives},
            {"target": target, "metric": "roc_auc_vital_score", "value": roc_auc(labels, score)},
            {"target": target, "metric": "average_precision_vital_score", "value": average_precision(labels, score)},
            {"target": target, "metric": "precision_at_top_decile", "value": precision_at_fraction(labels, score, 0.10)},
            {"target": target, "metric": "vital_red_precision", "value": red_tp / red_total if red_total else np.nan},
            {"target": target, "metric": "vital_red_precision_ci_low", "value": precision_low},
            {"target": target, "metric": "vital_red_precision_ci_high", "value": precision_high},
            {"target": target, "metric": "vital_red_recall", "value": red_tp / positives if positives else np.nan},
            {"target": target, "metric": "vital_red_recall_ci_low", "value": recall_low},
            {"target": target, "metric": "vital_red_recall_ci_high", "value": recall_high},
            {"target": target, "metric": "vital_red_count", "value": red_total},
            {"target": target, "metric": "confidence_interval_method", "value": "Wilson 95% for vital_red precision and recall"},
        ]

    metrics = pd.DataFrame(
        [row for target_name, labels in targets.items() for row in metric_rows(target_name, labels)]
    )
    prediction_columns = [
        "variation_id",
        "gene_baseline",
        "gene_score",
        "variant_key_baseline",
        "variant_key_score",
        "clinvar_id",
        "name",
        "clinical_group",
        "followup_clinical_group",
        "followup_clinical_significance",
        "review_status_baseline",
        "review_status_score",
        "review_strength",
        "review_score",
        "submitter_count",
        "match_category",
        "frequency_evidence_status",
        "strict_revised_to_b_or_vus",
        "broad_revised_or_destabilized",
        "clinical_significance_changed",
        "review_status_changed",
        "expanded_revised_or_review_changed",
        "global_af",
        "global_ac",
        "popmax_af",
        "popmax_ac",
        "frequency_signal_ac_ge_20",
        "vital_score",
        "vital_band",
        "vital_red_flag",
        "max_frequency_signal",
        "qualifying_frequency_ac",
        "vital_signal_reason",
    ]
    available = [column for column in prediction_columns if column in evaluated.columns]
    predictions = evaluated.loc[:, available].sort_values("vital_score", ascending=False)
    method_comparison, curve_points = make_historical_method_comparison(evaluated)
    threshold_sweep = vital_score_threshold_sweep(evaluated)
    enrichment = make_reclassification_enrichment(evaluated)
    ac_threshold_sensitivity = make_historical_ac_threshold_sensitivity(evaluated)
    save_table(metrics, data_path(output_prefix, "vital_historical_validation.csv"))
    save_table(predictions, data_path(output_prefix, "vital_historical_predictions.csv"))
    save_table(method_comparison, data_path(output_prefix, "vital_historical_method_comparison.csv"))
    save_table(curve_points, data_path(output_prefix, "vital_historical_curve_points.csv"))
    save_table(threshold_sweep, data_path(output_prefix, "vital_historical_threshold_sweep.csv"))
    save_table(enrichment, data_path(output_prefix, "vital_historical_enrichment.csv"))
    save_table(ac_threshold_sensitivity, data_path(output_prefix, "vital_historical_ac_threshold_sensitivity.csv"))
    plot_historical_curves(
        curve_points,
        figure_path(output_prefix, "vital_historical_curves.png"),
    )
    return metrics, predictions


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate VITAL reclassification predictions across two ClinVar "
            "variant_summary snapshots. The score table should be generated "
            "from the baseline snapshot."
        )
    )
    parser.add_argument("--score-table", type=Path, default=data_path("arrhythmia", "vital_scores.csv"))
    parser.add_argument("--baseline-variant-summary", type=Path, required=True)
    parser.add_argument("--followup-variant-summary", type=Path, required=True)
    parser.add_argument("--output-prefix", default="arrhythmia")
    parser.add_argument("--genes", nargs="+", default=ARRHYTHMIA_GENES)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    metrics, _ = evaluate_reclassification(
        score_table=args.score_table,
        baseline_summary=args.baseline_variant_summary,
        followup_summary=args.followup_variant_summary,
        output_prefix=args.output_prefix,
        genes=args.genes,
    )
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
