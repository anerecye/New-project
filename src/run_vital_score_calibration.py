from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
FIGURE_DIR = BASE_DIR / "figures"
SUPPLEMENTARY_DIR = BASE_DIR / "supplementary_tables"
INPUT = DATA_DIR / "vital_cross_disease_3000_2023_01_to_current_vital_historical_predictions.csv"
SCORE_INPUT = DATA_DIR / "vital_cross_disease_3000_2023_01_vital_scores.csv"


TARGETS = {
    "strict": "strict_revised_to_b_or_vus",
    "broad": "broad_revised_or_destabilized",
    "expanded": "expanded_revised_or_review_changed",
}

COMPONENT_COLUMNS = [
    "frequency_pressure_score",
    "ac_reliability_score",
    "popmax_enrichment_score",
    "variant_type_tension_score",
    "technical_detectability_score",
    "gene_constraint_score",
    "review_fragility_score",
]

REVIEW_FEATURES = [
    "review_fragility_score",
    "single_submitter_feature",
    "weak_review_feature",
]

TENSION_FEATURES = [
    "frequency_pressure_score",
    "ac_reliability_score",
    "popmax_enrichment_score",
    "variant_type_tension_score",
    "technical_detectability_score",
    "gene_constraint_score",
]

RESTRICTED_MODELS = {
    "review_only": REVIEW_FEATURES,
    "tension_only": TENSION_FEATURES,
    "review_plus_tension": REVIEW_FEATURES + TENSION_FEATURES,
}

GRID_GROUPS = {
    "af_pressure_multiplier": ["frequency_pressure_score"],
    "ac_reliability_multiplier": ["ac_reliability_score"],
    "review_fragility_multiplier": ["review_fragility_score"],
    "context_multiplier": [
        "popmax_enrichment_score",
        "variant_type_tension_score",
        "technical_detectability_score",
        "gene_constraint_score",
    ],
}


def bool_series(series: pd.Series) -> pd.Series:
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def pava_fit(x: np.ndarray, y: np.ndarray) -> pd.DataFrame:
    order = np.argsort(x, kind="mergesort")
    x_sorted = x[order]
    y_sorted = y[order]
    blocks: list[dict[str, float]] = []
    for xi, yi in zip(x_sorted, y_sorted):
        blocks.append({"lower": float(xi), "upper": float(xi), "sum": float(yi), "weight": 1.0})
        while len(blocks) >= 2:
            left = blocks[-2]
            right = blocks[-1]
            if left["sum"] / left["weight"] <= right["sum"] / right["weight"]:
                break
            merged = {
                "lower": left["lower"],
                "upper": right["upper"],
                "sum": left["sum"] + right["sum"],
                "weight": left["weight"] + right["weight"],
            }
            blocks[-2:] = [merged]
    rows = []
    for block in blocks:
        rows.append(
            {
                "score_min": block["lower"],
                "score_max": block["upper"],
                "n": int(block["weight"]),
                "calibrated_probability": block["sum"] / block["weight"],
            }
        )
    return pd.DataFrame(rows)


def pava_predict(model: pd.DataFrame, x: np.ndarray) -> np.ndarray:
    preds = np.zeros(len(x), dtype=float)
    for i, value in enumerate(x):
        match = model[(model["score_min"] <= value) & (value <= model["score_max"])]
        if match.empty:
            if value < model["score_min"].min():
                preds[i] = float(model.iloc[0]["calibrated_probability"])
            else:
                preds[i] = float(model.iloc[-1]["calibrated_probability"])
        else:
            preds[i] = float(match.iloc[0]["calibrated_probability"])
    return preds


def platt_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float, float]:
    mean = float(np.mean(x))
    sd = float(np.std(x)) or 1.0
    z = (x - mean) / sd

    def objective(params: np.ndarray) -> float:
        a, b = params
        logits = np.clip(a * z + b, -35, 35)
        probs = 1 / (1 + np.exp(-logits))
        eps = 1e-9
        return float(-np.sum(y * np.log(probs + eps) + (1 - y) * np.log(1 - probs + eps)))

    result = minimize(objective, np.array([1.0, np.log((y.mean() + 1e-6) / (1 - y.mean() + 1e-6))]), method="BFGS")
    a, b = result.x
    return float(a), float(b), mean, sd


def platt_predict(params: tuple[float, float, float, float], x: np.ndarray) -> np.ndarray:
    a, b, mean, sd = params
    logits = np.clip(a * ((x - mean) / sd) + b, -35, 35)
    return 1 / (1 + np.exp(-logits))


def stratified_folds(y: np.ndarray, n_folds: int = 5) -> np.ndarray:
    rng = np.random.default_rng(42)
    folds = np.zeros(len(y), dtype=int)
    for label in [0, 1]:
        idx = np.where(y == label)[0]
        rng.shuffle(idx)
        for fold, part in enumerate(np.array_split(idx, n_folds)):
            folds[part] = fold
    return folds


def brier_score(y: np.ndarray, pred: np.ndarray) -> float:
    return float(np.mean((pred - y) ** 2))


def expected_calibration_error(y: np.ndarray, pred: np.ndarray, n_bins: int = 10) -> float:
    edges = np.linspace(0, 1, n_bins + 1)
    total = len(y)
    ece = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (pred >= lo) & (pred <= hi if hi == 1 else pred < hi)
        if not mask.any():
            continue
        ece += (mask.sum() / total) * abs(float(y[mask].mean()) - float(pred[mask].mean()))
    return float(ece)


def roc_auc(y: np.ndarray, score: np.ndarray) -> float:
    y_bool = y.astype(bool)
    n_pos = int(y_bool.sum())
    n_neg = int((~y_bool).sum())
    if n_pos == 0 or n_neg == 0:
        return np.nan
    ranks = pd.Series(score).rank(method="average").to_numpy()
    rank_sum_pos = float(ranks[y_bool].sum())
    return float((rank_sum_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def average_precision(y: np.ndarray, score: np.ndarray) -> float:
    ordered = pd.DataFrame({"label": y.astype(bool), "score": score}).sort_values("score", ascending=False)
    n_pos = int(ordered["label"].sum())
    if n_pos == 0:
        return np.nan
    ordered["cum_pos"] = ordered["label"].cumsum()
    ordered["rank"] = np.arange(1, len(ordered) + 1)
    return float((ordered.loc[ordered["label"], "cum_pos"] / ordered.loc[ordered["label"], "rank"]).sum() / n_pos)


def spearman_rho(a: np.ndarray, b: np.ndarray) -> float:
    rank_a = pd.Series(a).rank(method="average").to_numpy(dtype=float)
    rank_b = pd.Series(b).rank(method="average").to_numpy(dtype=float)
    if np.std(rank_a) == 0 or np.std(rank_b) == 0:
        return np.nan
    return float(np.corrcoef(rank_a, rank_b)[0, 1])


def calibrate_target(df: pd.DataFrame, target_name: str, target_col: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    x = pd.to_numeric(df["vital_score"], errors="coerce").fillna(0).to_numpy(dtype=float)
    y = bool_series(df[target_col]).to_numpy(dtype=int)
    folds = stratified_folds(y)
    iso_oof = np.zeros(len(df), dtype=float)
    platt_oof = np.zeros(len(df), dtype=float)
    for fold in range(5):
        train = folds != fold
        test = folds == fold
        iso_model = pava_fit(x[train], y[train])
        iso_oof[test] = pava_predict(iso_model, x[test])
        platt_params = platt_fit(x[train], y[train])
        platt_oof[test] = platt_predict(platt_params, x[test])

    full_iso = pava_fit(x, y)
    full_platt = platt_fit(x, y)
    score_points = np.array([40, 50, 60, 70, 74, 80, 90], dtype=float)
    mapping = pd.DataFrame(
        {
            "target": target_name,
            "vital_score": score_points,
            "isotonic_probability": pava_predict(full_iso, score_points),
            "platt_probability": platt_predict(full_platt, score_points),
            "calibration_note": "probabilities_are_for_database_instability_endpoint_not_pathogenicity",
        }
    )

    metrics = pd.DataFrame(
        [
            {
                "target": target_name,
                "n": len(y),
                "event_count": int(y.sum()),
                "event_rate": float(y.mean()),
                "isotonic_brier_oof": brier_score(y, iso_oof),
                "platt_brier_oof": brier_score(y, platt_oof),
                "uncalibrated_brier_oof": brier_score(y, np.repeat(y.mean(), len(y))),
                "isotonic_ece_oof": expected_calibration_error(y, iso_oof),
                "platt_ece_oof": expected_calibration_error(y, platt_oof),
                "calibration_method": "5-fold out-of-fold isotonic PAVA and Platt logistic scaling on VITAL score only",
            }
        ]
    )

    bins = pd.cut(
        df["vital_score"],
        bins=[-0.001, 20, 40, 60, 70, 80, 100],
        labels=["0-20", "20-40", "40-60", "60-70", "70-80", "80-100"],
        include_lowest=True,
    )
    bin_rows = []
    for band, sub_index in pd.Series(df.index, index=bins).groupby(level=0, observed=True):
        idx = sub_index.to_numpy()
        if len(idx) == 0:
            continue
        bin_rows.append(
            {
                "target": target_name,
                "score_bin": str(band),
                "n": len(idx),
                "event_count": int(y[idx].sum()),
                "observed_event_rate": float(y[idx].mean()),
                "mean_isotonic_oof_probability": float(iso_oof[idx].mean()),
                "mean_platt_oof_probability": float(platt_oof[idx].mean()),
            }
        )
    return metrics, mapping, pd.DataFrame(bin_rows)


def build_empirical_frame(historical: pd.DataFrame) -> pd.DataFrame:
    scores = pd.read_csv(SCORE_INPUT)
    keep = ["variation_id", "weak_review_signal", *COMPONENT_COLUMNS]
    missing = [column for column in keep if column not in scores.columns]
    if missing:
        raise ValueError(f"Missing columns in {SCORE_INPUT}: {missing}")
    component_data = scores.loc[:, keep].copy()
    component_data["variation_id"] = component_data["variation_id"].astype(str)
    component_data = component_data.drop_duplicates("variation_id", keep="first")

    frame = historical.copy()
    frame["variation_id"] = frame["variation_id"].astype(str)
    frame = frame.merge(component_data, on="variation_id", how="left", validate="one_to_one")
    for column in COMPONENT_COLUMNS:
        frame[column] = pd.to_numeric(frame[column], errors="coerce").fillna(0.0)
    if "weak_review_signal" not in frame.columns or frame["weak_review_signal"].isna().all():
        frame["weak_review_signal"] = frame["review_strength"].isin(["single_submitter", "weak_or_no_assertion"])
    frame["weak_review_signal"] = bool_series(frame["weak_review_signal"])
    frame["weak_review_feature"] = frame["weak_review_signal"].astype(float)
    frame["single_submitter_feature"] = (
        pd.to_numeric(frame.get("submitter_count"), errors="coerce").fillna(99).le(1).astype(float)
    )
    return frame


def standardize_matrix(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = matrix.mean(axis=0)
    sd = matrix.std(axis=0)
    sd[sd == 0] = 1.0
    return (matrix - mean) / sd, mean, sd


def fit_nonnegative_logistic(x: np.ndarray, y: np.ndarray, l2: float = 0.25) -> np.ndarray:
    event_rate = np.clip(y.mean(), 1e-6, 1 - 1e-6)
    initial = np.zeros(x.shape[1] + 1, dtype=float)
    initial[0] = np.log(event_rate / (1 - event_rate))

    def objective(params: np.ndarray) -> float:
        intercept = params[0]
        coefficients = params[1:]
        logits = np.clip(intercept + x @ coefficients, -35, 35)
        probs = 1 / (1 + np.exp(-logits))
        eps = 1e-9
        loss = -np.sum(y * np.log(probs + eps) + (1 - y) * np.log(1 - probs + eps))
        penalty = 0.5 * l2 * float(np.sum(coefficients**2))
        return float(loss + penalty)

    bounds = [(None, None)] + [(0.0, None) for _ in range(x.shape[1])]
    result = minimize(objective, initial, method="L-BFGS-B", bounds=bounds)
    if not result.success:
        print(f"Warning: empirical logistic calibration did not fully converge: {result.message}")
    return result.x.astype(float)


def logistic_predict(params: np.ndarray, x: np.ndarray) -> np.ndarray:
    logits = np.clip(params[0] + x @ params[1:], -35, 35)
    return 1 / (1 + np.exp(-logits))


def empirical_logistic_calibration(df: pd.DataFrame, calibration_scope: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    y = bool_series(df["expanded_revised_or_review_changed"]).to_numpy(dtype=int)
    x_raw = df[COMPONENT_COLUMNS].to_numpy(dtype=float)
    x, _, _ = standardize_matrix(x_raw)
    folds = stratified_folds(y)
    oof = np.zeros(len(df), dtype=float)
    for fold in range(5):
        train = folds != fold
        test = folds == fold
        x_train, train_mean, train_sd = standardize_matrix(x_raw[train])
        x_test = (x_raw[test] - train_mean) / train_sd
        params = fit_nonnegative_logistic(x_train, y[train])
        oof[test] = logistic_predict(params, x_test)

    x_full, _, _ = standardize_matrix(x_raw)
    full_params = fit_nonnegative_logistic(x_full, y)
    full_prob = logistic_predict(full_params, x_full)
    expert_score = pd.to_numeric(df["vital_score"], errors="coerce").fillna(0).to_numpy(dtype=float)
    red = bool_series(df["vital_red_flag"])
    n_red = int(red.sum())
    top_n = max(n_red, 1)
    top_logistic_idx = np.argsort(-full_prob)[:top_n]
    expert_red_ids = set(df.loc[red, "variation_id"].astype(str))
    logistic_top_ids = set(df.iloc[top_logistic_idx]["variation_id"].astype(str))

    metrics = pd.DataFrame(
        [
            {
                "model": "expert_vital_score",
                "calibration_scope": calibration_scope,
                "target": "expanded",
                "n": len(df),
                "event_count": int(y.sum()),
                "roc_auc": roc_auc(y, expert_score),
                "average_precision": average_precision(y, expert_score),
                "brier_score": brier_score(y, np.clip(expert_score / 100, 0, 1)),
                "spearman_vs_expert_score": 1.0,
                "method_note": "Expert-specified VITAL score; not fitted to historical endpoint",
            },
            {
                "model": "component_logistic_oof",
                "calibration_scope": calibration_scope,
                "target": "expanded",
                "n": len(df),
                "event_count": int(y.sum()),
                "roc_auc": roc_auc(y, oof),
                "average_precision": average_precision(y, oof),
                "brier_score": brier_score(y, oof),
                "spearman_vs_expert_score": spearman_rho(expert_score, oof),
                "method_note": "5-fold out-of-fold non-negative logistic regression using VITAL component scores",
            },
            {
                "model": "component_logistic_full_fit_for_ranking_audit",
                "calibration_scope": calibration_scope,
                "target": "expanded",
                "n": len(df),
                "event_count": int(y.sum()),
                "roc_auc": roc_auc(y, full_prob),
                "average_precision": average_precision(y, full_prob),
                "brier_score": brier_score(y, full_prob),
                "spearman_vs_expert_score": spearman_rho(expert_score, full_prob),
                "method_note": "Full-data fit used only to compare component weighting and top-N queue overlap",
            },
        ]
    )

    coefficients = pd.DataFrame(
        {
            "calibration_scope": calibration_scope,
            "component_column": COMPONENT_COLUMNS,
            "standardized_nonnegative_logistic_coefficient": full_params[1:],
        }
    )
    coefficient_sum = coefficients["standardized_nonnegative_logistic_coefficient"].sum()
    coefficients["normalized_empirical_weight_fraction"] = (
        coefficients["standardized_nonnegative_logistic_coefficient"] / coefficient_sum
        if coefficient_sum > 0
        else np.nan
    )
    coefficients["expert_max_points"] = [45, 20, 10, 7, 8, 10, 10]
    coefficients["expert_weight_fraction"] = coefficients["expert_max_points"] / coefficients["expert_max_points"].sum()
    coefficients["calibration_note"] = (
        "Component weights are design-prespecified; this audit checks consistency and does not refit VITAL"
    )

    overlap = pd.DataFrame(
        [
            {
                "comparison": "expert_red_vs_topN_component_logistic",
                "calibration_scope": calibration_scope,
                "n_expert_red": n_red,
                "n_top_logistic": top_n,
                "overlap_count": len(expert_red_ids & logistic_top_ids),
                "overlap_fraction_of_expert_red": len(expert_red_ids & logistic_top_ids) / n_red if n_red else np.nan,
                "expert_red_expanded_events": int(bool_series(df.loc[red, "expanded_revised_or_review_changed"]).sum()),
                "top_logistic_expanded_events": int(y[top_logistic_idx].sum()),
                "topN_rule": "Top N by full-data empirical logistic probability, with N equal to expert red queue size",
            }
        ]
    )
    return metrics, coefficients, overlap


def grid_search_weight_multipliers(df: pd.DataFrame, calibration_scope: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    y = bool_series(df["expanded_revised_or_review_changed"]).to_numpy(dtype=int)
    expert_score = pd.to_numeric(df["vital_score"], errors="coerce").fillna(0).to_numpy(dtype=float)
    red = bool_series(df["vital_red_flag"])
    expert_red_ids = set(df.loc[red, "variation_id"].astype(str))
    n_red = max(int(red.sum()), 1)
    multipliers = [0.5, 0.75, 1.0, 1.25, 1.5]
    rows = []
    for af_multiplier in multipliers:
        for ac_multiplier in multipliers:
            for review_multiplier in multipliers:
                for context_multiplier in multipliers:
                    score = np.zeros(len(df), dtype=float)
                    multiplier_map = {
                        "af_pressure_multiplier": af_multiplier,
                        "ac_reliability_multiplier": ac_multiplier,
                        "review_fragility_multiplier": review_multiplier,
                        "context_multiplier": context_multiplier,
                    }
                    for group_name, columns in GRID_GROUPS.items():
                        score += multiplier_map[group_name] * df[columns].sum(axis=1).to_numpy(dtype=float)
                    score = np.clip(score, 0, 100)
                    top_idx = np.argsort(-score)[:n_red]
                    top_ids = set(df.iloc[top_idx]["variation_id"].astype(str))
                    rows.append(
                        {
                            "calibration_scope": calibration_scope,
                            **multiplier_map,
                            "roc_auc": roc_auc(y, score),
                            "average_precision": average_precision(y, score),
                            "spearman_vs_expert_score": spearman_rho(expert_score, score),
                            "topN_overlap_with_expert_red": len(top_ids & expert_red_ids),
                            "topN_expanded_events": int(y[top_idx].sum()),
                            "is_expert_weight_profile": bool(
                                af_multiplier == ac_multiplier == review_multiplier == context_multiplier == 1.0
                            ),
                        }
                    )
    grid = pd.DataFrame(rows).sort_values(
        ["average_precision", "roc_auc", "spearman_vs_expert_score"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    expert_row = grid[grid["is_expert_weight_profile"]].copy()
    if not expert_row.empty:
        expert_row["rank_by_average_precision"] = int(expert_row.index[0] + 1)
        expert_row["rank_denominator"] = len(grid)
    return grid, expert_row


def fit_oof_logistic(df: pd.DataFrame, features: list[str], target_col: str) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    y = bool_series(df[target_col]).to_numpy(dtype=int)
    x_raw = df[features].to_numpy(dtype=float)
    folds = stratified_folds(y)
    oof = np.zeros(len(df), dtype=float)
    for fold in range(5):
        train = folds != fold
        test = folds == fold
        x_train, train_mean, train_sd = standardize_matrix(x_raw[train])
        x_test = (x_raw[test] - train_mean) / train_sd
        params = fit_nonnegative_logistic(x_train, y[train])
        oof[test] = logistic_predict(params, x_test)
    x_full, _, _ = standardize_matrix(x_raw)
    params = fit_nonnegative_logistic(x_full, y)
    full_pred = logistic_predict(params, x_full)
    coefficients = pd.DataFrame(
        {
            "feature": features,
            "standardized_nonnegative_logistic_coefficient": params[1:],
        }
    )
    return oof, full_pred, coefficients


def bootstrap_metric_ci(y: np.ndarray, pred: np.ndarray, n_iter: int = 1000) -> dict[str, float]:
    rng = np.random.default_rng(42)
    values = {"roc_auc": [], "average_precision": [], "brier_score": []}
    indices = np.arange(len(y))
    for _ in range(n_iter):
        sample = rng.choice(indices, size=len(indices), replace=True)
        y_sample = y[sample]
        pred_sample = pred[sample]
        if len(np.unique(y_sample)) < 2:
            continue
        values["roc_auc"].append(roc_auc(y_sample, pred_sample))
        values["average_precision"].append(average_precision(y_sample, pred_sample))
        values["brier_score"].append(brier_score(y_sample, pred_sample))
    ci: dict[str, float] = {}
    for metric, metric_values in values.items():
        array = np.asarray(metric_values, dtype=float)
        ci[f"{metric}_ci_low"] = float(np.nanpercentile(array, 2.5)) if len(array) else np.nan
        ci[f"{metric}_ci_high"] = float(np.nanpercentile(array, 97.5)) if len(array) else np.nan
    return ci


def calibration_curve_table(y: np.ndarray, pred: np.ndarray, model_name: str, n_bins: int = 5) -> pd.DataFrame:
    table = pd.DataFrame({"y": y.astype(int), "pred": pred})
    table["bin"] = pd.qcut(table["pred"].rank(method="first"), q=n_bins, labels=False, duplicates="drop")
    rows = []
    for bin_id, sub in table.groupby("bin", dropna=False):
        rows.append(
            {
                "model": model_name,
                "calibration_bin": int(bin_id) + 1 if pd.notna(bin_id) else np.nan,
                "n": len(sub),
                "observed_event_rate": float(sub["y"].mean()),
                "mean_predicted_probability": float(sub["pred"].mean()),
                "min_predicted_probability": float(sub["pred"].min()),
                "max_predicted_probability": float(sub["pred"].max()),
            }
        )
    return pd.DataFrame(rows)


def restricted_feature_group_calibration(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    y = bool_series(df["expanded_revised_or_review_changed"]).to_numpy(dtype=int)
    expert_score = pd.to_numeric(df["vital_score"], errors="coerce").fillna(0).to_numpy(dtype=float)
    model_predictions: dict[str, np.ndarray] = {}
    metrics_rows = []
    coefficient_rows = []
    calibration_rows = []
    for model_name, features in RESTRICTED_MODELS.items():
        oof, full_pred, coefficients = fit_oof_logistic(df, features, "expanded_revised_or_review_changed")
        model_predictions[model_name] = full_pred
        metric_ci = bootstrap_metric_ci(y, oof)
        metrics_rows.append(
            {
                "model": model_name,
                "calibration_scope": "frequency_positive_records",
                "target": "expanded",
                "n": len(df),
                "event_count": int(y.sum()),
                "roc_auc_oof": roc_auc(y, oof),
                "average_precision_oof": average_precision(y, oof),
                "brier_score_oof": brier_score(y, oof),
                "spearman_vs_expert_vital_score": spearman_rho(expert_score, full_pred),
                **metric_ci,
                "model_note": "5-fold out-of-fold non-negative logistic regression restricted to frequency-positive variants",
            }
        )
        coefficients["model"] = model_name
        coefficients["calibration_scope"] = "frequency_positive_records"
        coefficients["feature_group"] = np.where(coefficients["feature"].isin(REVIEW_FEATURES), "review", "tension")
        coefficient_rows.append(coefficients)
        calibration_rows.append(calibration_curve_table(y, oof, model_name))

    coefficients_df = pd.concat(coefficient_rows, ignore_index=True)
    metrics_df = pd.DataFrame(metrics_rows)
    calibration_df = pd.concat(calibration_rows, ignore_index=True)

    dominance_rows = []
    for model_name, sub in coefficients_df.groupby("model"):
        review_sum = float(
            sub.loc[sub["feature_group"].eq("review"), "standardized_nonnegative_logistic_coefficient"].sum()
        )
        tension_sum = float(
            sub.loc[sub["feature_group"].eq("tension"), "standardized_nonnegative_logistic_coefficient"].sum()
        )
        dominance_rows.append(
            {
                "model": model_name,
                "review_coefficient_sum": review_sum,
                "tension_coefficient_sum": tension_sum,
                "tension_to_review_ratio": tension_sum / review_sum if review_sum else np.inf,
                "dominant_group": "tension" if tension_sum > review_sum else "review",
                "interpretation_note": "Dominance is coefficient magnitude after standardization within frequency-positive variants",
            }
        )
    dominance_df = pd.DataFrame(dominance_rows)

    case_rows = []
    high_tension = df["vital_band"].isin(["red_reclassification_priority", "orange_high_tension"]).fillna(False)
    for _, row in df.loc[high_tension].iterrows():
        record = {
            "variation_id": row.get("variation_id"),
            "gene": row.get("gene_baseline", row.get("gene_score", "")),
            "vital_band": row.get("vital_band"),
            "vital_score": row.get("vital_score"),
            "expanded_instability": bool_series(pd.Series([row.get("expanded_revised_or_review_changed")])).iloc[0],
        }
        row_position = row.name
        for model_name, pred in model_predictions.items():
            ranks = pd.Series(pred).rank(method="min", ascending=False).to_numpy(dtype=int)
            record[f"{model_name}_probability"] = float(pred[row_position])
            record[f"{model_name}_rank"] = int(ranks[row_position])
        case_rows.append(record)
    case_df = pd.DataFrame(case_rows).sort_values(["vital_band", "vital_score"], ascending=[False, False])
    return metrics_df, coefficients_df, calibration_df, dominance_df, case_df


def ac_gate_audit(df: pd.DataFrame) -> pd.DataFrame:
    signal = (
        pd.to_numeric(df["max_frequency_signal"], errors="coerce").fillna(0).gt(1e-5)
        & df["frequency_evidence_status"].eq("frequency_observed")
    )
    ac = pd.to_numeric(df["qualifying_frequency_ac"], errors="coerce").fillna(0)
    expanded = bool_series(df["expanded_revised_or_review_changed"])
    broad = bool_series(df["broad_revised_or_destabilized"])
    strict = bool_series(df["strict_revised_to_b_or_vus"])
    rows = []
    for threshold in [1, 2, 5, 10, 20, 50, 100]:
        mask = signal & ac.ge(threshold)
        weak = mask & df["review_strength"].isin(["single_submitter", "weak_or_no_assertion"])
        rows.append(
            {
                "ac_threshold": threshold,
                "ac_supported_frequency_flags": int(mask.sum()),
                "weak_review_ac_supported_flags": int(weak.sum()),
                "strict_event_rate": float(strict[mask].mean()) if mask.any() else np.nan,
                "broad_event_rate": float(broad[mask].mean()) if mask.any() else np.nan,
                "expanded_event_rate": float(expanded[mask].mean()) if mask.any() else np.nan,
                "red_priority_like_count_score_ge70_weak": int((weak & pd.to_numeric(df["vital_score"], errors="coerce").fillna(0).ge(70)).sum()),
                "threshold_note": "AC threshold is audited as an actionability gate, not learned into the score",
            }
        )
    return pd.DataFrame(rows)


def plot_calibration(bins: pd.DataFrame, mapping: pd.DataFrame) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    expanded_bins = bins[bins["target"].eq("expanded")].copy()
    strict_bins = bins[bins["target"].eq("strict")].copy()
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.0))

    for ax, table, title in [
        (axes[0], expanded_bins, "Expanded instability calibration"),
        (axes[1], strict_bins, "Strict downgrade calibration"),
    ]:
        x = np.arange(len(table))
        ax.bar(x - 0.18, table["observed_event_rate"], width=0.36, label="Observed", color="#5b8db8")
        ax.bar(x + 0.18, table["mean_isotonic_oof_probability"], width=0.36, label="Isotonic OOF", color="#f4a261")
        ax.set_xticks(x)
        ax.set_xticklabels(table["score_bin"], rotation=0)
        max_value = 0.0
        if len(table):
            max_value = float(
                np.nanmax(
                    [
                        table["observed_event_rate"].max(),
                        table["mean_isotonic_oof_probability"].max(),
                    ]
                )
            )
        ax.set_ylim(0, min(1.0, max(0.08, max_value * 1.25)))
        ax.set_title(title, weight="bold")
        ax.set_xlabel("VITAL score bin")
        ax.set_ylabel("Event probability")
        ax.grid(axis="y", alpha=0.25)
    axes[0].legend(frameon=False)
    fig.suptitle("VITAL score calibration on the 3,000-variant historical cross-disease set", weight="bold")
    fig.tight_layout()
    path = FIGURE_DIR / "vital_cross_disease_3000_score_calibration.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_empirical_weight_calibration(coefficients: pd.DataFrame, grid: pd.DataFrame) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2))

    coeff = coefficients.copy()
    coeff["component_label"] = coeff["component_column"].str.replace("_score", "", regex=False).str.replace("_", " ")
    x = np.arange(len(coeff))
    axes[0].bar(x - 0.18, coeff["expert_weight_fraction"], width=0.36, label="Expert max-point fraction", color="#4d7c8a")
    axes[0].bar(
        x + 0.18,
        coeff["normalized_empirical_weight_fraction"].fillna(0),
        width=0.36,
        label="Empirical logistic fraction",
        color="#d9822b",
    )
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(coeff["component_label"], rotation=35, ha="right")
    axes[0].set_ylabel("Relative component weight")
    axes[0].set_title("Component-weight comparison", weight="bold")
    axes[0].grid(axis="y", alpha=0.25)
    axes[0].legend(frameon=False)

    grid_plot = grid.head(20).copy().iloc[::-1]
    labels = [
        f"AF {row.af_pressure_multiplier:g}, AC {row.ac_reliability_multiplier:g}, "
        f"Review {row.review_fragility_multiplier:g}, Ctx {row.context_multiplier:g}"
        for row in grid_plot.itertuples()
    ]
    colors = ["#d9822b" if row.is_expert_weight_profile else "#6b8fb3" for row in grid_plot.itertuples()]
    axes[1].barh(np.arange(len(grid_plot)), grid_plot["average_precision"], color=colors)
    axes[1].set_yticks(np.arange(len(grid_plot)))
    axes[1].set_yticklabels(labels, fontsize=7)
    axes[1].set_xlabel("Average precision for expanded instability")
    axes[1].set_title("Top grid-search profiles", weight="bold")
    axes[1].grid(axis="x", alpha=0.25)

    fig.suptitle("VITAL weight-consistency audit", weight="bold")
    fig.tight_layout()
    path = FIGURE_DIR / "vital_cross_disease_3000_empirical_weight_calibration.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    print(f"Saved {path}")


def plot_restricted_calibration_models(
    metrics: pd.DataFrame,
    calibration: pd.DataFrame,
    dominance: pd.DataFrame,
) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))

    metric_plot = metrics.set_index("model").loc[["review_only", "tension_only", "review_plus_tension"]]
    x = np.arange(len(metric_plot))
    axes[0].bar(x - 0.2, metric_plot["roc_auc_oof"], width=0.2, label="ROC AUC", color="#4d7c8a")
    axes[0].bar(x, metric_plot["average_precision_oof"], width=0.2, label="PR AUC", color="#d9822b")
    axes[0].bar(x + 0.2, 1 - metric_plot["brier_score_oof"], width=0.2, label="1 - Brier", color="#6a994e")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(["Review", "Tension", "All"], rotation=0)
    axes[0].set_ylim(0, 1)
    axes[0].set_title("Restricted model performance", weight="bold")
    axes[0].grid(axis="y", alpha=0.25)
    axes[0].legend(frameon=False, fontsize=8)

    for model_name, sub in calibration.groupby("model"):
        label = {"review_only": "Review", "tension_only": "Tension", "review_plus_tension": "All"}[model_name]
        axes[1].plot(
            sub["mean_predicted_probability"],
            sub["observed_event_rate"],
            marker="o",
            label=label,
        )
    axes[1].plot([0, 1], [0, 1], linestyle="--", color="#888888", linewidth=1)
    axes[1].set_xlabel("Mean predicted probability")
    axes[1].set_ylabel("Observed expanded instability")
    axes[1].set_title("Calibration curve", weight="bold")
    axes[1].grid(alpha=0.25)
    axes[1].legend(frameon=False, fontsize=8)

    dom = dominance.set_index("model").loc[["review_only", "tension_only", "review_plus_tension"]]
    x = np.arange(len(dom))
    axes[2].bar(x - 0.18, dom["review_coefficient_sum"], width=0.36, label="Review", color="#8d6e63")
    axes[2].bar(x + 0.18, dom["tension_coefficient_sum"], width=0.36, label="Frequency tension", color="#3f7cac")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(["Review", "Tension", "All"])
    axes[2].set_ylabel("Standardized coefficient sum")
    axes[2].set_title("Feature dominance", weight="bold")
    axes[2].grid(axis="y", alpha=0.25)
    axes[2].legend(frameon=False, fontsize=8)

    fig.suptitle("Restricted VITAL weight-consistency audit within frequency-positive variants", weight="bold")
    fig.tight_layout()
    path = FIGURE_DIR / "vital_cross_disease_3000_restricted_calibration_models.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    print(f"Saved {path}")


def main() -> None:
    df = pd.read_csv(INPUT)
    empirical_df = build_empirical_frame(df)
    frequency_positive = (
        pd.to_numeric(empirical_df["max_frequency_signal"], errors="coerce").fillna(0).gt(1e-5)
        & empirical_df["frequency_evidence_status"].eq("frequency_observed")
    )
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    SUPPLEMENTARY_DIR.mkdir(parents=True, exist_ok=True)
    metrics_parts = []
    mapping_parts = []
    bin_parts = []
    for target_name, target_col in TARGETS.items():
        metrics, mapping, bins = calibrate_target(df, target_name, target_col)
        metrics_parts.append(metrics)
        mapping_parts.append(mapping)
        bin_parts.append(bins)
    metrics_df = pd.concat(metrics_parts, ignore_index=True)
    mapping_df = pd.concat(mapping_parts, ignore_index=True)
    bins_df = pd.concat(bin_parts, ignore_index=True)
    ac_df = ac_gate_audit(df)
    empirical_scopes = [
        ("all_records", empirical_df),
        ("frequency_positive_records", empirical_df.loc[frequency_positive].copy().reset_index(drop=True)),
    ]
    empirical_metrics_parts = []
    empirical_coefficient_parts = []
    empirical_overlap_parts = []
    grid_parts = []
    expert_grid_parts = []
    for scope_name, scope_df in empirical_scopes:
        scope_metrics, scope_coefficients, scope_overlap = empirical_logistic_calibration(scope_df, scope_name)
        scope_grid, scope_expert_grid_row = grid_search_weight_multipliers(scope_df, scope_name)
        empirical_metrics_parts.append(scope_metrics)
        empirical_coefficient_parts.append(scope_coefficients)
        empirical_overlap_parts.append(scope_overlap)
        grid_parts.append(scope_grid)
        expert_grid_parts.append(scope_expert_grid_row)
    empirical_metrics = pd.concat(empirical_metrics_parts, ignore_index=True)
    empirical_coefficients = pd.concat(empirical_coefficient_parts, ignore_index=True)
    empirical_overlap = pd.concat(empirical_overlap_parts, ignore_index=True)
    grid_df = pd.concat(grid_parts, ignore_index=True)
    expert_grid_row = pd.concat(expert_grid_parts, ignore_index=True)
    grid_top = grid_df.groupby("calibration_scope", group_keys=False).head(50).copy()
    restricted_df = empirical_df.loc[frequency_positive].copy().reset_index(drop=True)
    (
        restricted_metrics,
        restricted_coefficients,
        restricted_curve,
        restricted_dominance,
        restricted_cases,
    ) = restricted_feature_group_calibration(restricted_df)

    metrics_df.to_csv(DATA_DIR / "vital_cross_disease_3000_score_calibration_metrics.csv", index=False)
    mapping_df.to_csv(DATA_DIR / "vital_cross_disease_3000_score_calibration_mapping.csv", index=False)
    bins_df.to_csv(DATA_DIR / "vital_cross_disease_3000_score_calibration_bins.csv", index=False)
    ac_df.to_csv(DATA_DIR / "vital_cross_disease_3000_ac_gate_calibration.csv", index=False)
    empirical_metrics.to_csv(DATA_DIR / "vital_cross_disease_3000_empirical_weight_calibration_metrics.csv", index=False)
    empirical_coefficients.to_csv(DATA_DIR / "vital_cross_disease_3000_empirical_weight_coefficients.csv", index=False)
    empirical_overlap.to_csv(DATA_DIR / "vital_cross_disease_3000_empirical_weight_red_overlap.csv", index=False)
    grid_top.to_csv(DATA_DIR / "vital_cross_disease_3000_empirical_weight_grid_search_top.csv", index=False)
    expert_grid_row.to_csv(DATA_DIR / "vital_cross_disease_3000_empirical_weight_grid_search_expert_row.csv", index=False)
    restricted_metrics.to_csv(DATA_DIR / "vital_cross_disease_3000_restricted_calibration_model_metrics.csv", index=False)
    restricted_coefficients.to_csv(DATA_DIR / "vital_cross_disease_3000_restricted_calibration_coefficients.csv", index=False)
    restricted_curve.to_csv(DATA_DIR / "vital_cross_disease_3000_restricted_calibration_curve.csv", index=False)
    restricted_dominance.to_csv(DATA_DIR / "vital_cross_disease_3000_restricted_calibration_feature_dominance.csv", index=False)
    restricted_cases.to_csv(DATA_DIR / "vital_cross_disease_3000_restricted_calibration_high_tension_cases.csv", index=False)
    mapping_df.to_csv(SUPPLEMENTARY_DIR / "Supplementary_Table_S17_VITAL_score_calibration.tsv", sep="\t", index=False)
    s21 = pd.concat(
        [
            empirical_metrics.assign(table_section="model_metrics"),
            empirical_coefficients.assign(table_section="component_coefficients"),
            empirical_overlap.assign(table_section="red_queue_overlap"),
            expert_grid_row.assign(table_section="expert_grid_profile"),
            restricted_metrics.assign(table_section="restricted_model_metrics"),
            restricted_coefficients.assign(table_section="restricted_coefficients"),
            restricted_dominance.assign(table_section="restricted_feature_dominance"),
        ],
        ignore_index=True,
        sort=False,
    )
    s21.fillna("NA").to_csv(
        SUPPLEMENTARY_DIR / "Supplementary_Table_S21_VITAL_empirical_weight_calibration.tsv",
        sep="\t",
        index=False,
    )
    plot_calibration(bins_df, mapping_df)
    plot_empirical_weight_calibration(
        empirical_coefficients[empirical_coefficients["calibration_scope"].eq("frequency_positive_records")],
        grid_df[grid_df["calibration_scope"].eq("frequency_positive_records")].head(50),
    )
    plot_restricted_calibration_models(restricted_metrics, restricted_curve, restricted_dominance)
    print("Saved VITAL score calibration outputs")


if __name__ == "__main__":
    main()
