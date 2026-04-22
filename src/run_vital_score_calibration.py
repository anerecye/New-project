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


TARGETS = {
    "strict": "strict_revised_to_b_or_vus",
    "broad": "broad_revised_or_destabilized",
    "expanded": "expanded_revised_or_review_changed",
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
    print(f"Saved {path}")


def main() -> None:
    df = pd.read_csv(INPUT)
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

    metrics_df.to_csv(DATA_DIR / "vital_cross_disease_3000_score_calibration_metrics.csv", index=False)
    mapping_df.to_csv(DATA_DIR / "vital_cross_disease_3000_score_calibration_mapping.csv", index=False)
    bins_df.to_csv(DATA_DIR / "vital_cross_disease_3000_score_calibration_bins.csv", index=False)
    ac_df.to_csv(DATA_DIR / "vital_cross_disease_3000_ac_gate_calibration.csv", index=False)
    mapping_df.to_csv(SUPPLEMENTARY_DIR / "Supplementary_Table_S17_VITAL_score_calibration.tsv", sep="\t", index=False)
    plot_calibration(bins_df, mapping_df)
    print("Saved VITAL score calibration outputs")


if __name__ == "__main__":
    main()
