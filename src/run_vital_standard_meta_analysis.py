from __future__ import annotations

from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires pandas, numpy, and matplotlib. Install them with: "
        "python -m pip install -r requirements.txt"
    ) from exc

from vital_standard import BASE_DIR, DATA_DIR, build_gene_problem_matrix, build_problem_intersections, load_processed_prediction_tables, summarize_standard_view


FIGURE_DIR = BASE_DIR / "figures"

META_SUMMARY_OUT = DATA_DIR / "vital_standard_meta_summary.csv"
META_OVERALL_OUT = DATA_DIR / "vital_standard_meta_overall.csv"
GENE_MATRIX_OUT = DATA_DIR / "vital_standard_gene_problem_matrix.csv"
INTERSECTIONS_OUT = DATA_DIR / "vital_standard_gene_problem_intersections.csv"
UPSET_OUT = FIGURE_DIR / "vital_standard_gene_problem_upset.png"

CATEGORY_COLUMNS = [
    ("evaluation_limited", "Evaluation-limited", "#8b98a5"),
    ("popmax_review", "Popmax review", "#d18f2f"),
    ("model_conflict", "Model conflict", "#b94a48"),
    ("recessive_route", "Recessive route", "#4b6cb7"),
    ("sv_required", "SV/CNV layer", "#7b5ea7"),
]


def save_table(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, lineterminator="\n")
    print(f"Saved {path.relative_to(BASE_DIR)} ({len(df)} rows)")


def compute_overall_rows(summary: pd.DataFrame) -> pd.DataFrame:
    if summary.empty:
        return pd.DataFrame(
            columns=[
                "scope",
                "variant_count",
                "pass_rate_mean",
                "nonpass_rate_mean",
                "nonpass_rate_min",
                "nonpass_rate_max",
                "evaluation_limited_rate_mean",
                "popmax_review_rate_mean",
                "model_conflict_rate_mean",
                "recessive_route_rate_mean",
            ]
        )

    rate_cols = {
        "evaluation_limited_rate_mean": "evaluation_limited_count",
        "popmax_review_rate_mean": "popmax_review_count",
        "model_conflict_rate_mean": "model_conflict_count",
        "recessive_route_rate_mean": "recessive_route_count",
    }

    rows = []
    for scope, frame in {
        "all_meta_datasets": summary,
        "external_disease_panels_only": summary.loc[summary["dataset"].ne("arrhythmia")].copy(),
    }.items():
        if frame.empty:
            continue
        row = {
            "scope": scope,
            "variant_count": int(frame["variant_count"].sum()),
            "pass_rate_mean": float(frame["pass_rate"].mean()),
            "nonpass_rate_mean": float(frame["nonpass_rate"].mean()),
            "nonpass_rate_min": float(frame["nonpass_rate"].min()),
            "nonpass_rate_max": float(frame["nonpass_rate"].max()),
        }
        for out_col, count_col in rate_cols.items():
            row[out_col] = float((frame[count_col] / frame["variant_count"]).mean())
        rows.append(row)
    return pd.DataFrame(rows)


def plot_upset(gene_matrix: pd.DataFrame, intersections: pd.DataFrame, output_path: Path) -> None:
    if gene_matrix.empty or intersections.empty:
        return

    plot_df = intersections.loc[intersections["gene_count"].gt(0)].head(12).copy()
    if plot_df.empty:
        return

    fig = plt.figure(figsize=(12, 7))
    grid = fig.add_gridspec(2, 2, height_ratios=[3.0, 2.4], width_ratios=[1.6, 4.8], hspace=0.06, wspace=0.06)
    ax_blank = fig.add_subplot(grid[0, 0])
    ax_blank.axis("off")
    ax_top = fig.add_subplot(grid[0, 1])
    ax_left = fig.add_subplot(grid[1, 0])
    ax_matrix = fig.add_subplot(grid[1, 1], sharex=ax_top)

    x_positions = np.arange(len(plot_df))
    ax_top.bar(x_positions, plot_df["gene_count"], color="#2f4858", width=0.72)
    ax_top.set_ylabel("Genes")
    ax_top.set_title("Problem intersections across VITAL-flagged genes", fontsize=13, pad=10)
    ax_top.set_xticks([])
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)
    for xpos, value in zip(x_positions, plot_df["gene_count"]):
        ax_top.text(xpos, value + 0.35, str(int(value)), ha="center", va="bottom", fontsize=8)

    category_labels = [label for _, label, _ in CATEGORY_COLUMNS]
    category_counts = [int(gene_matrix[column].sum()) for column, _, _ in CATEGORY_COLUMNS]
    category_colors = [color for _, _, color in CATEGORY_COLUMNS]
    y_positions = np.arange(len(CATEGORY_COLUMNS))
    ax_left.barh(y_positions, category_counts, color=category_colors, alpha=0.9)
    ax_left.set_yticks(y_positions)
    ax_left.set_yticklabels(category_labels, fontsize=9)
    ax_left.invert_yaxis()
    ax_left.set_xlabel("Genes")
    ax_left.spines["top"].set_visible(False)
    ax_left.spines["right"].set_visible(False)

    ax_matrix.set_yticks(y_positions)
    ax_matrix.set_yticklabels(category_labels, fontsize=9)
    ax_matrix.invert_yaxis()
    ax_matrix.set_xlim(-0.6, len(plot_df) - 0.4)
    ax_matrix.set_xticks(x_positions)
    ax_matrix.set_xticklabels([f"I{i + 1}" for i in range(len(plot_df))], fontsize=9)
    ax_matrix.set_xlabel("Top intersections")
    ax_matrix.spines["top"].set_visible(False)
    ax_matrix.spines["right"].set_visible(False)

    for xpos, (_, row) in enumerate(plot_df.iterrows()):
        active_rows: list[int] = []
        for y, (column, _, color) in enumerate(CATEGORY_COLUMNS):
            ax_matrix.scatter(xpos, y, s=55, color="#d9dde3", zorder=2)
            if bool(row[column]):
                ax_matrix.scatter(xpos, y, s=75, color=color, edgecolor="#1f2933", linewidth=0.6, zorder=3)
                active_rows.append(y)
        if len(active_rows) >= 2:
            ax_matrix.plot([xpos, xpos], [min(active_rows), max(active_rows)], color="#64748b", linewidth=1.2, zorder=1)

    note = "I1..In sorted by gene count; categories show evaluation, population, model, recessive-route, and SV/CNV burden."
    fig.text(0.5, 0.01, note, ha="center", va="bottom", fontsize=9, color="#475569")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output_path.relative_to(BASE_DIR)}")


def main() -> None:
    combined = load_processed_prediction_tables(include_auxiliary=False)
    meta = combined.loc[combined["include_in_meta"]].copy()
    summary = summarize_standard_view(meta)
    gene_matrix = build_gene_problem_matrix(meta)
    intersections = build_problem_intersections(gene_matrix)
    overall = compute_overall_rows(summary)

    save_table(summary, META_SUMMARY_OUT)
    save_table(overall, META_OVERALL_OUT)
    save_table(gene_matrix, GENE_MATRIX_OUT)
    save_table(intersections, INTERSECTIONS_OUT)
    plot_upset(gene_matrix, intersections, UPSET_OUT)


if __name__ == "__main__":
    main()
