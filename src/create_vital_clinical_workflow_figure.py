from __future__ import annotations

from pathlib import Path
from textwrap import fill

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


BASE_DIR = Path(__file__).resolve().parents[1]
FIGURE_PATH = BASE_DIR / "figures" / "vital_clinical_workflow.png"


def add_box(ax, x, y, w, h, title, body, face, edge="#263238", title_color="#101820"):
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.018,rounding_size=0.025",
        linewidth=1.4,
        edgecolor=edge,
        facecolor=face,
    )
    ax.add_patch(box)
    ax.text(
        x + 0.04,
        y + h - 0.12,
        title,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
        color=title_color,
    )
    ax.text(
        x + 0.04,
        y + h - 0.30,
        fill(body, 28),
        ha="left",
        va="top",
        fontsize=10.2,
        color="#263238",
        linespacing=1.2,
    )


def add_arrow(ax, start, end, color="#455a64"):
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=18,
        linewidth=1.8,
        color=color,
        shrinkA=4,
        shrinkB=4,
    )
    ax.add_patch(arrow)


def main() -> None:
    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(14, 7.2))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 7.2)
    ax.axis("off")

    ax.text(
        0.35,
        6.82,
        "VITAL clinical reclassification risk prioritization workflow",
        fontsize=18,
        fontweight="bold",
        ha="left",
        va="top",
        color="#101820",
    )
    ax.text(
        0.35,
        6.48,
        "Designed to reduce false-positive variant re-evaluation burden while preserving human clinical decision-making.",
        fontsize=11.5,
        ha="left",
        va="top",
        color="#455a64",
    )

    y = 4.15
    w = 2.45
    h = 1.75
    add_box(
        ax,
        0.35,
        y,
        w,
        h,
        "ClinVar P/LP",
        "1,731 arrhythmia assertions carried forward, including no-frequency-evidence variants.",
        "#e8f1fb",
    )
    add_box(
        ax,
        3.15,
        y,
        w,
        h,
        "Naive AF screen",
        "115 variants flagged by popmax/global AF >1e-5.",
        "#fff3d6",
    )
    add_box(
        ax,
        5.95,
        y,
        w,
        h,
        "VITAL evidence filter",
        "AF, popmax, AC, variant type, gene constraint, detectability, and review fragility.",
        "#e8f6ef",
    )
    add_box(
        ax,
        8.75,
        y,
        w,
        h,
        "Review queue",
        "3 red-priority variants; 112 naive alerts withheld from urgent review.",
        "#fde7e3",
    )
    add_box(
        ax,
        11.55,
        y,
        2.1,
        h,
        "Clinical review",
        "Expert ACMG/AMP review determines retain P/LP, VUS, B/LB, or validation need.",
        "#f1edf8",
    )

    for sx in [2.8, 5.6, 8.4, 11.2]:
        add_arrow(ax, (sx, y + h / 2), (sx + 0.32, y + h / 2))

    add_box(
        ax,
        0.35,
        1.74,
        3.95,
        1.55,
        "Burden reduction",
        "97.4% fewer urgent action flags than naive AF screening. Operational benchmark: 53 false positives reduced to 0.",
        "#eef3f7",
    )
    add_box(
        ax,
        4.75,
        1.74,
        4.1,
        1.55,
        "Workflow time saved",
        "At 30-60 min per first-pass review, suppressing 112 naive alerts avoids about 56-112 reviewer-hours per audit.",
        "#eef3f7",
    )
    add_box(
        ax,
        9.3,
        1.74,
        4.35,
        1.55,
        "Human decision-making preserved",
        "VITAL prioritizes records for review; it does not assign benignity or replace expert interpretation.",
        "#eef3f7",
    )

    ax.text(
        0.35,
        0.82,
        "Gray no-frequency-evidence variants remain visible and are not treated as AF=0.",
        fontsize=11,
        ha="left",
        va="center",
        color="#455a64",
    )
    ax.text(
        13.65,
        0.82,
        "Output: auditable review-priority queue",
        fontsize=11,
        ha="right",
        va="center",
        color="#455a64",
    )

    fig.tight_layout(pad=0.4)
    fig.savefig(FIGURE_PATH, dpi=220)
    plt.close(fig)
    print(f"Saved {FIGURE_PATH}")


if __name__ == "__main__":
    main()
