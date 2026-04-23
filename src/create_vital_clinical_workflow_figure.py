from __future__ import annotations

from pathlib import Path
from textwrap import fill

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


BASE_DIR = Path(__file__).resolve().parents[1]
FIGURE_PATH = BASE_DIR / "figures" / "vital_clinical_workflow.png"


def add_box(
    ax,
    x,
    y,
    w,
    h,
    title,
    body,
    face,
    edge="#263238",
    title_color="#101820",
    title_size=12.5,
    body_size=9.4,
    wrap=34,
):
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
        fontsize=title_size,
        fontweight="bold",
        color=title_color,
    )
    ax.text(
        x + 0.04,
        y + h - 0.30,
        fill(body, wrap),
        ha="left",
        va="top",
        fontsize=body_size,
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
    fig, ax = plt.subplots(figsize=(15.5, 9.2))
    ax.set_xlim(0, 15.5)
    ax.set_ylim(0, 9.2)
    ax.axis("off")

    ax.text(
        0.4,
        8.86,
        "Laboratory decision tree for public-label review",
        fontsize=19,
        fontweight="bold",
        ha="left",
        va="top",
        color="#101820",
    )
    ax.text(
        0.4,
        8.50,
        "Batch decision support for clinical laboratories with bioinformatics support; final classification remains expert-led.",
        fontsize=11.2,
        ha="left",
        va="top",
        color="#455a64",
    )

    add_box(
        ax,
        0.4,
        6.78,
        2.35,
        1.18,
        "1. Batch score",
        "ClinVar P/LP records plus cached gnomAD AF, popmax, AC, variant type, detectability, gene context, and review status.",
        "#e8f1fb",
        wrap=29,
    )
    add_box(
        ax,
        3.25,
        6.78,
        2.35,
        1.18,
        "2. Assign band",
        "Red-priority, orange high-tension, yellow watchlist, gray no-frequency evidence, or routine green/blue.",
        "#fff3d6",
        wrap=30,
    )
    add_box(
        ax,
        6.1,
        6.78,
        2.55,
        1.18,
        "3. Route action",
        "Route by band plus inheritance mechanism, phenotype relevance, and local laboratory policy.",
        "#e8f6ef",
        wrap=31,
    )
    add_box(
        ax,
        9.25,
        6.78,
        2.55,
        1.18,
        "4. Expert review",
        "Variant scientist, laboratory director, and disease expert verify ClinVar, gnomAD, literature, phenotype, and mechanism.",
        "#fde7e3",
        wrap=31,
    )
    add_box(
        ax,
        12.35,
        6.78,
        2.75,
        1.18,
        "5. Record outcome",
        "Document retain P/LP, VUS, LB/B, context-specific assertion, validation need, or deferred monitoring.",
        "#f1edf8",
        wrap=32,
    )

    for sx in [2.8, 5.65, 8.7, 11.85]:
        add_arrow(ax, (sx, 7.37), (sx + 0.36, 7.37))

    add_box(
        ax,
        0.4,
        4.72,
        3.4,
        1.38,
        "Red-priority",
        "Urgent re-review, ideally within the next variant-review cycle or 30 days for patient-facing assertions. Confirm live ClinVar/gnomAD, phenotype, inheritance, literature, and local evidence.",
        "#f8d7da",
        wrap=39,
    )
    add_box(
        ax,
        4.15,
        4.72,
        3.4,
        1.38,
        "Orange high-tension",
        "Scheduled batch review. Escalate if AC-supported, high-actionability gene, weak review, or phenotype/literature mismatch is present.",
        "#ffe5cc",
        wrap=39,
    )
    add_box(
        ax,
        7.9,
        4.72,
        3.4,
        1.38,
        "Yellow watchlist",
        "No urgent action. Re-query at routine ClinVar/gnomAD release updates and review in aggregate if the band, AC, or review status changes.",
        "#fff4bf",
        wrap=39,
    )
    add_box(
        ax,
        11.65,
        4.72,
        3.45,
        1.38,
        "Gray no-frequency evidence",
        "Do not treat as AF=0. Perform representation-aware checks; use orthogonal validation if an indel, duplication, or complex allele could affect a patient-facing decision.",
        "#e6e6e6",
        wrap=40,
    )

    add_box(
        ax,
        0.4,
        2.78,
        4.55,
        1.12,
        "Who reviews red-priority calls?",
        "Primary variant scientist prepares evidence summary; laboratory director signs classification action; disease expert or panel resolves mechanism-sensitive cases.",
        "#eef3f7",
        wrap=49,
    )
    add_box(
        ax,
        5.45,
        2.78,
        4.55,
        1.12,
        "What is documented?",
        "Component breakdown, live database check date, evidence reviewed, ACMG/ClinGen criteria affected, final disposition, and next recheck date.",
        "#eef3f7",
        wrap=49,
    )
    add_box(
        ax,
        10.55,
        2.78,
        4.55,
        1.12,
        "What is not automated?",
        "Pathogenicity, benignity, penetrance, phase, phenotype match, segregation, and functional interpretation remain human-review decisions.",
        "#eef3f7",
        wrap=49,
    )

    ax.text(
        0.4,
        1.52,
        "Operational examples from the arrhythmia audit: 115 naive AF alerts -> 3 red-priority calls; 112 urgent triggers withheld from immediate review.",
        fontsize=10.4,
        ha="left",
        va="center",
        color="#455a64",
    )
    ax.text(0.4, 1.18, "Output: auditable queue, not automatic reclassification.", fontsize=10.4, ha="left", va="center", color="#455a64")

    fig.tight_layout(pad=0.4)
    fig.savefig(FIGURE_PATH, dpi=220)
    plt.close(fig)
    print(f"Saved {FIGURE_PATH}")


if __name__ == "__main__":
    main()
