from __future__ import annotations

from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires matplotlib. Install dependencies or use the bundled workspace Python."
    ) from exc


BASE_DIR = Path(__file__).resolve().parents[1]
FIGURE_PATH = BASE_DIR / "figures" / "vital_gray_zone_workflow.png"


def add_box(
    ax,
    xy: tuple[float, float],
    width: float,
    height: float,
    text: str,
    fill: str,
    edge: str = "#2b3440",
    fontsize: int = 8,
) -> None:
    box = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        linewidth=1.1,
        edgecolor=edge,
        facecolor=fill,
    )
    ax.add_patch(box)
    ax.text(
        xy[0] + width / 2,
        xy[1] + height / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        color="#111827",
        wrap=True,
    )


def add_arrow(ax, start: tuple[float, float], end: tuple[float, float]) -> None:
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=10,
        linewidth=1.0,
        color="#374151",
        shrinkA=4,
        shrinkB=4,
    )
    ax.add_patch(arrow)


def main() -> None:
    fig, ax = plt.subplots(figsize=(12.2, 7.4))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    add_box(
        ax,
        (0.36, 0.83),
        0.28,
        0.09,
        "Gray no-frequency-evidence variant\n(no exact usable AF; not AF=0)",
        "#e5e7eb",
        fontsize=9,
    )

    add_box(
        ax,
        (0.08, 0.64),
        0.25,
        0.12,
        "SNV / MNV path\nrepresentation likely simpler",
        "#d8eef0",
        fontsize=9,
    )
    add_box(
        ax,
        (0.38, 0.64),
        0.24,
        0.12,
        "Indel / duplication /\ncomplex allele path",
        "#f2e6c9",
        fontsize=9,
    )
    add_box(
        ax,
        (0.69, 0.64),
        0.24,
        0.12,
        "Exact match without AF\nor query ambiguity",
        "#e2d9f3",
        fontsize=9,
    )

    add_box(
        ax,
        (0.05, 0.42),
        0.29,
        0.14,
        "Prioritize now if:\nLOF/splice OR critical domain\nAND dominant/high-actionability gene\nOR LOEUF < 0.5 if available",
        "#bcdde3",
        fontsize=8,
    )
    add_box(
        ax,
        (0.37, 0.42),
        0.27,
        0.14,
        "Representation QC first:\nleft-normalize REF/ALT, check transcript,\nrepeat/GC context, breakpoint plausibility",
        "#ead39d",
        fontsize=8,
    )
    add_box(
        ax,
        (0.68, 0.42),
        0.27,
        0.14,
        "Re-query / cache audit:\nconfirm current gnomAD status,\ncheck genome/exome blocks,\nrecord API/cache status",
        "#d7c9ee",
        fontsize=8,
    )

    add_box(
        ax,
        (0.07, 0.21),
        0.26,
        0.13,
        "Clinical evidence review:\nphenotype match, inheritance,\nsegregation, functional evidence,\nClinVar review strength",
        "#c9dfc6",
        fontsize=8,
    )
    add_box(
        ax,
        (0.38, 0.21),
        0.25,
        0.13,
        "Orthogonal validation if\npatient-facing decision:\nPCR/Sanger, MLPA/array,\nlong-read or lab evidence review",
        "#f0c7b9",
        fontsize=8,
    )
    add_box(
        ax,
        (0.69, 0.21),
        0.24,
        0.13,
        "Deferred gray queue:\nlow-impact or well-supported records\nwithout immediate clinical consequence",
        "#d6d9de",
        fontsize=8,
    )

    add_box(
        ax,
        (0.28, 0.055),
        0.44,
        0.10,
        "Queue outcome: expedited expert review, batch review,\northogonal validation, or scheduled re-query at next release",
        "#dbeafe",
        fontsize=9,
    )

    add_arrow(ax, (0.45, 0.83), (0.20, 0.76))
    add_arrow(ax, (0.50, 0.83), (0.50, 0.76))
    add_arrow(ax, (0.56, 0.83), (0.81, 0.76))
    add_arrow(ax, (0.20, 0.64), (0.20, 0.56))
    add_arrow(ax, (0.50, 0.64), (0.50, 0.56))
    add_arrow(ax, (0.81, 0.64), (0.81, 0.56))
    add_arrow(ax, (0.20, 0.42), (0.20, 0.34))
    add_arrow(ax, (0.50, 0.42), (0.50, 0.34))
    add_arrow(ax, (0.81, 0.42), (0.81, 0.34))
    add_arrow(ax, (0.20, 0.21), (0.41, 0.16))
    add_arrow(ax, (0.50, 0.21), (0.50, 0.16))
    add_arrow(ax, (0.81, 0.21), (0.61, 0.16))

    ax.text(
        0.5,
        0.985,
        "VITAL gray-zone workflow: absence is a managed evidence state, not rarity",
        ha="center",
        va="top",
        fontsize=13,
        fontweight="bold",
        color="#111827",
    )
    ax.text(
        0.5,
        0.015,
        "Thresholds such as LOEUF < 0.5 are operational examples when external constraint is available; final priority remains phenotype- and inheritance-aware.",
        ha="center",
        va="bottom",
        fontsize=7.5,
        color="#4b5563",
    )

    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_PATH, dpi=220, bbox_inches="tight", pad_inches=0.16)
    plt.close(fig)
    print(f"Saved {FIGURE_PATH}")


if __name__ == "__main__":
    main()
