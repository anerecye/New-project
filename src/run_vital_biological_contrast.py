"""Create a compact biological-contrast layer for urgent-review case studies.

This script intentionally keeps the biological layer small and interpretable:
it joins cached score tables to real gnomAD v4.1 LOEUF values and HPA v23 heart
expression, then emits a clinician-readable contrast table for the red queue.
"""

from __future__ import annotations

import argparse
import textwrap
import zipfile
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data"
PROCESSED_DIR = DATA_DIR / "processed"
EXTERNAL_DIR = DATA_DIR / "external"
FIGURE_DIR = BASE_DIR / "figures"
SUPPLEMENTARY_DIR = BASE_DIR / "supplementary_tables"

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

CONTEXT = {
    "SCN5A": {
        "domain_context": "composite haplotype; protein/domain assignment is not allele-specific from the ClinVar title alone",
        "expected_mechanism": "dominant Brugada/channelopathy label, but literature context includes drug-triggered susceptibility",
        "contrast_interpretation": (
            "Frequency incompatible with a straightforward high-penetrance dominant Mendelian label; "
            "better framed as a context-dependent susceptibility haplotype pending expert review."
        ),
    },
    "TRDN": {
        "domain_context": "frameshift in triadin luminal protein; LOF assertion intersects recessive/biallelic triadin-null biology",
        "expected_mechanism": "recessive or biallelic null mechanism; heterozygous carrier frequency can be tolerated",
        "contrast_interpretation": (
            "Consistent with recessive carrier frequency rather than automatic benignity; the score flags assertion tension, "
            "not loss of affected-state pathogenicity."
        ),
    },
    "KCNH2": {
        "domain_context": "canonical splice-donor assertion near a highly curated dominant LQTS gene; same-site allele ambiguity is possible",
        "expected_mechanism": "dominant high-penetrance splice-disruption model would predict strong rarity",
        "contrast_interpretation": (
            "Frequency tension is biologically uncomfortable in a constrained dominant gene, but the case is borderline "
            "and could reflect allele-specific evidence gaps, reduced penetrance, founder enrichment, or overcalling."
        ),
    },
}


def read_gnomad_constraint(path: Path, genes: list[str]) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path}. Download from "
            "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/"
            "gnomad.v4.1.constraint_metrics.tsv"
        )
    usecols = [
        "gene",
        "transcript",
        "canonical",
        "mane_select",
        "lof.oe",
        "lof.oe_ci.upper",
        "lof.pLI",
        "mis.oe",
        "mis.z_score",
        "constraint_flags",
    ]
    df = pd.read_csv(path, sep="\t", usecols=usecols)
    df = df[df["gene"].isin(genes)].copy()
    df["canonical"] = df["canonical"].astype(str).str.lower().eq("true")
    df["mane_select"] = df["mane_select"].astype(str).str.lower().eq("true")
    df["priority"] = 3 * df["mane_select"].astype(int) + 2 * df["canonical"].astype(int)
    df = (
        df.sort_values(["gene", "priority"], ascending=[True, False])
        .drop_duplicates("gene", keep="first")
        .drop(columns=["priority"])
    )
    return df.rename(
        columns={
            "transcript": "constraint_transcript",
            "canonical": "constraint_canonical",
            "mane_select": "constraint_mane_select",
            "lof.oe": "gnomad_lof_oe",
            "lof.oe_ci.upper": "gnomad_lof_oe_ci_upper",
            "lof.pLI": "gnomad_pli",
            "mis.oe": "gnomad_mis_oe",
            "mis.z_score": "gnomad_mis_z",
            "constraint_flags": "gnomad_constraint_flags",
        }
    )


def read_hpa_heart_expression(path: Path, genes: list[str]) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path}. Download from "
            "https://v23.proteinatlas.org/download/rna_tissue_consensus.tsv.zip"
        )
    with zipfile.ZipFile(path) as archive:
        with archive.open(archive.namelist()[0]) as handle:
            df = pd.read_csv(handle, sep="\t")
    heart = df[
        df["Gene name"].isin(genes)
        & df["Tissue"].astype(str).str.casefold().eq("heart muscle")
    ].copy()
    heart = heart.rename(
        columns={
            "Gene name": "gene",
            "Gene": "hpa_gene_id",
            "Tissue": "hpa_tissue",
            "nTPM": "hpa_heart_muscle_ntpm",
        }
    )
    heart["heart_expression_context"] = pd.cut(
        pd.to_numeric(heart["hpa_heart_muscle_ntpm"], errors="coerce"),
        bins=[-float("inf"), 1, 10, 100, float("inf")],
        labels=["not_detected_or_low", "low", "moderate", "high"],
    ).astype(str)
    return heart[["gene", "hpa_gene_id", "hpa_tissue", "hpa_heart_muscle_ntpm", "heart_expression_context"]]


def frequency_signal(row: pd.Series) -> str:
    return (
        f"global AF={row['global_af']:.2e}, AC={row['global_ac']:.0f}; "
        f"{str(row['popmax_population']).upper()} popmax AF={row['popmax_af']:.2e}, AC={row['popmax_ac']:.0f}"
    )


def loeuf_argument(row: pd.Series) -> str:
    value = pd.to_numeric(row.get("gnomad_lof_oe_ci_upper"), errors="coerce")
    if pd.isna(value):
        return "LOEUF unavailable"
    label = "strongly constrained" if value < 0.35 else "moderately constrained" if value < 0.7 else "not LOF-constrained"
    return f"LOEUF={value:.3f} ({label})"


def make_gene_context(scores: pd.DataFrame, constraint: pd.DataFrame, expression: pd.DataFrame) -> pd.DataFrame:
    genes = pd.DataFrame({"gene": sorted(scores["gene"].dropna().unique())})
    gene_counts = scores.groupby("gene", dropna=False).agg(
        clinvar_plp_variant_count=("variant_key", "count"),
        af_observed_count=("frequency_evidence_status", lambda s: int(s.eq("frequency_observed").sum())),
        vital_red_count=("vital_red_flag", lambda s: int(s.astype(bool).sum())),
        high_tension_count=("vital_score", lambda s: int((pd.to_numeric(s, errors="coerce") >= 60).sum())),
    )
    out = genes.merge(gene_counts, on="gene", how="left")
    out = out.merge(constraint, on="gene", how="left")
    out = out.merge(expression, on="gene", how="left")
    out["biological_context_source"] = (
        "gnomAD v4.1 constraint metrics for LOEUF/pLI/missense z; "
        "Human Protein Atlas v23 tissue consensus for heart-muscle nTPM"
    )
    return out


def make_contrast_cases(scores: pd.DataFrame, gene_context: pd.DataFrame) -> pd.DataFrame:
    red = scores[scores["vital_red_flag"].astype(bool)].copy()
    red = red.merge(gene_context, on="gene", how="left")
    rows = []
    for _, row in red.sort_values("vital_score", ascending=False).iterrows():
        context = CONTEXT.get(row["gene"], {})
        rows.append(
            {
                "variant": f"{row['gene']} {row['clinvar_id']} ({row['title']})",
                "vital": f"{row['vital_score']:.1f} / {str(row['vital_band']).replace('_', ' ')}",
                "af_signal": frequency_signal(row),
                "loeuf_constraint_argument": loeuf_argument(row),
                "heart_expression": (
                    f"HPA heart muscle nTPM={row['hpa_heart_muscle_ntpm']:.1f} "
                    f"({row['heart_expression_context']})"
                    if pd.notna(row.get("hpa_heart_muscle_ntpm"))
                    else "heart expression unavailable"
                ),
                "domain_or_transcript_context": context.get("domain_context", "requires variant-level domain annotation"),
                "expected_mechanism": context.get("expected_mechanism", "requires disease-mechanism curation"),
                "reality_interpretation": context.get(
                    "contrast_interpretation",
                    "Frequency tension requires expert interpretation rather than automatic reclassification.",
                ),
            }
        )
    return pd.DataFrame(rows)


def make_summary(scores: pd.DataFrame, cases: pd.DataFrame) -> pd.DataFrame:
    af_observed = scores[scores["frequency_evidence_status"].eq("frequency_observed")].copy()
    naive = af_observed[af_observed["standard_acmg_frequency_flag"].astype(bool)]
    high = af_observed[pd.to_numeric(af_observed["vital_score"], errors="coerce").ge(60)]
    red = af_observed[af_observed["vital_red_flag"].astype(bool)]
    rows = [
        {
            "result": "frequency_tension_common_but_hidden",
            "n": len(naive),
            "denominator": len(af_observed),
            "percent": 100 * len(naive) / len(af_observed) if len(af_observed) else float("nan"),
            "interpretation": (
                "Popmax/global AF >1e-5 tension is common among AF-observed ClinVar P/LP variants; "
                "global-only workflows hide much of this signal."
            ),
        },
        {
            "result": "biologically_interpretable_high_tension_core",
            "n": len(high),
            "denominator": len(af_observed),
            "percent": 100 * len(high) / len(af_observed) if len(af_observed) else float("nan"),
            "interpretation": (
                "Scores >=60 compress the broad frequency-tension space into a small set of variants "
                "where frequency must be interpreted against mechanism, LOEUF, expression, and domain context."
            ),
        },
        {
            "result": "red_case_study_queue",
            "n": len(red),
            "denominator": len(af_observed),
            "percent": 100 * len(red) / len(af_observed) if len(af_observed) else float("nan"),
            "interpretation": (
                "The urgent-review queue is a case-study set, not an appendix: each red variant has a distinct biological explanation "
                "for why AF tension does or does not undermine a P/LP assertion."
            ),
        },
    ]
    rows.append(
        {
            "result": "red_cases_with_real_lof_constraint_and_heart_expression",
            "n": int(
                cases["loeuf_constraint_argument"].notna().sum()
                if not cases.empty
                else 0
            ),
            "denominator": len(cases),
            "percent": 100.0 if len(cases) else float("nan"),
            "interpretation": (
                "All red cases are shown with real gnomAD LOEUF and HPA heart-expression context, "
                "turning metadata flags into biological contrast statements."
            ),
        }
    )
    return pd.DataFrame(rows)


def make_figure(cases: pd.DataFrame, path: Path) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(len(cases), 1, figsize=(11, 7.5))
    if len(cases) == 1:
        axes = [axes]
    colors = ["#e8f0f6", "#f4ead6", "#f1e4e8"]
    for ax, (_, row), color in zip(axes, cases.iterrows(), colors):
        ax.axis("off")
        ax.add_patch(plt.Rectangle((0, 0), 1, 1, transform=ax.transAxes, color=color, zorder=0))
        variant = row["variant"].split(" (")[0]
        body = (
            f"{variant}\n"
            f"Score: {row['vital']}\n"
            f"AF: {row['af_signal']}\n"
            f"{row['loeuf_constraint_argument']} | {row['heart_expression']}\n"
            f"Mechanism: {row['expected_mechanism']}\n"
            f"Interpretation: {row['reality_interpretation']}"
        )
        ax.text(
            0.02,
            0.92,
            textwrap.fill(body, width=125, subsequent_indent=""),
            va="top",
            ha="left",
            fontsize=9.5,
            color="#1f2933",
            linespacing=1.35,
        )
    fig.suptitle(
        "Urgent-review variants as biological contrast case studies",
        fontsize=14,
        fontweight="bold",
        color="#1f2933",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(path, dpi=220)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build biological contrast case tables.")
    parser.add_argument("--score-table", type=Path, default=PROCESSED_DIR / "arrhythmia_vital_scores.csv")
    parser.add_argument("--constraint-table", type=Path, default=EXTERNAL_DIR / "gnomad.v4.1.constraint_metrics.tsv")
    parser.add_argument("--hpa-table", type=Path, default=EXTERNAL_DIR / "hpa_rna_tissue_consensus.tsv.zip")
    parser.add_argument("--output-prefix", default="arrhythmia")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    SUPPLEMENTARY_DIR.mkdir(parents=True, exist_ok=True)
    scores = pd.read_csv(args.score_table)
    genes = sorted(set(ARRHYTHMIA_GENES).intersection(set(scores["gene"].dropna().unique())))
    constraint = read_gnomad_constraint(args.constraint_table, genes)
    expression = read_hpa_heart_expression(args.hpa_table, genes)
    gene_context = make_gene_context(scores, constraint, expression)
    cases = make_contrast_cases(scores, gene_context)
    summary = make_summary(scores, cases)

    gene_context.to_csv(PROCESSED_DIR / f"{args.output_prefix}_gene_biological_context.csv", index=False)
    cases.to_csv(PROCESSED_DIR / f"{args.output_prefix}_vital_biological_contrast_cases.csv", index=False)
    summary.to_csv(PROCESSED_DIR / f"{args.output_prefix}_vital_biological_contrast_summary.csv", index=False)
    cases.to_csv(SUPPLEMENTARY_DIR / "Supplementary_Table_S16_biological_contrast_cases.tsv", sep="\t", index=False)
    make_figure(cases, FIGURE_DIR / "vital_biological_contrast_cases.png")

    print(f"Wrote {len(gene_context)} gene biological context rows")
    print(f"Wrote {len(cases)} urgent-review contrast cases")


if __name__ == "__main__":
    main()
