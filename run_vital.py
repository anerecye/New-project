from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
DEFAULT_SCORES = ROOT / "data" / "processed" / "arrhythmia_vital_scores.csv"

COMPONENTS = [
    ("AF pressure", "frequency_pressure_score", 45),
    ("AC reliability", "ac_reliability_score", 20),
    ("Popmax enrichment", "popmax_enrichment_score", 10),
    ("Variant type tension", "variant_type_tension_score", 6),
    ("Technical detectability", "technical_detectability_score", 8),
    ("Gene constraint", "gene_constraint_score", 10),
    ("Review fragility", "review_fragility_score", 10),
]

BAND_LABELS = {
    "red_reclassification_priority": "RED (urgent review)",
    "orange_high_tension": "ORANGE (high tension)",
    "yellow_watchlist": "YELLOW (watchlist)",
    "blue_low_support_signal": "BLUE (low-support signal)",
    "green_frequency_consistent": "GREEN (frequency-consistent)",
    "gray_no_frequency_evidence": "GRAY (no frequency evidence)",
}


def truthy(value: object) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def fmt_num(value: object, digits: int = 3) -> str:
    if value is None or pd.isna(value):
        return "NA"
    value = float(value)
    if digits == 0:
        return f"{value:,.0f}"
    if value == 0:
        return "0"
    if abs(value) < 0.001 or abs(value) >= 1000:
        return f"{value:.{digits}e}"
    return f"{value:.{digits}g}"


def fmt_score(value: object) -> str:
    if value is None or pd.isna(value):
        return "NA"
    return f"{float(value):.1f}"


def normalize_vcv(value: object) -> str:
    raw = str(value).strip()
    if not raw:
        return raw
    upper = raw.upper()
    if upper.startswith("VCV"):
        digits = "".join(ch for ch in upper[3:] if ch.isdigit())
        return f"VCV{digits.zfill(9)}" if digits else upper
    if raw.isdigit():
        return f"VCV{raw.zfill(9)}"
    return upper


def load_scores(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise SystemExit(f"Score table not found: {path}")
    return pd.read_csv(
        path,
        dtype={
            "clinvar_id": "string",
            "variation_id": "string",
            "gene": "string",
            "title": "string",
        },
    )


def find_variant(scores: pd.DataFrame, query: object) -> pd.Series | None:
    vcv = normalize_vcv(query)
    digits = vcv[3:].lstrip("0") if vcv.startswith("VCV") else str(query).strip()

    mask = scores["clinvar_id"].astype("string").str.upper().eq(vcv)
    if not mask.any() and "variation_id" in scores.columns:
        mask = scores["variation_id"].astype("string").str.replace(r"\.0$", "", regex=True).eq(digits)
    if not mask.any():
        return None
    return scores.loc[mask].iloc[0]


def band_label(row: pd.Series) -> str:
    band = str(row.get("vital_band", "NA"))
    return BAND_LABELS.get(band, band)


def max_af(row: pd.Series) -> float | None:
    values = []
    for col in ("global_af", "popmax_af", "max_frequency_signal"):
        value = row.get(col)
        if value is not None and not pd.isna(value):
            values.append(float(value))
    return max(values) if values else None


def key_signals(row: pd.Series) -> list[str]:
    signals: list[str] = []

    evidence = str(row.get("frequency_evidence_status", ""))
    if evidence and evidence != "frequency_observed":
        signals.append("no usable exact frequency evidence; keep as gray uncertainty, not AF=0")

    max_frequency = max_af(row)
    if max_frequency is not None:
        if max_frequency > 1e-4:
            signals.append(f"high population frequency signal (max AF {fmt_num(max_frequency)})")
        elif max_frequency > 1e-5:
            signals.append(f"population frequency above rare-disease review threshold (max AF {fmt_num(max_frequency)})")

    if truthy(row.get("frequency_signal_ac_ge_20")):
        signals.append(f"AC-supported evidence (qualifying AC {fmt_num(row.get('qualifying_frequency_ac'), 0)})")
    elif row.get("qualifying_frequency_ac") is not None and not pd.isna(row.get("qualifying_frequency_ac")):
        if float(row.get("qualifying_frequency_ac")) > 0:
            signals.append(f"low-AC frequency signal (qualifying AC {fmt_num(row.get('qualifying_frequency_ac'), 0)})")

    ratio = row.get("popmax_global_ratio")
    if ratio is not None and not pd.isna(ratio) and float(ratio) >= 10:
        pop = row.get("popmax_population")
        signals.append(f"population-specific enrichment (popmax/global ratio {fmt_num(ratio)}; popmax={pop})")

    if truthy(row.get("weak_review_signal")):
        signals.append("weak or single-submitter ClinVar review support")

    variant_type = str(row.get("variant_type", "")).lower()
    if variant_type in {"deletion", "insertion", "duplication"}:
        signals.append(f"{variant_type} representation may be technically harder than SNV")

    functional_class = str(row.get("functional_class", "")).lower()
    if functional_class in {"lof", "splice_or_intronic"}:
        signals.append(f"functional annotation: {functional_class}")

    if truthy(row.get("vital_red_flag")):
        signals.append("manual expert re-review recommended; not automatic benign classification")

    return signals or ["no major frequency-tension signal in cached score table"]


def component_lines(row: pd.Series) -> list[str]:
    lines = []
    for label, column, maximum in COMPONENTS:
        if column in row:
            lines.append(f"  {label:<24} {fmt_score(row.get(column)):>5}/{maximum}")
    return lines


def print_single(row: pd.Series, query: str) -> None:
    print("=" * 64)
    print("Cached review summary")
    print("=" * 64)
    print(f"Query:        {query}")
    print(f"ClinVar ID:   {row.get('clinvar_id', 'NA')}")
    print(f"Gene:         {row.get('gene', 'NA')}")
    print(f"Variant:      {row.get('title', 'NA')}")
    print()
    print(f"Score:        {fmt_score(row.get('vital_score'))}")
    print(f"Band:         {band_label(row)}")
    print(f"Red flag:     {truthy(row.get('vital_red_flag'))}")
    print()
    print("Frequency:")
    print(f"  Global AF / AC: {fmt_num(row.get('global_af'))} / {fmt_num(row.get('global_ac'), 0)}")
    print(
        "  Popmax AF / AC: "
        f"{fmt_num(row.get('popmax_af'))} / {fmt_num(row.get('popmax_ac'), 0)}"
        f" ({row.get('popmax_population', 'NA')})"
    )
    print()
    print("Review:")
    print(f"  ClinVar review: {row.get('review_status', 'NA')}")
    print(f"  Review strength: {row.get('review_strength', 'NA')}")
    print()
    print("Key signals:")
    for signal in key_signals(row):
        print(f"- {signal}")
    print()
    print("Component scores:")
    for line in component_lines(row):
        print(line)
    print("=" * 64)


def input_id_column(df: pd.DataFrame) -> str:
    candidates = ["vcv", "VCV", "VCV_ID", "vcv_id", "clinvar_id", "ClinVarID", "variation_id"]
    for col in candidates:
        if col in df.columns:
            return col
    raise SystemExit(f"Input CSV needs one of these columns: {', '.join(candidates)}")


def batch_results(scores: pd.DataFrame, input_csv: Path) -> pd.DataFrame:
    if not input_csv.exists():
        raise SystemExit(f"Input CSV not found: {input_csv}")
    queries = pd.read_csv(input_csv, dtype="string")
    id_col = input_id_column(queries)

    rows = []
    for query in queries[id_col].fillna(""):
        hit = find_variant(scores, query)
        if hit is None:
            rows.append(
                {
                    "query": query,
                    "status": "not_found_in_cached_scores",
                    "gene": "NA",
                    "clinvar_id": normalize_vcv(query),
                    "vital_score": "NA",
                    "band": "NA",
                    "red_flag": False,
                    "max_af": "NA",
                    "qualifying_ac": "NA",
                    "key_signals": "not found in offline demo score table",
                }
            )
            continue
        signals = "; ".join(key_signals(hit)[:4])
        rows.append(
            {
                "query": query,
                "status": "found",
                "gene": hit.get("gene", "NA"),
                "clinvar_id": hit.get("clinvar_id", "NA"),
                "vital_score": fmt_score(hit.get("vital_score")),
                "band": band_label(hit),
                "red_flag": truthy(hit.get("vital_red_flag")),
                "max_af": fmt_num(max_af(hit)),
                "qualifying_ac": fmt_num(hit.get("qualifying_frequency_ac"), 0),
                "key_signals": signals,
            }
        )
    return pd.DataFrame(rows)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fast offline cached-score demo for ClinVar VCV records.",
    )
    parser.add_argument("--vcv", help="ClinVar VCV ID, for example VCV000440850")
    parser.add_argument("--input", type=Path, help="CSV with a vcv/clinvar_id/variation_id column")
    parser.add_argument("--output", type=Path, help="Optional output CSV for batch mode")
    parser.add_argument("--scores", type=Path, default=DEFAULT_SCORES, help="Cached score CSV")
    parser.add_argument("--examples", action="store_true", help="Print demo VCV IDs")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])

    if args.examples:
        print("Try these cached demo variants:")
        print("  VCV000440850  SCN5A  RED anchor")
        print("  VCV001325231  TRDN   RED near-stable")
        print("  VCV004535537  KCNH2  RED borderline")
        return 0

    if not args.vcv and not args.input:
        raise SystemExit("Use --vcv VCV000440850 or --input data/examples/sample_variants.csv")

    scores = load_scores(args.scores)

    if args.vcv:
        hit = find_variant(scores, args.vcv)
        if hit is None:
            print(f"Not found in cached score table: {args.vcv}")
            print("Tip: this offline demo currently uses arrhythmia_vital_scores.csv.")
            return 2
        print_single(hit, args.vcv)

    if args.input:
        results = batch_results(scores, args.input)
        print(results.to_string(index=False))
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            results.to_csv(args.output, index=False)
            print(f"\nSaved: {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
