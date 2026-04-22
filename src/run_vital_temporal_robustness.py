from __future__ import annotations

import argparse
import time
from itertools import combinations
from pathlib import Path

try:
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import requests
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires pandas, numpy, matplotlib, and requests. Install them with: "
        "python -m pip install -r requirements.txt"
    ) from exc

import advanced_variant_analyses as ava
from score_vital_from_variant_summary import load_plp_snapshot
from validate_vital_reclassification import ARRHYTHMIA_GENES


BASE_DIR = Path(__file__).resolve().parents[1]
RAW_DIR = BASE_DIR / "data" / "raw"
DATA_DIR = BASE_DIR / "data" / "processed"
FIGURE_DIR = BASE_DIR / "figures"
SUPPLEMENT_DIR = BASE_DIR / "supplementary_tables"

RAW_SNAPSHOTS = [
    ("2023-01", RAW_DIR / "variant_summary_2023-01.txt.gz"),
    ("2024-01", RAW_DIR / "variant_summary_2024-01.txt.gz"),
]
SNAPSHOT_URLS = {
    "2023-01": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/2023/variant_summary_2023-01.txt.gz",
    "2024-01": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/2024/variant_summary_2024-01.txt.gz",
}
CURRENT_SNAPSHOT_NAME = "2026-04-current"
CURRENT_CANONICAL_SCORES = DATA_DIR / "arrhythmia_vital_scores.csv"

GNOMAD_MATCH_CACHE_FILES = [
    DATA_DIR / "arrhythmia_2023_01_gnomad_matched.csv",
    DATA_DIR / "arrhythmia_gnomad_matched.csv",
    DATA_DIR / "arrhythmia_intermediate_only_2024_2025_gnomad_matched.csv",
]

POPULATION_CACHE_FILES = [
    DATA_DIR / "arrhythmia_2023_01_population_af.csv",
    DATA_DIR / "arrhythmia_population_af.csv",
    DATA_DIR / "arrhythmia_intermediate_only_2024_2025_population_af.csv",
]

AC_THRESHOLDS = [5, 10, 20, 50]
RECESSIVE_CONTEXT_GENES = {"CASQ2", "TRDN"}


def save_csv(df: pd.DataFrame, output_path: Path, sep: str = ",") -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = output_path.with_name(f".{output_path.name}.tmp")
    table = df.copy()
    if sep == "\t":
        object_columns = table.select_dtypes(include=["object", "string"]).columns
        table.loc[:, object_columns] = table.loc[:, object_columns].replace("", "NA")
    last_error: PermissionError | None = None
    for attempt in range(1, 7):
        try:
            table.to_csv(temp_path, sep=sep, index=False, na_rep="NA", lineterminator="\n")
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


def download_snapshot(snapshot_name: str, output_path: Path, retries: int = 20) -> None:
    url = SNAPSHOT_URLS[snapshot_name]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = output_path.with_suffix(output_path.suffix + ".download")
    existing = temp_path.stat().st_size if temp_path.exists() else 0
    headers = {"Range": f"bytes={existing}-"} if existing else {}
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            mode = "ab" if existing else "wb"
            print(f"Downloading {snapshot_name} attempt {attempt}/{retries}: {url}")
            with requests.get(url, stream=True, timeout=60, headers=headers) as response:
                if response.status_code == 416:
                    temp_path.replace(output_path)
                    return
                response.raise_for_status()
                with temp_path.open(mode) as handle:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            handle.write(chunk)
            temp_path.replace(output_path)
            print(f"Downloaded {output_path}")
            return
        except Exception as exc:  # noqa: BLE001 - download retry should keep original error.
            last_error = exc
            delay = min(5 * attempt, 60)
            print(f"Download retry for {snapshot_name} in {delay}s: {exc}")
            time.sleep(delay)
            existing = temp_path.stat().st_size if temp_path.exists() else 0
            headers = {"Range": f"bytes={existing}-"} if existing else {}
    if last_error is not None:
        raise last_error


def ensure_raw_snapshots(download_missing: bool) -> None:
    missing = [(snapshot, path) for snapshot, path in RAW_SNAPSHOTS if not path.exists()]
    if not missing:
        return
    if not download_missing:
        missing_paths = ", ".join(str(path) for _, path in missing)
        raise FileNotFoundError(f"Missing ClinVar snapshots: {missing_paths}")
    for snapshot, path in missing:
        download_snapshot(snapshot, path)


def truthy(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def normalize_variation_id(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.endswith(".0"):
        text = text[:-2]
    return text


def load_cache_table(paths: list[Path]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for path in paths:
        if not path.exists():
            print(f"Skipping missing cache: {path}")
            continue
        frame = pd.read_csv(path, low_memory=False)
        if "variant_key" not in frame.columns:
            continue
        frame["cache_source_file"] = path.name
        frames.append(frame)
    if not frames:
        raise FileNotFoundError("No usable cached gnomAD tables were found.")
    return pd.concat(frames, ignore_index=True)


def build_gnomad_match_cache() -> pd.DataFrame:
    cache = load_cache_table(GNOMAD_MATCH_CACHE_FILES)
    for column in ["gnomad_af", "gnomad_ac", "gnomad_an"]:
        cache[column] = pd.to_numeric(cache.get(column), errors="coerce")
    cache["match_category"] = cache.get("match_category", "").fillna("").astype(str)
    cache["cache_rank"] = np.select(
        [
            cache["match_category"].eq("exact_match") & cache["gnomad_af"].notna(),
            cache["match_category"].eq("exact_match"),
            cache["match_category"].eq("allele_discordance"),
            cache["match_category"].eq("no_gnomad_record"),
        ],
        [0, 1, 2, 3],
        default=4,
    )
    keep = [
        "variant_key",
        "match_category",
        "gnomad_af",
        "gnomad_ac",
        "gnomad_an",
        "gnomad_status",
        "cache_source_file",
        "cache_rank",
    ]
    cache = cache.loc[:, [column for column in keep if column in cache.columns]].copy()
    cache = (
        cache.sort_values(["variant_key", "cache_rank"])
        .drop_duplicates("variant_key", keep="first")
        .drop(columns=["cache_rank"])
        .reset_index(drop=True)
    )
    cache = cache.rename(
        columns={
            "match_category": "cached_match_category",
            "gnomad_af": "cached_gnomad_af",
            "gnomad_ac": "cached_gnomad_ac",
            "gnomad_an": "cached_gnomad_an",
            "gnomad_status": "cached_gnomad_status",
        }
    )
    return cache


def build_population_cache() -> pd.DataFrame:
    cache = load_cache_table(POPULATION_CACHE_FILES)
    for column in ["global_af", "global_ac", "global_an", "popmax_af", "popmax_ac"]:
        if column in cache.columns:
            cache[column] = pd.to_numeric(cache[column], errors="coerce")
    cache["population_cache_rank"] = np.where(cache.get("global_af").notna(), 0, 1)
    cache = (
        cache.sort_values(["variant_key", "population_cache_rank"])
        .drop_duplicates("variant_key", keep="first")
        .drop(columns=["population_cache_rank"])
        .reset_index(drop=True)
    )
    return cache


def annotate_gnomad(variants_df: pd.DataFrame, match_cache: pd.DataFrame) -> pd.DataFrame:
    annotated = ava.prepare_annotations(variants_df)
    annotated = annotated.merge(match_cache, on="variant_key", how="left")
    annotated["gnomad_cache_status"] = np.where(
        annotated["cached_match_category"].notna(),
        "from_local_cache",
        "not_in_local_cache_no_frequency_evidence",
    )
    annotated["match_category"] = annotated["cached_match_category"].fillna("gnomad_query_error")
    annotated["gnomad_af"] = pd.to_numeric(annotated["cached_gnomad_af"], errors="coerce")
    annotated["gnomad_ac"] = pd.to_numeric(annotated["cached_gnomad_ac"], errors="coerce")
    annotated["gnomad_an"] = pd.to_numeric(annotated["cached_gnomad_an"], errors="coerce")
    annotated["gnomad_status"] = annotated.get("cached_gnomad_status", "").fillna(
        annotated["gnomad_cache_status"]
    )
    annotated = annotated.drop(
        columns=[
            column
            for column in [
                "cached_match_category",
                "cached_gnomad_af",
                "cached_gnomad_ac",
                "cached_gnomad_an",
                "cached_gnomad_status",
            ]
            if column in annotated.columns
        ]
    )
    return ava.prepare_annotations(annotated)


def compute_vital(
    snapshot_df: pd.DataFrame,
    match_cache: pd.DataFrame,
    population_cache: pd.DataFrame,
) -> pd.DataFrame:
    annotated = annotate_gnomad(snapshot_df, match_cache)
    population_df = population_cache[
        population_cache["variant_key"].isin(annotated["variant_key"])
    ].copy()
    scores, _, _, _, _ = ava.make_vital_score_tables(annotated, population_df)
    scores["inheritance_flag"] = np.where(
        scores["gene"].isin(RECESSIVE_CONTEXT_GENES),
        "recessive_context_required",
        "standard",
    )
    scores["variation_id"] = scores["variation_id"].map(normalize_variation_id)
    scores["clinvar_id"] = scores["clinvar_id"].fillna("")
    scores["snapshot"] = snapshot_df.attrs.get("snapshot", "")
    return add_red_threshold_columns(scores)


def load_current_canonical_scores() -> pd.DataFrame:
    if not CURRENT_CANONICAL_SCORES.exists():
        raise FileNotFoundError(f"Missing current canonical VITAL scores: {CURRENT_CANONICAL_SCORES}")
    scores = pd.read_csv(CURRENT_CANONICAL_SCORES, low_memory=False)
    scores["snapshot"] = CURRENT_SNAPSHOT_NAME
    scores["variation_id"] = scores["variation_id"].map(normalize_variation_id)
    scores["inheritance_flag"] = np.where(
        scores["gene"].isin(RECESSIVE_CONTEXT_GENES),
        "recessive_context_required",
        "standard",
    )
    scores["gnomad_cache_status"] = "canonical_current_score_table"
    return add_red_threshold_columns(scores)


def add_red_threshold_columns(scores: pd.DataFrame) -> pd.DataFrame:
    result = scores.copy()
    result["vital_score"] = pd.to_numeric(result.get("vital_score"), errors="coerce")
    result["qualifying_frequency_ac"] = pd.to_numeric(
        result.get("qualifying_frequency_ac"), errors="coerce"
    )
    frequency_signal = truthy(result.get("standard_acmg_frequency_flag", pd.Series(False, index=result.index)))
    weak_review = truthy(result.get("weak_review_signal", pd.Series(False, index=result.index)))
    for cutoff in AC_THRESHOLDS:
        result[f"ac_supported_ge_{cutoff}"] = frequency_signal & (
            result["qualifying_frequency_ac"].fillna(0) >= cutoff
        )
        result[f"red_ac_ge_{cutoff}"] = (
            result["vital_score"].fillna(-1).ge(ava.VITAL_ACTION_THRESHOLD)
            & weak_review
            & result[f"ac_supported_ge_{cutoff}"]
        )
    result["vital_red_flag"] = result["red_ac_ge_20"]
    result["is_red"] = result["red_ac_ge_20"]
    return result


def snapshot_summary(scored: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for snapshot, sub in scored.groupby("snapshot", sort=False):
        frequency_observed = sub["frequency_evidence_status"].eq("frequency_observed")
        af_flags = truthy(sub["standard_acmg_frequency_flag"])
        row = {
            "snapshot": snapshot,
            "total": len(sub),
            "frequency_observed": int(frequency_observed.sum()),
            "af_flags_popmax_or_global_gt_1e_5": int(af_flags.sum()),
            "ac_supported_ge_20": int(sub["ac_supported_ge_20"].sum()),
            "red_ac_ge_20": int(sub["red_ac_ge_20"].sum()),
            "red_standard_context": int(
                (sub["red_ac_ge_20"] & sub["inheritance_flag"].eq("standard")).sum()
            ),
            "red_recessive_context_required": int(
                (sub["red_ac_ge_20"] & sub["inheritance_flag"].eq("recessive_context_required")).sum()
            ),
            "median_vital_observed": float(sub.loc[frequency_observed, "vital_score"].median()),
            "max_vital": float(sub["vital_score"].max(skipna=True)),
        }
        rows.append(row)
    return pd.DataFrame(rows)


def ac_threshold_sensitivity(scored: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for snapshot, sub in scored.groupby("snapshot", sort=False):
        for cutoff in AC_THRESHOLDS:
            red = sub[f"red_ac_ge_{cutoff}"]
            ac_supported = sub[f"ac_supported_ge_{cutoff}"]
            red_ids = sub.loc[red, "clinvar_id"].fillna("").astype(str).tolist()
            rows.append(
                {
                    "snapshot": snapshot,
                    "ac_cutoff": cutoff,
                    "ac_supported_frequency_flags": int(ac_supported.sum()),
                    "red_count": int(red.sum()),
                    "red_standard_context": int((red & sub["inheritance_flag"].eq("standard")).sum()),
                    "red_recessive_context_required": int(
                        (red & sub["inheritance_flag"].eq("recessive_context_required")).sum()
                    ),
                    "red_clinvar_ids": "|".join(red_ids),
                    "red_genes": "|".join(sub.loc[red, "gene"].fillna("").astype(str).tolist()),
                }
            )
    return pd.DataFrame(rows)


def red_membership(scored: pd.DataFrame) -> pd.DataFrame:
    red_any = scored[[f"red_ac_ge_{cutoff}" for cutoff in AC_THRESHOLDS]].any(axis=1)
    columns = [
        "snapshot",
        "gene",
        "clinvar_id",
        "variation_id",
        "variant_key",
        "title",
        "review_strength",
        "inheritance_flag",
        "max_frequency_signal",
        "qualifying_frequency_ac",
        "vital_score",
        "vital_band",
        "vital_signal_reason",
        *[f"red_ac_ge_{cutoff}" for cutoff in AC_THRESHOLDS],
    ]
    return scored.loc[red_any, [column for column in columns if column in scored.columns]].sort_values(
        ["snapshot", "vital_score"], ascending=[True, False]
    )


def red_overlap(scored: pd.DataFrame) -> pd.DataFrame:
    sets = {
        snapshot: set(sub.loc[sub["red_ac_ge_20"], "variant_key"].astype(str))
        for snapshot, sub in scored.groupby("snapshot", sort=False)
    }
    rows: list[dict[str, object]] = []
    for left, right in combinations(sets, 2):
        left_set = sets[left]
        right_set = sets[right]
        overlap = left_set & right_set
        union = left_set | right_set
        rows.append(
            {
                "left_snapshot": left,
                "right_snapshot": right,
                "left_red_count": len(left_set),
                "right_red_count": len(right_set),
                "overlap_count": len(overlap),
                "jaccard": len(overlap) / len(union) if union else np.nan,
                "overlap_variant_keys": "|".join(sorted(overlap)),
                "left_only_variant_keys": "|".join(sorted(left_set - right_set)),
                "right_only_variant_keys": "|".join(sorted(right_set - left_set)),
            }
        )
    return pd.DataFrame(rows)


def inheritance_summary(scored: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (snapshot, inheritance_flag), sub in scored.groupby(["snapshot", "inheritance_flag"], sort=False):
        rows.append(
            {
                "snapshot": snapshot,
                "inheritance_flag": inheritance_flag,
                "total": len(sub),
                "af_flags_popmax_or_global_gt_1e_5": int(truthy(sub["standard_acmg_frequency_flag"]).sum()),
                "ac_supported_ge_20": int(sub["ac_supported_ge_20"].sum()),
                "red_ac_ge_20": int(sub["red_ac_ge_20"].sum()),
                "median_vital_observed": float(
                    sub.loc[sub["frequency_evidence_status"].eq("frequency_observed"), "vital_score"].median()
                ),
            }
        )
    return pd.DataFrame(rows)


def af_flags_table(scored: pd.DataFrame) -> pd.DataFrame:
    mask = truthy(scored["standard_acmg_frequency_flag"])
    columns = [
        "snapshot",
        "gene",
        "clinvar_id",
        "variation_id",
        "variant_key",
        "title",
        "review_strength",
        "inheritance_flag",
        "global_af",
        "global_ac",
        "popmax_af",
        "popmax_ac",
        "popmax_population",
        "max_frequency_signal",
        "qualifying_frequency_ac",
        "frequency_signal_ac_ge_20",
        "vital_score",
        "vital_band",
        "vital_red_flag",
    ]
    return scored.loc[mask, [column for column in columns if column in scored.columns]].sort_values(
        ["snapshot", "max_frequency_signal"], ascending=[True, False]
    )


def make_distribution_plot(scored: pd.DataFrame, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    colors = {"2023-01": "#355C7D", "2024-01": "#D95F02", "2026-04": "#1B9E77"}
    bins = np.linspace(0, 100, 41)
    for snapshot, sub in scored.groupby("snapshot", sort=False):
        values = pd.to_numeric(sub["vital_score"], errors="coerce").dropna()
        if values.empty:
            continue
        ax.hist(
            values,
            bins=bins,
            histtype="step",
            linewidth=2,
            density=True,
            label=snapshot,
            color=colors.get(snapshot),
        )
    ax.axvline(70, color="#9B1C31", linestyle="--", linewidth=1.5, label="red gate")
    ax.set_xlabel("VITAL score")
    ax.set_ylabel("Density")
    ax.set_title("Temporal VITAL score distribution")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=220)
    plt.close(fig)
    print(f"Saved {output_path}")


def run(output_prefix: str, download_missing: bool) -> dict[str, Path]:
    ensure_raw_snapshots(download_missing=download_missing)
    match_cache = build_gnomad_match_cache()
    population_cache = build_population_cache()

    scored_frames: list[pd.DataFrame] = []
    for snapshot_name, path in RAW_SNAPSHOTS:
        if not path.exists():
            raise FileNotFoundError(f"Missing ClinVar snapshot: {path}")
        snapshot = load_plp_snapshot(path, genes=ARRHYTHMIA_GENES)
        snapshot.attrs["snapshot"] = snapshot_name
        scored = compute_vital(snapshot, match_cache, population_cache)
        scored_frames.append(scored)
    scored_frames.append(load_current_canonical_scores())

    scored_all = pd.concat(scored_frames, ignore_index=True)
    outputs = {
        "scored_variants": DATA_DIR / f"{output_prefix}_vital_scored_variants.tsv",
        "snapshot_summary": DATA_DIR / f"{output_prefix}_snapshot_summary.tsv",
        "ac_threshold_sensitivity": DATA_DIR / f"{output_prefix}_ac_threshold_sensitivity.tsv",
        "red_threshold_membership": DATA_DIR / f"{output_prefix}_red_threshold_membership.tsv",
        "red_queue": DATA_DIR / f"{output_prefix}_red_queue.tsv",
        "red_overlap": DATA_DIR / f"{output_prefix}_red_overlap.tsv",
        "inheritance_summary": DATA_DIR / f"{output_prefix}_inheritance_summary.tsv",
        "af_flags": DATA_DIR / f"{output_prefix}_af_flags.tsv",
        "distribution_figure": FIGURE_DIR / f"{output_prefix}_vital_score_distribution.png",
        "supplement_scored_variants": SUPPLEMENT_DIR / "Supplementary_Table_S25_temporal_robustness_scored_variants.tsv",
        "supplement_summary": SUPPLEMENT_DIR / "Supplementary_Table_S26_temporal_robustness_summary.tsv",
        "supplement_ac_sensitivity": SUPPLEMENT_DIR / "Supplementary_Table_S27_temporal_ac_threshold_sensitivity.tsv",
    }

    summary = snapshot_summary(scored_all)
    ac_sensitivity = ac_threshold_sensitivity(scored_all)
    membership = red_membership(scored_all)
    red_queue = scored_all.loc[scored_all["red_ac_ge_20"]].copy().sort_values(
        ["snapshot", "vital_score"], ascending=[True, False]
    )
    overlap = red_overlap(scored_all)
    inherit = inheritance_summary(scored_all)
    af_flags = af_flags_table(scored_all)

    save_csv(scored_all, outputs["scored_variants"], sep="\t")
    save_csv(summary, outputs["snapshot_summary"], sep="\t")
    save_csv(ac_sensitivity, outputs["ac_threshold_sensitivity"], sep="\t")
    save_csv(membership, outputs["red_threshold_membership"], sep="\t")
    save_csv(red_queue, outputs["red_queue"], sep="\t")
    save_csv(overlap, outputs["red_overlap"], sep="\t")
    save_csv(inherit, outputs["inheritance_summary"], sep="\t")
    save_csv(af_flags, outputs["af_flags"], sep="\t")
    save_csv(scored_all, outputs["supplement_scored_variants"], sep="\t")
    save_csv(summary, outputs["supplement_summary"], sep="\t")
    save_csv(ac_sensitivity, outputs["supplement_ac_sensitivity"], sep="\t")
    make_distribution_plot(scored_all, outputs["distribution_figure"])

    print("\nTemporal robustness summary")
    print(summary.to_string(index=False))
    print("\nAC threshold sensitivity")
    print(ac_sensitivity.to_string(index=False))
    print("\nRed overlap")
    print(overlap.to_string(index=False))
    return outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Temporal robustness analysis for arrhythmia ClinVar VITAL scoring."
    )
    parser.add_argument(
        "--output-prefix",
        default="arrhythmia_temporal_robustness",
        help="Prefix for processed outputs and figure.",
    )
    parser.add_argument(
        "--no-download",
        action="store_true",
        help="Require local data/raw ClinVar snapshots instead of downloading missing archives.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = run(args.output_prefix, download_missing=not args.no_download)
    print("\nOutputs")
    for name, path in outputs.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
