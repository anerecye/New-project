from __future__ import annotations

import argparse
import logging
import os
import sys
import time
from pathlib import Path
from xml.etree import ElementTree as ET

try:
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import requests
    import seaborn as sns
    from scipy.stats import fisher_exact, mannwhitneyu
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires the full project dependencies. Install them with: "
        "python -m pip install -r requirements.txt"
    ) from exc


SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parents[0]
sys.path.insert(0, str(SCRIPT_DIR))

from advanced_variant_analyses import prepare_annotations, sanitize_output_prefix


logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

DATA_DIR = BASE_DIR / "data" / "processed"
FIGURE_DIR = BASE_DIR / "figures"


def prefixed_name(prefix: str, name: str) -> str:
    return f"{prefix}_{name}" if prefix else name


def data_path(prefix: str, name: str) -> Path:
    return DATA_DIR / prefixed_name(prefix, name)


def figure_path(prefix: str, name: str) -> Path:
    return FIGURE_DIR / prefixed_name(prefix, name)


def resolve_prefix(value: str | None) -> str:
    cleaned = sanitize_output_prefix(value)
    if cleaned:
        return cleaned
    if (DATA_DIR / "arrhythmia_gnomad_matched.csv").exists():
        return "arrhythmia"
    return ""


def save_table(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    log.info("Saved: %s (%s rows)", path, len(df))


def load_matched(prefix: str) -> pd.DataFrame:
    path = data_path(prefix, "gnomad_matched.csv")
    if not path.exists():
        raise FileNotFoundError(f"Matched gnomAD table not found: {path}")
    return pd.read_csv(path)


def add_common_flags(df: pd.DataFrame) -> pd.DataFrame:
    result = prepare_annotations(df)
    result["is_kcnh2"] = result["gene"].eq("KCNH2")
    result["non_overlap"] = result["match_category"].isin(
        ["allele_discordance", "no_gnomad_record"]
    )
    result["overlap_group"] = np.where(result["non_overlap"], "non_overlap", "exact")
    result["broad_type"] = np.where(result["is_indel"], "indel", "SNV_or_MNV")
    result["kb_bin"] = (result["pos"] // 1000).astype("Int64")
    result["indel_size"] = (
        result["ref"].fillna("").astype(str).str.len()
        - result["alt"].fillna("").astype(str).str.len()
    ).abs()
    return result


def fetch_ucsc_sequence(chrom: str, start0: int, end0: int) -> tuple[str, str]:
    url = (
        "https://api.genome.ucsc.edu/getData/sequence"
        f"?genome=hg38;chrom={chrom};start={start0};end={end0}"
    )
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    payload = response.json()
    return str(payload["dna"]).upper(), url


def fetch_ucsc_track(chrom: str, start0: int, end0: int, track: str) -> tuple[list[dict[str, object]], str]:
    url = (
        "https://api.genome.ucsc.edu/getData/track"
        f"?genome=hg38;track={track};chrom={chrom};start={start0};end={end0}"
    )
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    payload = response.json()
    return list(payload.get(track, [])), url


def gc_fraction(sequence: str) -> float:
    valid = [base for base in sequence.upper() if base in "ACGT"]
    if not valid:
        return np.nan
    return (valid.count("G") + valid.count("C")) / len(valid)


def add_reference_gc(
    variants: pd.DataFrame,
    chrom: str,
    radii: tuple[int, ...],
) -> tuple[pd.DataFrame, str]:
    if variants.empty:
        return variants.copy(), ""

    max_radius = max(radii)
    min_pos = int(variants["pos"].min())
    max_pos = int(variants["pos"].max())
    start0 = max(0, min_pos - max_radius - 1)
    end0 = max_pos + max_radius
    sequence, url = fetch_ucsc_sequence(chrom, start0, end0)

    result = variants.copy()

    def window_gc(pos: object, radius: int) -> float:
        position = int(pos)
        window_start0 = position - radius - 1
        window_end0 = position + radius
        i0 = max(0, window_start0 - start0)
        i1 = min(len(sequence), window_end0 - start0)
        return gc_fraction(sequence[i0:i1])

    for radius in radii:
        result[f"local_gc_{radius}bp"] = result["pos"].map(
            lambda pos, radius=radius: window_gc(pos, radius)
        )
    result["reference_sequence_source"] = url
    return result, url


def add_kcnh2_reference_gc(kcnh2: pd.DataFrame, radii: tuple[int, ...]) -> tuple[pd.DataFrame, str]:
    return add_reference_gc(kcnh2, "chr7", radii)


def variant_reference_interval(row: pd.Series, padding: int = 0) -> tuple[int, int]:
    position0 = int(row["pos"]) - 1
    span = max(1, len(str(row.get("ref", ""))), len(str(row.get("alt", ""))))
    start0 = max(0, position0 - padding)
    end0 = position0 + span + padding
    return start0, end0


def track_interval(item: dict[str, object], track: str) -> tuple[int, int]:
    if track == "rmsk":
        return int(item["genoStart"]), int(item["genoEnd"])
    return int(item["chromStart"]), int(item["chromEnd"])


def overlap_length(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def annotate_track_overlap(
    variants: pd.DataFrame,
    track_items: list[dict[str, object]],
    track: str,
    padding: int = 10,
) -> pd.DataFrame:
    result = variants.copy()
    exact_flags: list[bool] = []
    window_flags: list[bool] = []
    names: list[str] = []
    classes: list[str] = []
    families: list[str] = []
    max_overlaps: list[int] = []

    for _, row in result.iterrows():
        exact_start, exact_end = variant_reference_interval(row, padding=0)
        window_start, window_end = variant_reference_interval(row, padding=padding)
        hit_names: list[str] = []
        hit_classes: list[str] = []
        hit_families: list[str] = []
        exact_hit = False
        window_hit = False
        max_overlap = 0

        for item in track_items:
            item_start, item_end = track_interval(item, track)
            exact_overlap = overlap_length(exact_start, exact_end, item_start, item_end)
            window_overlap = overlap_length(window_start, window_end, item_start, item_end)
            if exact_overlap > 0:
                exact_hit = True
            if window_overlap > 0:
                window_hit = True
                max_overlap = max(max_overlap, window_overlap)
                if track == "rmsk":
                    hit_names.append(str(item.get("repName", "")))
                    hit_classes.append(str(item.get("repClass", "")))
                    hit_families.append(str(item.get("repFamily", "")))
                else:
                    hit_names.append(str(item.get("sequence", item.get("name", ""))))
                    hit_classes.append("Simple_repeat")
                    hit_families.append(str(item.get("period", "")))

        exact_flags.append(exact_hit)
        window_flags.append(window_hit)
        names.append(";".join(sorted({value for value in hit_names if value})))
        classes.append(";".join(sorted({value for value in hit_classes if value})))
        families.append(";".join(sorted({value for value in hit_families if value})))
        max_overlaps.append(max_overlap)

    result[f"{track}_overlap_variant"] = exact_flags
    result[f"{track}_overlap_10bp_window"] = window_flags
    result[f"{track}_names_10bp_window"] = names
    result[f"{track}_classes_10bp_window"] = classes
    result[f"{track}_families_10bp_window"] = families
    result[f"{track}_max_overlap_bp_10bp_window"] = max_overlaps
    return result


def add_repeatmasker_annotations(kcnh2: pd.DataFrame) -> tuple[pd.DataFrame, str, str]:
    if kcnh2.empty:
        return kcnh2.copy(), "", ""
    start0 = max(0, int(kcnh2["pos"].min()) - 250)
    end0 = int(kcnh2["pos"].max()) + 250
    rmsk_items, rmsk_url = fetch_ucsc_track("chr7", start0, end0, "rmsk")
    simple_items, simple_url = fetch_ucsc_track("chr7", start0, end0, "simpleRepeat")
    result = annotate_track_overlap(kcnh2, rmsk_items, "rmsk", padding=10)
    result = annotate_track_overlap(result, simple_items, "simpleRepeat", padding=10)
    result["repeatmasker_source"] = rmsk_url
    result["simple_repeat_source"] = simple_url
    return result, rmsk_url, simple_url


def duplication_size_bucket(size: object) -> str:
    try:
        value = int(size)
    except (TypeError, ValueError):
        return "unknown"
    if value <= 10:
        return "short_1_10bp"
    if value <= 50:
        return "medium_11_50bp"
    return "long_gt_50bp"


def make_overall_summary(df: pd.DataFrame, kcnh2: pd.DataFrame, source_url: str) -> pd.DataFrame:
    global_rate = float(df["non_overlap"].mean())
    observed = int(kcnh2["non_overlap"].sum())
    total = int(len(kcnh2))
    expected = total * global_rate

    table = pd.crosstab(df["is_kcnh2"], df["non_overlap"])
    for index in [False, True]:
        if index not in table.index:
            table.loc[index, :] = 0
    for column in [False, True]:
        if column not in table.columns:
            table.loc[:, column] = 0
    table = table.sort_index().sort_index(axis=1)
    odds, p_value = fisher_exact(
        [
            [table.loc[True, True], table.loc[True, False]],
            [table.loc[False, True], table.loc[False, False]],
        ],
        alternative="greater",
    )

    rows = [
        {
            "metric": "global_non_overlap_fraction",
            "value": global_rate,
            "note": "All genes in the matched table.",
        },
        {
            "metric": "kcnh2_total_variants",
            "value": total,
            "note": "",
        },
        {
            "metric": "kcnh2_non_overlap_count",
            "value": observed,
            "note": "",
        },
        {
            "metric": "kcnh2_non_overlap_fraction",
            "value": observed / total if total else np.nan,
            "note": "",
        },
        {
            "metric": "kcnh2_expected_non_overlap_count_at_global_rate",
            "value": expected,
            "note": "",
        },
        {
            "metric": "kcnh2_excess_non_overlap_count",
            "value": observed - expected,
            "note": "Observed minus expected at the all-gene non-overlap rate.",
        },
        {
            "metric": "kcnh2_non_overlap_fisher_or_vs_other_genes",
            "value": odds,
            "note": "",
        },
        {
            "metric": "kcnh2_non_overlap_fisher_p_greater",
            "value": p_value,
            "note": "",
        },
    ]

    if "local_gc_100bp" in kcnh2.columns:
        exact_gc = kcnh2.loc[~kcnh2["non_overlap"], "local_gc_100bp"].dropna()
        non_gc = kcnh2.loc[kcnh2["non_overlap"], "local_gc_100bp"].dropna()
        rows.extend(
            [
                {
                    "metric": "kcnh2_exact_median_local_gc_100bp",
                    "value": float(exact_gc.median()) if not exact_gc.empty else np.nan,
                    "note": source_url,
                },
                {
                    "metric": "kcnh2_non_overlap_median_local_gc_100bp",
                    "value": float(non_gc.median()) if not non_gc.empty else np.nan,
                    "note": source_url,
                },
                {
                    "metric": "kcnh2_local_gc_100bp_mannwhitney_p_two_sided",
                    "value": (
                        mannwhitneyu(non_gc, exact_gc, alternative="two-sided").pvalue
                        if len(exact_gc) and len(non_gc)
                        else np.nan
                    ),
                    "note": "Tests non-overlap vs exact local GC distribution.",
                },
            ]
        )

    return pd.DataFrame(rows)


def parse_clinvar_date(value: object) -> pd.Timestamp:
    return pd.to_datetime(value, errors="coerce")


def fetch_clinvar_vcv_metadata(
    variation_ids: pd.Series,
    email: str | None = None,
    batch_size: int = 100,
    pause: float = 0.35,
) -> pd.DataFrame:
    resolved_email = email or os.getenv("ENTREZ_EMAIL") or ""
    ids = [
        str(int(float(value)))
        for value in variation_ids.dropna().astype(str)
        if str(value).strip()
    ]
    ids = sorted(set(ids), key=lambda value: int(value))
    if not ids:
        return pd.DataFrame()

    rows: list[dict[str, object]] = []
    for start in range(0, len(ids), batch_size):
        batch = ids[start : start + batch_size]
        params = {
            "db": "clinvar",
            "id": ",".join(batch),
            "rettype": "vcv",
            "retmode": "xml",
            "is_variationid": "true",
            "tool": "codex",
        }
        if resolved_email:
            params["email"] = resolved_email
        response = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            params=params,
            timeout=60,
        )
        response.raise_for_status()
        root = ET.fromstring(response.text)
        for archive in root.findall(".//VariationArchive"):
            rows.append(
                {
                    "variation_id": archive.attrib.get("VariationID", ""),
                    "clinvar_id": archive.attrib.get("Accession", ""),
                    "variation_type_xml": archive.attrib.get("VariationType", ""),
                    "number_of_submissions_xml": archive.attrib.get("NumberOfSubmissions", ""),
                    "number_of_submitters_xml": archive.attrib.get("NumberOfSubmitters", ""),
                    "date_created": archive.attrib.get("DateCreated", ""),
                    "date_last_updated": archive.attrib.get("DateLastUpdated", ""),
                    "most_recent_submission": archive.attrib.get("MostRecentSubmission", ""),
                    "record_type": archive.attrib.get("RecordType", ""),
                }
            )
        time.sleep(pause)

    metadata = pd.DataFrame(rows).drop_duplicates(subset="variation_id")
    for column in [
        "number_of_submissions_xml",
        "number_of_submitters_xml",
    ]:
        if column in metadata.columns:
            metadata[column] = pd.to_numeric(metadata[column], errors="coerce")
    for column in ["date_created", "date_last_updated", "most_recent_submission"]:
        if column in metadata.columns:
            metadata[f"{column}_parsed"] = metadata[column].map(parse_clinvar_date)
            metadata[f"{column}_year"] = metadata[f"{column}_parsed"].dt.year
    return metadata


def load_or_fetch_clinvar_metadata(
    kcnh2: pd.DataFrame,
    cache_path: Path,
    email: str | None,
    fetch: bool,
) -> pd.DataFrame:
    if cache_path.exists() and not fetch:
        log.info("Loaded existing ClinVar VCV metadata cache: %s", cache_path)
        return pd.read_csv(cache_path)
    if not fetch:
        log.warning("ClinVar VCV metadata fetch disabled and no cache exists: %s", cache_path)
        return pd.DataFrame()
    metadata = fetch_clinvar_vcv_metadata(kcnh2["variation_id"], email=email)
    save_table(metadata, cache_path)
    return metadata


def make_submission_date_summary(kcnh2: pd.DataFrame) -> pd.DataFrame:
    if "date_created_year" not in kcnh2.columns:
        return pd.DataFrame(
            columns=[
                "date_field",
                "year",
                "variant_count",
                "non_overlap_count",
                "non_overlap_fraction",
            ]
        )

    rows: list[dict[str, object]] = []
    for field in ["date_created_year", "most_recent_submission_year", "date_last_updated_year"]:
        if field not in kcnh2.columns:
            continue
        by_year = (
            kcnh2.dropna(subset=[field])
            .groupby(field)
            .agg(
                variant_count=("variant_key", "count"),
                non_overlap_count=("non_overlap", "sum"),
                duplication_count=("variant_type", lambda values: int((values == "duplication").sum())),
            )
            .reset_index()
        )
        for _, row in by_year.iterrows():
            rows.append(
                {
                    "date_field": field,
                    "year": int(row[field]),
                    "variant_count": int(row["variant_count"]),
                    "non_overlap_count": int(row["non_overlap_count"]),
                    "non_overlap_fraction": (
                        row["non_overlap_count"] / row["variant_count"]
                        if row["variant_count"]
                        else np.nan
                    ),
                    "duplication_count": int(row["duplication_count"]),
                }
            )

    if "date_created_year" in kcnh2.columns:
        dated = kcnh2.dropna(subset=["date_created_year"]).copy()
        if not dated.empty:
            median_year = dated["date_created_year"].median()
            dated["created_period"] = np.where(
                dated["date_created_year"] <= median_year,
                f"created_<=_{int(median_year)}",
                f"created_>_{int(median_year)}",
            )
            by_period = (
                dated.groupby("created_period")
                .agg(
                    variant_count=("variant_key", "count"),
                    non_overlap_count=("non_overlap", "sum"),
                )
                .reset_index()
            )
            for _, row in by_period.iterrows():
                rows.append(
                    {
                        "date_field": "date_created_median_split",
                        "year": row["created_period"],
                        "variant_count": int(row["variant_count"]),
                        "non_overlap_count": int(row["non_overlap_count"]),
                        "non_overlap_fraction": (
                            row["non_overlap_count"] / row["variant_count"]
                            if row["variant_count"]
                            else np.nan
                        ),
                        "duplication_count": np.nan,
                    }
                )
            if len(by_period) == 2:
                old = dated[dated["date_created_year"] <= median_year]
                new = dated[dated["date_created_year"] > median_year]
                odds, p_value = fisher_exact(
                    [
                        [int(old["non_overlap"].sum()), int((~old["non_overlap"]).sum())],
                        [int(new["non_overlap"].sum()), int((~new["non_overlap"]).sum())],
                    ],
                    alternative="greater",
                )
                rows.append(
                    {
                        "date_field": "date_created_median_split_test",
                        "year": f"created_<=_{int(median_year)}_vs_>_{int(median_year)}",
                        "variant_count": len(dated),
                        "non_overlap_count": int(dated["non_overlap"].sum()),
                        "non_overlap_fraction": dated["non_overlap"].mean(),
                        "duplication_count": np.nan,
                        "odds_ratio_older_vs_newer": odds,
                        "fisher_p_greater_older_vs_newer": p_value,
                    }
                )

    return pd.DataFrame(rows)


def make_stratified_tests(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for column in ["broad_type", "variant_type", "functional_class"]:
        for stratum, group in df.groupby(column, dropna=False):
            if len(group) < 10:
                continue
            table = pd.crosstab(group["is_kcnh2"], group["non_overlap"])
            for index in [False, True]:
                if index not in table.index:
                    table.loc[index, :] = 0
            for value in [False, True]:
                if value not in table.columns:
                    table.loc[:, value] = 0
            table = table.sort_index().sort_index(axis=1)
            odds, p_value = fisher_exact(
                [
                    [table.loc[True, True], table.loc[True, False]],
                    [table.loc[False, True], table.loc[False, False]],
                ],
                alternative="greater",
            )
            kcnh2_total = int(table.loc[True].sum())
            other_total = int(table.loc[False].sum())
            rows.append(
                {
                    "stratification": column,
                    "stratum": str(stratum),
                    "kcnh2_total": kcnh2_total,
                    "kcnh2_non_overlap": int(table.loc[True, True]),
                    "kcnh2_non_overlap_fraction": (
                        table.loc[True, True] / kcnh2_total if kcnh2_total else np.nan
                    ),
                    "other_total": other_total,
                    "other_non_overlap": int(table.loc[False, True]),
                    "other_non_overlap_fraction": (
                        table.loc[False, True] / other_total if other_total else np.nan
                    ),
                    "odds_ratio": odds,
                    "fisher_p_greater": p_value,
                }
            )
    return pd.DataFrame(rows).sort_values(
        ["stratification", "kcnh2_non_overlap"], ascending=[True, False]
    )


def make_type_contribution(df: pd.DataFrame) -> pd.DataFrame:
    kcnh2 = df[df["is_kcnh2"]].copy()
    rows: list[dict[str, object]] = []
    for variant_type, group in kcnh2.groupby("variant_type", dropna=False):
        other = df[(~df["is_kcnh2"]) & (df["variant_type"].eq(variant_type))]
        if other.empty:
            continue
        other_rate = float(other["non_overlap"].mean())
        observed = int(group["non_overlap"].sum())
        expected = len(group) * other_rate
        rows.append(
            {
                "variant_type": variant_type,
                "kcnh2_total": len(group),
                "kcnh2_observed_non_overlap": observed,
                "other_gene_non_overlap_rate_same_type": other_rate,
                "kcnh2_expected_non_overlap_at_other_rate": expected,
                "kcnh2_excess_non_overlap_same_type_adjusted": observed - expected,
            }
        )
    return pd.DataFrame(rows).sort_values(
        "kcnh2_excess_non_overlap_same_type_adjusted", ascending=False
    )


def make_duplication_size_summary(kcnh2: pd.DataFrame) -> pd.DataFrame:
    duplications = kcnh2[kcnh2["variant_type"].eq("duplication")].copy()
    if duplications.empty:
        return pd.DataFrame()

    duplications["duplication_size_bucket"] = duplications["indel_size"].map(
        duplication_size_bucket
    )
    rows: list[dict[str, object]] = []
    for bucket, group in duplications.groupby("duplication_size_bucket", dropna=False):
        repeat_columns = [
            column
            for column in [
                "rmsk_overlap_10bp_window",
                "simpleRepeat_overlap_10bp_window",
            ]
            if column in group.columns
        ]
        row: dict[str, object] = {
            "duplication_size_bucket": bucket,
            "variant_count": len(group),
            "non_overlap_count": int(group["non_overlap"].sum()),
            "non_overlap_fraction": group["non_overlap"].mean(),
            "allele_discordance_count": int((group["match_category"] == "allele_discordance").sum()),
            "no_gnomad_record_count": int((group["match_category"] == "no_gnomad_record").sum()),
            "exact_match_count": int((group["match_category"] == "exact_match").sum()),
            "median_duplication_size_bp": group["indel_size"].median(),
            "max_duplication_size_bp": group["indel_size"].max(),
        }
        for column in repeat_columns:
            row[f"{column}_fraction"] = group[column].fillna(False).astype(bool).mean()
        rows.append(row)

    return pd.DataFrame(rows).sort_values(
        ["duplication_size_bucket", "variant_count"], ascending=[True, False]
    )


def make_duplication_sensitivity(kcnh2: pd.DataFrame, global_non_overlap_rate: float) -> pd.DataFrame:
    filters = [
        ("none_removed", pd.Series(False, index=kcnh2.index)),
        (
            "remove_long_duplications_gt_50bp",
            kcnh2["variant_type"].eq("duplication") & (kcnh2["indel_size"] > 50),
        ),
        ("remove_all_duplications", kcnh2["variant_type"].eq("duplication")),
    ]
    if "rmsk_overlap_10bp_window" in kcnh2.columns:
        filters.append(
            (
                "remove_repeatmasker_window_duplications",
                kcnh2["variant_type"].eq("duplication")
                & kcnh2["rmsk_overlap_10bp_window"].fillna(False).astype(bool),
            )
        )

    rows: list[dict[str, object]] = []
    for scenario, remove_mask in filters:
        remaining = kcnh2.loc[~remove_mask].copy()
        observed = int(remaining["non_overlap"].sum())
        expected = len(remaining) * global_non_overlap_rate
        rows.append(
            {
                "scenario": scenario,
                "removed_variant_count": int(remove_mask.sum()),
                "remaining_variant_count": len(remaining),
                "remaining_non_overlap_count": observed,
                "remaining_non_overlap_fraction": (
                    observed / len(remaining) if len(remaining) else np.nan
                ),
                "expected_non_overlap_count_at_global_rate": expected,
                "excess_non_overlap_count": observed - expected,
            }
        )
    return pd.DataFrame(rows)


def make_repeat_overlap_summary(kcnh2: pd.DataFrame) -> pd.DataFrame:
    repeat_columns = [
        "rmsk_overlap_variant",
        "rmsk_overlap_10bp_window",
        "simpleRepeat_overlap_variant",
        "simpleRepeat_overlap_10bp_window",
    ]
    rows: list[dict[str, object]] = []
    for column in repeat_columns:
        if column not in kcnh2.columns:
            continue
        for subset_name, subset in [
            ("all_kcnh2", kcnh2),
            ("duplications", kcnh2[kcnh2["variant_type"].eq("duplication")]),
            ("long_duplications_gt_50bp", kcnh2[
                kcnh2["variant_type"].eq("duplication") & (kcnh2["indel_size"] > 50)
            ]),
        ]:
            if subset.empty:
                continue
            exact = subset[~subset["non_overlap"]]
            non = subset[subset["non_overlap"]]
            exact_pos = int(exact[column].fillna(False).astype(bool).sum())
            non_pos = int(non[column].fillna(False).astype(bool).sum())
            table = [
                [non_pos, len(non) - non_pos],
                [exact_pos, len(exact) - exact_pos],
            ]
            odds, p_value = (
                fisher_exact(table, alternative="two-sided")
                if len(exact) and len(non)
                else (np.nan, np.nan)
            )
            rows.append(
                {
                    "subset": subset_name,
                    "repeat_feature": column,
                    "exact_positive": exact_pos,
                    "exact_total": len(exact),
                    "exact_fraction": exact_pos / len(exact) if len(exact) else np.nan,
                    "non_overlap_positive": non_pos,
                    "non_overlap_total": len(non),
                    "non_overlap_fraction": non_pos / len(non) if len(non) else np.nan,
                    "odds_ratio_non_overlap_vs_exact": odds,
                    "fisher_p_two_sided": p_value,
                }
            )
    return pd.DataFrame(rows)


def make_gc_control_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for gene in ["KCNH2", "SCN5A"]:
        group = df[df["gene"].eq(gene)].copy()
        if group.empty:
            continue
        exact = group[~group["non_overlap"]]
        non = group[group["non_overlap"]]
        row: dict[str, object] = {
            "gene": gene,
            "variant_count": len(group),
            "exact_match_count": len(exact),
            "non_overlap_count": len(non),
            "non_overlap_fraction": len(non) / len(group) if len(group) else np.nan,
            "allele_discordance_count": int((group["match_category"] == "allele_discordance").sum()),
            "no_gnomad_record_count": int((group["match_category"] == "no_gnomad_record").sum()),
            "indel_fraction": group["is_indel"].mean(),
            "duplication_fraction": group["variant_type"].eq("duplication").mean(),
            "median_indel_size": group["indel_size"].median(),
        }
        if "local_gc_100bp" in group.columns:
            row["median_local_gc_100bp"] = group["local_gc_100bp"].median()
            row["exact_median_local_gc_100bp"] = exact["local_gc_100bp"].median()
            row["non_overlap_median_local_gc_100bp"] = non["local_gc_100bp"].median()
        rows.append(row)

    result = pd.DataFrame(rows)
    if len(result) == 2:
        k = df[df["gene"].eq("KCNH2")]
        s = df[df["gene"].eq("SCN5A")]
        table = [
            [int(k["non_overlap"].sum()), int((~k["non_overlap"]).sum())],
            [int(s["non_overlap"].sum()), int((~s["non_overlap"]).sum())],
        ]
        odds, p_value = fisher_exact(table, alternative="two-sided")
        result["kcnh2_vs_scn5a_non_overlap_or"] = odds
        result["kcnh2_vs_scn5a_non_overlap_fisher_p"] = p_value
    return result


def make_kb_bins(kcnh2: pd.DataFrame) -> pd.DataFrame:
    agg_map = {
        "variant_count": ("variant_key", "count"),
        "non_overlap_count": ("non_overlap", "sum"),
        "allele_discordance_count": (
            "match_category",
            lambda values: int((values == "allele_discordance").sum()),
        ),
        "no_gnomad_record_count": (
            "match_category",
            lambda values: int((values == "no_gnomad_record").sum()),
        ),
        "exact_match_count": (
            "match_category",
            lambda values: int((values == "exact_match").sum()),
        ),
        "indel_fraction": ("is_indel", "mean"),
        "duplication_count": (
            "variant_type",
            lambda values: int((values == "duplication").sum()),
        ),
        "median_indel_size": ("indel_size", "median"),
    }
    if "local_gc_100bp" in kcnh2.columns:
        agg_map["median_local_gc_100bp"] = ("local_gc_100bp", "median")
    result = kcnh2.groupby("kb_bin").agg(**agg_map).reset_index()
    result["non_overlap_fraction"] = result["non_overlap_count"] / result["variant_count"]
    return result.sort_values("variant_count", ascending=False)


def make_variant_context(kcnh2: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "variant_key",
        "gene",
        "chrom",
        "pos",
        "ref",
        "alt",
        "match_category",
        "overlap_group",
        "title",
        "variant_type",
        "functional_class",
        "is_indel",
        "indel_size",
        "is_cpg_transition_proxy",
        "is_hotspot_position",
        "is_local_hotspot_5bp",
        "position_variant_count",
        "nearby_5bp_variant_count",
        "kb_bin",
        "gnomad_af",
        "gnomad_ac",
        "gnomad_an",
        "review_status",
        "date_created",
        "date_created_year",
        "date_last_updated",
        "date_last_updated_year",
        "most_recent_submission",
        "most_recent_submission_year",
        "number_of_submitters_xml",
        "number_of_submissions_xml",
        "rmsk_overlap_variant",
        "rmsk_overlap_10bp_window",
        "rmsk_names_10bp_window",
        "rmsk_classes_10bp_window",
        "rmsk_families_10bp_window",
        "simpleRepeat_overlap_variant",
        "simpleRepeat_overlap_10bp_window",
        "simpleRepeat_names_10bp_window",
        "simpleRepeat_families_10bp_window",
    ]
    columns.extend(
        column for column in kcnh2.columns if column.startswith("local_gc_")
    )
    if "reference_sequence_source" in kcnh2.columns:
        columns.append("reference_sequence_source")
    return kcnh2.loc[:, [column for column in columns if column in kcnh2.columns]].sort_values(
        ["match_category", "pos", "variant_key"]
    )


def plot_diagnostics(
    kcnh2: pd.DataFrame,
    kb_bins: pd.DataFrame,
    duplication_summary: pd.DataFrame,
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    counts = kcnh2["match_category"].value_counts().reindex(
        ["exact_match", "allele_discordance", "no_gnomad_record"], fill_value=0
    )
    sns.barplot(x=counts.index, y=counts.values, color="#2563EB", ax=axes[0, 0])
    axes[0, 0].set_title("KCNH2 match categories")
    axes[0, 0].set_xlabel("")
    axes[0, 0].set_ylabel("Variant count")
    axes[0, 0].tick_params(axis="x", rotation=20)

    type_counts = (
        kcnh2.groupby(["overlap_group", "variant_type"])
        .size()
        .rename("count")
        .reset_index()
    )
    type_pivot = type_counts.pivot_table(
        index="overlap_group", columns="variant_type", values="count", fill_value=0
    )
    type_pivot = type_pivot.div(type_pivot.sum(axis=1), axis=0)
    type_pivot.plot(kind="bar", stacked=True, ax=axes[0, 1], colormap="tab20")
    axes[0, 1].set_title("Variant-type composition")
    axes[0, 1].set_xlabel("")
    axes[0, 1].set_ylabel("Fraction")
    axes[0, 1].tick_params(axis="x", rotation=0)
    axes[0, 1].legend(title="Variant type", frameon=False, fontsize=8)

    if "local_gc_100bp" in kcnh2.columns:
        sns.boxplot(
            data=kcnh2,
            x="match_category",
            y="local_gc_100bp",
            order=["exact_match", "allele_discordance", "no_gnomad_record"],
            color="#93C5FD",
            ax=axes[1, 0],
        )
        axes[1, 0].set_title("Local hg38 GC (+/-100 bp)")
        axes[1, 0].set_xlabel("")
        axes[1, 0].set_ylabel("GC fraction")
        axes[1, 0].tick_params(axis="x", rotation=20)
    else:
        axes[1, 0].axis("off")

    plot_bins = kb_bins.sort_values("kb_bin").copy()
    width = 0.8
    axes[1, 1].bar(
        plot_bins["kb_bin"].astype(str),
        plot_bins["allele_discordance_count"],
        width=width,
        label="allele_discordance",
        color="#F59E0B",
    )
    axes[1, 1].bar(
        plot_bins["kb_bin"].astype(str),
        plot_bins["no_gnomad_record_count"],
        width=width,
        bottom=plot_bins["allele_discordance_count"],
        label="no_gnomad_record",
        color="#9CA3AF",
    )
    axes[1, 1].bar(
        plot_bins["kb_bin"].astype(str),
        plot_bins["exact_match_count"],
        width=width,
        bottom=plot_bins["allele_discordance_count"] + plot_bins["no_gnomad_record_count"],
        label="exact_match",
        color="#16A34A",
    )
    axes[1, 1].set_title("KCNH2 1 kb bins")
    axes[1, 1].set_xlabel("chr7 position / 1 kb")
    axes[1, 1].set_ylabel("Variant count")
    axes[1, 1].tick_params(axis="x", rotation=70)
    axes[1, 1].legend(frameon=False, fontsize=8)

    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    log.info("Saved: %s", output_path)

    if not duplication_summary.empty:
        dup_path = output_path.with_name(output_path.stem.replace("diagnostics", "duplication_sizes") + output_path.suffix)
        fig, ax = plt.subplots(figsize=(8, 4.5))
        plot_df = duplication_summary.copy()
        ax.bar(
            plot_df["duplication_size_bucket"],
            plot_df["non_overlap_count"],
            label="non_overlap",
            color="#F59E0B",
        )
        ax.bar(
            plot_df["duplication_size_bucket"],
            plot_df["exact_match_count"],
            bottom=plot_df["non_overlap_count"],
            label="exact_match",
            color="#16A34A",
        )
        ax.set_title("KCNH2 duplication size buckets")
        ax.set_xlabel("")
        ax.set_ylabel("Variant count")
        ax.tick_params(axis="x", rotation=20)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(dup_path, dpi=180)
        plt.close(fig)
        log.info("Saved: %s", dup_path)


def run(
    output_prefix: str | None = None,
    fetch_reference_sequence: bool = True,
    fetch_repeatmasker: bool = True,
    use_clinvar_dates: bool = True,
    fetch_clinvar_dates: bool = False,
    email: str | None = None,
    make_plots: bool = True,
) -> dict[str, Path]:
    prefix = resolve_prefix(output_prefix)
    annotated = add_common_flags(load_matched(prefix))
    kcnh2 = annotated[annotated["is_kcnh2"]].copy()

    source_url = ""
    if fetch_reference_sequence and not kcnh2.empty:
        kcnh2, source_url = add_kcnh2_reference_gc(kcnh2, radii=(25, 50, 100, 200))
        scn5a_mask = annotated["gene"].eq("SCN5A")
        if scn5a_mask.any():
            scn5a, _ = add_reference_gc(
                annotated.loc[scn5a_mask].copy(),
                "chr3",
                radii=(100,),
            )
            for column in [column for column in scn5a.columns if column.startswith("local_gc_")]:
                annotated.loc[scn5a.index, column] = scn5a[column]
            annotated.loc[scn5a.index, "reference_sequence_source"] = scn5a[
                "reference_sequence_source"
            ]

    if fetch_repeatmasker and not kcnh2.empty:
        kcnh2, _, _ = add_repeatmasker_annotations(kcnh2)

    if use_clinvar_dates:
        metadata_cache = data_path(prefix, "kcnh2_clinvar_vcv_metadata.csv")
        metadata = load_or_fetch_clinvar_metadata(
            kcnh2,
            metadata_cache,
            email=email,
            fetch=fetch_clinvar_dates or not metadata_cache.exists(),
        )
        if not metadata.empty:
            kcnh2["variation_id"] = kcnh2["variation_id"].astype(str)
            kcnh2["clinvar_id"] = kcnh2["clinvar_id"].astype(str)
            metadata["variation_id"] = metadata["variation_id"].astype(str)
            metadata["clinvar_id"] = metadata["clinvar_id"].astype(str)
            overlap = [
                column
                for column in metadata.columns
                if column in kcnh2.columns and column not in {"variation_id", "clinvar_id"}
            ]
            if overlap:
                kcnh2 = kcnh2.drop(columns=overlap)
            kcnh2 = kcnh2.merge(metadata, on=["variation_id", "clinvar_id"], how="left")

    paths = {
        "summary": data_path(prefix, "kcnh2_non_overlap_summary.csv"),
        "stratified_tests": data_path(prefix, "kcnh2_non_overlap_stratified_tests.csv"),
        "type_contribution": data_path(prefix, "kcnh2_non_overlap_type_contribution.csv"),
        "duplication_size_summary": data_path(prefix, "kcnh2_duplication_size_summary.csv"),
        "duplication_sensitivity": data_path(prefix, "kcnh2_duplication_sensitivity.csv"),
        "repeat_overlap_summary": data_path(prefix, "kcnh2_repeat_overlap_summary.csv"),
        "gc_control_summary": data_path(prefix, "kcnh2_gc_control_summary.csv"),
        "submission_date_summary": data_path(prefix, "kcnh2_submission_date_summary.csv"),
        "kb_bins": data_path(prefix, "kcnh2_non_overlap_kb_bins.csv"),
        "variant_context": data_path(prefix, "kcnh2_non_overlap_variant_context.csv"),
    }

    summary = make_overall_summary(annotated, kcnh2, source_url)
    stratified = make_stratified_tests(annotated)
    type_contribution = make_type_contribution(annotated)
    duplication_summary = make_duplication_size_summary(kcnh2)
    duplication_sensitivity = make_duplication_sensitivity(
        kcnh2,
        global_non_overlap_rate=float(annotated["non_overlap"].mean()),
    )
    repeat_summary = make_repeat_overlap_summary(kcnh2)
    gc_control = make_gc_control_summary(
        pd.concat(
            [
                kcnh2,
                annotated[annotated["gene"].eq("SCN5A")],
            ],
            ignore_index=True,
            sort=False,
        )
    )
    submission_dates = make_submission_date_summary(kcnh2)
    kb_bins = make_kb_bins(kcnh2)
    variant_context = make_variant_context(kcnh2)

    save_table(summary, paths["summary"])
    save_table(stratified, paths["stratified_tests"])
    save_table(type_contribution, paths["type_contribution"])
    save_table(duplication_summary, paths["duplication_size_summary"])
    save_table(duplication_sensitivity, paths["duplication_sensitivity"])
    save_table(repeat_summary, paths["repeat_overlap_summary"])
    save_table(gc_control, paths["gc_control_summary"])
    save_table(submission_dates, paths["submission_date_summary"])
    save_table(kb_bins, paths["kb_bins"])
    save_table(variant_context, paths["variant_context"])

    if make_plots:
        plot_diagnostics(
            kcnh2,
            kb_bins,
            duplication_summary,
            figure_path(prefix, "kcnh2_non_overlap_diagnostics.png"),
        )

    return paths


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Investigate why KCNH2 has the strongest non-overlap signal."
    )
    parser.add_argument("--output-prefix", help="Prefix used by the matched gnomAD CSV.")
    parser.add_argument(
        "--skip-reference-sequence",
        action="store_true",
        help="Do not fetch hg38 sequence from the UCSC API for local GC windows.",
    )
    parser.add_argument(
        "--skip-repeatmasker",
        action="store_true",
        help="Do not fetch UCSC RepeatMasker and simpleRepeat tracks.",
    )
    parser.add_argument(
        "--skip-clinvar-dates",
        action="store_true",
        help="Do not fetch or create ClinVar VCV date metadata.",
    )
    parser.add_argument(
        "--refresh-clinvar-dates",
        action="store_true",
        help="Refresh ClinVar VCV date metadata even if the cache exists.",
    )
    parser.add_argument("--email", help="Optional email for NCBI E-utilities.")
    parser.add_argument("--no-plots", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = run(
        output_prefix=args.output_prefix,
        fetch_reference_sequence=not args.skip_reference_sequence,
        fetch_repeatmasker=not args.skip_repeatmasker,
        use_clinvar_dates=not args.skip_clinvar_dates,
        fetch_clinvar_dates=not args.skip_clinvar_dates and args.refresh_clinvar_dates,
        email=args.email,
        make_plots=not args.no_plots,
    )
    print("KCNH2 diagnostic outputs")
    for label, path in outputs.items():
        print(f"{label}: {path}")


if __name__ == "__main__":
    main()
