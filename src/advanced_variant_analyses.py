from __future__ import annotations

import argparse
import json
import logging
import math
import re
import time
from pathlib import Path
from typing import Iterable

try:
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import requests
    import seaborn as sns
    from scipy.stats import fisher_exact
    from statsmodels.stats.multitest import multipletests
    from tqdm import tqdm
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This script requires the full project dependencies. Install them with: "
        "python -m pip install -r requirements.txt"
    ) from exc


logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
FIGURE_DIR = BASE_DIR / "figures"

GNOMAD_GRAPHQL = "https://gnomad.broadinstitute.org/api"
GNOMAD_DATASET = "gnomad_r4"
AF_ULTRA_RARE = 1e-5
AF_RARE = 1e-4
MIN_RECLASSIFICATION_AC = 20
LIFT_PSEUDOCOUNT = 1e-12

POPULATIONS = ("afr", "amr", "asj", "eas", "fin", "mid", "nfe", "sas", "remaining")
POPULATION_LABELS = {
    "afr": "AFR",
    "amr": "AMR",
    "asj": "ASJ",
    "eas": "EAS",
    "fin": "FIN",
    "mid": "MID",
    "nfe": "NFE",
    "sas": "SAS",
    "remaining": "Remaining",
}

GNOMAD_POPULATION_QUERY = """
query VariantPopulationQuery($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variantId
    exome {
      af
      ac
      an
      faf95 {
        popmax
        popmax_population
      }
      faf99 {
        popmax
        popmax_population
      }
      populations {
        id
        ac
        an
      }
    }
  }
}
"""

GNOMAD_EXOME_GENOME_QUERY = """
query VariantExomeGenomeQuery($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variantId
    exome {
      af
      ac
      an
    }
    genome {
      af
      ac
      an
    }
  }
}
"""


def sanitize_output_prefix(value: str | None) -> str:
    if value is None:
        return ""
    cleaned = "".join(
        char if char.isalnum() or char in {"-", "_"} else "_" for char in value.strip()
    ).strip("_")
    return cleaned


def resolve_output_prefix(value: str | None) -> str:
    cleaned = sanitize_output_prefix(value)
    if cleaned:
        return cleaned
    if (DATA_DIR / "arrhythmia_gnomad_matched.csv").exists():
        return "arrhythmia"
    return ""


def prefixed_name(prefix: str, name: str) -> str:
    return f"{prefix}_{name}" if prefix else name


def data_path(prefix: str, name: str) -> Path:
    return DATA_DIR / prefixed_name(prefix, name)


def figure_path(prefix: str, name: str) -> Path:
    return FIGURE_DIR / prefixed_name(prefix, name)


def save_table(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    log.info("Saved: %s (%s rows)", output_path, len(df))


def normalize_chrom(value: object) -> str:
    chrom = str(value).strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return chrom


def gnomad_variant_id(row: pd.Series) -> str:
    return f"{normalize_chrom(row['chrom'])}-{int(row['pos'])}-{row['ref']}-{row['alt']}"


def to_float(value: object) -> float:
    try:
        if pd.isna(value):
            return math.nan
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def to_int(value: object) -> int | None:
    try:
        if pd.isna(value):
            return None
        return int(float(value))
    except (TypeError, ValueError):
        return None


def safe_af(ac: object, an: object) -> float:
    ac_value = to_float(ac)
    an_value = to_float(an)
    if math.isnan(ac_value) or math.isnan(an_value) or an_value <= 0:
        return math.nan
    return ac_value / an_value


def parse_bool_series(series: pd.Series) -> pd.Series:
    normalized = series.fillna("").astype(str).str.strip().str.lower()
    return normalized.isin({"true", "1", "yes"})


def ensure_variant_key(df: pd.DataFrame) -> pd.DataFrame:
    result = df.copy()
    if "variant_key" not in result.columns:
        result["variant_key"] = (
            result["chrom"].astype(str)
            + ":"
            + result["pos"].astype(str)
            + ":"
            + result["ref"].astype(str)
            + ":"
            + result["alt"].astype(str)
        )
    return result


def first_available_column(df: pd.DataFrame, candidates: Iterable[str]) -> pd.Series:
    for column in candidates:
        if column in df.columns:
            return df[column]
    return pd.Series("", index=df.index, dtype="object")


def infer_submitter_count_from_review(status: object) -> float:
    text = str(status or "").strip().lower()
    if "single submitter" in text:
        return 1.0
    if "multiple submitters" in text:
        return 2.0
    return math.nan


def classify_review_strength(status: object) -> tuple[str, int]:
    text = str(status or "").strip().lower()
    if not text:
        return "missing_review", 0
    if "practice guideline" in text:
        return "practice_guideline", 4
    if "reviewed by expert panel" in text:
        return "expert_panel", 3
    if "multiple submitters" in text and "no conflicts" in text:
        return "multiple_submitters_no_conflicts", 2
    if "single submitter" in text:
        return "single_submitter", 1
    if "no assertion" in text or "no criteria" in text:
        return "weak_or_no_assertion", 0
    return "other_review", 1


def clean_existing_text(value: object) -> str:
    text = " ".join(str(value or "").strip().split())
    if text.lower() in {"nan", "none", "na", "n/a"}:
        return ""
    return text


def infer_variant_type(row: pd.Series) -> str:
    existing = clean_existing_text(row.get("variant_type", ""))
    if existing:
        return normalize_variant_type(existing)

    title = str(row.get("title", row.get("Name", ""))).lower()
    ref = str(row.get("ref", ""))
    alt = str(row.get("alt", ""))

    if "deletion" in title or "del" in title and len(ref) > len(alt):
        return "deletion"
    if "duplication" in title or "dup" in title and len(alt) > len(ref):
        return "duplication"
    if "insertion" in title or "ins" in title and len(alt) > len(ref):
        return "insertion"
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV_or_substitution"
    if len(ref) > len(alt):
        return "deletion"
    if len(alt) > len(ref):
        return "insertion"
    return "other"


def normalize_variant_type(value: object) -> str:
    text = clean_existing_text(value).lower()
    if not text:
        return "other"
    if "single nucleotide" in text or text in {"snv", "snp"}:
        return "SNV"
    if "deletion" in text:
        return "deletion"
    if "duplication" in text:
        return "duplication"
    if "insertion" in text:
        return "insertion"
    if "indel" in text:
        return "indel"
    if "substitution" in text:
        return "SNV"
    if "mnv" in text:
        return "MNV_or_substitution"
    return clean_existing_text(value)


def classify_functional_effect(row: pd.Series) -> str:
    consequence = str(row.get("molecular_consequence", ""))
    title = str(row.get("title", row.get("Name", "")))
    protein = str(row.get("protein_change", ""))
    variant_type = str(row.get("variant_type", ""))
    text = f"{consequence} {title} {protein} {variant_type}".lower()

    if any(token in text for token in ["frameshift", "stop gained", "nonsense"]):
        return "LOF"
    if re.search(r"p\.[a-z]{3}\d+.*fs", text):
        return "LOF"
    if "ter" in text or "*" in text:
        return "LOF"
    if any(token in text for token in ["splice donor", "splice acceptor"]):
        return "LOF"
    if "start lost" in text or "initiator codon" in text:
        return "LOF"
    if any(token in text for token in ["splice region", "intron variant", "intron"]):
        return "splice_or_intronic"
    if "missense" in text:
        return "missense"
    if pd.notna(row.get("title", pd.NA)):
        # ClinVar titles often carry protein changes even when XML consequence is unavailable.
        if re.search(r"p\.[A-Za-z]{3}\d+[A-Za-z]{3}", title):
            return "missense"
        if re.search(r"c\.[0-9_*]+[+-][12][A-Z]?>", title):
            return "LOF"
        if re.search(r"c\.[0-9_*]+[+-]\d+", title):
            return "splice_or_intronic"
    if "synonymous" in text:
        return "synonymous"
    if any(token in text for token in ["inframe", "in-frame"]):
        return "inframe_indel"
    return "other_or_unresolved"


def add_nearby_cluster_counts(df: pd.DataFrame, window_bp: int = 5) -> pd.Series:
    counts = pd.Series(1, index=df.index, dtype="int64")
    if df.empty:
        return counts

    grouped = df.groupby(["chrom", "gene"], dropna=False)
    for _, group in grouped:
        positions = pd.to_numeric(group["pos"], errors="coerce")
        valid = positions.dropna().astype(int)
        if valid.empty:
            continue
        ordered_index = valid.sort_values().index
        ordered_positions = valid.loc[ordered_index].to_numpy()
        left = np.searchsorted(ordered_positions, ordered_positions - window_bp, side="left")
        right = np.searchsorted(ordered_positions, ordered_positions + window_bp, side="right")
        counts.loc[ordered_index] = right - left
    return counts


def prepare_annotations(df: pd.DataFrame) -> pd.DataFrame:
    result = ensure_variant_key(df)
    result["chrom"] = first_available_column(result, ["chrom", "Chromosome", "CHROM"]).astype(str)
    result["pos"] = pd.to_numeric(first_available_column(result, ["pos", "Start", "POS"]), errors="coerce")
    result["ref"] = first_available_column(result, ["ref", "ReferenceAllele", "REF"]).astype(str)
    result["alt"] = first_available_column(result, ["alt", "AlternateAllele", "ALT"]).astype(str)
    result["gene"] = first_available_column(result, ["gene", "GeneSymbol"]).astype(str)
    result["review_status"] = first_available_column(
        result, ["review_status", "best_review_status", "ReviewStatus_values", "ReviewStatus"]
    ).astype(str)
    result["title"] = first_available_column(result, ["title", "Name"]).astype(str)
    result["gnomad_af"] = pd.to_numeric(result.get("gnomad_af", result.get("AF")), errors="coerce")
    result["gnomad_ac"] = pd.to_numeric(result.get("gnomad_ac", result.get("AC")), errors="coerce")
    result["gnomad_an"] = pd.to_numeric(result.get("gnomad_an", result.get("AN")), errors="coerce")

    if "match_category" not in result.columns:
        if "gnomad_match" in result.columns:
            match = parse_bool_series(result["gnomad_match"])
            result["match_category"] = np.where(match, "exact_match", "no_gnomad_record")
        else:
            result["match_category"] = "unknown"

    observed_submitters = pd.to_numeric(
        result.get("submitter_count", pd.Series(np.nan, index=result.index)),
        errors="coerce",
    )
    inferred_submitters = result["review_status"].map(infer_submitter_count_from_review)
    result["submitter_count_observed"] = observed_submitters
    result["submitter_count"] = observed_submitters.combine_first(inferred_submitters)
    result["submitter_count_source"] = np.select(
        [
            observed_submitters.notna(),
            inferred_submitters.notna(),
        ],
        [
            "clinvar_xml",
            "review_status_lower_bound",
        ],
        default="unavailable",
    )

    review = result["review_status"].map(classify_review_strength)
    result["review_strength"] = [item[0] for item in review]
    result["review_score"] = [item[1] for item in review]

    if "variant_type" not in result.columns:
        result["variant_type"] = ""
    if "molecular_consequence" not in result.columns:
        result["molecular_consequence"] = ""
    if "protein_change" not in result.columns:
        result["protein_change"] = ""

    result["variant_type"] = result.apply(infer_variant_type, axis=1)
    result["functional_class"] = result.apply(classify_functional_effect, axis=1)
    result["is_snv"] = (result["ref"].str.len() == 1) & (result["alt"].str.len() == 1)
    result["is_indel"] = result["ref"].str.len() != result["alt"].str.len()
    result["is_cpg_transition_proxy"] = result.apply(
        lambda row: bool(
            row["is_snv"]
            and (
                (str(row["ref"]).upper(), str(row["alt"]).upper()) in {("C", "T"), ("G", "A")}
            )
        ),
        axis=1,
    )
    result["position_key"] = result["chrom"].map(normalize_chrom) + ":" + result["pos"].astype("Int64").astype(str)
    result["position_variant_count"] = result.groupby("position_key")["variant_key"].transform("count")
    result["nearby_5bp_variant_count"] = add_nearby_cluster_counts(result, window_bp=5)
    result["is_hotspot_position"] = result["position_variant_count"] >= 2
    result["is_local_hotspot_5bp"] = result["nearby_5bp_variant_count"] >= 3
    return result


def post_gnomad_graphql(
    session: requests.Session,
    query: str,
    variables: dict[str, object],
    retries: int = 4,
    base_delay: float = 1.0,
) -> dict[str, object]:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            response = session.post(
                GNOMAD_GRAPHQL,
                json={"query": query, "variables": variables},
                timeout=30,
            )
            payload = response.json()
            if payload.get("errors"):
                raise RuntimeError(json.dumps(payload["errors"], ensure_ascii=True))
            response.raise_for_status()
            return payload
        except (requests.RequestException, RuntimeError, ValueError) as exc:
            last_error = exc
            if attempt == retries:
                break
            delay = base_delay * (2 ** (attempt - 1))
            log.warning("gnomAD query failed on attempt %s/%s: %s", attempt, retries, exc)
            time.sleep(delay)
    raise RuntimeError(f"gnomAD query failed after retries: {last_error}")


def fetch_population_af_record(
    session: requests.Session,
    row: pd.Series,
    dataset: str = GNOMAD_DATASET,
) -> dict[str, object]:
    variant_id = gnomad_variant_id(row)
    base: dict[str, object] = {
        "variant_key": row["variant_key"],
        "variant_id": variant_id,
        "population_query_status": "ok",
    }

    try:
        payload = post_gnomad_graphql(
            session,
            GNOMAD_POPULATION_QUERY,
            {"variantId": variant_id, "dataset": dataset},
        )
    except Exception as exc:
        base["population_query_status"] = f"error:{exc}"
        return base

    variant = payload.get("data", {}).get("variant")
    exome = variant.get("exome") if variant else None
    if not exome:
        base["population_query_status"] = "missing_exome"
        return base

    base["global_af"] = to_float(exome.get("af"))
    base["global_ac"] = to_int(exome.get("ac"))
    base["global_an"] = to_int(exome.get("an"))
    faf95 = exome.get("faf95") or {}
    faf99 = exome.get("faf99") or {}
    base["faf95_popmax"] = to_float(faf95.get("popmax"))
    base["faf95_popmax_population"] = faf95.get("popmax_population") or ""
    base["faf99_popmax"] = to_float(faf99.get("popmax"))
    base["faf99_popmax_population"] = faf99.get("popmax_population") or ""

    populations = {
        str(item.get("id", "")).lower(): item
        for item in (exome.get("populations") or [])
    }
    popmax_population = ""
    popmax_af = math.nan
    popmax_ac = math.nan
    for population in POPULATIONS:
        item = populations.get(population, {})
        ac = to_int(item.get("ac"))
        an = to_int(item.get("an"))
        af = safe_af(ac, an)
        base[f"{population}_ac"] = ac
        base[f"{population}_an"] = an
        base[f"{population}_af"] = af
        if not math.isnan(af) and (math.isnan(popmax_af) or af > popmax_af):
            popmax_population = population
            popmax_af = af
            popmax_ac = float(ac or 0)

    if not math.isnan(popmax_af) and popmax_af == 0:
        popmax_population = ""
    base["popmax_population"] = popmax_population
    base["popmax_af"] = popmax_af
    base["popmax_ac"] = popmax_ac
    return base


def load_or_fetch_population_af(
    annotated: pd.DataFrame,
    output_path: Path,
    dataset: str = GNOMAD_DATASET,
    pause: float = 0.25,
    fetch: bool = True,
    force: bool = False,
) -> pd.DataFrame:
    exact = annotated[
        (annotated["match_category"] == "exact_match") & annotated["gnomad_af"].notna()
    ].copy()
    existing = pd.DataFrame()
    if output_path.exists() and not force:
        existing = pd.read_csv(output_path)
        log.info("Loaded existing population AF cache: %s", output_path)

    has_population_columns = any(f"{population}_af" in existing.columns for population in POPULATIONS)
    if not has_population_columns:
        existing_keys: set[str] = set()
    else:
        usable_existing = existing.copy()
        if "population_query_status" in usable_existing.columns:
            usable_existing = usable_existing[
                usable_existing["population_query_status"].astype(str).str.startswith("ok")
                | usable_existing["population_query_status"].isna()
            ]
        existing_keys = set(usable_existing.get("variant_key", pd.Series(dtype=str)).astype(str))
    missing = exact.loc[~exact["variant_key"].astype(str).isin(existing_keys)].copy()

    records: list[dict[str, object]] = []
    if fetch and not missing.empty:
        with requests.Session() as session:
            session.headers.update({"User-Agent": "advanced-variant-analyses/1.0"})
            for _, row in tqdm(missing.iterrows(), total=len(missing), desc="gnomAD population AF"):
                records.append(fetch_population_af_record(session, row, dataset=dataset))
                time.sleep(pause)
    elif not fetch and existing.empty:
        log.warning("Population AF fetch disabled and no cache exists; population analyses will be sparse.")
        records = [
            {
                "variant_key": row["variant_key"],
                "variant_id": gnomad_variant_id(row),
                "population_query_status": "not_fetched",
                "global_af": row["gnomad_af"],
                "global_ac": row["gnomad_ac"],
                "global_an": row["gnomad_an"],
            }
            for _, row in exact.iterrows()
        ]

    fetched = pd.DataFrame(records)
    if existing.empty:
        combined = fetched
    elif fetched.empty:
        combined = existing
    else:
        combined = pd.concat([existing, fetched], ignore_index=True)

    if not combined.empty:
        combined = combined.drop_duplicates(subset="variant_key", keep="last").reset_index(drop=True)
    save_table(combined, output_path)
    return combined


def fetch_exome_genome_af_record(
    session: requests.Session,
    row: pd.Series,
    dataset: str = GNOMAD_DATASET,
) -> dict[str, object]:
    variant_id = gnomad_variant_id(row)
    base: dict[str, object] = {
        "variant_key": row["variant_key"],
        "variant_id": variant_id,
        "gene": row.get("gene", ""),
        "clinvar_id": row.get("clinvar_id", ""),
        "variation_id": row.get("variation_id", ""),
        "title": row.get("title", ""),
        "review_status": row.get("review_status", ""),
        "variant_type": row.get("variant_type", ""),
        "functional_class": row.get("functional_class", ""),
        "exome_genome_query_status": "ok",
    }
    try:
        payload = post_gnomad_graphql(
            session,
            GNOMAD_EXOME_GENOME_QUERY,
            {"variantId": variant_id, "dataset": dataset},
        )
        variant = payload.get("data", {}).get("variant")
    except Exception as exc:
        base["exome_genome_query_status"] = f"error:{exc}"
        variant = None

    if not variant:
        base["exome_genome_query_status"] = (
            base["exome_genome_query_status"]
            if str(base["exome_genome_query_status"]).startswith("error:")
            else "no_gnomad_record"
        )
        return base

    for source in ("exome", "genome"):
        block = variant.get(source) or {}
        base[f"{source}_af"] = to_float(block.get("af"))
        base[f"{source}_ac"] = to_int(block.get("ac"))
        base[f"{source}_an"] = to_int(block.get("an"))
    return base


def load_or_fetch_exome_genome_af(
    annotated: pd.DataFrame,
    output_path: Path,
    dataset: str = GNOMAD_DATASET,
    pause: float = 0.25,
    fetch: bool = True,
    force: bool = False,
) -> pd.DataFrame:
    exact = annotated[annotated["match_category"] == "exact_match"].copy()
    existing = pd.DataFrame()
    if output_path.exists() and not force:
        existing = pd.read_csv(output_path)
        log.info("Loaded existing exome/genome AF cache: %s", output_path)

    records: list[dict[str, object]] = []
    if fetch:
        usable_existing = existing.copy()
        if "exome_genome_query_status" in usable_existing.columns:
            usable_existing = usable_existing[
                usable_existing["exome_genome_query_status"].fillna("").astype(str).eq("ok")
            ]
        existing_keys = set(usable_existing.get("variant_key", pd.Series(dtype="object")).astype(str))
        missing = exact.loc[~exact["variant_key"].astype(str).isin(existing_keys)].copy()
        if not missing.empty:
            with requests.Session() as session:
                session.headers.update({"User-Agent": "advanced-variant-analyses/1.0"})
                for _, row in tqdm(missing.iterrows(), total=len(missing), desc="gnomAD exome/genome AF"):
                    records.append(fetch_exome_genome_af_record(session, row, dataset=dataset))
                    time.sleep(pause)
    elif existing.empty:
        log.warning("Exome/genome AF fetch disabled and no cache exists; genome sensitivity will be sparse.")
        records = [
            {
                "variant_key": row["variant_key"],
                "variant_id": gnomad_variant_id(row),
                "gene": row.get("gene", ""),
                "clinvar_id": row.get("clinvar_id", ""),
                "variation_id": row.get("variation_id", ""),
                "title": row.get("title", ""),
                "review_status": row.get("review_status", ""),
                "variant_type": row.get("variant_type", ""),
                "functional_class": row.get("functional_class", ""),
                "exome_genome_query_status": "not_fetched",
                "exome_af": row.get("gnomad_af", np.nan),
                "exome_ac": row.get("gnomad_ac", np.nan),
                "exome_an": row.get("gnomad_an", np.nan),
            }
            for _, row in exact.iterrows()
        ]

    fetched = pd.DataFrame(records)
    if existing.empty:
        combined = fetched
    elif fetched.empty:
        combined = existing
    else:
        combined = pd.concat([existing, fetched], ignore_index=True)

    if not combined.empty:
        combined = combined.drop_duplicates(subset="variant_key", keep="last").reset_index(drop=True)
        for source in ("exome", "genome"):
            for suffix in ("af", "ac", "an"):
                column = f"{source}_{suffix}"
                if column not in combined.columns:
                    combined[column] = np.nan
                combined[column] = pd.to_numeric(combined[column], errors="coerce")
    save_table(combined, output_path)
    return combined


def make_exome_genome_af_summary(comparison: pd.DataFrame) -> pd.DataFrame:
    if comparison.empty:
        return pd.DataFrame(
            columns=[
                "group",
                "label",
                "variant_count",
                "query_ok_count",
                "query_error_count",
                "exome_ac_positive",
                "genome_ac_positive",
                "genome_only_ac_positive",
                "exome_only_ac_positive",
                "both_ac_positive",
                "neither_ac_positive",
                "max_exome_af",
                "max_genome_af",
                "genome_af_gt_1e_5",
                "genome_af_gt_1e_4",
            ]
        )

    df = comparison.copy()
    for column in ["exome_af", "exome_ac", "genome_af", "genome_ac"]:
        if column not in df.columns:
            df[column] = np.nan
        df[column] = pd.to_numeric(df[column], errors="coerce")

    rows: list[dict[str, object]] = []

    def add_summary(group: str, label: str, group_df: pd.DataFrame) -> None:
        exome_positive = group_df["exome_ac"].fillna(0) > 0
        genome_positive = group_df["genome_ac"].fillna(0) > 0
        query_status = group_df.get(
            "exome_genome_query_status",
            pd.Series("", index=group_df.index, dtype="object"),
        ).fillna("").astype(str)
        query_ok = query_status.eq("ok")
        rows.append(
            {
                "group": group,
                "label": label,
                "variant_count": len(group_df),
                "query_ok_count": int(query_ok.sum()),
                "query_error_count": int((~query_ok).sum()),
                "exome_ac_positive": int(exome_positive.sum()),
                "genome_ac_positive": int(genome_positive.sum()),
                "genome_only_ac_positive": int((genome_positive & ~exome_positive).sum()),
                "exome_only_ac_positive": int((exome_positive & ~genome_positive).sum()),
                "both_ac_positive": int((genome_positive & exome_positive).sum()),
                "neither_ac_positive": int((~genome_positive & ~exome_positive).sum()),
                "max_exome_af": float(group_df["exome_af"].max()) if group_df["exome_af"].notna().any() else np.nan,
                "max_genome_af": float(group_df["genome_af"].max()) if group_df["genome_af"].notna().any() else np.nan,
                "genome_af_gt_1e_5": int((group_df["genome_af"] > AF_ULTRA_RARE).sum()),
                "genome_af_gt_1e_4": int((group_df["genome_af"] > AF_RARE).sum()),
            }
        )

    add_summary("overall", "all_exact_matches", df)
    for column in ["variant_type", "functional_class"]:
        if column not in df.columns:
            continue
        for label, group_df in df.groupby(column, dropna=False):
            add_summary(column, str(label), group_df)

    return pd.DataFrame(rows)


def merge_population_fields(annotated: pd.DataFrame, population_df: pd.DataFrame) -> pd.DataFrame:
    if population_df.empty:
        result = annotated.copy()
        result["global_af"] = result["gnomad_af"]
        result["popmax_af"] = np.nan
        result["popmax_population"] = ""
        return result

    allowed_columns = {
        "variant_key",
        "population_query_status",
        "global_af",
        "global_ac",
        "global_an",
        "faf95_popmax",
        "faf95_popmax_population",
        "faf99_popmax",
        "faf99_popmax_population",
        "popmax_population",
        "popmax_af",
        "popmax_ac",
    }
    allowed_columns.update(
        f"{population}_{suffix}"
        for population in POPULATIONS
        for suffix in ("af", "ac", "an")
    )
    pop_columns = [column for column in population_df.columns if column in allowed_columns]
    result = annotated.merge(population_df.loc[:, pop_columns], on="variant_key", how="left")
    result["global_af"] = pd.to_numeric(result.get("global_af"), errors="coerce")
    result["global_af"] = result["global_af"].combine_first(result["gnomad_af"])
    result["popmax_af"] = pd.to_numeric(result.get("popmax_af"), errors="coerce")
    result["popmax_population"] = result.get(
        "popmax_population", pd.Series("", index=result.index, dtype="object")
    ).fillna("").astype(str)
    return result


def add_popmax_ac_from_population_columns(df: pd.DataFrame) -> pd.DataFrame:
    result = df.copy()
    if "popmax_ac" not in result.columns:
        result["popmax_ac"] = np.nan
    result["popmax_ac"] = pd.to_numeric(result["popmax_ac"], errors="coerce")

    missing = result["popmax_ac"].isna()
    if missing.any() and "popmax_population" in result.columns:
        for population in POPULATIONS:
            column = f"{population}_ac"
            if column not in result.columns:
                continue
            mask = missing & result["popmax_population"].fillna("").astype(str).str.lower().eq(population)
            result.loc[mask, "popmax_ac"] = pd.to_numeric(result.loc[mask, column], errors="coerce")
    return result


def make_population_tables(annotated: pd.DataFrame, population_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    wide = merge_population_fields(
        annotated[(annotated["match_category"] == "exact_match") & annotated["gnomad_af"].notna()],
        population_df,
    )
    wide = add_popmax_ac_from_population_columns(wide)
    for population in POPULATIONS:
        for suffix in ("af", "ac", "an"):
            column = f"{population}_{suffix}"
            if column not in wide.columns:
                wide[column] = np.nan
        wide[f"{population}_af"] = pd.to_numeric(wide[f"{population}_af"], errors="coerce")

    identity_columns = [
        "variant_key",
        "gene",
        "clinvar_id",
        "variation_id",
        "title",
        "review_status",
        "review_strength",
        "submitter_count",
        "submitter_count_source",
        "variant_type",
        "functional_class",
        "global_af",
        "global_ac",
        "global_an",
        "popmax_population",
        "popmax_af",
        "popmax_ac",
        "faf95_popmax",
        "faf95_popmax_population",
        "faf99_popmax",
        "faf99_popmax_population",
    ]
    identity_columns = [column for column in identity_columns if column in wide.columns]
    pop_columns = [
        f"{population}_{suffix}"
        for population in POPULATIONS
        for suffix in ("af", "ac", "an")
        if f"{population}_{suffix}" in wide.columns
    ]
    wide_table = wide.loc[:, [*identity_columns, *pop_columns]].copy()

    long_records: list[dict[str, object]] = []
    for _, row in wide.iterrows():
        for population in POPULATIONS:
            long_records.append(
                {
                    "variant_key": row["variant_key"],
                    "gene": row.get("gene", ""),
                    "population": population,
                    "population_label": POPULATION_LABELS.get(population, population.upper()),
                    "population_af": row.get(f"{population}_af", np.nan),
                    "population_ac": row.get(f"{population}_ac", np.nan),
                    "population_an": row.get(f"{population}_an", np.nan),
                    "global_af": row.get("global_af", np.nan),
                    "global_ac": row.get("global_ac", np.nan),
                    "popmax_population": row.get("popmax_population", ""),
                    "popmax_af": row.get("popmax_af", np.nan),
                    "popmax_ac": row.get("popmax_ac", np.nan),
                    "functional_class": row.get("functional_class", ""),
                    "review_strength": row.get("review_strength", ""),
                }
            )
    long_table = pd.DataFrame(long_records)
    if not long_table.empty:
        long_table["population_af"] = pd.to_numeric(long_table["population_af"], errors="coerce")
        long_table["global_af"] = pd.to_numeric(long_table["global_af"], errors="coerce")

    summary_rows: list[dict[str, object]] = []
    for population in POPULATIONS:
        af = pd.to_numeric(wide[f"{population}_af"], errors="coerce")
        ac = pd.to_numeric(wide[f"{population}_ac"], errors="coerce")
        an = pd.to_numeric(wide[f"{population}_an"], errors="coerce")
        nonzero = af[af > 0]
        summary_rows.append(
            {
                "population": population,
                "population_label": POPULATION_LABELS.get(population, population.upper()),
                "n_variants_with_an": int(an.notna().sum()),
                "n_variants_ac_positive": int((ac.fillna(0) > 0).sum()),
                "fraction_ac_positive": float((ac.fillna(0) > 0).mean()) if len(ac) else np.nan,
                "median_af_nonzero": float(nonzero.median()) if not nonzero.empty else np.nan,
                "max_af": float(af.max()) if af.notna().any() else np.nan,
                "n_af_gt_1e_5": int((af > AF_ULTRA_RARE).sum()),
                "n_af_gt_1e_4": int((af > AF_RARE).sum()),
            }
        )
    summary = pd.DataFrame(summary_rows)

    outlier_source = wide.copy()
    outlier_source["global_af"] = pd.to_numeric(outlier_source["global_af"], errors="coerce")
    outlier_source["popmax_af"] = pd.to_numeric(outlier_source["popmax_af"], errors="coerce")
    outlier_source["popmax_global_ratio"] = (
        (outlier_source["popmax_af"] + LIFT_PSEUDOCOUNT)
        / (outlier_source["global_af"] + LIFT_PSEUDOCOUNT)
    )
    outlier_source["log10_popmax_lift"] = np.log10(outlier_source["popmax_global_ratio"])
    candidate_lifts = outlier_source.loc[
        (outlier_source["global_af"] > 0) & (outlier_source["popmax_af"] > 0),
        "log10_popmax_lift",
    ].dropna()
    if len(candidate_lifts) >= 4:
        q1 = candidate_lifts.quantile(0.25)
        q3 = candidate_lifts.quantile(0.75)
        lift_threshold = max(math.log10(3), q3 + 1.5 * (q3 - q1))
    else:
        lift_threshold = math.log10(3)
    outlier_source["is_robust_lift_outlier"] = outlier_source["log10_popmax_lift"] >= lift_threshold
    outlier_source["rare_global_pop_enriched"] = (
        (outlier_source["global_af"] <= AF_RARE)
        & (outlier_source["popmax_af"] >= AF_ULTRA_RARE)
        & (outlier_source["popmax_global_ratio"] >= 3)
    )
    outliers = outlier_source.loc[
        outlier_source["rare_global_pop_enriched"]
        | (
            (outlier_source["global_af"] <= AF_ULTRA_RARE)
            & (outlier_source["popmax_af"] > AF_ULTRA_RARE)
        )
    ].copy()
    outliers["popmax_population_label"] = outliers["popmax_population"].map(
        lambda value: POPULATION_LABELS.get(str(value), str(value).upper())
    )
    outlier_columns = [
        "variant_key",
        "gene",
        "title",
        "review_status",
        "submitter_count",
        "variant_type",
        "functional_class",
        "global_af",
        "global_ac",
        "popmax_population",
        "popmax_population_label",
        "popmax_af",
        "popmax_ac",
        "popmax_global_ratio",
        "log10_popmax_lift",
        "is_robust_lift_outlier",
        "rare_global_pop_enriched",
    ]
    outliers = outliers.loc[:, [column for column in outlier_columns if column in outliers.columns]]
    outliers = outliers.sort_values(
        ["popmax_af", "popmax_global_ratio"], ascending=[False, False]
    ).reset_index(drop=True)
    return wide_table, long_table, summary, outliers


def make_reclassification_risk_table(annotated: pd.DataFrame, population_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    exact = merge_population_fields(
        annotated[(annotated["match_category"] == "exact_match") & annotated["gnomad_af"].notna()],
        population_df,
    )
    exact = add_popmax_ac_from_population_columns(exact)
    for column, fallback_column in [("global_ac", "gnomad_ac"), ("global_an", "gnomad_an")]:
        if column not in exact.columns:
            exact[column] = exact.get(fallback_column, np.nan)
        exact[column] = pd.to_numeric(exact[column], errors="coerce")
    exact["popmax_ac"] = pd.to_numeric(exact["popmax_ac"], errors="coerce")
    exact["global_af"] = pd.to_numeric(exact["global_af"], errors="coerce")
    exact["popmax_af"] = pd.to_numeric(exact["popmax_af"], errors="coerce")
    exact["max_frequency_signal"] = exact[["global_af", "popmax_af"]].max(axis=1, skipna=True)
    popmax_is_max = exact["popmax_af"].notna() & (
        exact["global_af"].isna() | exact["popmax_af"].ge(exact["global_af"])
    )
    exact["max_frequency_source"] = np.where(popmax_is_max, "popmax", "global")
    exact.loc[exact["max_frequency_signal"].isna(), "max_frequency_source"] = ""
    exact["max_frequency_ac"] = exact["global_ac"]
    exact.loc[popmax_is_max, "max_frequency_ac"] = exact.loc[popmax_is_max, "popmax_ac"]
    exact["global_frequency_signal_ac_ge_20"] = (
        (exact["global_af"] > AF_ULTRA_RARE) & (exact["global_ac"] >= MIN_RECLASSIFICATION_AC)
    )
    exact["popmax_frequency_signal_ac_ge_20"] = (
        (exact["popmax_af"] > AF_ULTRA_RARE) & (exact["popmax_ac"] >= MIN_RECLASSIFICATION_AC)
    )
    exact["frequency_signal_ac_ge_20"] = (
        exact["global_frequency_signal_ac_ge_20"] | exact["popmax_frequency_signal_ac_ge_20"]
    )
    exact["frequency_signal_ac_filter"] = np.where(
        exact["frequency_signal_ac_ge_20"],
        f"passes_AC_ge_{MIN_RECLASSIFICATION_AC}",
        f"below_AC_{MIN_RECLASSIFICATION_AC}_or_missing",
    )
    exact["weak_review_signal"] = (exact["review_score"] <= 1) | (exact["submitter_count"].fillna(99) <= 1)
    exact["frequency_threshold"] = np.select(
        [
            exact["max_frequency_signal"] > AF_RARE,
            exact["max_frequency_signal"] > AF_ULTRA_RARE,
        ],
        [
            "AF_gt_1e_4",
            "AF_gt_1e_5",
        ],
        default="below_1e_5",
    )
    exact["potential_misclassification_signal"] = (
        (exact["max_frequency_signal"] > AF_ULTRA_RARE)
        & exact["weak_review_signal"]
        & exact["frequency_signal_ac_ge_20"]
    )
    exact["risk_level"] = np.select(
        [
            (exact["max_frequency_signal"] > AF_RARE) & exact["weak_review_signal"],
            (exact["max_frequency_signal"] > AF_ULTRA_RARE) & exact["weak_review_signal"],
            exact["max_frequency_signal"] > AF_RARE,
            exact["max_frequency_signal"] > AF_ULTRA_RARE,
        ],
        [
            "high_frequency_weak_review",
            "moderate_frequency_weak_review",
            "high_frequency_stronger_review",
            "moderate_frequency_stronger_review",
        ],
        default="low_frequency",
    )
    exact["signal_reason"] = exact.apply(
        lambda row: "; ".join(
            reason
            for reason in [
                row["frequency_threshold"] if row["frequency_threshold"] != "below_1e_5" else "",
                "weak_review" if row["weak_review_signal"] else "",
                (
                    f"popmax_{row['popmax_population']}"
                    if pd.notna(row.get("popmax_af")) and row.get("popmax_population")
                    else ""
                ),
                row["frequency_signal_ac_filter"],
            ]
            if reason
        ),
        axis=1,
    )
    table = exact.loc[exact["max_frequency_signal"] > AF_ULTRA_RARE].copy()
    risk_order = {
        "high_frequency_weak_review": 0,
        "moderate_frequency_weak_review": 1,
        "high_frequency_stronger_review": 2,
        "moderate_frequency_stronger_review": 3,
        "low_frequency": 4,
    }
    table["_risk_order"] = table["risk_level"].map(risk_order).fillna(9)
    columns = [
        "variant_key",
        "gene",
        "clinvar_id",
        "variation_id",
        "title",
        "clinsig",
        "review_status",
        "review_strength",
        "review_score",
        "submitter_count",
        "submitter_count_source",
        "variant_type",
        "functional_class",
        "global_af",
        "global_ac",
        "global_an",
        "popmax_af",
        "popmax_ac",
        "popmax_population",
        "max_frequency_signal",
        "max_frequency_source",
        "max_frequency_ac",
        "frequency_threshold",
        "frequency_signal_ac_filter",
        "frequency_signal_ac_ge_20",
        "weak_review_signal",
        "potential_misclassification_signal",
        "risk_level",
        "signal_reason",
    ]
    table = table.loc[:, [column for column in columns if column in table.columns]]
    table = table.sort_values(
        ["risk_level", "max_frequency_signal"],
        key=lambda series: series.map(risk_order).fillna(series) if series.name == "risk_level" else series,
        ascending=[True, False],
    ).reset_index(drop=True)

    summary = (
        table.groupby(["risk_level", "review_strength", "frequency_signal_ac_filter"], dropna=False)
        .size()
        .rename("variant_count")
        .reset_index()
        .sort_values(["risk_level", "frequency_signal_ac_filter", "variant_count"], ascending=[True, True, False])
    )
    return table, summary


def make_gene_variant_type_tables(annotated: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = (
        annotated.groupby(["gene", "functional_class"], dropna=False)
        .size()
        .rename("variant_count")
        .reset_index()
    )
    totals = counts.groupby("gene")["variant_count"].transform("sum")
    counts["fraction_within_gene"] = counts["variant_count"] / totals.where(totals > 0, other=1)
    counts = counts.sort_values(["gene", "variant_count"], ascending=[True, False]).reset_index(drop=True)

    af_df = annotated[
        (annotated["match_category"] == "exact_match") & annotated["gnomad_af"].notna()
    ].copy()
    if af_df.empty:
        af_summary = pd.DataFrame(
            columns=[
                "gene",
                "functional_class",
                "n_af_covered",
                "median_af",
                "max_af",
                "n_af_gt_1e_5",
                "n_af_gt_1e_4",
            ]
        )
    else:
        af_summary = (
            af_df.groupby(["gene", "functional_class"], dropna=False)
            .agg(
                n_af_covered=("variant_key", "count"),
                median_af=("gnomad_af", "median"),
                max_af=("gnomad_af", "max"),
                n_af_gt_1e_5=("gnomad_af", lambda values: int((values > AF_ULTRA_RARE).sum())),
                n_af_gt_1e_4=("gnomad_af", lambda values: int((values > AF_RARE).sum())),
            )
            .reset_index()
            .sort_values(["gene", "n_af_covered"], ascending=[True, False])
        )
    return counts, af_summary


def fisher_feature_row(df: pd.DataFrame, feature: str) -> dict[str, object]:
    exact = df[df["overlap_group"] == "exact_match"]
    non = df[df["overlap_group"] == "non_overlap"]
    exact_pos = int(exact[feature].fillna(False).astype(bool).sum())
    non_pos = int(non[feature].fillna(False).astype(bool).sum())
    exact_total = len(exact)
    non_total = len(non)
    table = [[non_pos, non_total - non_pos], [exact_pos, exact_total - exact_pos]]
    odds, p_value = fisher_exact(table, alternative="two-sided") if exact_total and non_total else (np.nan, np.nan)
    return {
        "feature": feature,
        "exact_positive": exact_pos,
        "exact_total": exact_total,
        "exact_fraction": exact_pos / exact_total if exact_total else np.nan,
        "non_overlap_positive": non_pos,
        "non_overlap_total": non_total,
        "non_overlap_fraction": non_pos / non_total if non_total else np.nan,
        "odds_ratio_non_overlap_vs_exact": odds,
        "fisher_p_two_sided": p_value,
    }


def make_non_overlap_tables(
    annotated: pd.DataFrame,
    coverage_df: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    analysis = annotated[
        annotated["match_category"].isin(["exact_match", "allele_discordance", "no_gnomad_record"])
    ].copy()
    analysis["overlap_group"] = np.where(
        analysis["match_category"] == "exact_match",
        "exact_match",
        "non_overlap",
    )
    analysis["feature_LOF"] = analysis["functional_class"].eq("LOF")
    analysis["feature_missense"] = analysis["functional_class"].eq("missense")
    analysis["feature_splice_or_intronic"] = analysis["functional_class"].eq("splice_or_intronic")

    feature_columns = [
        "is_snv",
        "is_indel",
        "is_cpg_transition_proxy",
        "is_hotspot_position",
        "is_local_hotspot_5bp",
        "feature_LOF",
        "feature_missense",
        "feature_splice_or_intronic",
    ]
    feature_summary = pd.DataFrame([fisher_feature_row(analysis, feature) for feature in feature_columns])
    if feature_summary["fisher_p_two_sided"].notna().any():
        valid = feature_summary["fisher_p_two_sided"].notna()
        _, q_values, _, _ = multipletests(
            feature_summary.loc[valid, "fisher_p_two_sided"],
            method="fdr_bh",
        )
        feature_summary.loc[valid, "bh_q"] = q_values
    else:
        feature_summary["bh_q"] = np.nan

    gene_rows: list[dict[str, object]] = []
    exact_total = int((analysis["overlap_group"] == "exact_match").sum())
    non_total = int((analysis["overlap_group"] == "non_overlap").sum())
    for gene, group in analysis.groupby("gene", dropna=False):
        exact_gene = int(((group["overlap_group"] == "exact_match")).sum())
        non_gene = int(((group["overlap_group"] == "non_overlap")).sum())
        contingency = [
            [non_gene, non_total - non_gene],
            [exact_gene, exact_total - exact_gene],
        ]
        odds, p_value = (
            fisher_exact(contingency, alternative="greater")
            if exact_total and non_total
            else (np.nan, np.nan)
        )
        allele_discordance = int((group["match_category"] == "allele_discordance").sum())
        no_record = int((group["match_category"] == "no_gnomad_record").sum())
        gene_rows.append(
            {
                "gene": gene,
                "exact_match_count": exact_gene,
                "non_overlap_count": non_gene,
                "allele_discordance_count": allele_discordance,
                "no_gnomad_record_count": no_record,
                "non_overlap_fraction_within_gene": non_gene / len(group) if len(group) else np.nan,
                "odds_ratio_non_overlap_enrichment": odds,
                "fisher_p_greater": p_value,
            }
        )
    gene_enrichment = pd.DataFrame(gene_rows)
    if not gene_enrichment.empty and gene_enrichment["fisher_p_greater"].notna().any():
        valid = gene_enrichment["fisher_p_greater"].notna()
        _, q_values, _, _ = multipletests(
            gene_enrichment.loc[valid, "fisher_p_greater"],
            method="fdr_bh",
        )
        gene_enrichment.loc[valid, "bh_q"] = q_values
    gene_enrichment = gene_enrichment.sort_values(
        ["non_overlap_count", "non_overlap_fraction_within_gene"],
        ascending=[False, False],
    ).reset_index(drop=True)

    gene_distribution = (
        analysis.groupby(["gene", "match_category"], dropna=False)
        .size()
        .rename("variant_count")
        .reset_index()
        .sort_values(["gene", "variant_count"], ascending=[True, False])
    )
    variant_columns = [
        "variant_key",
        "gene",
        "match_category",
        "title",
        "review_status",
        "variant_type",
        "functional_class",
        "is_snv",
        "is_cpg_transition_proxy",
        "is_hotspot_position",
        "is_local_hotspot_5bp",
        "position_variant_count",
        "nearby_5bp_variant_count",
        "gnomad_variants_at_pos_for_match",
        "gnomad_region_query_status",
    ]
    variant_level = analysis.loc[
        analysis["overlap_group"] == "non_overlap",
        [column for column in variant_columns if column in analysis.columns],
    ].copy()

    if coverage_df is None or coverage_df.empty:
        coverage_summary = pd.DataFrame(
            columns=["match_category", "has_any_gnomad", "count", "fraction"]
        )
    else:
        cov = coverage_df.copy()
        if "has_any_gnomad" in cov.columns:
            cov["has_any_gnomad"] = parse_bool_series(cov["has_any_gnomad"])
            coverage_summary = (
                cov.groupby(["match_category", "has_any_gnomad"], dropna=False)
                .size()
                .rename("count")
                .reset_index()
            )
            totals = coverage_summary.groupby("match_category")["count"].transform("sum")
            coverage_summary["fraction"] = coverage_summary["count"] / totals.where(totals > 0, other=1)
        else:
            coverage_summary = pd.DataFrame(
                columns=["match_category", "has_any_gnomad", "count", "fraction"]
            )

    return feature_summary, gene_enrichment, gene_distribution, variant_level, coverage_summary


def plot_population_distribution(long_table: pd.DataFrame, output_path: Path) -> None:
    plot_df = long_table[long_table["population_af"] > 0].copy()
    if plot_df.empty:
        log.warning("Skipping population AF figure because all population AF values are zero/missing.")
        return
    plot_df["log10_population_af"] = np.log10(plot_df["population_af"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 5))
    sns.boxplot(
        data=plot_df,
        x="population_label",
        y="log10_population_af",
        order=[POPULATION_LABELS[p] for p in POPULATIONS],
        color="#93C5FD",
        fliersize=0,
        ax=ax,
    )
    sns.stripplot(
        data=plot_df,
        x="population_label",
        y="log10_population_af",
        order=[POPULATION_LABELS[p] for p in POPULATIONS],
        color="#1D4ED8",
        size=2,
        alpha=0.35,
        ax=ax,
    )
    ax.axhline(np.log10(AF_ULTRA_RARE), color="#DC2626", linestyle="--", linewidth=1)
    ax.axhline(np.log10(AF_RARE), color="#F97316", linestyle="--", linewidth=1)
    ax.set_xlabel("gnomAD genetic ancestry group")
    ax.set_ylabel("log10(population AF), nonzero values")
    ax.set_title("Population-specific AF distribution for exact gnomAD matches")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    log.info("Saved: %s", output_path)


def plot_population_outliers(outliers: pd.DataFrame, output_path: Path) -> None:
    if outliers.empty:
        log.warning("Skipping population outlier figure because no candidates were found.")
        return
    plot_df = outliers.head(15).copy()
    plot_df["label"] = plot_df["gene"].astype(str) + " " + plot_df["variant_key"].astype(str)
    plot_df["popmax_population_label"] = plot_df["popmax_population"].map(
        lambda value: POPULATION_LABELS.get(str(value), str(value).upper())
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, max(4.5, 0.35 * len(plot_df))))
    sns.barplot(
        data=plot_df,
        y="label",
        x="popmax_af",
        hue="popmax_population_label",
        dodge=False,
        ax=ax,
    )
    ax.axvline(AF_ULTRA_RARE, color="#DC2626", linestyle="--", linewidth=1)
    ax.set_xscale("log")
    ax.set_xlabel("Popmax AF")
    ax.set_ylabel("")
    ax.set_title("Globally rare variants enriched in one gnomAD population")
    ax.legend(title="Popmax", frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    log.info("Saved: %s", output_path)


def plot_reclassification_risk(table: pd.DataFrame, output_path: Path) -> None:
    if table.empty:
        log.warning("Skipping reclassification-risk figure because table is empty.")
        return
    counts = table["risk_level"].value_counts().reset_index()
    counts.columns = ["risk_level", "count"]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    sns.barplot(data=counts, x="count", y="risk_level", color="#EF4444", ax=ax)
    ax.set_xlabel("Variant count")
    ax.set_ylabel("")
    ax.set_title("High-frequency ClinVar P/LP variants by review-risk signal")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    log.info("Saved: %s", output_path)


def plot_variant_type_by_gene(counts: pd.DataFrame, output_path: Path) -> None:
    if counts.empty:
        log.warning("Skipping variant-type figure because counts are empty.")
        return
    top_genes = (
        counts.groupby("gene")["variant_count"].sum().sort_values(ascending=False).head(20).index
    )
    plot_df = counts[counts["gene"].isin(top_genes)].copy()
    pivot = (
        plot_df.pivot_table(
            index="gene",
            columns="functional_class",
            values="variant_count",
            aggfunc="sum",
            fill_value=0,
        )
        .loc[top_genes]
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ax = pivot.plot(
        kind="barh",
        stacked=True,
        figsize=(9, max(5, 0.35 * len(pivot))),
        colormap="tab20",
        edgecolor="white",
    )
    ax.set_xlabel("ClinVar P/LP variant count")
    ax.set_ylabel("")
    ax.set_title("Functional class composition by gene")
    ax.invert_yaxis()
    ax.legend(title="Functional class", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close()
    log.info("Saved: %s", output_path)


def plot_non_overlap_features(feature_summary: pd.DataFrame, output_path: Path) -> None:
    if feature_summary.empty:
        log.warning("Skipping non-overlap feature figure because summary is empty.")
        return
    plot_df = feature_summary.melt(
        id_vars=["feature"],
        value_vars=["exact_fraction", "non_overlap_fraction"],
        var_name="group",
        value_name="fraction",
    )
    plot_df["group"] = plot_df["group"].map(
        {"exact_fraction": "Exact match", "non_overlap_fraction": "Non-overlap"}
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 5))
    sns.barplot(data=plot_df, y="feature", x="fraction", hue="group", ax=ax)
    ax.set_xlabel("Fraction of variants")
    ax.set_ylabel("")
    ax.set_title("Feature structure in exact matches vs non-overlap variants")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    log.info("Saved: %s", output_path)


def plot_non_overlap_genes(gene_enrichment: pd.DataFrame, output_path: Path) -> None:
    if gene_enrichment.empty:
        log.warning("Skipping non-overlap gene figure because enrichment table is empty.")
        return
    plot_df = gene_enrichment.head(15).copy()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, max(4.5, 0.32 * len(plot_df))))
    sns.barplot(data=plot_df, y="gene", x="non_overlap_count", color="#F59E0B", ax=ax)
    ax.set_xlabel("Non-overlap variant count")
    ax.set_ylabel("")
    ax.set_title("Genes contributing most non-overlap variants")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    log.info("Saved: %s", output_path)


def load_optional_coverage(prefix: str) -> pd.DataFrame:
    path = data_path(prefix, "coverage_analysis.csv")
    if not path.exists():
        log.warning("Coverage cache not found: %s", path)
        return pd.DataFrame()
    return pd.read_csv(path)


def run_advanced_analyses(
    output_prefix: str | None = None,
    dataset: str = GNOMAD_DATASET,
    gnomad_pause: float = 0.25,
    fetch_population_af: bool = True,
    force_population_fetch: bool = False,
    fetch_exome_genome_af: bool = True,
    force_exome_genome_fetch: bool = False,
    make_plots: bool = True,
) -> dict[str, Path]:
    prefix = resolve_output_prefix(output_prefix)
    matched_path = data_path(prefix, "gnomad_matched.csv")
    if not matched_path.exists():
        raise FileNotFoundError(f"Matched gnomAD file not found: {matched_path}")

    matched = pd.read_csv(matched_path)
    annotated = prepare_annotations(matched)
    coverage = load_optional_coverage(prefix)

    output_paths = {
        "population_af": data_path(prefix, "population_af.csv"),
        "population_af_long": data_path(prefix, "population_af_long.csv"),
        "population_af_summary": data_path(prefix, "population_af_summary.csv"),
        "population_af_outliers": data_path(prefix, "population_af_outliers.csv"),
        "exome_genome_af_comparison": data_path(prefix, "exome_genome_af_comparison.csv"),
        "exome_genome_af_summary": data_path(prefix, "exome_genome_af_summary.csv"),
        "reclassification_risk": data_path(prefix, "reclassification_risk.csv"),
        "reclassification_risk_summary": data_path(prefix, "reclassification_risk_summary.csv"),
        "gene_variant_type_summary": data_path(prefix, "gene_variant_type_summary.csv"),
        "gene_lof_missense_af_summary": data_path(prefix, "gene_lof_missense_af_summary.csv"),
        "non_overlap_feature_summary": data_path(prefix, "non_overlap_feature_summary.csv"),
        "non_overlap_gene_enrichment": data_path(prefix, "non_overlap_gene_enrichment.csv"),
        "non_overlap_gene_distribution": data_path(prefix, "non_overlap_gene_distribution.csv"),
        "non_overlap_variant_level": data_path(prefix, "non_overlap_variant_level.csv"),
        "non_overlap_coverage_summary": data_path(prefix, "non_overlap_coverage_summary.csv"),
    }

    exome_genome_df = load_or_fetch_exome_genome_af(
        annotated,
        output_paths["exome_genome_af_comparison"],
        dataset=dataset,
        pause=gnomad_pause,
        fetch=fetch_exome_genome_af,
        force=force_exome_genome_fetch,
    )
    exome_genome_summary = make_exome_genome_af_summary(exome_genome_df)
    save_table(exome_genome_df, output_paths["exome_genome_af_comparison"])
    save_table(exome_genome_summary, output_paths["exome_genome_af_summary"])

    population_df = load_or_fetch_population_af(
        annotated,
        output_paths["population_af"],
        dataset=dataset,
        pause=gnomad_pause,
        fetch=fetch_population_af,
        force=force_population_fetch,
    )
    population_wide, population_long, population_summary, population_outliers = make_population_tables(
        annotated, population_df
    )
    save_table(population_wide, output_paths["population_af"])
    save_table(population_long, output_paths["population_af_long"])
    save_table(population_summary, output_paths["population_af_summary"])
    save_table(population_outliers, output_paths["population_af_outliers"])

    risk_table, risk_summary = make_reclassification_risk_table(annotated, population_df)
    save_table(risk_table, output_paths["reclassification_risk"])
    save_table(risk_summary, output_paths["reclassification_risk_summary"])

    gene_type_counts, gene_af_summary = make_gene_variant_type_tables(annotated)
    save_table(gene_type_counts, output_paths["gene_variant_type_summary"])
    save_table(gene_af_summary, output_paths["gene_lof_missense_af_summary"])

    (
        non_overlap_features,
        non_overlap_genes,
        non_overlap_gene_distribution,
        non_overlap_variant_level,
        non_overlap_coverage,
    ) = make_non_overlap_tables(annotated, coverage)
    save_table(non_overlap_features, output_paths["non_overlap_feature_summary"])
    save_table(non_overlap_genes, output_paths["non_overlap_gene_enrichment"])
    save_table(non_overlap_gene_distribution, output_paths["non_overlap_gene_distribution"])
    save_table(non_overlap_variant_level, output_paths["non_overlap_variant_level"])
    save_table(non_overlap_coverage, output_paths["non_overlap_coverage_summary"])

    if make_plots:
        plot_population_distribution(
            population_long,
            figure_path(prefix, "population_af_distribution.png"),
        )
        plot_population_outliers(
            population_outliers,
            figure_path(prefix, "population_af_outliers.png"),
        )
        plot_reclassification_risk(
            risk_table,
            figure_path(prefix, "reclassification_risk.png"),
        )
        plot_variant_type_by_gene(
            gene_type_counts,
            figure_path(prefix, "variant_type_by_gene.png"),
        )
        plot_non_overlap_features(
            non_overlap_features,
            figure_path(prefix, "non_overlap_features.png"),
        )
        plot_non_overlap_genes(
            non_overlap_genes,
            figure_path(prefix, "non_overlap_genes.png"),
        )

    return output_paths


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run ancestry-stratified, reclassification-risk, gene-feature, "
            "and non-overlap analyses from the API pipeline CSV outputs."
        )
    )
    parser.add_argument(
        "--output-prefix",
        help=(
            "Prefix used by arrhythmia_variant_pipeline.py outputs. If omitted, "
            "auto-detects arrhythmia when arrhythmia_gnomad_matched.csv exists."
        ),
    )
    parser.add_argument("--dataset", default=GNOMAD_DATASET)
    parser.add_argument("--gnomad-pause", type=float, default=0.25)
    parser.add_argument(
        "--no-fetch-population-af",
        action="store_true",
        help="Use an existing population AF cache and do not query gnomAD.",
    )
    parser.add_argument(
        "--force-population-fetch",
        action="store_true",
        help="Ignore the existing population AF cache and query gnomAD again.",
    )
    parser.add_argument(
        "--no-fetch-exome-genome-af",
        action="store_true",
        help="Use an existing exome/genome AF cache and do not query gnomAD.",
    )
    parser.add_argument(
        "--force-exome-genome-fetch",
        action="store_true",
        help="Ignore the existing exome/genome AF cache and query gnomAD again.",
    )
    parser.add_argument("--no-plots", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = run_advanced_analyses(
        output_prefix=args.output_prefix,
        dataset=args.dataset,
        gnomad_pause=args.gnomad_pause,
        fetch_population_af=not args.no_fetch_population_af,
        force_population_fetch=args.force_population_fetch,
        fetch_exome_genome_af=not args.no_fetch_exome_genome_af,
        force_exome_genome_fetch=args.force_exome_genome_fetch,
        make_plots=not args.no_plots,
    )
    print("Advanced analysis outputs")
    for name, path in outputs.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
