from __future__ import annotations

import argparse
import json
import logging
import os
import time
from datetime import datetime, timezone
from http.client import IncompleteRead
from pathlib import Path
from typing import Callable, Iterable
from xml.etree import ElementTree as ET

try:
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import numpy as np
    import pandas as pd
    import requests
    import seaborn as sns
    from Bio import Entrez
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

AF_ULTRA_RARE = 1e-5
AF_RARE = 1e-4

GNOMAD_GRAPHQL = "https://gnomad.broadinstitute.org/api"
GNOMAD_DATASET = "gnomad_r4"
GNOMAD_PAUSE = 0.7
GNOMAD_RETRIES = 4
GNOMAD_RETRY_BASE_DELAY = 1.0
GNOMAD_REGION_RETRY_BASE_DELAY = 1.5
CLINVAR_PAUSE = 0.35
CLINVAR_BATCH_SIZE = 50
CLINVAR_EFETCH_RETRIES = 3
CLINVAR_EFETCH_RETRY_BASE_DELAY = 1.0
CACHE_METADATA_VERSION = 1
PIPELINE_VERSION = "arrhythmia-api-v1"
LOW_AF_COVERED_WARNING = 10
TRDN_AF_COVERED_WARNING = 15
ARTICLE_ALLELE_DISCORDANCE = 849
ARTICLE_NO_GNOMAD_RECORD = 962
MATCH_RATIO_WARNING_MIN_UNMATCHED = 50

DEFAULT_CLINVAR_CACHE = DATA_DIR / "clinvar_variants.csv"
DEFAULT_GNOMAD_CACHE = DATA_DIR / "gnomad_matched.csv"
DEFAULT_COVERAGE_CACHE = DATA_DIR / "coverage_analysis.csv"
DEFAULT_FISHER_PATH = DATA_DIR / "per_gene_fisher.csv"
DEFAULT_QC_SUMMARY_PATH = DATA_DIR / "run_qc_summary.csv"

DEFAULT_AF_FIGURE = FIGURE_DIR / "af_distribution.png"
DEFAULT_MATCH_FIGURE = FIGURE_DIR / "match_categories.png"
DEFAULT_HEATMAP_FIGURE = FIGURE_DIR / "gene_heatmap.png"
DEFAULT_COVERAGE_FIGURE = FIGURE_DIR / "coverage_analysis.png"


GNOMAD_QUERY = """
query VariantQuery($variantId: String!, $dataset: DatasetId!) {
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

GNOMAD_REGION_QUERY = """
query RegionVariants(
  $chrom: String!
  $start: Int!
  $stop: Int!
  $dataset: DatasetId!
  $referenceGenome: ReferenceGenomeId!
) {
  region(
    chrom: $chrom
    start: $start
    stop: $stop
    reference_genome: $referenceGenome
  ) {
    variants(dataset: $dataset) {
      variantId
      pos
      ref
      alt
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
}
"""

PATHOGENIC_CLASSIFICATIONS = {
    "pathogenic",
    "likely pathogenic",
    "pathogenic/likely pathogenic",
    "likely pathogenic/pathogenic",
}


def set_entrez_email(email: str | None) -> None:
    resolved_email = email or os.getenv("ENTREZ_EMAIL")
    if not resolved_email:
        raise SystemExit(
            "ClinVar fetch requires an Entrez email. Pass --email or set ENTREZ_EMAIL."
        )
    Entrez.email = resolved_email


def now_utc_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def canonical_genes(genes: Iterable[str]) -> list[str]:
    return sorted({str(gene).strip().upper() for gene in genes if str(gene).strip()})


def sanitize_output_prefix(value: str | None) -> str:
    if not value:
        return ""
    cleaned = "".join(
        char if char.isalnum() or char in {"-", "_"} else "_" for char in value.strip()
    ).strip("_")
    if not cleaned:
        raise SystemExit("--output-prefix must contain at least one letter or digit.")
    return cleaned


def prefixed_path(path: Path, prefix: str | None) -> Path:
    cleaned = sanitize_output_prefix(prefix)
    if not cleaned:
        return path
    return path.with_name(f"{cleaned}_{path.name}")


def cache_metadata_path(path: Path) -> Path:
    return path.with_name(path.name + ".meta.json")


def make_cache_metadata(
    cache_type: str,
    genes: Iterable[str],
    row_count: int,
    dataset: str | None = None,
    extra: dict[str, object] | None = None,
) -> dict[str, object]:
    metadata: dict[str, object] = {
        "cache_metadata_version": CACHE_METADATA_VERSION,
        "pipeline_version": PIPELINE_VERSION,
        "cache_type": cache_type,
        "generated_at_utc": now_utc_iso(),
        "genes": canonical_genes(genes),
        "row_count": int(row_count),
        "af_ultra_rare": AF_ULTRA_RARE,
        "af_rare": AF_RARE,
    }
    if dataset is not None:
        metadata["gnomad_dataset"] = dataset
    if extra:
        metadata.update(extra)
    return metadata


def write_cache_metadata(path: Path, metadata: dict[str, object]) -> None:
    metadata_path = cache_metadata_path(path)
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")


def load_cache_metadata(path: Path) -> dict[str, object] | None:
    metadata_path = cache_metadata_path(path)
    if not metadata_path.exists():
        return None
    try:
        return json.loads(metadata_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        log.warning("Could not parse cache metadata %s: %s", metadata_path, exc)
        return None


def validate_cache_metadata(
    path: Path,
    expected: dict[str, object],
    strict: bool = False,
) -> None:
    metadata = load_cache_metadata(path)
    if metadata is None:
        message = f"Cache metadata is missing for {path}"
        if strict:
            raise SystemExit(message)
        log.warning(message)
        return

    mismatches = []
    for key, expected_value in expected.items():
        actual_value = metadata.get(key)
        if actual_value != expected_value:
            mismatches.append(f"{key}: expected {expected_value!r}, found {actual_value!r}")

    if not mismatches:
        return

    message = f"Cache metadata mismatch for {path}: " + "; ".join(mismatches)
    if strict:
        raise SystemExit(message)
    log.warning(message)


def normalize_chrom(chrom: object) -> str:
    value = str(chrom).strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    return value


def normalize_classification(value: str | None) -> str:
    return " ".join(str(value or "").strip().lower().replace("_", " ").split())


def is_pathogenic_or_likely_pathogenic(value: str | None) -> bool:
    normalized = normalize_classification(value)
    return normalized in PATHOGENIC_CLASSIFICATIONS


def first_text(root: ET.Element, xpath: str) -> str:
    element = root.find(xpath)
    if element is None or element.text is None:
        return ""
    return element.text.strip()


def extract_clinvar_classification(archive: ET.Element) -> tuple[str, str]:
    description_paths = [
        ".//InterpretedRecord/Classifications/GermlineClassification/Description",
        ".//ClassifiedRecord/Classifications/GermlineClassification/Description",
        ".//GermlineClassification/Description",
    ]
    review_paths = [
        ".//InterpretedRecord/Classifications/GermlineClassification/ReviewStatus",
        ".//ClassifiedRecord/Classifications/GermlineClassification/ReviewStatus",
        ".//GermlineClassification/ReviewStatus",
    ]

    description = ""
    for path in description_paths:
        description = first_text(archive, path)
        if description:
            break

    review_status = ""
    for path in review_paths:
        review_status = first_text(archive, path)
        if review_status:
            break

    return description, review_status


def ordered_unique(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        cleaned = " ".join(str(value or "").strip().split())
        if not cleaned or cleaned in seen:
            continue
        seen.add(cleaned)
        ordered.append(cleaned)
    return ordered


def parse_int_attribute(element: ET.Element, attribute: str) -> int | None:
    value = element.attrib.get(attribute, "")
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def extract_variant_type(archive: ET.Element) -> str:
    archive_type = archive.attrib.get("VariationType", "")
    candidates = [archive_type]
    candidates.extend(item.text or "" for item in archive.findall(".//VariantType"))
    for candidate in ordered_unique(candidates):
        if candidate.lower() != "variation":
            return candidate
    return ordered_unique(candidates)[0] if ordered_unique(candidates) else ""


def extract_molecular_consequences(archive: ET.Element) -> str:
    consequences = [
        item.attrib.get("Type", "")
        for item in archive.findall(".//MolecularConsequence")
    ]
    return ";".join(ordered_unique(consequences))


def extract_protein_changes(archive: ET.Element) -> str:
    changes = [
        item.attrib.get("change", "")
        for item in archive.findall(".//ProteinExpression")
    ]
    return ";".join(ordered_unique(changes))


def extract_grch38_vcf_location(archive: ET.Element) -> dict[str, object] | None:
    locations = [
        location
        for location in archive.findall(".//SequenceLocation")
        if location.attrib.get("Assembly") == "GRCh38"
        and location.attrib.get("positionVCF")
        and location.attrib.get("referenceAlleleVCF")
        and location.attrib.get("alternateAlleleVCF")
    ]
    if not locations:
        return None

    location = next(
        (item for item in locations if item.attrib.get("forDisplay") == "true"),
        locations[0],
    )

    try:
        return {
            "chrom": normalize_chrom(location.attrib["Chr"]),
            "pos": int(location.attrib["positionVCF"]),
            "ref": location.attrib["referenceAlleleVCF"],
            "alt": location.attrib["alternateAlleleVCF"],
        }
    except (KeyError, ValueError):
        return None


def increment_counter(counter: dict[str, int] | None, key: str) -> None:
    if counter is not None:
        counter[key] = counter.get(key, 0) + 1


def parse_clinvar_vcv_xml(
    raw_xml: str, gene: str, stats: dict[str, int] | None = None
) -> list[dict[str, object]]:
    root = ET.fromstring(raw_xml)
    records: list[dict[str, object]] = []

    for archive in root.findall(".//VariationArchive"):
        increment_counter(stats, "vcv_records")
        if first_text(archive, "./Species") != "Homo sapiens":
            increment_counter(stats, "non_human")
            continue

        clinsig, review_status = extract_clinvar_classification(archive)
        if not is_pathogenic_or_likely_pathogenic(clinsig):
            increment_counter(stats, "non_plp")
            continue

        location = extract_grch38_vcf_location(archive)
        if location is None:
            increment_counter(stats, "missing_grch38_vcf_location")
            continue

        increment_counter(stats, "kept_plp_grch38")
        records.append(
            {
                "gene": gene,
                "clinvar_id": archive.attrib.get("Accession", ""),
                "variation_id": archive.attrib.get("VariationID", ""),
                "variant_type": extract_variant_type(archive),
                "molecular_consequence": extract_molecular_consequences(archive),
                "protein_change": extract_protein_changes(archive),
                "submitter_count": parse_int_attribute(archive, "NumberOfSubmitters"),
                "submission_count": parse_int_attribute(archive, "NumberOfSubmissions"),
                "date_last_updated": archive.attrib.get("DateLastUpdated", ""),
                "date_created": archive.attrib.get("DateCreated", ""),
                "most_recent_submission": archive.attrib.get("MostRecentSubmission", ""),
                "chrom": location["chrom"],
                "pos": location["pos"],
                "ref": location["ref"],
                "alt": location["alt"],
                "clinsig": clinsig,
                "review_status": review_status,
                "title": archive.attrib.get("VariationName", ""),
            }
        )

    return records


def fetch_clinvar_ids(
    gene: str,
    max_records: int | None = None,
    page_size: int = 1000,
    pause: float = CLINVAR_PAUSE,
) -> list[str]:
    query = f"{gene}[GENE]"
    with Entrez.esearch(db="clinvar", term=query, retmax=0) as handle:
        initial = Entrez.read(handle)

    total = int(initial.get("Count", 0))
    if total == 0:
        return []

    limit = min(total, max_records) if max_records is not None else total
    ids: list[str] = []
    for retstart in range(0, limit, page_size):
        retmax = min(page_size, limit - retstart)
        with Entrez.esearch(
            db="clinvar",
            term=query,
            retstart=retstart,
            retmax=retmax,
        ) as handle:
            result = Entrez.read(handle)
        ids.extend(result.get("IdList", []))
        time.sleep(pause)

    if max_records is not None and total > max_records:
        log.warning(
            "ClinVar gene %s has %s records; limited to %s by --clinvar-retmax",
            gene,
            total,
            max_records,
        )
    return ids


def fetch_clinvar_vcv_xml_batch(
    batch: list[str],
    gene: str,
    retries: int = CLINVAR_EFETCH_RETRIES,
    base_delay: float = CLINVAR_EFETCH_RETRY_BASE_DELAY,
) -> str | None:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with Entrez.efetch(
                db="clinvar",
                id=",".join(batch),
                rettype="vcv",
                retmode="xml",
                is_variationid="true",
            ) as handle:
                return handle.read()
        except (IncompleteRead, OSError, RuntimeError) as exc:
            last_error = exc

        if attempt == retries:
            break

        delay = base_delay * (2 ** (attempt - 1))
        log.warning(
            "ClinVar efetch failed for gene %s batch size %s on attempt %s/%s: %s; retrying in %.1fs",
            gene,
            len(batch),
            attempt,
            retries,
            last_error,
            delay,
        )
        time.sleep(delay)

    log.warning(
        "ClinVar efetch failed for gene %s batch size %s after %s attempts: %s",
        gene,
        len(batch),
        retries,
        last_error,
    )
    return None


def fetch_clinvar_batch_records(
    gene: str,
    batch: list[str],
    stats: dict[str, int],
    pause: float = CLINVAR_PAUSE,
) -> list[dict[str, object]]:
    raw = fetch_clinvar_vcv_xml_batch(batch, gene)
    if raw is not None:
        return parse_clinvar_vcv_xml(raw, gene, stats=stats)

    if len(batch) == 1:
        log.warning("Skipping ClinVar VariationID %s for gene %s after failed efetch", batch[0], gene)
        return []

    midpoint = len(batch) // 2
    log.warning(
        "Splitting failed ClinVar batch for gene %s from %s into %s and %s records",
        gene,
        len(batch),
        midpoint,
        len(batch) - midpoint,
    )
    records: list[dict[str, object]] = []
    records.extend(fetch_clinvar_batch_records(gene, batch[:midpoint], stats, pause=pause))
    time.sleep(pause)
    records.extend(fetch_clinvar_batch_records(gene, batch[midpoint:], stats, pause=pause))
    return records


def fetch_clinvar_variants(
    genes: list[str],
    retmax: int | None = None,
    batch_size: int = CLINVAR_BATCH_SIZE,
    pause: float = CLINVAR_PAUSE,
) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    total_stats: dict[str, int] = {}

    for gene in tqdm(genes, desc="ClinVar fetch"):
        gene_stats: dict[str, int] = {}
        try:
            ids = fetch_clinvar_ids(gene, max_records=retmax, pause=pause)
            if not ids:
                log.warning("ClinVar returned no records for gene %s", gene)
                continue

            for batch_start in range(0, len(ids), batch_size):
                batch = ids[batch_start : batch_start + batch_size]
                records.extend(fetch_clinvar_batch_records(gene, batch, gene_stats, pause=pause))
                time.sleep(pause)

            for key, value in gene_stats.items():
                total_stats[key] = total_stats.get(key, 0) + value

            missing_grch38 = gene_stats.get("missing_grch38_vcf_location", 0)
            if missing_grch38:
                log.warning(
                    "ClinVar gene %s: discarded %s P/LP records without GRCh38 VCF coordinates",
                    gene,
                    missing_grch38,
                )

        except Exception as exc:
            log.warning("ClinVar error for gene %s: %s", gene, exc)

    df = pd.DataFrame(records)
    if total_stats:
        log.info(
            "ClinVar parsing summary: %s VCV records, %s P/LP GRCh38 kept, "
            "%s P/LP discarded without GRCh38 VCF coordinates",
            total_stats.get("vcv_records", 0),
            total_stats.get("kept_plp_grch38", 0),
            total_stats.get("missing_grch38_vcf_location", 0),
        )
    if df.empty:
        return df

    df = df.dropna(subset=["chrom", "pos", "ref", "alt"]).copy()
    df = df[df["pos"].astype(int) > 0]
    df["variant_key"] = (
        df["chrom"].astype(str)
        + ":"
        + df["pos"].astype(str)
        + ":"
        + df["ref"].astype(str)
        + ":"
        + df["alt"].astype(str)
    )
    df = df.drop_duplicates(subset="variant_key").reset_index(drop=True)
    log.info("ClinVar: %s unique P/LP variants across %s genes", len(df), len(genes))
    return df


def gnomad_variant_id(row: pd.Series) -> str:
    return f"{normalize_chrom(row['chrom'])}-{int(row['pos'])}-{row['ref']}-{row['alt']}"


def is_retryable_http_error(exc: requests.HTTPError) -> bool:
    response = exc.response
    if response is None:
        return True
    return response.status_code in {429, 500, 502, 503, 504}


def _with_retry(
    fn: Callable[[], dict[str, object]],
    retries: int = GNOMAD_RETRIES,
    base_delay: float = GNOMAD_RETRY_BASE_DELAY,
    context: str = "request",
) -> dict[str, object]:
    last_error: Exception | None = None

    for attempt in range(1, retries + 1):
        try:
            return fn()
        except requests.HTTPError as exc:
            if not is_retryable_http_error(exc):
                raise
            last_error = exc
        except (requests.Timeout, requests.ConnectionError, requests.RequestException, ValueError) as exc:
            last_error = exc

        if attempt == retries:
            break

        delay = base_delay * (2 ** (attempt - 1))
        log.warning(
            "%s failed on attempt %s/%s: %s; retrying in %.1fs",
            context,
            attempt,
            retries,
            last_error,
            delay,
        )
        time.sleep(delay)

    if last_error is not None:
        raise last_error
    raise RuntimeError(f"{context} failed without an exception")


def post_gnomad_graphql(
    session: requests.Session,
    query: str,
    variables: dict[str, object],
    retries: int = GNOMAD_RETRIES,
    base_delay: float = GNOMAD_RETRY_BASE_DELAY,
    context: str = "gnomAD GraphQL request",
) -> dict[str, object]:
    def request_once() -> dict[str, object]:
        response = session.post(
            GNOMAD_GRAPHQL,
            json={"query": query, "variables": variables},
            timeout=20,
        )
        payload = response.json()
        if payload.get("errors"):
            raise RuntimeError(json.dumps(payload["errors"], ensure_ascii=True))
        response.raise_for_status()
        return payload

    payload = _with_retry(
        request_once,
        retries=retries,
        base_delay=base_delay,
        context=context,
    )
    return payload


def fetch_gnomad_af(
    session: requests.Session, variant_id: str, dataset: str = GNOMAD_DATASET
) -> dict[str, object]:
    try:
        payload = post_gnomad_graphql(
            session,
            GNOMAD_QUERY,
            {"variantId": variant_id, "dataset": dataset},
            context=f"gnomAD variant {variant_id}",
        )
        variant = payload.get("data", {}).get("variant")
        if variant is None:
            return {
                "matched": False,
                "af": None,
                "ac": None,
                "an": None,
                "genome_af": None,
                "genome_ac": None,
                "genome_an": None,
                "status": "no_gnomad_record",
            }

        exome = variant.get("exome")
        genome = variant.get("genome") or {}
        if exome is None:
            return {
                "matched": True,
                "af": None,
                "ac": None,
                "an": None,
                "genome_af": genome.get("af"),
                "genome_ac": genome.get("ac"),
                "genome_an": genome.get("an"),
                "status": "no_exome_data",
            }

        return {
            "matched": True,
            "af": exome.get("af"),
            "ac": exome.get("ac"),
            "an": exome.get("an"),
            "genome_af": genome.get("af"),
            "genome_ac": genome.get("ac"),
            "genome_an": genome.get("an"),
            "status": "matched",
        }
    except RuntimeError as exc:
        if "Variant not found" in str(exc):
            return {
                "matched": False,
                "af": None,
                "ac": None,
                "an": None,
                "genome_af": None,
                "genome_ac": None,
                "genome_an": None,
                "status": "no_gnomad_record",
            }
        return {
            "matched": False,
            "af": None,
            "ac": None,
            "an": None,
            "genome_af": None,
            "genome_ac": None,
            "genome_an": None,
            "status": f"error:{exc}",
        }
    except Exception as exc:
        return {
            "matched": False,
            "af": None,
            "ac": None,
            "an": None,
            "genome_af": None,
            "genome_ac": None,
            "genome_an": None,
            "status": f"error:{exc}",
        }


def fetch_gnomad_region(
    session: requests.Session,
    chrom: str,
    start: int,
    stop: int,
    dataset: str = GNOMAD_DATASET,
) -> list[dict[str, object]] | None:
    try:
        payload = post_gnomad_graphql(
            session,
            GNOMAD_REGION_QUERY,
            {
                "chrom": normalize_chrom(chrom),
                "start": max(1, int(start)),
                "stop": max(1, int(stop)),
                "dataset": dataset,
                "referenceGenome": "GRCh38",
            },
            base_delay=GNOMAD_REGION_RETRY_BASE_DELAY,
            context=f"gnomAD region {normalize_chrom(chrom)}:{start}-{stop}",
        )
        region = payload.get("data", {}).get("region")
        return region.get("variants", []) if region else []
    except Exception as exc:
        log.warning(
            "gnomAD region query failed for %s:%s-%s: %s",
            normalize_chrom(chrom),
            start,
            stop,
            exc,
        )
        return None


def has_discordant_allele_at_position(
    variants: list[dict[str, object]], row: pd.Series
) -> bool:
    pos = int(row["pos"])
    ref = str(row["ref"])
    alt = str(row["alt"])
    for variant in variants:
        if int(variant.get("pos", -1)) != pos:
            continue
        if str(variant.get("ref", "")) != ref or str(variant.get("alt", "")) != alt:
            return True
    return False


def match_gnomad(
    df: pd.DataFrame,
    pause: float = GNOMAD_PAUSE,
    dataset: str = GNOMAD_DATASET,
    classify_unmatched_with_region: bool = True,
) -> pd.DataFrame:
    statuses: list[str] = []
    afs: list[object] = []
    acs: list[object] = []
    ans: list[object] = []
    genome_afs: list[object] = []
    genome_acs: list[object] = []
    genome_ans: list[object] = []
    match_categories: list[str] = []
    variants_at_position_counts: list[object] = []
    region_statuses: list[str] = []

    with requests.Session() as session:
        session.headers.update({"User-Agent": "arrhythmia-variant-pipeline/1.0"})

        for _, row in tqdm(df.iterrows(), total=len(df), desc="gnomAD match"):
            result = fetch_gnomad_af(session, gnomad_variant_id(row), dataset=dataset)
            status = str(result["status"])
            region_count: object = 0
            region_status = "not_queried"

            if status in {"matched", "no_exome_data"}:
                category = "exact_match"
            elif status == "no_gnomad_record" and classify_unmatched_with_region:
                region_variants = fetch_gnomad_region(
                    session,
                    chrom=str(row["chrom"]),
                    start=int(row["pos"]),
                    stop=int(row["pos"]),
                    dataset=dataset,
                )
                if region_variants is None:
                    region_count = pd.NA
                    region_status = "error"
                    category = "query_error"
                else:
                    region_status = "queried"
                    region_count = len(region_variants)
                    category = (
                        "allele_discordance"
                        if has_discordant_allele_at_position(region_variants, row)
                        else "no_gnomad_record"
                    )
            elif status == "no_gnomad_record":
                category = "no_gnomad_record"
            elif status.startswith("error:"):
                category = "query_error"
            else:
                category = "allele_discordance"

            statuses.append(status)
            afs.append(result["af"])
            acs.append(result["ac"])
            ans.append(result["an"])
            genome_afs.append(result["genome_af"])
            genome_acs.append(result["genome_ac"])
            genome_ans.append(result["genome_an"])
            match_categories.append(category)
            variants_at_position_counts.append(region_count)
            region_statuses.append(region_status)
            time.sleep(pause)

    matched_df = df.copy()
    matched_df["gnomad_status"] = statuses
    matched_df["gnomad_af"] = pd.to_numeric(pd.Series(afs), errors="coerce").to_numpy()
    matched_df["gnomad_ac"] = pd.to_numeric(pd.Series(acs), errors="coerce").to_numpy()
    matched_df["gnomad_an"] = pd.to_numeric(pd.Series(ans), errors="coerce").to_numpy()
    matched_df["gnomad_genome_af"] = pd.to_numeric(pd.Series(genome_afs), errors="coerce").to_numpy()
    matched_df["gnomad_genome_ac"] = pd.to_numeric(pd.Series(genome_acs), errors="coerce").to_numpy()
    matched_df["gnomad_genome_an"] = pd.to_numeric(pd.Series(genome_ans), errors="coerce").to_numpy()
    matched_df["match_category"] = match_categories
    matched_df["gnomad_variants_at_pos_for_match"] = variants_at_position_counts
    matched_df["gnomad_region_query_status"] = region_statuses

    log.info("gnomAD match categories:\n%s", matched_df["match_category"].value_counts())
    return matched_df


def classify_frequencies(df: pd.DataFrame) -> pd.DataFrame:
    def classify_af(af: object) -> str:
        if pd.isna(af):
            return "no_af_data"
        if float(af) <= AF_ULTRA_RARE:
            return "ultra_rare"
        if float(af) <= AF_RARE:
            return "rare"
        return "higher_freq"

    result = df.copy()
    result["freq_class"] = None
    exact_mask = result["match_category"] == "exact_match"
    result.loc[exact_mask, "freq_class"] = result.loc[exact_mask, "gnomad_af"].apply(
        classify_af
    )
    return result


def per_gene_fisher(df: pd.DataFrame, genes: Iterable[str] | None = None) -> pd.DataFrame:
    matched = df[
        (df["match_category"] == "exact_match") & df["gnomad_af"].notna()
    ].copy()
    matched["is_outlier"] = matched["gnomad_af"].astype(float) > AF_ULTRA_RARE

    gene_order = list(genes) if genes is not None else sorted(matched["gene"].unique())
    rows: list[dict[str, object]] = []

    for gene in gene_order:
        gene_mask = matched["gene"] == gene
        n_gene = int(gene_mask.sum())
        n_other = int((~gene_mask).sum())

        if 0 < n_gene < LOW_AF_COVERED_WARNING:
            log.warning(
                "Gene %s has only %s AF-covered exact matches; Fisher p-value may be unstable.",
                gene,
                n_gene,
            )
        if str(gene).upper() == "TRDN" and n_gene < TRDN_AF_COVERED_WARNING:
            log.warning(
                "TRDN has only %s AF-covered exact matches; compare cautiously with the article.",
                n_gene,
            )

        if n_gene == 0 or n_other == 0:
            rows.append(
                {
                    "gene": gene,
                    "n_af_covered": n_gene,
                    "n_outlier": 0,
                    "pct_outlier": 0.0,
                    "odds_ratio": np.nan,
                    "fisher_p": np.nan,
                }
            )
            continue

        a = int(matched.loc[gene_mask, "is_outlier"].sum())
        b = int(n_gene - a)
        c = int(matched.loc[~gene_mask, "is_outlier"].sum())
        d = int(n_other - c)

        odds, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append(
            {
                "gene": gene,
                "n_af_covered": n_gene,
                "n_outlier": a,
                "pct_outlier": round(100 * a / n_gene, 1) if n_gene else 0.0,
                "odds_ratio": round(float(odds), 3),
                "fisher_p": p_value,
            }
        )

    result = pd.DataFrame(rows)
    if result.empty:
        return result

    result["bh_q"] = np.nan
    result["significant"] = False
    valid_mask = result["fisher_p"].notna()
    if valid_mask.any():
        reject, q_values, _, _ = multipletests(
            result.loc[valid_mask, "fisher_p"], method="fdr_bh"
        )
        result.loc[valid_mask, "bh_q"] = q_values
        result.loc[valid_mask, "significant"] = reject

    return result.sort_values(["fisher_p", "gene"], na_position="last").reset_index(
        drop=True
    )


def coverage_analysis(
    df: pd.DataFrame,
    pause: float = GNOMAD_PAUSE,
    sample_n: int | None = 200,
    dataset: str = GNOMAD_DATASET,
) -> pd.DataFrame:
    unmatched = df[
        df["match_category"].isin(["no_gnomad_record", "allele_discordance"])
    ].copy()

    if sample_n is not None and sample_n > 0 and len(unmatched) > sample_n:
        unmatched = unmatched.sample(sample_n, random_state=42)

    results: list[dict[str, object]] = []
    with requests.Session() as session:
        session.headers.update({"User-Agent": "arrhythmia-variant-pipeline/1.0"})

        for _, row in tqdm(
            unmatched.iterrows(), total=len(unmatched), desc="Coverage analysis"
        ):
            variants_at_pos = fetch_gnomad_region(
                session,
                chrom=str(row["chrom"]),
                start=int(row["pos"]) - 1,
                stop=int(row["pos"]) + 1,
                dataset=dataset,
            )
            query_status = "error" if variants_at_pos is None else "ok"
            variants_for_summary = variants_at_pos or []
            alt_alleles = sorted(
                {
                    str(variant.get("alt", ""))
                    for variant in variants_for_summary
                    if variant.get("alt")
                }
            )
            variant_ids = sorted(
                {
                    str(variant.get("variantId", ""))
                    for variant in variants_for_summary
                    if variant.get("variantId")
                }
            )
            results.append(
                {
                    "variant_key": row["variant_key"],
                    "gene": row["gene"],
                    "match_category": row["match_category"],
                    "coverage_query_status": query_status,
                    "gnomad_variants_at_pos": (
                        pd.NA if variants_at_pos is None else len(variants_at_pos)
                    ),
                    "has_any_gnomad": (
                        pd.NA if variants_at_pos is None else len(variants_at_pos) > 0
                    ),
                    "alt_alleles_at_pos": ";".join(alt_alleles),
                    "gnomad_variant_ids_at_pos": ";".join(variant_ids),
                }
            )
            time.sleep(pause)

    cov_df = pd.DataFrame(results)
    if not cov_df.empty:
        log.info(
            "Coverage query status:\n%s",
            cov_df["coverage_query_status"].value_counts(dropna=False),
        )
        ok_cov_df = cov_df[cov_df["coverage_query_status"] == "ok"]
        if not ok_cov_df.empty:
            log.info(
                "Coverage analysis summary:\n%s",
                ok_cov_df.groupby(["match_category", "has_any_gnomad"]).size(),
            )
    return cov_df


def plot_af_distribution(df: pd.DataFrame, output_path: Path = DEFAULT_AF_FIGURE) -> None:
    matched = df[(df["match_category"] == "exact_match") & df["gnomad_af"].notna()].copy()
    if matched.empty:
        log.warning("Skipping AF distribution figure because no matched AF values exist.")
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    matched["log10_af"] = np.log10(matched["gnomad_af"].astype(float).clip(lower=1e-9))

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(
        matched["log10_af"],
        bins=40,
        color="#2563EB",
        alpha=0.85,
        edgecolor="white",
        linewidth=0.4,
    )
    ax.axvline(
        np.log10(AF_ULTRA_RARE),
        color="#DC2626",
        linestyle="--",
        linewidth=1.5,
        label=f"AF = {AF_ULTRA_RARE:.0e}",
    )
    ax.axvline(
        np.log10(AF_RARE),
        color="#F97316",
        linestyle="--",
        linewidth=1.5,
        label=f"AF = {AF_RARE:.0e}",
    )
    ax.set_xlabel("log10(Allele Frequency)")
    ax.set_ylabel("Number of variants")
    ax.set_title("Allele-frequency distribution of ClinVar-matched variants")
    ax.legend(frameon=False)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"10^{{{x:.0f}}}"))
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    log.info("Saved: %s", output_path)


def plot_match_categories(
    df: pd.DataFrame, output_path: Path = DEFAULT_MATCH_FIGURE
) -> None:
    if df.empty or "match_category" not in df:
        log.warning("Skipping match-category figure because no match categories exist.")
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    counts = df["match_category"].value_counts()
    labels = {
        "exact_match": "Exact gnomAD match",
        "allele_discordance": "Positional match,\nallele discordance",
        "no_gnomad_record": "No gnomAD record",
        "query_error": "Query error",
    }
    display_labels = [labels.get(index, index) for index in counts.index]

    colors = ["#16A34A", "#F59E0B", "#9CA3AF", "#EF4444"]
    fig, ax = plt.subplots(figsize=(8, 4))
    bars = ax.barh(display_labels, counts.values, color=colors[: len(counts)], edgecolor="white")
    for bar, value in zip(bars, counts.values):
        ax.text(
            bar.get_width() + counts.max() * 0.01,
            bar.get_y() + bar.get_height() / 2,
            f"{value} ({100 * value / len(df):.1f}%)",
            va="center",
            fontsize=10,
        )
    ax.set_xlabel("Number of variants")
    ax.set_title("ClinVar P/LP variant match categories")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    log.info("Saved: %s", output_path)


def plot_gene_heatmap(
    fisher_df: pd.DataFrame, output_path: Path = DEFAULT_HEATMAP_FIGURE
) -> None:
    if fisher_df.empty:
        log.warning("Skipping gene heatmap because Fisher results are empty.")
        return

    plot_df = fisher_df[fisher_df["n_af_covered"] >= 3].copy()
    if plot_df.empty:
        log.warning("Skipping gene heatmap because fewer than 3 AF-covered variants exist.")
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plot_df = plot_df.sort_values("pct_outlier", ascending=False)

    fig, ax = plt.subplots(figsize=(7, max(4, len(plot_df) * 0.45)))
    pivot = plot_df.set_index("gene")[["pct_outlier"]]
    sns.heatmap(
        pivot,
        annot=True,
        fmt=".1f",
        cmap="YlOrRd",
        linewidths=0.5,
        ax=ax,
        cbar_kws={"label": "% outlier variants"},
        annot_kws={"size": 9},
    )

    for index, (_, row) in enumerate(plot_df.iterrows()):
        if bool(row.get("significant", False)):
            ax.text(
                1.02,
                index + 0.5,
                "*",
                transform=ax.get_yaxis_transform(),
                va="center",
                fontsize=11,
                color="#DC2626",
            )

    ax.set_title("Per-gene outlier burden (AF > 1e-5)\n* = BH q < 0.05")
    ax.set_ylabel("")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved: %s", output_path)


def plot_coverage_summary(
    cov_df: pd.DataFrame, output_path: Path = DEFAULT_COVERAGE_FIGURE
) -> None:
    if cov_df.empty:
        log.warning("Skipping coverage figure because coverage results are empty.")
        return

    plot_df = cov_df.copy()
    if "coverage_query_status" in plot_df.columns:
        failed_queries = int((plot_df["coverage_query_status"] != "ok").sum())
        if failed_queries:
            log.warning(
                "Coverage figure excludes %s rows with failed region queries.",
                failed_queries,
            )
        plot_df = plot_df[plot_df["coverage_query_status"] == "ok"].copy()

    plot_df = plot_df.dropna(subset=["has_any_gnomad"]).copy()
    if plot_df.empty:
        log.warning("Skipping coverage figure because no successful coverage rows exist.")
        return

    plot_df["has_any_gnomad"] = (
        plot_df["has_any_gnomad"].astype(str).str.strip().str.lower().isin({"true", "1"})
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary = (
        plot_df.groupby(["match_category", "has_any_gnomad"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=[False, True], fill_value=0)
    )
    summary.columns = ["No gnomAD at position", "Other gnomAD variant at position"]
    summary.index = [str(index).replace("_", " ").title() for index in summary.index]

    ax = summary.plot(
        kind="bar",
        stacked=True,
        figsize=(7, 5),
        color=["#6B7280", "#3B82F6"],
        edgecolor="white",
    )
    ax.set_title("gnomAD coverage at unmatched ClinVar positions")
    ax.set_ylabel("Number of variants sampled")
    ax.set_xlabel("")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=180)
    plt.close()
    log.info("Saved: %s", output_path)


def save_table(
    df: pd.DataFrame,
    output_path: Path,
    metadata: dict[str, object] | None = None,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    if metadata is not None:
        write_cache_metadata(output_path, metadata)


def load_table(
    path: Path,
    expected_metadata: dict[str, object] | None = None,
    strict_metadata: bool = False,
) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    if expected_metadata is not None:
        validate_cache_metadata(path, expected_metadata, strict=strict_metadata)
    return pd.read_csv(path)


def log_match_summary(matched_df: pd.DataFrame) -> None:
    total = len(matched_df)
    exact = matched_df[matched_df["match_category"] == "exact_match"]
    exact_pct = 100 * len(exact) / total if total else 0
    allele_discordance = int((matched_df["match_category"] == "allele_discordance").sum())
    no_record = int((matched_df["match_category"] == "no_gnomad_record").sum())

    log.info("Match summary:")
    log.info("  Exact matches:      %s (%.1f%%)", len(exact), exact_pct)
    log.info("  Allele discordance: %s", allele_discordance)
    log.info("  No gnomAD record:   %s", no_record)

    if not exact.empty:
        log.info("Frequency distribution for exact matches:")
        log.info("\n%s", exact["freq_class"].value_counts(dropna=False).to_string())
        af_values = exact["gnomad_af"].dropna()
        if not af_values.empty:
            log.info("  Median AF: %.2e", af_values.median())
            log.info("  Max AF:    %.2e", af_values.max())


def make_run_qc_summary(matched_df: pd.DataFrame, cov_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    total = len(matched_df)
    category_counts = (
        matched_df["match_category"].value_counts(dropna=False)
        if "match_category" in matched_df
        else pd.Series(dtype=int)
    )

    for category in ["exact_match", "allele_discordance", "no_gnomad_record", "query_error"]:
        count = int(category_counts.get(category, 0))
        rows.append(
            {
                "section": "match_category",
                "metric": category,
                "count": count,
                "fraction": count / total if total else np.nan,
                "value": "",
            }
        )

    allele_discordance = int(category_counts.get("allele_discordance", 0))
    no_record = int(category_counts.get("no_gnomad_record", 0))
    unmatched_total = allele_discordance + no_record
    observed_ratio = (
        allele_discordance / no_record if no_record > 0 else np.nan
    )
    article_ratio = ARTICLE_ALLELE_DISCORDANCE / ARTICLE_NO_GNOMAD_RECORD
    rows.extend(
        [
            {
                "section": "match_qc",
                "metric": "allele_discordance_to_no_record_ratio",
                "count": "",
                "fraction": "",
                "value": observed_ratio,
            },
            {
                "section": "match_qc",
                "metric": "article_ratio_849_to_962",
                "count": "",
                "fraction": "",
                "value": article_ratio,
            },
        ]
    )

    if "gnomad_genome_af" in matched_df.columns:
        exact_mask = matched_df["match_category"] == "exact_match"
        exact_total = int(exact_mask.sum())
        genome_af = pd.to_numeric(matched_df.loc[exact_mask, "gnomad_genome_af"], errors="coerce")
        genome_ac = pd.to_numeric(
            matched_df.loc[exact_mask, "gnomad_genome_ac"]
            if "gnomad_genome_ac" in matched_df.columns
            else pd.Series(dtype=float),
            errors="coerce",
        )
        no_exome_mask = exact_mask & pd.to_numeric(matched_df["gnomad_af"], errors="coerce").isna()
        genome_data_count = int(genome_af.notna().sum())
        genome_positive_count = int((genome_ac.fillna(0) > 0).sum())
        no_exome_with_genome_count = int(
            (
                no_exome_mask
                & pd.to_numeric(matched_df["gnomad_genome_af"], errors="coerce").notna()
            ).sum()
        )
        rows.extend(
            [
                {
                    "section": "genome_qc",
                    "metric": "exact_matches_with_genome_af",
                    "count": genome_data_count,
                    "fraction": genome_data_count / exact_total if exact_total else np.nan,
                    "value": "",
                },
                {
                    "section": "genome_qc",
                    "metric": "exact_matches_with_genome_ac_positive",
                    "count": genome_positive_count,
                    "fraction": genome_positive_count / exact_total if exact_total else np.nan,
                    "value": "",
                },
                {
                    "section": "genome_qc",
                    "metric": "no_exome_af_exact_matches_with_genome_af",
                    "count": no_exome_with_genome_count,
                    "fraction": (
                        no_exome_with_genome_count / int(no_exome_mask.sum())
                        if int(no_exome_mask.sum())
                        else np.nan
                    ),
                    "value": "",
                },
            ]
        )

    if unmatched_total > 0:
        log.info(
            "Unmatched split: allele_discordance=%s, no_gnomad_record=%s, ratio=%.3g "
            "(article ratio %.3g from 849/962).",
            allele_discordance,
            no_record,
            observed_ratio,
            article_ratio,
        )
        if (
            unmatched_total >= MATCH_RATIO_WARNING_MIN_UNMATCHED
            and pd.notna(observed_ratio)
            and (observed_ratio < 0.5 or observed_ratio > 2.0)
        ):
            log.warning(
                "allele_discordance:no_gnomad_record ratio %.3g is far from the "
                "near-balanced article split; inspect ClinVar filtering and "
                "coverage_query_status.",
                observed_ratio,
            )

    if not cov_df.empty and "coverage_query_status" in cov_df:
        coverage_total = len(cov_df)
        status_counts = cov_df["coverage_query_status"].value_counts(dropna=False)
        for status, count_value in status_counts.items():
            count = int(count_value)
            rows.append(
                {
                    "section": "coverage_query_status",
                    "metric": str(status),
                    "count": count,
                    "fraction": count / coverage_total if coverage_total else np.nan,
                    "value": "",
                }
            )

        error_count = int(status_counts.drop(labels=["ok"], errors="ignore").sum())
        if error_count:
            log.warning(
                "Coverage analysis has %s failed/non-ok region queries; inspect %s rows.",
                error_count,
                "coverage_query_status",
            )

    return pd.DataFrame(rows)


def run_pipeline(
    genes: list[str] = ARRHYTHMIA_GENES,
    email: str | None = None,
    output_prefix: str | None = None,
    clinvar_retmax: int | None = None,
    clinvar_batch_size: int = CLINVAR_BATCH_SIZE,
    skip_clinvar: bool = False,
    clinvar_cache: Path = DEFAULT_CLINVAR_CACHE,
    skip_gnomad: bool = False,
    gnomad_cache: Path = DEFAULT_GNOMAD_CACHE,
    skip_coverage: bool = False,
    coverage_cache: Path = DEFAULT_COVERAGE_CACHE,
    fisher_path: Path = DEFAULT_FISHER_PATH,
    qc_summary_path: Path = DEFAULT_QC_SUMMARY_PATH,
    af_figure_path: Path = DEFAULT_AF_FIGURE,
    match_figure_path: Path = DEFAULT_MATCH_FIGURE,
    heatmap_figure_path: Path = DEFAULT_HEATMAP_FIGURE,
    coverage_figure_path: Path = DEFAULT_COVERAGE_FIGURE,
    coverage_sample_n: int | None = 200,
    gnomad_pause: float = GNOMAD_PAUSE,
    dataset: str = GNOMAD_DATASET,
    make_plots: bool = True,
    classify_unmatched_with_region: bool = True,
    strict_cache_metadata: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if clinvar_batch_size < 1:
        raise SystemExit("--clinvar-batch-size must be at least 1.")

    if not skip_clinvar:
        set_entrez_email(email)

    clinvar_cache = prefixed_path(clinvar_cache, output_prefix)
    gnomad_cache = prefixed_path(gnomad_cache, output_prefix)
    coverage_cache = prefixed_path(coverage_cache, output_prefix)
    fisher_path = prefixed_path(fisher_path, output_prefix)
    qc_summary_path = prefixed_path(qc_summary_path, output_prefix)
    af_figure_path = prefixed_path(af_figure_path, output_prefix)
    match_figure_path = prefixed_path(match_figure_path, output_prefix)
    heatmap_figure_path = prefixed_path(heatmap_figure_path, output_prefix)
    coverage_figure_path = prefixed_path(coverage_figure_path, output_prefix)

    gene_signature = canonical_genes(genes)
    clinvar_expected_metadata = {
        "cache_type": "clinvar_variants",
        "genes": gene_signature,
        "pipeline_version": PIPELINE_VERSION,
        "clinvar_retmax": clinvar_retmax,
        "clinvar_batch_size": clinvar_batch_size,
    }
    gnomad_expected_metadata = {
        "cache_type": "gnomad_matched",
        "genes": gene_signature,
        "gnomad_dataset": dataset,
        "pipeline_version": PIPELINE_VERSION,
        "region_match_qc": classify_unmatched_with_region,
    }
    coverage_expected_metadata = {
        "cache_type": "coverage_analysis",
        "genes": gene_signature,
        "gnomad_dataset": dataset,
        "pipeline_version": PIPELINE_VERSION,
        "coverage_sample_n": coverage_sample_n,
    }

    if skip_clinvar:
        log.info("Loading ClinVar cache: %s", clinvar_cache)
        clinvar_df = load_table(
            clinvar_cache,
            expected_metadata=clinvar_expected_metadata,
            strict_metadata=strict_cache_metadata,
        )
    else:
        clinvar_df = fetch_clinvar_variants(
            genes,
            retmax=clinvar_retmax,
            batch_size=clinvar_batch_size,
        )
        save_table(
            clinvar_df,
            clinvar_cache,
            metadata=make_cache_metadata(
                "clinvar_variants",
                genes,
                len(clinvar_df),
                extra={
                    "clinvar_retmax": clinvar_retmax,
                    "clinvar_batch_size": clinvar_batch_size,
                },
            ),
        )
        log.info("ClinVar variants saved: %s", clinvar_cache)

    log.info("Total unique ClinVar variants: %s", len(clinvar_df))

    if skip_gnomad:
        log.info("Loading gnomAD cache: %s", gnomad_cache)
        matched_df = load_table(
            gnomad_cache,
            expected_metadata=gnomad_expected_metadata,
            strict_metadata=strict_cache_metadata,
        )
    else:
        matched_df = match_gnomad(
            clinvar_df,
            pause=gnomad_pause,
            dataset=dataset,
            classify_unmatched_with_region=classify_unmatched_with_region,
        )
        save_table(
            matched_df,
            gnomad_cache,
            metadata=make_cache_metadata(
                "gnomad_matched",
                genes,
                len(matched_df),
                dataset=dataset,
                extra={"region_match_qc": classify_unmatched_with_region},
            ),
        )
        log.info("gnomAD matched data saved: %s", gnomad_cache)

    matched_df = classify_frequencies(matched_df)
    if not skip_gnomad:
        save_table(
            matched_df,
            gnomad_cache,
            metadata=make_cache_metadata(
                "gnomad_matched",
                genes,
                len(matched_df),
                dataset=dataset,
                extra={"region_match_qc": classify_unmatched_with_region},
            ),
        )
    log_match_summary(matched_df)

    fisher_df = per_gene_fisher(matched_df, genes=genes)
    save_table(
        fisher_df,
        fisher_path,
        metadata=make_cache_metadata(
            "per_gene_fisher",
            genes,
            len(fisher_df),
            dataset=dataset,
        ),
    )
    log.info("Per-gene results saved: %s", fisher_path)

    if skip_coverage:
        log.info("Loading coverage cache: %s", coverage_cache)
        cov_df = load_table(
            coverage_cache,
            expected_metadata=coverage_expected_metadata,
            strict_metadata=strict_cache_metadata,
        )
    else:
        cov_df = coverage_analysis(
            matched_df,
            pause=gnomad_pause,
            sample_n=coverage_sample_n,
            dataset=dataset,
        )
        save_table(
            cov_df,
            coverage_cache,
            metadata=make_cache_metadata(
                "coverage_analysis",
                genes,
                len(cov_df),
                dataset=dataset,
                extra={"coverage_sample_n": coverage_sample_n},
            ),
        )
        log.info("Coverage analysis saved: %s", coverage_cache)

    qc_summary_df = make_run_qc_summary(matched_df, cov_df)
    save_table(
        qc_summary_df,
        qc_summary_path,
        metadata=make_cache_metadata(
            "run_qc_summary",
            genes,
            len(qc_summary_df),
            dataset=dataset,
            extra={"coverage_sample_n": coverage_sample_n},
        ),
    )
    log.info("Run QC summary saved: %s", qc_summary_path)

    if make_plots:
        plot_af_distribution(matched_df, af_figure_path)
        plot_match_categories(matched_df, match_figure_path)
        plot_gene_heatmap(fisher_df, heatmap_figure_path)
        plot_coverage_summary(cov_df, coverage_figure_path)

    log.info("Pipeline complete.")
    return matched_df, fisher_df, cov_df


def parse_gene_list(values: list[str] | None) -> list[str]:
    if not values:
        return ARRHYTHMIA_GENES

    tokens: list[str] = []
    for value in values:
        tokens.extend(value.replace(";", ",").split(","))
    genes = [item.strip().upper() for item in tokens]
    return [gene for gene in genes if gene]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch ClinVar P/LP arrhythmia variants and match them to gnomAD r4."
    )
    parser.add_argument(
        "--genes",
        nargs="+",
        help="Gene list. Accepts space-, comma-, or semicolon-separated values.",
    )
    parser.add_argument("--email", help="Entrez email. Can also use ENTREZ_EMAIL.")
    parser.add_argument(
        "--output-prefix",
        help="Prefix all generated CSV, metadata, and PNG output filenames.",
    )
    parser.add_argument(
        "--clinvar-retmax",
        type=int,
        help="Optional maximum ClinVar records to inspect per gene.",
    )
    parser.add_argument(
        "--clinvar-batch-size",
        type=int,
        default=CLINVAR_BATCH_SIZE,
        help="ClinVar VCV records per efetch request. Lower this for large genes.",
    )
    parser.add_argument("--skip-clinvar", action="store_true", help="Load ClinVar cache.")
    parser.add_argument("--skip-gnomad", action="store_true", help="Load gnomAD cache.")
    parser.add_argument("--skip-coverage", action="store_true", help="Load coverage cache.")
    parser.add_argument(
        "--strict-cache-metadata",
        action="store_true",
        help="Fail instead of warning when cache metadata is missing or stale.",
    )
    parser.add_argument("--no-plots", action="store_true", help="Skip figure generation.")
    parser.add_argument(
        "--no-region-match-qc",
        action="store_true",
        help="Do not query same-position gnomAD records for unmatched variants.",
    )
    parser.add_argument("--clinvar-cache", type=Path, default=DEFAULT_CLINVAR_CACHE)
    parser.add_argument("--gnomad-cache", type=Path, default=DEFAULT_GNOMAD_CACHE)
    parser.add_argument("--coverage-cache", type=Path, default=DEFAULT_COVERAGE_CACHE)
    parser.add_argument("--fisher-output", type=Path, default=DEFAULT_FISHER_PATH)
    parser.add_argument("--qc-summary-output", type=Path, default=DEFAULT_QC_SUMMARY_PATH)
    parser.add_argument("--coverage-sample-n", type=int, default=200)
    parser.add_argument("--gnomad-pause", type=float, default=GNOMAD_PAUSE)
    parser.add_argument("--dataset", default=GNOMAD_DATASET)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_pipeline(
        genes=parse_gene_list(args.genes),
        email=args.email,
        output_prefix=args.output_prefix,
        clinvar_retmax=args.clinvar_retmax,
        clinvar_batch_size=args.clinvar_batch_size,
        skip_clinvar=args.skip_clinvar,
        clinvar_cache=args.clinvar_cache,
        skip_gnomad=args.skip_gnomad,
        gnomad_cache=args.gnomad_cache,
        skip_coverage=args.skip_coverage,
        coverage_cache=args.coverage_cache,
        fisher_path=args.fisher_output,
        qc_summary_path=args.qc_summary_output,
        coverage_sample_n=args.coverage_sample_n,
        gnomad_pause=args.gnomad_pause,
        dataset=args.dataset,
        make_plots=not args.no_plots,
        classify_unmatched_with_region=not args.no_region_match_qc,
        strict_cache_metadata=args.strict_cache_metadata,
    )


if __name__ == "__main__":
    main()
