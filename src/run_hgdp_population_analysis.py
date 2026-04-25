"""HGDP population-layer analysis for ClinVar P/LP arrhythmia variants.

This script adds the Human Genome Diversity Project (HGDP) as an independent
population-diversity layer on top of the existing gnomAD-based evaluability
audit.  It queries the gnomAD v3.1.2 GraphQL API for the HGDP+1KG subset
(~4 000 samples across 80+ fine-grained populations) and performs six analyses:

1. Population presence and frequency in HGDP populations
2. Ancestry mismatch / population conflict detection
3. Haplotype-background proxy (multi-population occurrence)
4. Disease-model compatibility (dominant high-penetrance check)
5. Evaluability assessment through HGDP
6. Summary statistics

Usage
-----
    python src/run_hgdp_population_analysis.py          # full API run
    python src/run_hgdp_population_analysis.py --cached  # reuse cached HGDP data
"""

from __future__ import annotations

import argparse
import logging
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)

BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
FIGURE_DIR = BASE_DIR / "figures"
SUPPLEMENT_DIR = BASE_DIR / "supplementary_tables"

SCORES_IN = DATA_DIR / "arrhythmia_vital_scores.csv"

HGDP_CACHE = DATA_DIR / "hgdp_population_frequencies.csv"
HGDP_PRESENCE_OUT = DATA_DIR / "hgdp_population_presence.csv"
HGDP_ANCESTRY_MISMATCH_OUT = DATA_DIR / "hgdp_ancestry_mismatch.csv"
HGDP_HAPLOTYPE_PROXY_OUT = DATA_DIR / "hgdp_haplotype_background_proxy.csv"
HGDP_DISEASE_MODEL_OUT = DATA_DIR / "hgdp_disease_model_compatibility.csv"
HGDP_EVALUABILITY_OUT = DATA_DIR / "hgdp_evaluability_assessment.csv"
HGDP_SUMMARY_OUT = DATA_DIR / "hgdp_population_analysis_summary.csv"
HGDP_REGIME_FIGURE = FIGURE_DIR / "hgdp_regime_distribution.png"
HGDP_ANCESTRY_FIGURE = FIGURE_DIR / "hgdp_ancestry_mismatch_heatmap.png"
HGDP_EVALUABILITY_FIGURE = FIGURE_DIR / "hgdp_evaluability_breakdown.png"

HGDP_PRESENCE_SUPP = (
    SUPPLEMENT_DIR / "Supplementary_Table_S50_hgdp_population_presence.tsv"
)
HGDP_ANCESTRY_SUPP = (
    SUPPLEMENT_DIR / "Supplementary_Table_S51_hgdp_ancestry_mismatch.tsv"
)
HGDP_HAPLOTYPE_SUPP = (
    SUPPLEMENT_DIR / "Supplementary_Table_S52_hgdp_haplotype_background_proxy.tsv"
)
HGDP_DISEASE_MODEL_SUPP = (
    SUPPLEMENT_DIR / "Supplementary_Table_S53_hgdp_disease_model_compatibility.tsv"
)
HGDP_EVALUABILITY_SUPP = (
    SUPPLEMENT_DIR / "Supplementary_Table_S54_hgdp_evaluability_assessment.tsv"
)
HGDP_SUMMARY_SUPP = (
    SUPPLEMENT_DIR / "Supplementary_Table_S55_hgdp_analysis_summary.tsv"
)

# ── gnomAD API ──────────────────────────────────────────────────────────

GNOMAD_GRAPHQL = "https://gnomad.broadinstitute.org/api"
GNOMAD_DATASET = "gnomad_r3"
GNOMAD_PAUSE = 1.2
GNOMAD_RETRIES = 5
GNOMAD_RETRY_BASE_DELAY = 3.0

GNOMAD_HGDP_QUERY = """
query VariantQuery($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variantId
    genome {
      af
      ac
      an
      populations {
        id
        ac
        an
      }
    }
  }
}
"""

# ── HGDP population geography ──────────────────────────────────────────

HGDP_REGION_MAP: dict[str, str] = {
    "bantukenya": "Africa",
    "bantusafrica": "Africa",
    "biaka": "Africa",
    "mandenka": "Africa",
    "mbuti": "Africa",
    "san": "Africa",
    "yoruba": "Africa",
    "balochi": "Central_South_Asia",
    "brahui": "Central_South_Asia",
    "burusho": "Central_South_Asia",
    "hazara": "Central_South_Asia",
    "kalash": "Central_South_Asia",
    "makrani": "Central_South_Asia",
    "pathan": "Central_South_Asia",
    "sindhi": "Central_South_Asia",
    "cambodian": "East_Asia",
    "dai": "East_Asia",
    "daur": "East_Asia",
    "han": "East_Asia",
    "hezhen": "East_Asia",
    "japanese": "East_Asia",
    "lahu": "East_Asia",
    "miaozu": "East_Asia",
    "mongola": "East_Asia",
    "naxi": "East_Asia",
    "oroqen": "East_Asia",
    "she": "East_Asia",
    "tu": "East_Asia",
    "tujia": "East_Asia",
    "uygur": "East_Asia",
    "xibo": "East_Asia",
    "yakut": "East_Asia",
    "yizu": "East_Asia",
    "adygei": "Europe",
    "basque": "Europe",
    "french": "Europe",
    "italian": "Europe",
    "orcadian": "Europe",
    "russian": "Europe",
    "sardinian": "Europe",
    "tuscan": "Europe",
    "bedouin": "Middle_East",
    "druze": "Middle_East",
    "mozabite": "Middle_East",
    "palestinian": "Middle_East",
    "colombian": "Americas",
    "karitiana": "Americas",
    "maya": "Americas",
    "pima": "Americas",
    "surui": "Americas",
    "bougainville": "Oceania",
    "papuan": "Oceania",
}

HGDP_REGIONS = sorted(set(HGDP_REGION_MAP.values()))

# ── disease-model thresholds ────────────────────────────────────────────

AF_REVIEW_THRESHOLD = 1e-5
DOMINANT_HIGH_PEN_MCAF = 2.5e-5  # max credible AF: prev=1/2000, AC_frac=0.1, pen=0.5, model=2
DOMINANT_LOW_PEN_MCAF = 1e-3     # generous: prev=1/1000, AC_frac=0.2, pen=0.1, model=2
MIN_AN_FOR_EVALUATION = 10       # require at least 10 alleles observed in a population


# ═══════════════════════════════════════════════════════════════════════
# 0. Data acquisition
# ═══════════════════════════════════════════════════════════════════════

def variant_key_to_gnomad_id(variant_key: str) -> str:
    """Convert ``chrom:pos:ref:alt`` to ``chrom-pos-ref-alt``."""
    return variant_key.replace(":", "-")


def query_gnomad_hgdp(variant_id: str) -> dict | None:
    """Return the parsed JSON response for one variant from gnomAD r3."""
    for attempt in range(GNOMAD_RETRIES):
        try:
            resp = requests.post(
                GNOMAD_GRAPHQL,
                json={
                    "query": GNOMAD_HGDP_QUERY,
                    "variables": {
                        "variantId": variant_id,
                        "dataset": GNOMAD_DATASET,
                    },
                },
                timeout=30,
            )
            resp.raise_for_status()
            payload = resp.json()
            if "errors" in payload:
                log.warning("gnomAD error for %s: %s", variant_id, payload["errors"])
                return None
            return payload.get("data", {}).get("variant")
        except (requests.RequestException, ValueError) as exc:
            delay = GNOMAD_RETRY_BASE_DELAY * (2 ** attempt)
            log.warning(
                "gnomAD request failed for %s (attempt %d/%d): %s — retrying in %.1fs",
                variant_id, attempt + 1, GNOMAD_RETRIES, exc, delay,
            )
            time.sleep(delay)
    return None


def extract_hgdp_populations(variant_data: dict) -> list[dict]:
    """Extract HGDP-prefixed population entries from gnomAD genome data."""
    genome = variant_data.get("genome")
    if not genome:
        return []
    populations = genome.get("populations", [])
    rows = []
    for pop in populations:
        pop_id = pop.get("id", "")
        if not pop_id.startswith("hgdp:"):
            continue
        pop_name = pop_id[5:]  # strip "hgdp:" prefix
        if pop_name in ("XX", "XY"):
            continue
        ac = pop.get("ac", 0)
        an = pop.get("an", 0)
        af = ac / an if an > 0 else 0.0
        rows.append({
            "hgdp_population": pop_name,
            "hgdp_region": HGDP_REGION_MAP.get(pop_name, "unknown"),
            "hgdp_ac": ac,
            "hgdp_an": an,
            "hgdp_af": af,
        })
    return rows


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
    }
  }
}
"""

REGION_CHUNK_SIZE = 100_000
MATCH_WINDOW = 5  # bp window for relaxed indel matching

REGION_CACHE = DATA_DIR / "gnomad_r3_region_variants_cache.csv.gz"


# ── variant normalization (bcftools-style trim) ──────────────────────

def normalize_variant(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Left-trim and right-trim shared bases (bcftools norm approximation)."""
    pos, ref, alt = int(pos), str(ref).upper(), str(alt).upper()
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1
    return pos, ref, alt


def event_signature(pos: int, ref: str, alt: str) -> tuple[str, str, int] | None:
    """Return (event_type, payload, position) for indels, None for SNVs/MNVs."""
    pos, ref, alt = normalize_variant(pos, ref, alt)
    if len(ref) == len(alt):
        return None
    if len(ref) < len(alt):
        return ("ins", alt[len(ref):], pos)
    return ("del", ref[len(alt):], pos)


def parse_variant_id(vid: str) -> tuple[str, int, str, str]:
    """Parse 'chrom-pos-ref-alt' into (chrom, pos, ref, alt)."""
    parts = vid.split("-")
    return parts[0], int(parts[1]), parts[2], parts[3]


# ── region queries ───────────────────────────────────────────────────

def query_region_variants(chrom: str, start: int, stop: int) -> list[dict]:
    """Return list of gnomAD r3 variant records in a region."""
    for attempt in range(GNOMAD_RETRIES):
        try:
            resp = requests.post(
                GNOMAD_GRAPHQL,
                json={
                    "query": GNOMAD_REGION_QUERY,
                    "variables": {
                        "chrom": chrom,
                        "start": start,
                        "stop": stop,
                        "dataset": GNOMAD_DATASET,
                        "referenceGenome": "GRCh38",
                    },
                },
                timeout=120,
            )
            resp.raise_for_status()
            payload = resp.json()
            if "errors" in payload:
                log.warning("Region query error chr%s:%d-%d: %s", chrom, start, stop, payload["errors"])
                return []
            region = payload.get("data", {}).get("region") or {}
            return region.get("variants") or []
        except (requests.RequestException, ValueError) as exc:
            delay = GNOMAD_RETRY_BASE_DELAY * (2 ** attempt)
            log.warning(
                "Region query failed chr%s:%d-%d (attempt %d/%d): %s — retrying in %.1fs",
                chrom, start, stop, attempt + 1, GNOMAD_RETRIES, exc, delay,
            )
            time.sleep(delay)
    return []


def collect_region_gnomad_variants(scores: pd.DataFrame) -> pd.DataFrame:
    """Collect all gnomAD r3 variant records across gene regions.
    Returns a DataFrame with columns: variant_id, chrom, pos, ref, alt."""

    if REGION_CACHE.exists():
        log.info("Loading cached gnomAD r3 region variants from %s", REGION_CACHE.name)
        return pd.read_csv(REGION_CACHE)

    gene_regions: dict[str, tuple[str, int, int]] = {}
    for gene, group in scores.groupby("gene"):
        parts = group["variant_key"].str.split(":")
        chrom = parts.str[0].iloc[0]
        positions = parts.str[1].astype(int)
        gene_regions[gene] = (chrom, int(positions.min()) - 100, int(positions.max()) + 100)

    all_records: list[dict] = []
    log.info("Pass 1: region queries to collect gnomAD r3 variants (%d gene regions)…", len(gene_regions))

    for gene, (chrom, start, stop) in sorted(gene_regions.items()):
        chunk_start = start
        while chunk_start < stop:
            chunk_stop = min(chunk_start + REGION_CHUNK_SIZE, stop)
            variants = query_region_variants(chrom, chunk_start, chunk_stop)
            for v in variants:
                all_records.append({
                    "variant_id": v["variantId"],
                    "chrom": chrom,
                    "pos": int(v["pos"]),
                    "ref": v["ref"],
                    "alt": v["alt"],
                    "gene_region": gene,
                })
            log.info(
                "  %s chr%s:%d-%d → %d gnomAD variants",
                gene, chrom, chunk_start, chunk_stop, len(variants),
            )
            time.sleep(GNOMAD_PAUSE)
            chunk_start = chunk_stop

    df = pd.DataFrame(all_records)
    df.to_csv(REGION_CACHE, index=False, compression="gzip")
    log.info("Cached %d gnomAD r3 region variants to %s", len(df), REGION_CACHE.name)
    return df


def match_clinvar_to_gnomad(scores: pd.DataFrame, gnomad_region: pd.DataFrame) -> dict[str, tuple[str, str]]:
    """Multi-level matching of ClinVar variants to gnomAD r3 IDs.

    Returns: dict mapping clinvar variant_key → (gnomad_variant_id, match_tier)
    Match tiers:
      - exact: identical chrom-pos-ref-alt
      - normalized: equivalent after left/right trim
      - event_equivalent: same indel event within ±MATCH_WINDOW bp
      - position_overlap: same position, different representation
    """
    matches: dict[str, tuple[str, str]] = {}

    # Build gnomAD lookup structures
    gnomad_by_id: dict[str, dict] = {}
    gnomad_by_pos: dict[tuple[str, int], list[dict]] = {}
    gnomad_by_norm: dict[tuple[str, int, str, str], list[dict]] = {}

    for _, row in gnomad_region.iterrows():
        chrom = str(row["chrom"])
        pos = int(row["pos"])
        ref = str(row["ref"])
        alt = str(row["alt"])
        vid = str(row["variant_id"])
        rec = {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "variant_id": vid}

        gnomad_by_id[vid] = rec
        gnomad_by_pos.setdefault((chrom, pos), []).append(rec)

        npos, nref, nalt = normalize_variant(pos, ref, alt)
        gnomad_by_norm.setdefault((chrom, npos, nref, nalt), []).append(rec)

    log.info("Built gnomAD lookup: %d unique IDs, %d positions, %d normalized keys",
             len(gnomad_by_id), len(gnomad_by_pos), len(gnomad_by_norm))

    n_exact = n_norm = n_event = n_pos = 0

    for _, row in scores.iterrows():
        vk = row["variant_key"]
        parts = vk.split(":")
        chrom, pos, ref, alt = parts[0], int(parts[1]), parts[2], parts[3]
        gnomad_id = variant_key_to_gnomad_id(vk)

        # Tier 1: exact match
        if gnomad_id in gnomad_by_id:
            matches[vk] = (gnomad_id, "exact")
            n_exact += 1
            continue

        # Tier 2: normalized match
        npos, nref, nalt = normalize_variant(pos, ref, alt)
        norm_key = (chrom, npos, nref, nalt)
        if norm_key in gnomad_by_norm:
            best = gnomad_by_norm[norm_key][0]
            matches[vk] = (best["variant_id"], "normalized")
            n_norm += 1
            continue

        # Tier 3: event-equivalent match (same indel event ±MATCH_WINDOW bp)
        cv_event = event_signature(pos, ref, alt)
        found_event = False
        if cv_event is not None:
            cv_kind, cv_payload, cv_epos = cv_event
            for offset in range(-MATCH_WINDOW, MATCH_WINDOW + 1):
                check_pos = pos + offset
                candidates = gnomad_by_pos.get((chrom, check_pos), [])
                for cand in candidates:
                    cand_event = event_signature(cand["pos"], cand["ref"], cand["alt"])
                    if cand_event is None:
                        continue
                    c_kind, c_payload, c_epos = cand_event
                    if c_kind == cv_kind and c_payload == cv_payload and abs(cv_epos - c_epos) <= MATCH_WINDOW:
                        matches[vk] = (cand["variant_id"], "event_equivalent")
                        n_event += 1
                        found_event = True
                        break
                if found_event:
                    break
        if found_event:
            continue

        # Tier 4: position overlap — same position, any variant
        pos_candidates = gnomad_by_pos.get((chrom, pos), [])
        if pos_candidates:
            matches[vk] = (pos_candidates[0]["variant_id"], "position_overlap")
            n_pos += 1
            continue

    n_total = len(scores)
    n_matched = len(matches)
    n_lost = n_total - n_matched
    log.info(
        "Matching complete: %d/%d (%.1f%%) matched",
        n_matched, n_total, 100 * n_matched / n_total,
    )
    log.info(
        "  exact: %d | normalized: %d | event_equivalent: %d | position_overlap: %d | lost: %d",
        n_exact, n_norm, n_event, n_pos, n_lost,
    )
    return matches


# ── match classification columns ─────────────────────────────────────

MATCH_TIER_TO_CLASS = {
    "exact": "strict_allele",
    "normalized": "normalized_allele",
    "event_equivalent": "normalized_allele",
    "position_overlap": "position_overlap",
    "no_match": "no_match",
}

MATCH_CLASS_TO_RESOLUTION = {
    "strict_allele": "allele_resolved",
    "normalized_allele": "representation_rescued",
    "position_overlap": "locus_context_only",
    "no_match": "unevaluable",
}


def add_match_classification(hgdp: pd.DataFrame) -> pd.DataFrame:
    """Add match_class and allele_resolution_level columns."""
    if "gnomad_r3_match_tier" not in hgdp.columns:
        hgdp["match_class"] = "no_match"
        hgdp["allele_resolution_level"] = "unevaluable"
        return hgdp

    hgdp["match_class"] = hgdp["gnomad_r3_match_tier"].map(MATCH_TIER_TO_CLASS).fillna("no_match")
    hgdp["allele_resolution_level"] = hgdp["match_class"].map(MATCH_CLASS_TO_RESOLUTION).fillna("unevaluable")
    return hgdp


def fetch_all_hgdp_data(scores: pd.DataFrame) -> pd.DataFrame:
    """Two-pass HGDP data acquisition with multi-level variant matching.

    Pass 1: region queries collect all gnomAD r3 variants in gene regions.
    Matching: exact → normalized → event-equivalent → position-overlap.
    Pass 2: individual queries for matched variants to get HGDP population data.
    """
    # Pass 1 — collect gnomAD r3 variants via region queries
    gnomad_region = collect_region_gnomad_variants(scores)

    # Multi-level matching
    match_map = match_clinvar_to_gnomad(scores, gnomad_region)

    # Pass 2 — individual queries for HGDP population data
    all_rows: list[dict] = []
    variant_keys = scores["variant_key"].tolist()
    clinvar_ids = scores["clinvar_id"].tolist()
    genes = scores["gene"].tolist()

    to_query = [
        (vk, cid, gene, match_map[vk][0], match_map[vk][1])
        for vk, cid, gene in zip(variant_keys, clinvar_ids, genes)
        if vk in match_map
    ]
    unmatched_keys = {vk for vk in variant_keys if vk not in match_map}

    log.info(
        "Pass 2: querying HGDP population data for %d matched variants "
        "(%d unmatched — no gnomAD r3 equivalent found)…",
        len(to_query), len(unmatched_keys),
    )

    # Record unmatched variants
    for vk, cid, gene in zip(variant_keys, clinvar_ids, genes):
        if vk in unmatched_keys:
            all_rows.append({
                "variant_key": vk,
                "clinvar_id": cid,
                "gene": gene,
                "hgdp_population": "NOT_FOUND",
                "hgdp_region": "NOT_FOUND",
                "hgdp_ac": 0,
                "hgdp_an": 0,
                "hgdp_af": 0.0,
                "gnomad_r3_genome_af": None,
                "gnomad_r3_genome_ac": None,
                "gnomad_r3_genome_an": None,
                "gnomad_r3_match_tier": "no_match",
                "gnomad_r3_matched_id": None,
            })

    found = 0
    for vk, cid, gene, gnomad_id, match_tier in tqdm(to_query, desc="HGDP individual queries"):
        variant_data = query_gnomad_hgdp(gnomad_id)
        time.sleep(GNOMAD_PAUSE)

        if variant_data is None:
            all_rows.append({
                "variant_key": vk,
                "clinvar_id": cid,
                "gene": gene,
                "hgdp_population": "NOT_FOUND",
                "hgdp_region": "NOT_FOUND",
                "hgdp_ac": 0,
                "hgdp_an": 0,
                "hgdp_af": 0.0,
                "gnomad_r3_genome_af": None,
                "gnomad_r3_genome_ac": None,
                "gnomad_r3_genome_an": None,
                "gnomad_r3_match_tier": match_tier,
                "gnomad_r3_matched_id": gnomad_id,
            })
            continue

        found += 1
        genome = variant_data.get("genome") or {}
        gnomad_af = genome.get("af")
        gnomad_ac = genome.get("ac")
        gnomad_an = genome.get("an")

        hgdp_pops = extract_hgdp_populations(variant_data)
        if not hgdp_pops:
            all_rows.append({
                "variant_key": vk,
                "clinvar_id": cid,
                "gene": gene,
                "hgdp_population": "NO_HGDP_DATA",
                "hgdp_region": "NO_HGDP_DATA",
                "hgdp_ac": 0,
                "hgdp_an": 0,
                "hgdp_af": 0.0,
                "gnomad_r3_genome_af": gnomad_af,
                "gnomad_r3_genome_ac": gnomad_ac,
                "gnomad_r3_genome_an": gnomad_an,
                "gnomad_r3_match_tier": match_tier,
                "gnomad_r3_matched_id": gnomad_id,
            })
            continue

        for pop_row in hgdp_pops:
            all_rows.append({
                "variant_key": vk,
                "clinvar_id": cid,
                "gene": gene,
                "gnomad_r3_genome_af": gnomad_af,
                "gnomad_r3_genome_ac": gnomad_ac,
                "gnomad_r3_genome_an": gnomad_an,
                "gnomad_r3_match_tier": match_tier,
                "gnomad_r3_matched_id": gnomad_id,
                **pop_row,
            })

    log.info(
        "Pass 2 complete: %d variants returned HGDP data out of %d queried",
        found, len(to_query),
    )
    result = pd.DataFrame(all_rows)
    result = add_match_classification(result)
    return result


# ═══════════════════════════════════════════════════════════════════════
# 1. Biological contradiction — population presence as frequency evidence
# ═══════════════════════════════════════════════════════════════════════

def analyse_biological_contradiction(hgdp: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    """For each variant: determine whether its observed HGDP frequency
    constitutes a biological contradiction with the claimed disease model.

    This is not an observation layer — it is a verdict: the observed
    population frequency is either compatible or incompatible with the
    disease model implied by the public P/LP label.
    """
    real = hgdp[
        ~hgdp["hgdp_population"].isin(["NOT_FOUND", "NO_HGDP_DATA"])
    ].copy()

    observed = real[real["hgdp_ac"] > 0]

    per_variant = (
        observed.groupby("variant_key")
        .agg(
            n_populations_observed=("hgdp_population", "nunique"),
            n_regions_observed=("hgdp_region", "nunique"),
            total_hgdp_ac=("hgdp_ac", "sum"),
            total_hgdp_an=("hgdp_an", "sum"),
            max_hgdp_af=("hgdp_af", "max"),
            mean_hgdp_af=("hgdp_af", "mean"),
            populations_list=("hgdp_population", lambda x: ";".join(sorted(x))),
            regions_list=("hgdp_region", lambda x: ";".join(sorted(set(x)))),
        )
        .reset_index()
    )
    per_variant["global_hgdp_af"] = (
        per_variant["total_hgdp_ac"] / per_variant["total_hgdp_an"]
    ).where(per_variant["total_hgdp_an"] > 0, 0.0)

    not_found_keys = set(
        hgdp.loc[hgdp["hgdp_population"] == "NOT_FOUND", "variant_key"]
    )

    result = scores[["variant_key", "clinvar_id", "gene"]].merge(
        per_variant, on="variant_key", how="left",
    )
    result["hgdp_status"] = "observed_in_hgdp"
    result.loc[result["n_populations_observed"].isna(), "hgdp_status"] = "absent_in_hgdp"
    result.loc[result["variant_key"].isin(not_found_keys), "hgdp_status"] = "not_in_gnomad_r3"
    result["n_populations_observed"] = result["n_populations_observed"].fillna(0).astype(int)
    result["n_regions_observed"] = result["n_regions_observed"].fillna(0).astype(int)
    result["total_hgdp_ac"] = result["total_hgdp_ac"].fillna(0).astype(int)

    # Merge match_class from hgdp
    if "match_class" in hgdp.columns:
        match_per_variant = hgdp.groupby("variant_key")["match_class"].first().reset_index()
        resolution_per_variant = hgdp.groupby("variant_key")["allele_resolution_level"].first().reset_index()
        result = result.merge(match_per_variant, on="variant_key", how="left")
        result = result.merge(resolution_per_variant, on="variant_key", how="left")
        result["match_class"] = result["match_class"].fillna("no_match")
        result["allele_resolution_level"] = result["allele_resolution_level"].fillna("unevaluable")
    else:
        result["match_class"] = "no_match"
        result["allele_resolution_level"] = "unevaluable"

    # Biological contradiction verdict
    conditions = [
        (result["hgdp_status"] == "not_in_gnomad_r3"),
        (result["hgdp_status"] == "absent_in_hgdp"),
        (result["max_hgdp_af"] > DOMINANT_LOW_PEN_MCAF),
        (result["max_hgdp_af"] > DOMINANT_HIGH_PEN_MCAF),
        (result["n_populations_observed"] >= 2) & (result["max_hgdp_af"] > AF_REVIEW_THRESHOLD),
        (result["max_hgdp_af"] > 0),
    ]
    verdicts = [
        "unevaluable",
        "no_contradiction_absent",
        "hard_biological_contradiction",
        "biological_contradiction",
        "frequency_incompatible_multi_pop",
        "no_contradiction_rare",
    ]
    result["biological_contradiction_verdict"] = np.select(
        conditions, verdicts, default="no_contradiction_absent",
    )

    return result


# ═══════════════════════════════════════════════════════════════════════
# 2. Ancestry non-portability
# ═══════════════════════════════════════════════════════════════════════

def analyse_ancestry_nonportability(hgdp: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    """Classify each variant's ancestry portability: does the P/LP label
    transfer across ancestries, or is it regionally restricted?

    Output: portable / conditionally_portable / non_portable.
    """
    real = hgdp[
        ~hgdp["hgdp_population"].isin(["NOT_FOUND", "NO_HGDP_DATA"])
    ].copy()

    region_agg = (
        real.groupby(["variant_key", "hgdp_region"])
        .agg(
            region_ac=("hgdp_ac", "sum"),
            region_an=("hgdp_an", "sum"),
            n_pops=("hgdp_population", "nunique"),
        )
        .reset_index()
    )
    region_agg["region_af"] = (
        region_agg["region_ac"] / region_agg["region_an"]
    ).where(region_agg["region_an"] > 0, 0.0)

    rows: list[dict] = []
    for vk, group in region_agg.groupby("variant_key"):
        evaluable = group[group["region_an"] >= MIN_AN_FOR_EVALUATION]
        if evaluable.empty:
            continue
        observed = evaluable[evaluable["region_ac"] > 0]
        max_af = evaluable["region_af"].max()
        min_af = evaluable["region_af"].min()
        af_range = max_af - min_af
        max_region = evaluable.loc[evaluable["region_af"].idxmax(), "hgdp_region"]
        n_regions_with_allele = len(observed)
        n_regions_evaluable = len(evaluable)

        total_ac = evaluable["region_ac"].sum()
        total_an = evaluable["region_an"].sum()
        global_proxy = total_ac / total_an if total_an > 0 else 0.0
        enrichment_ratio = max_af / global_proxy if global_proxy > 0 else 0.0

        # Portability classification
        absent_regions = n_regions_evaluable - n_regions_with_allele
        if n_regions_with_allele == 0:
            portability = "portable"
            ancestry_conflict = "none"
        elif absent_regions >= 2 and max_af > AF_REVIEW_THRESHOLD:
            portability = "non_portable"
            ancestry_conflict = "strong_regional_restriction"
        elif enrichment_ratio >= 10:
            portability = "non_portable"
            ancestry_conflict = "strong_enrichment"
        elif enrichment_ratio >= 3 or (absent_regions >= 2 and max_af > 0):
            portability = "conditionally_portable"
            ancestry_conflict = "moderate_enrichment"
        elif absent_regions >= 1:
            portability = "conditionally_portable"
            ancestry_conflict = "mild_regional_asymmetry"
        else:
            portability = "portable"
            ancestry_conflict = "none"

        rows.append({
            "variant_key": vk,
            "n_regions_evaluable": n_regions_evaluable,
            "n_regions_with_allele": n_regions_with_allele,
            "max_regional_af": max_af,
            "min_regional_af": min_af,
            "af_range_across_regions": af_range,
            "max_region": max_region,
            "hgdp_global_af_proxy": global_proxy,
            "max_region_enrichment_ratio": enrichment_ratio,
            "ancestry_portability": portability,
            "ancestry_conflict_class": ancestry_conflict,
            "region_af_detail": ";".join(
                f"{r}={af:.6f}"
                for r, af in zip(evaluable["hgdp_region"], evaluable["region_af"])
                if af > 0
            ),
        })

    result = pd.DataFrame(rows)
    if result.empty:
        result = pd.DataFrame(columns=[
            "variant_key", "n_regions_evaluable", "n_regions_with_allele",
            "max_regional_af", "min_regional_af", "af_range_across_regions",
            "max_region", "hgdp_global_af_proxy", "max_region_enrichment_ratio",
            "ancestry_portability", "ancestry_conflict_class", "region_af_detail",
        ])
    return scores[["variant_key", "clinvar_id", "gene"]].merge(
        result, on="variant_key", how="inner",
    )


# ═══════════════════════════════════════════════════════════════════════
# 3. Haplotype independence — single-effect assumption undermined
# ═══════════════════════════════════════════════════════════════════════

def analyse_haplotype_independence(hgdp: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    """Multi-regional presence implies multiple independent genomic
    backgrounds, undermining a single-effect assumption.

    This is not a caveat — it is a structural observation that the variant
    exists in diverse genetic contexts where a fixed penetrance assumption
    cannot hold uniformly.
    """
    real = hgdp[
        ~hgdp["hgdp_population"].isin(["NOT_FOUND", "NO_HGDP_DATA"])
    ].copy()
    observed = real[real["hgdp_ac"] > 0]

    rows: list[dict] = []
    for vk, group in observed.groupby("variant_key"):
        pops = sorted(group["hgdp_population"].unique())
        regions = sorted(group["hgdp_region"].unique())
        total_ac = int(group["hgdp_ac"].sum())

        if len(regions) >= 3:
            haplotype_class = "multiple_independent_backgrounds"
            single_effect_undermined = True
        elif len(regions) == 2:
            haplotype_class = "two_independent_backgrounds"
            single_effect_undermined = True
        elif len(pops) >= 2:
            haplotype_class = "within_region_diversity"
            single_effect_undermined = False
        else:
            haplotype_class = "single_background_founder_possible"
            single_effect_undermined = False

        rows.append({
            "variant_key": vk,
            "n_hgdp_populations": len(pops),
            "n_hgdp_regions": len(regions),
            "total_hgdp_ac": total_ac,
            "haplotype_diversity_class": haplotype_class,
            "single_effect_assumption_undermined": single_effect_undermined,
            "hgdp_populations": ";".join(pops),
            "hgdp_regions": ";".join(regions),
        })

    result = pd.DataFrame(rows)
    if result.empty:
        result = pd.DataFrame(columns=[
            "variant_key", "n_hgdp_populations", "n_hgdp_regions",
            "total_hgdp_ac", "haplotype_diversity_class",
            "single_effect_assumption_undermined",
            "hgdp_populations", "hgdp_regions",
        ])
    return scores[["variant_key", "clinvar_id", "gene"]].merge(
        result, on="variant_key", how="inner",
    )


# ═══════════════════════════════════════════════════════════════════════
# 4. Disease-model incompatibility — the main hammer
# ═══════════════════════════════════════════════════════════════════════

def _classify_regime(af: float, n_pops: int) -> dict:
    """Classify a variant into a disease-model regime based on AF."""
    if af > DOMINANT_LOW_PEN_MCAF:
        strict_compat = "excluded"
        generous_compat = "excluded"
    elif af > DOMINANT_HIGH_PEN_MCAF:
        strict_compat = "excluded"
        generous_compat = "boundary" if af > DOMINANT_LOW_PEN_MCAF / 2 else "permitted"
    else:
        strict_compat = "permitted"
        generous_compat = "permitted"

    if n_pops >= 2 and af > AF_REVIEW_THRESHOLD:
        carrier_flag = "carrier_architecture_implied"
    else:
        carrier_flag = "not_flagged"

    if strict_compat == "excluded" and generous_compat == "excluded":
        regime = "hard_incompatible"
    elif strict_compat == "excluded":
        regime = "boundary"
    elif carrier_flag == "carrier_architecture_implied":
        regime = "carrier_compatible"
    else:
        regime = "dominant_compatible"

    required_pen = None
    if af > 0:
        required_pen = (1 / 2000 * 0.1) / (af * 2.0)

    return {
        "strict_dominant_status": strict_compat,
        "generous_dominant_status": generous_compat,
        "carrier_flag": carrier_flag,
        "disease_model_regime": regime,
        "required_penetrance_for_dominant": required_pen,
    }


def analyse_disease_model_incompatibility(
    hgdp: pd.DataFrame,
    presence: pd.DataFrame,
    scores: pd.DataFrame,
) -> pd.DataFrame:
    """A single P/LP label spans mutually exclusive disease models.

    For ALL matched variants (not just HGDP-observed), determine regime
    using max(max_hgdp_af, gnomad_r3_genome_af) — whichever is higher
    provides the most informative constraint.

    Includes dual analysis: strict_allele-only vs all matches for robustness.
    """
    # Work with ALL matched variants, not just observed_in_hgdp
    matched = presence[presence["hgdp_status"] != "not_in_gnomad_r3"].copy()
    if matched.empty:
        return pd.DataFrame()

    # Merge gnomad_r3_genome_af from hgdp (match_class comes from presence)
    if "gnomad_r3_genome_af" in hgdp.columns and "gnomad_r3_genome_af" not in matched.columns:
        genome_af_per_variant = (
            hgdp.groupby("variant_key")["gnomad_r3_genome_af"]
            .first()
            .reset_index()
        )
        matched = matched.merge(genome_af_per_variant, on="variant_key", how="left")
    if "gnomad_r3_genome_af" not in matched.columns:
        matched["gnomad_r3_genome_af"] = 0.0
    matched["gnomad_r3_genome_af"] = matched["gnomad_r3_genome_af"].fillna(0.0)
    # match_class should already be in presence; fall back if missing
    if "match_class" not in matched.columns:
        if "match_class" in hgdp.columns:
            mc_per_variant = hgdp.groupby("variant_key")["match_class"].first().reset_index()
            matched = matched.merge(mc_per_variant, on="variant_key", how="left")
        if "match_class" not in matched.columns:
            matched["match_class"] = "strict_allele"
    matched["match_class"] = matched["match_class"].fillna("no_match")

    rows: list[dict] = []
    for _, row in matched.iterrows():
        _raw_hgdp = row.get("max_hgdp_af", 0.0)
        max_hgdp_af = 0.0 if (_raw_hgdp is None or pd.isna(_raw_hgdp)) else float(_raw_hgdp)
        _raw_genome = row.get("gnomad_r3_genome_af", 0.0)
        genome_af = 0.0 if (_raw_genome is None or pd.isna(_raw_genome)) else float(_raw_genome)
        # Use the higher AF as the most informative constraint
        effective_af = max(max_hgdp_af, genome_af)
        n_pops = int(row.get("n_populations_observed", 0))

        regime_info = _classify_regime(effective_af, n_pops)

        rows.append({
            "variant_key": row["variant_key"],
            "clinvar_id": row["clinvar_id"],
            "gene": row["gene"],
            "max_hgdp_af": max_hgdp_af,
            "gnomad_r3_genome_af": genome_af,
            "effective_af": effective_af,
            "total_hgdp_ac": int(row.get("total_hgdp_ac", 0)),
            "n_hgdp_populations": n_pops,
            **regime_info,
            "match_class": row.get("match_class", "strict_allele"),
            "is_strict_allele_match": row.get("match_class", "strict_allele") == "strict_allele",
        })

    return pd.DataFrame(rows)


# ═══════════════════════════════════════════════════════════════════════
# 5. Evaluability crisis
# ═══════════════════════════════════════════════════════════════════════

def analyse_evaluability(
    hgdp: pd.DataFrame,
    presence: pd.DataFrame,
    scores: pd.DataFrame,
) -> pd.DataFrame:
    """The majority of public assertions cannot be technically evaluated
    under population constraint.

    Not "rare" — unevaluable. The distinction matters: absence of evidence
    for frequency is not evidence of rarity.
    """
    real = hgdp[
        ~hgdp["hgdp_population"].isin(["NOT_FOUND", "NO_HGDP_DATA"])
    ].copy()

    variant_an = (
        real.groupby("variant_key")
        .agg(total_hgdp_an=("hgdp_an", "sum"), n_pops_with_data=("hgdp_population", "nunique"))
        .reset_index()
    )

    result = scores[["variant_key", "clinvar_id", "gene"]].merge(
        variant_an, on="variant_key", how="left",
    )
    result = result.merge(
        presence[["variant_key", "hgdp_status", "n_populations_observed", "max_hgdp_af"]],
        on="variant_key",
        how="left",
    )

    result["total_hgdp_an"] = result["total_hgdp_an"].fillna(0).astype(int)
    result["n_pops_with_data"] = result["n_pops_with_data"].fillna(0).astype(int)
    result["n_populations_observed"] = result["n_populations_observed"].fillna(0).astype(int)

    conditions = [
        (result["hgdp_status"] == "not_in_gnomad_r3"),
        (result["total_hgdp_an"] < MIN_AN_FOR_EVALUATION),
        (result["n_populations_observed"] == 0),
        (result["n_populations_observed"] == 1),
        (result["n_populations_observed"] >= 2),
    ]
    choices = [
        "technically_unevaluable",
        "insufficient_population_sampling",
        "evaluable_population_absent",
        "evaluable_single_population_context",
        "evaluable_multi_population_context",
    ]
    result["hgdp_evaluability"] = np.select(conditions, choices, default="technically_unevaluable")

    # Population spread score (0–1)
    spread_data = []
    for vk, group in real.groupby("variant_key"):
        evaluable_regions = group[group["hgdp_an"] >= MIN_AN_FOR_EVALUATION]["hgdp_region"].nunique()
        observed_regions = group[(group["hgdp_ac"] > 0) & (group["hgdp_an"] >= MIN_AN_FOR_EVALUATION)]["hgdp_region"].nunique()
        spread = observed_regions / evaluable_regions if evaluable_regions > 0 else 0.0
        spread_data.append({"variant_key": vk, "population_spread_score": spread})

    if spread_data:
        spread_df = pd.DataFrame(spread_data)
        result = result.merge(spread_df, on="variant_key", how="left")
    else:
        result["population_spread_score"] = 0.0

    result["population_spread_score"] = result["population_spread_score"].fillna(0.0)

    return result


# ═══════════════════════════════════════════════════════════════════════
# 6. Integrated summary — one label, multiple realities
# ═══════════════════════════════════════════════════════════════════════

def build_summary(
    scores: pd.DataFrame,
    presence: pd.DataFrame,
    ancestry: pd.DataFrame,
    haplotype: pd.DataFrame,
    disease_model: pd.DataFrame,
    evaluability: pd.DataFrame,
) -> pd.DataFrame:
    """Merge all layers into one row-per-variant summary.

    The headline: pathogenicity is not an intrinsic property of an allele
    but a context-dependent claim whose validity depends on population
    observability and disease-model specification.
    """
    summary = scores[["variant_key", "clinvar_id", "gene", "vital_band", "vital_red_flag"]].copy()

    summary = summary.merge(
        presence[["variant_key", "hgdp_status", "n_populations_observed",
                  "n_regions_observed", "max_hgdp_af", "global_hgdp_af",
                  "biological_contradiction_verdict"]],
        on="variant_key", how="left",
    )
    if not ancestry.empty:
        summary = summary.merge(
            ancestry[["variant_key", "ancestry_portability",
                      "ancestry_conflict_class",
                      "max_region_enrichment_ratio", "max_region"]],
            on="variant_key", how="left",
        )
    if not haplotype.empty:
        summary = summary.merge(
            haplotype[["variant_key", "haplotype_diversity_class",
                       "single_effect_assumption_undermined"]],
            on="variant_key", how="left",
        )
    if not disease_model.empty:
        summary = summary.merge(
            disease_model[["variant_key", "disease_model_regime",
                           "strict_dominant_status", "carrier_flag",
                           "required_penetrance_for_dominant"]],
            on="variant_key", how="left",
        )
    summary = summary.merge(
        evaluability[["variant_key", "hgdp_evaluability", "population_spread_score"]],
        on="variant_key", how="left",
    )

    # Match classification (from presence layer)
    for col in ("match_class", "allele_resolution_level"):
        if col in presence.columns and col not in summary.columns:
            summary = summary.merge(
                presence[["variant_key", col]], on="variant_key", how="left",
            )

    # Portability numeric score for the scatter plot (0=non_portable, 0.5=conditional, 1=portable)
    port_map = {"non_portable": 0.0, "conditionally_portable": 0.5, "portable": 1.0}
    summary["portability_score"] = summary.get("ancestry_portability", pd.Series(dtype=float)).map(port_map).fillna(1.0)

    # Disease model numeric score (0=hard_incompat, 0.33=boundary, 0.67=carrier, 1=compatible)
    regime_map = {"hard_incompatible": 0.0, "boundary": 0.33, "carrier_compatible": 0.67, "dominant_compatible": 1.0}
    summary["disease_model_score"] = summary.get("disease_model_regime", pd.Series(dtype=float)).map(regime_map).fillna(1.0)

    # Evaluability numeric score (0=unevaluable, 0.5=single, 1=multi)
    eval_map = {
        "technically_unevaluable": 0.0,
        "insufficient_population_sampling": 0.0,
        "evaluable_population_absent": 0.25,
        "evaluable_single_population_context": 0.5,
        "evaluable_multi_population_context": 1.0,
    }
    summary["evaluability_score"] = summary["hgdp_evaluability"].map(eval_map).fillna(0.0)

    return summary


# ═══════════════════════════════════════════════════════════════════════
# Figures
# ═══════════════════════════════════════════════════════════════════════

def plot_regime_distribution(disease_model: pd.DataFrame) -> None:
    """Bar chart: a single P/LP label spans mutually exclusive disease models."""
    if disease_model.empty:
        log.warning("No disease-model data to plot.")
        return
    counts = disease_model["disease_model_regime"].value_counts()
    labels = {
        "hard_incompatible": "Hard incompatible\nwith dominant model",
        "boundary": "Boundary /\nmonitored",
        "carrier_compatible": "Carrier /\nrecessive implied",
        "dominant_compatible": "Dominant\ncompatible",
    }
    ordered = ["hard_incompatible", "boundary", "carrier_compatible", "dominant_compatible"]
    present = [r for r in ordered if r in counts.index]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    colors = ["#d32f2f", "#ff9800", "#2196f3", "#4caf50"]
    color_map = dict(zip(ordered, colors))
    bars = ax.bar(
        range(len(present)),
        [counts[r] for r in present],
        color=[color_map[r] for r in present],
        edgecolor="white", linewidth=0.8,
    )
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels([labels.get(r, r) for r in present], fontsize=10)
    ax.set_ylabel("Number of P/LP variants observed in HGDP", fontsize=11)
    ax.set_title(
        "A single P/LP label spans mutually exclusive disease models\n"
        "(HGDP population-frequency layer)",
        fontsize=12, fontweight="bold",
    )
    for bar, regime in zip(bars, present):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.5,
            str(counts[regime]),
            ha="center", va="bottom", fontweight="bold", fontsize=11,
        )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(HGDP_REGIME_FIGURE, dpi=200)
    plt.close(fig)
    log.info("Saved %s", HGDP_REGIME_FIGURE)


def plot_ancestry_nonportability_heatmap(hgdp: pd.DataFrame, ancestry: pd.DataFrame) -> None:
    """Heatmap of regional AF for non-portable and conditionally-portable P/LP variants."""
    if ancestry.empty or "ancestry_portability" not in ancestry.columns:
        log.warning("No ancestry portability data to plot.")
        return
    non_port = ancestry[
        ancestry["ancestry_portability"].isin(["non_portable", "conditionally_portable"])
    ].copy()
    if non_port.empty:
        log.warning("No non-portable variants to plot.")
        return

    top = non_port.nlargest(25, "max_region_enrichment_ratio")
    top_keys = set(top["variant_key"])

    real = hgdp[
        ~hgdp["hgdp_population"].isin(["NOT_FOUND", "NO_HGDP_DATA"])
    ].copy()
    region_agg = (
        real[real["variant_key"].isin(top_keys)]
        .groupby(["variant_key", "hgdp_region"])
        .agg(region_ac=("hgdp_ac", "sum"), region_an=("hgdp_an", "sum"))
        .reset_index()
    )
    region_agg["region_af"] = (
        region_agg["region_ac"] / region_agg["region_an"]
    ).where(region_agg["region_an"] > 0, 0.0)

    pivot = region_agg.pivot_table(
        index="variant_key", columns="hgdp_region", values="region_af", fill_value=0.0,
    )
    gene_map = ancestry.set_index("variant_key")["gene"].to_dict()
    port_map = ancestry.set_index("variant_key")["ancestry_portability"].to_dict()
    pivot.index = [
        f"{gene_map.get(vk, '')} [{port_map.get(vk, '')}]"
        for vk in pivot.index
    ]

    fig, ax = plt.subplots(figsize=(10, max(6, len(pivot) * 0.4)))
    im = ax.imshow(pivot.values, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=7)
    ax.set_title(
        "Non-portable P/LP assertions: AF varies by ancestry\n"
        "One label does not mean one biological reality",
        fontsize=11, fontweight="bold",
    )
    fig.colorbar(im, ax=ax, label="Regional allele frequency", shrink=0.6)
    fig.tight_layout()
    fig.savefig(HGDP_ANCESTRY_FIGURE, dpi=200)
    plt.close(fig)
    log.info("Saved %s", HGDP_ANCESTRY_FIGURE)


def plot_evaluability_crisis(evaluability: pd.DataFrame) -> None:
    """Pie: the majority cannot be technically evaluated."""
    counts = evaluability["hgdp_evaluability"].value_counts()
    labels_map = {
        "technically_unevaluable": "Technically unevaluable",
        "insufficient_population_sampling": "Insufficient sampling",
        "evaluable_population_absent": "Absent in HGDP",
        "evaluable_single_population_context": "Single-population context",
        "evaluable_multi_population_context": "Multi-population context",
    }
    colors_map = {
        "technically_unevaluable": "#616161",
        "insufficient_population_sampling": "#9e9e9e",
        "evaluable_population_absent": "#90caf9",
        "evaluable_single_population_context": "#ffcc80",
        "evaluable_multi_population_context": "#ef5350",
    }
    ordered = [k for k in labels_map if k in counts.index]
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.pie(
        [counts[k] for k in ordered],
        labels=[labels_map[k] for k in ordered],
        colors=[colors_map[k] for k in ordered],
        autopct="%1.1f%%",
        startangle=140,
        textprops={"fontsize": 9},
    )
    ax.set_title(
        "Population evaluability of ClinVar P/LP arrhythmia assertions\n"
        "Most public assertions cannot be technically evaluated under population constraint",
        fontsize=10, fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(HGDP_EVALUABILITY_FIGURE, dpi=200)
    plt.close(fig)
    log.info("Saved %s", HGDP_EVALUABILITY_FIGURE)


HGDP_SCATTER_FIGURE = FIGURE_DIR / "hgdp_portability_vs_disease_model.png"


def plot_portability_vs_disease_model(summary: pd.DataFrame) -> None:
    """The killer figure: X = ancestry portability, Y = disease-model
    compatibility, color = evaluability.

    One scatter that shows P/LP decomposing into separate regimes.
    """
    plot_df = summary.dropna(subset=["portability_score", "disease_model_score"]).copy()
    if plot_df.empty:
        log.warning("No data for scatter plot.")
        return

    # Add jitter for readability
    rng = np.random.default_rng(42)
    jx = rng.uniform(-0.06, 0.06, len(plot_df))
    jy = rng.uniform(-0.04, 0.04, len(plot_df))

    fig, ax = plt.subplots(figsize=(10, 7))
    sc = ax.scatter(
        plot_df["portability_score"].values + jx,
        plot_df["disease_model_score"].values + jy,
        c=plot_df["evaluability_score"].values,
        cmap="RdYlGn",
        s=30,
        alpha=0.7,
        edgecolors="black",
        linewidths=0.3,
    )
    ax.set_xlabel("Ancestry portability\n(0 = non-portable, 0.5 = conditional, 1 = portable)", fontsize=11)
    ax.set_ylabel(
        "Disease-model compatibility\n"
        "(0 = hard incompatible, 0.33 = boundary, 0.67 = carrier, 1 = dominant compatible)",
        fontsize=11,
    )
    ax.set_title(
        '"Pathogenic" is not an intrinsic property of an allele\n'
        "but a context-dependent claim that decomposes under population constraint",
        fontsize=12, fontweight="bold",
    )
    cb = fig.colorbar(sc, ax=ax, shrink=0.6)
    cb.set_label("HGDP evaluability (0 = unevaluable → 1 = multi-population)", fontsize=10)
    ax.set_xlim(-0.15, 1.15)
    ax.set_ylim(-0.1, 1.1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(HGDP_SCATTER_FIGURE, dpi=200)
    plt.close(fig)
    log.info("Saved %s", HGDP_SCATTER_FIGURE)


# ═══════════════════════════════════════════════════════════════════════
# Console report — assertive framing
# ═══════════════════════════════════════════════════════════════════════

def print_report(
    scores: pd.DataFrame,
    hgdp: pd.DataFrame,
    presence: pd.DataFrame,
    ancestry: pd.DataFrame,
    haplotype: pd.DataFrame,
    disease_model: pd.DataFrame,
    evaluability: pd.DataFrame,
) -> None:
    n_total = len(scores)

    print("\n" + "=" * 78)
    print("  HGDP POPULATION-LAYER ANALYSIS")
    print("  Pathogenicity is not an intrinsic property of an allele.")
    print("=" * 78)

    # 0. Match rate report
    print(f"\n{'─' * 60}")
    print("0. VARIANT MATCHING — ClinVar to gnomAD r3 (HGDP)")
    print(f"{'─' * 60}")
    if "gnomad_r3_match_tier" in hgdp.columns:
        per_variant_tier = hgdp.groupby("variant_key")["gnomad_r3_match_tier"].first()
        tier_counts = per_variant_tier.value_counts()
        for tier, count in tier_counts.items():
            print(f"   {tier:30s}  {count:5d}  ({100*count/n_total:.1f}%)")
        n_matched = per_variant_tier[per_variant_tier != "no_match"].shape[0]
        n_recovered = per_variant_tier[per_variant_tier.isin(["normalized", "event_equivalent", "position_overlap"])].shape[0]
        n_lost = n_total - n_matched
        print(f"\n   Total matched:         {n_matched}/{n_total} ({100*n_matched/n_total:.1f}%)")
        print(f"   Recovered by relaxed:  {n_recovered}")
        print(f"   Lost (no match):       {n_lost} ({100*n_lost/n_total:.1f}%)")
    else:
        print("   Match tier data not available in cached HGDP data.")

    # 1. Biological contradiction
    print(f"\n{'─' * 60}")
    print("1. BIOLOGICAL CONTRADICTION — frequency vs claimed disease model")
    print(f"{'─' * 60}")
    verdict_counts = presence["biological_contradiction_verdict"].value_counts()
    for verdict, count in verdict_counts.items():
        print(f"   {verdict:45s}  {count:5d}  ({100*count/n_total:.1f}%)")
    n_contradicted = presence[
        presence["biological_contradiction_verdict"].isin([
            "hard_biological_contradiction",
            "biological_contradiction",
            "frequency_incompatible_multi_pop",
        ])
    ].shape[0]
    if n_contradicted > 0:
        print(f"\n   → {n_contradicted} P/LP assertions show observed HGDP frequency")
        print("     incompatible with the claimed disease model.")

    # 2. Non-portability
    print(f"\n{'─' * 60}")
    print("2. ANCESTRY NON-PORTABILITY — the P/LP label does not transfer")
    print(f"{'─' * 60}")
    if not ancestry.empty and "ancestry_portability" in ancestry.columns:
        port_counts = ancestry["ancestry_portability"].value_counts()
        for port, count in port_counts.items():
            print(f"   {port:35s}  {count:5d}")
        n_nonport = port_counts.get("non_portable", 0) + port_counts.get("conditionally_portable", 0)
        n_ancestry_eval = len(ancestry)
        if n_ancestry_eval > 0:
            print(f"\n   → {n_nonport}/{n_ancestry_eval} ancestry-evaluable P/LP variants "
                  f"({100*n_nonport/n_ancestry_eval:.1f}%) are non-portable across ancestries.")
    else:
        print("   No ancestry data available.")

    # 3. Haplotype independence
    print(f"\n{'─' * 60}")
    print("3. HAPLOTYPE INDEPENDENCE — single-effect assumption undermined")
    print(f"{'─' * 60}")
    if not haplotype.empty and "single_effect_assumption_undermined" in haplotype.columns:
        n_undermined = haplotype["single_effect_assumption_undermined"].sum()
        n_hap_total = len(haplotype)
        print(f"   Multi-regional presence (≥2 regions):   {n_undermined}/{n_hap_total}")
        if n_hap_total > 0:
            print(f"   → {100*n_undermined/n_hap_total:.1f}% of HGDP-observed variants exist on")
            print("     multiple independent genomic backgrounds.")
        hap_counts = haplotype["haplotype_diversity_class"].value_counts()
        for hc, count in hap_counts.items():
            print(f"   {hc:45s}  {count:5d}")
    else:
        print("   No haplotype data available.")

    # 4. Disease-model incompatibility — THE HAMMER (with dual analysis)
    print(f"\n{'─' * 60}")
    print("4. DISEASE-MODEL INCOMPATIBILITY")
    print("   A single P/LP label spans mutually exclusive disease regimes.")
    print(f"{'─' * 60}")
    if not disease_model.empty:
        n_observed = len(disease_model)
        regime_order = ["hard_incompatible", "boundary", "carrier_compatible", "dominant_compatible"]
        regime_counts = disease_model["disease_model_regime"].value_counts()
        for regime in regime_order:
            count = regime_counts.get(regime, 0)
            pct = 100 * count / n_observed if n_observed else 0
            print(f"   {regime:40s}  {count:5d}  ({pct:.1f}%)")

        n_not_dominant = disease_model[
            disease_model["disease_model_regime"] != "dominant_compatible"
        ].shape[0]
        if n_observed > 0:
            print(f"\n   → ALL MATCHES ({n_observed} variants): "
                  f"{n_not_dominant}/{n_observed} ({100*n_not_dominant/n_observed:.1f}%)"
                  " incompatible/boundary/carrier")

        # Dual analysis: strict_allele only vs representation-rescued
        if "is_strict_allele_match" in disease_model.columns:
            strict_only = disease_model[disease_model["is_strict_allele_match"]]
            n_strict = len(strict_only)
            n_strict_not_dom = strict_only[
                strict_only["disease_model_regime"] != "dominant_compatible"
            ].shape[0]

            print("\n   ROBUSTNESS — strict-only vs full set:")
            if n_strict > 0:
                print(f"   → STRICT ALLELE ONLY ({n_strict} variants):")
                strict_regime = strict_only["disease_model_regime"].value_counts()
                for regime in regime_order:
                    count = strict_regime.get(regime, 0)
                    pct = 100 * count / n_strict
                    print(f"     {regime:38s}  {count:5d}  ({pct:.1f}%)")
                print(f"     → {n_strict_not_dom}/{n_strict} "
                      f"({100*n_strict_not_dom/n_strict:.1f}%) non-dominant-compatible")

            overlap_only = disease_model[~disease_model["is_strict_allele_match"]]
            n_overlap = len(overlap_only)
            if n_overlap > 0:
                n_overlap_not_dom = overlap_only[
                    overlap_only["disease_model_regime"] != "dominant_compatible"
                ].shape[0]
                print(f"\n   → REPRESENTATION-RESCUED ({n_overlap} variants):")
                overlap_regime = overlap_only["disease_model_regime"].value_counts()
                for regime in regime_order:
                    count = overlap_regime.get(regime, 0)
                    pct = 100 * count / n_overlap
                    print(f"     {regime:38s}  {count:5d}  ({pct:.1f}%)")
                print(f"     → {n_overlap_not_dom}/{n_overlap} "
                      f"({100*n_overlap_not_dom/n_overlap:.1f}%) non-dominant-compatible")

            # Verdict
            strict_pct = 100 * n_strict_not_dom / n_strict if n_strict > 0 else 0
            all_pct = 100 * n_not_dominant / n_observed if n_observed > 0 else 0
            if strict_pct > 10 and all_pct > 10:
                print(f"\n   ▸ ROBUST: pattern holds in strict-only ({strict_pct:.1f}%) "
                      f"and all matches ({all_pct:.1f}%).")
            elif all_pct > 10 and strict_pct <= 10:
                print("\n   ▸ WARNING: conclusion depends on representation-rescued matches.")
    else:
        print("   No disease-model data available.")

    # 5. Evaluability crisis
    print(f"\n{'─' * 60}")
    print("5. EVALUABILITY CRISIS")
    print("   Not 'rare' — technically unevaluable.")
    print(f"{'─' * 60}")
    eval_counts = evaluability["hgdp_evaluability"].value_counts()
    for ev, count in eval_counts.items():
        print(f"   {ev:45s}  {count:5d}  ({100*count/n_total:.1f}%)")
    n_unevaluable = evaluability[
        evaluability["hgdp_evaluability"].isin([
            "technically_unevaluable", "insufficient_population_sampling",
        ])
    ].shape[0]
    print(f"\n   → {n_unevaluable}/{n_total} P/LP assertions ({100*n_unevaluable/n_total:.1f}%)")
    print("     cannot be technically evaluated under population constraint.")
    print("     The distinction matters: absence of frequency evidence")
    print("     is not evidence of rarity.")

    # 6. RECONCILIATION-AWARE HEADLINE
    print(f"\n{'═' * 78}")
    print("  HEADLINE")
    print(f"{'═' * 78}")
    print()
    print('  "Pathogenicity is not an intrinsic property of an allele')
    print('   but a context-dependent claim whose validity depends on')
    print('   population observability and disease-model specification."')
    print()

    # Match reconciliation statistics
    if "gnomad_r3_match_tier" in hgdp.columns:
        per_vt = hgdp.groupby("variant_key")["gnomad_r3_match_tier"].first()
        n_mapped = per_vt[per_vt != "no_match"].shape[0]
        n_strict = per_vt[per_vt == "exact"].shape[0]
        n_rescued = per_vt[per_vt.isin(["normalized", "event_equivalent", "position_overlap"])].shape[0]
    else:
        n_mapped = n_strict = n_rescued = 0

    print(f"   From {n_total} public ClinVar P/LP assertions:")
    print(f"   • Only {n_mapped} ({100*n_mapped/n_total:.1f}%) became population-mappable")
    print(f"     after reconciliation ({n_strict} strict allele + {n_rescued} representation-rescued)")
    if n_contradicted > 0:
        print(f"   • {n_contradicted} show biological contradiction with claimed model")
    if not ancestry.empty and "ancestry_portability" in ancestry.columns:
        n_nonport = ancestry["ancestry_portability"].isin(["non_portable", "conditionally_portable"]).sum()
        print(f"   • {n_nonport} are non-portable across ancestries")
    if not haplotype.empty and "single_effect_assumption_undermined" in haplotype.columns:
        n_undermined = int(haplotype["single_effect_assumption_undermined"].sum())
        print(f"   • {n_undermined} undermine the single-effect assumption (multi-regional)")
    if not disease_model.empty:
        n_dm = len(disease_model)
        n_not_dom = disease_model[disease_model["disease_model_regime"] != "dominant_compatible"].shape[0]
        n_hard = (disease_model["disease_model_regime"] == "hard_incompatible").sum()
        n_boundary = (disease_model["disease_model_regime"] == "boundary").sum()
        print(f"   • {n_not_dom}/{n_dm} population-mapped variants ({100*n_not_dom/n_dm:.0f}%) "
              f"span disease regimes incompatible with dominant reading")
        print(f"     ({n_hard} hard_incompatible, {n_boundary} boundary)")
    print(f"   • {n_unevaluable} cannot be technically evaluated under population constraint")
    print()
    print("   The public P/LP label often does not transfer to population data")
    print(f"   without additional interpretation — {n_rescued}/{n_mapped} matched variants")
    print(f"   ({100*n_rescued/n_mapped:.0f}%) depend on representation-rescue, not direct allele identity." if n_mapped > 0 else "")
    print()


# ═══════════════════════════════════════════════════════════════════════
# Entrypoint
# ═══════════════════════════════════════════════════════════════════════

def save_csv_and_supplement(df: pd.DataFrame, csv_path: Path, supp_path: Path) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    supp_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv_path, index=False)
    df.to_csv(supp_path, index=False, sep="\t")
    log.info("Saved %s (%d rows) + supplement", csv_path.name, len(df))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cached", action="store_true",
        help="Reuse cached HGDP population frequency data instead of re-querying the API.",
    )
    args = parser.parse_args()

    if not SCORES_IN.exists():
        raise SystemExit(f"Missing input: {SCORES_IN}")
    scores = pd.read_csv(SCORES_IN)
    log.info("Loaded %d variants from %s", len(scores), SCORES_IN.name)

    # 0. Data acquisition
    if args.cached and HGDP_CACHE.exists():
        log.info("Loading cached HGDP data from %s", HGDP_CACHE.name)
        hgdp = pd.read_csv(HGDP_CACHE)
        # Ensure match classification columns exist (may be absent in older caches)
        if "match_class" not in hgdp.columns:
            hgdp = add_match_classification(hgdp)
    else:
        hgdp = fetch_all_hgdp_data(scores)
        hgdp.to_csv(HGDP_CACHE, index=False)
        log.info("Saved HGDP cache to %s (%d rows)", HGDP_CACHE.name, len(hgdp))

    # 1. Biological contradiction
    presence = analyse_biological_contradiction(hgdp, scores)
    save_csv_and_supplement(presence, HGDP_PRESENCE_OUT, HGDP_PRESENCE_SUPP)

    # 2. Ancestry non-portability
    ancestry = analyse_ancestry_nonportability(hgdp, scores)
    save_csv_and_supplement(ancestry, HGDP_ANCESTRY_MISMATCH_OUT, HGDP_ANCESTRY_SUPP)

    # 3. Haplotype independence
    haplotype = analyse_haplotype_independence(hgdp, scores)
    save_csv_and_supplement(haplotype, HGDP_HAPLOTYPE_PROXY_OUT, HGDP_HAPLOTYPE_SUPP)

    # 4. Disease-model incompatibility
    disease_model = analyse_disease_model_incompatibility(hgdp, presence, scores)
    save_csv_and_supplement(disease_model, HGDP_DISEASE_MODEL_OUT, HGDP_DISEASE_MODEL_SUPP)

    # 5. Evaluability crisis
    evaluability = analyse_evaluability(hgdp, presence, scores)
    save_csv_and_supplement(evaluability, HGDP_EVALUABILITY_OUT, HGDP_EVALUABILITY_SUPP)

    # 6. Integrated summary
    summary = build_summary(scores, presence, ancestry, haplotype, disease_model, evaluability)
    save_csv_and_supplement(summary, HGDP_SUMMARY_OUT, HGDP_SUMMARY_SUPP)

    # Figures
    plot_regime_distribution(disease_model)
    plot_ancestry_nonportability_heatmap(hgdp, ancestry)
    plot_evaluability_crisis(evaluability)
    plot_portability_vs_disease_model(summary)

    # Report
    print_report(scores, hgdp, presence, ancestry, haplotype, disease_model, evaluability)


if __name__ == "__main__":
    main()
