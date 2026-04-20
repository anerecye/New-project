# ClinVar vs gnomAD Allele Frequency Analysis

This project compares allele frequencies from ClinVar arrhythmia-associated variants with population frequencies from gnomAD. The analysis focuses on inherited arrhythmia genes and asks whether clinically reported pathogenic variants occupy a frequency-constrained subset of population variation.

## Key Result

ClinVar pathogenic and likely pathogenic variants are strongly shifted toward rare allele frequencies relative to the queried gnomAD population variation. The main figure shows that clinical databases preferentially capture variants below common rare-disease frequency thresholds, consistent with a frequency-constrained bias in clinical variant databases.

## Repository Structure

- `src/`: analysis and data-preparation scripts
- `data/processed/`: small processed CSV/TSV outputs used by the analysis
- `figures/`: PNG figures generated from processed outputs
- `requirements.txt`: Python dependencies

Large raw inputs such as full ClinVar downloads, gnomAD VCFs, VCF indexes, local environments, and generated cache files are intentionally excluded.

## How to Run

Create an environment and install dependencies:

```bash
python -m venv .venv
python -m pip install -r requirements.txt
```

Regenerate the summary tables and gnomAD AF distribution figure from the included processed merge:

```bash
python src/create_metric_tables.py
python src/summarize_gnomad_merge.py
python src/make_publication_figure.py
```

To rerun the full pipeline from raw ClinVar data, place `variant_summary.txt.gz` in `data/raw/` and run:

```bash
python src/clean_clinvar.py
python src/qc_variant_keys.py
python src/collapse_variant_records.py
python src/export_gnomad_lookup_input.py
python src/make_bed_for_gnomad.py
bash src/split_bed_by_chromosome.sh
bash src/extract_gnomad_exomes.sh
python src/merge_gnomad.py --gnomad data/processed/gnomad_subset_with_header.tsv
python src/create_metric_tables.py
python src/summarize_gnomad_merge.py
python src/make_publication_figure.py
```

The gnomAD extraction scripts require `bcftools` and internet access to query public gnomAD VCFs.

## API-Based ClinVar and gnomAD Pipeline

An alternative end-to-end API runner is available for reproducing the ClinVar
P/LP to gnomAD r4 matching workflow without downloading full gnomAD VCFs:

```bash
python src/arrhythmia_variant_pipeline.py --email anerecye@gmail.com
```

This script:

- fetches ClinVar records for the 20 inherited-arrhythmia genes through NCBI
  Entrez and filters to germline pathogenic or likely pathogenic variants;
- logs ClinVar P/LP records discarded because no GRCh38 VCF coordinate is
  available in the VCV XML;
- matches variants to the gnomAD GraphQL API using exact `chrom-pos-ref-alt`
  variant IDs, with exponential backoff for transient GraphQL failures, and
  stores both exome and genome AF/AC/AN fields when available;
- labels unmatched variants as `allele_discordance` when another gnomAD allele
  exists at the same position, otherwise `no_gnomad_record`;
- classifies matched AF values by `1e-5` and `1e-4` thresholds;
- writes per-gene Fisher exact tests with BH correction and warns when
  per-gene AF-covered counts are small;
- writes `run_qc_summary.csv`, including the
  `allele_discordance:no_gnomad_record` split and coverage-query status counts;
- samples unmatched positions for a coverage check and writes four figures.

NCBI requires an email address for Entrez calls. You can pass it with `--email`
or set `ENTREZ_EMAIL`. Full gnomAD matching is intentionally rate-limited and
can take 30-60 minutes. CSV caches are accompanied by `.meta.json` files that
record gene list, gnomAD dataset, thresholds, and related settings. Reuse cached
outputs on repeated runs:

```bash
python src/arrhythmia_variant_pipeline.py --skip-clinvar --skip-gnomad --skip-coverage
```

Use `--strict-cache-metadata` to fail instead of warning if a cache was produced
with a different gene list, dataset, or cache mode.

Large ClinVar genes can produce interrupted XML downloads. The API runner uses
`--clinvar-batch-size 50` by default and will retry and split failed ClinVar
batches before skipping any individual VariationID.

The API runner also supports space-separated gene lists and output prefixes,
which is useful for running an arrhythmia set and a non-arrhythmia control set
side by side:

```bash
python src/arrhythmia_variant_pipeline.py --email anerecye@gmail.com --output-prefix arrhythmia --clinvar-batch-size 50 --genes KCNQ1 KCNH2 SCN5A KCNE1 KCNE2 RYR2 CASQ2 TRDN CALM1 CALM2 CALM3 ANK2 SCN4B KCNJ2 HCN4 CACNA1C CACNB2 CACNA2D1 AKAP9 SNTA1
python src/arrhythmia_variant_pipeline.py --email anerecye@gmail.com --output-prefix control --clinvar-batch-size 50 --genes BRCA1 BRCA2 MLH1 APC
```

With prefixes, outputs are written as files such as
`data/processed/control_gnomad_matched.csv`,
`data/processed/control_run_qc_summary.csv`, and
`figures/control_af_distribution.png`.

## Advanced Stratified and QC Analyses

Run the advanced analysis layer after the API pipeline has produced a
`*_gnomad_matched.csv` file:

```bash
python src/advanced_variant_analyses.py --output-prefix arrhythmia
```

This adds five manuscript-oriented analysis blocks:

- ancestry-stratified gnomAD exome AFs for AFR, AMR, ASJ, EAS, FIN, MID, NFE,
  SAS, and remaining groups, including computed `popmax_af` and
  globally-rare/population-enriched candidates;
- gnomAD v4.1 exome-vs-genome sensitivity for exact matches, including
  duplication-specific genome presence checks;
- reclassification-risk screening for ClinVar P/LP variants with AF above
  `1e-5` or `1e-4`, review strength, submitter-count signal, and a
  `potential_misclassification_signal` flag that requires at least one
  threshold-exceeding frequency source with `AC >= 20`;
- gene-level functional features, including LOF, missense,
  splice-or-intronic, in-frame, synonymous, and unresolved classes;
- non-overlap structure for `allele_discordance` and `no_gnomad_record`
  variants, including gene enrichment, indel/SNV balance, hotspot-position
  markers, coverage-query summaries, and a CpG-transition proxy.

Main CSV outputs include:

- `data/processed/arrhythmia_population_af.csv`
- `data/processed/arrhythmia_population_af_long.csv`
- `data/processed/arrhythmia_population_af_summary.csv`
- `data/processed/arrhythmia_population_af_outliers.csv`
- `data/processed/arrhythmia_exome_genome_af_comparison.csv`
- `data/processed/arrhythmia_exome_genome_af_summary.csv`
- `data/processed/arrhythmia_reclassification_risk.csv`
- `data/processed/arrhythmia_reclassification_risk_summary.csv`
- `data/processed/arrhythmia_gene_variant_type_summary.csv`
- `data/processed/arrhythmia_gene_lof_missense_af_summary.csv`
- `data/processed/arrhythmia_gene_non_overlap_summary.csv`
- `data/processed/arrhythmia_non_overlap_feature_summary.csv`
- `data/processed/arrhythmia_non_overlap_gene_enrichment.csv`
- `data/processed/arrhythmia_non_overlap_variant_level.csv`

Machine-readable supplementary tables are generated as TSV files:

```bash
python src/generate_supplementary_tables.py --output-prefix arrhythmia
```

This writes `supplementary_tables/Supplementary_Table_S2_reclassification_candidates.tsv`
and `supplementary_tables/Supplementary_Table_S5_gene_non_overlap_summary.tsv`.

Main figures include:

- `figures/arrhythmia_population_af_distribution.png`
- `figures/arrhythmia_population_af_outliers.png`
- `figures/arrhythmia_reclassification_risk.png`
- `figures/arrhythmia_variant_type_by_gene.png`
- `figures/arrhythmia_non_overlap_features.png`
- `figures/arrhythmia_non_overlap_genes.png`

Population AF queries are cached in `arrhythmia_population_af.csv`. Reuse the
cache without new gnomAD requests:

```bash
python src/advanced_variant_analyses.py --output-prefix arrhythmia --no-fetch-population-af
```

Exome/genome AF queries are cached in
`arrhythmia_exome_genome_af_comparison.csv`. Reuse this cache with
`--no-fetch-exome-genome-af`; use `--force-exome-genome-fetch` to refresh it.

For existing ClinVar caches that were generated before the XML parser included
`NumberOfSubmitters`, submitter counts in the risk table are lower-bound
inferences from `review_status` and are marked with
`submitter_count_source=review_status_lower_bound`. Fresh ClinVar fetches from
the updated API pipeline write exact `submitter_count` and `submission_count`
fields from the VCV XML. The CpG field is a C>T/G>A transition proxy because
the current processed inputs do not include flanking reference sequence.

## Reviewer-Facing Interpretation Notes

The BRCA1, BRCA2, MLH1, and APC comparator is retained only as a
non-arrhythmia cancer-predisposition comparator. It should not be interpreted
as a clean technical control for arrhythmia genes because founder effects and
reduced penetrance can allow some cancer predisposition P/LP variants to appear
in population data at non-zero frequency. The observed exact-match difference
between the arrhythmia set and this comparator (20.2% versus 9.3% in the
current cached run) is therefore descriptive; a cardiomyopathy-gene comparator
would be a stronger like-for-like control.

`allele_discordance` means that gnomAD has at least one record at the same
GRCh38 position as the ClinVar P/LP variant, but not the same REF/ALT allele.
This category is a representation and locus-context flag, not a pathogenicity
classification. It may reflect a different alternate allele at the same base,
left/right alignment around indels, repeat-context representation differences,
or local mapping artifacts. Variants without any same-position gnomAD record
remain labeled `no_gnomad_record`.

Reclassification candidates are screened by AF and explicitly annotated by
allele-count support. `Supplementary_Table_S2_reclassification_candidates.tsv`
keeps the full 115 AF-threshold candidates for auditability, but the
`potential_misclassification_signal` flag now requires weak review evidence and
at least one threshold-exceeding global or popmax frequency with `AC >= 20`.
In the current cache, 9/115 candidates pass this AC support filter and 3/115
also have weak review evidence.

The KCNH2 submission-date trend should be treated as descriptive only. The
median-split Fisher test is not conventionally significant in the current
cache (`p=0.081918`), so the manuscript should avoid claiming a definite
temporal improvement. The year-by-year distribution is provided in
`data/processed/arrhythmia_kcnh2_submission_date_summary.csv`.

TRDN counts depend on which pipeline/table is used. The legacy reviewer-QC
table reports 6/23 AF-covered TRDN variants above `1e-5` (26.1%,
BH `q=0.030765`), whereas the current API pipeline reports 7/31 (22.6%,
BH `q=0.00049`). Manuscript text and tables should not mix these denominators.

For the strongest observed non-overlap gene signal, run a focused KCNH2
diagnostic pass:

```bash
python src/investigate_kcnh2_non_overlap.py --output-prefix arrhythmia
```

This writes `arrhythmia_kcnh2_non_overlap_summary.csv`,
`arrhythmia_kcnh2_non_overlap_stratified_tests.csv`,
`arrhythmia_kcnh2_non_overlap_type_contribution.csv`,
`arrhythmia_kcnh2_duplication_size_summary.csv`,
`arrhythmia_kcnh2_duplication_sensitivity.csv`,
`arrhythmia_kcnh2_repeat_overlap_summary.csv`,
`arrhythmia_kcnh2_gc_control_summary.csv`,
`arrhythmia_kcnh2_submission_date_summary.csv`,
`arrhythmia_kcnh2_non_overlap_kb_bins.csv`,
`arrhythmia_kcnh2_non_overlap_variant_context.csv`, and
`figures/arrhythmia_kcnh2_non_overlap_diagnostics.png`. By default it fetches
the hg38 chr7 KCNH2 interval from the UCSC sequence API for local GC windows,
UCSC `rmsk` and `simpleRepeat` tracks for repeat overlap, and ClinVar VCV
metadata for submission dates. Pass `--skip-reference-sequence`,
`--skip-repeatmasker`, or `--skip-clinvar-dates` to disable those remote
queries.

## Figure

`figures/figure_final_publication.png` compares log10 allele-frequency distributions for gnomAD population variants and ClinVar clinical variants. Dashed vertical lines mark approximate rare-disease frequency thresholds at `1e-5` and `1e-4`. The clinical distribution is concentrated below these thresholds, illustrating the frequency-constrained nature of clinically cataloged pathogenic variation.

## Main Metric Outputs

- `data/processed/main_metrics.csv`: headline counts, outlier fractions, maximum AF, and median AF.
- `data/processed/af_category_counts.csv`: counts by ultra-rare, rare, and common AF categories.
- `data/processed/clinvar_vs_gnomad_af_summary.csv`: descriptive AF comparison for ClinVar-matched variants versus the queried gnomAD population subset.
- `data/processed/outlier_variants.csv`: all ClinVar variants with `AF > 1e-5`.
- `data/processed/gene_af_stats.csv`: per-gene AF summary statistics.
- `data/processed/outlier_counts_by_gene.csv`: per-gene `AF > 1e-5` outlier counts, denominators, within-gene fractions, and each gene's share of all outliers.
- `data/processed/top_variants.csv`: top 10 variants ranked by AF, not the full outlier list.
- `data/processed/exomes_genomes_af_comparison.csv`: exome/genome AF comparison when overlapping genome inputs are available. In the current small processed example, no overlapping exome/genome variant keys are available.
- `data/processed/match_qc_summary.csv`: exact-match and unmatched-category counts for ClinVar-to-gnomAD matching.
- `data/processed/match_qc_by_gene.csv`: per-gene exact-match and unmatched-category counts.
- `data/processed/gene_outlier_enrichment_tests.csv`: Fisher exact tests for per-gene enrichment of `AF > 1e-5` variants.
- `figures/af_distribution_updated.png`: threshold-marked AF histogram.
- `figures/af_categories.png`: AF category barplot.
