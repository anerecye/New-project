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
