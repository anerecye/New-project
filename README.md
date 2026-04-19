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
python src/summarize_gnomad_merge.py
python src/make_publication_figure.py
```

The gnomAD extraction scripts require `bcftools` and internet access to query public gnomAD VCFs.

## Figure

`figures/figure_final_publication.png` compares log10 allele-frequency distributions for gnomAD population variants and ClinVar clinical variants. Dashed vertical lines mark approximate rare-disease frequency thresholds at `1e-5` and `1e-4`. The clinical distribution is concentrated below these thresholds, illustrating the frequency-constrained nature of clinically cataloged pathogenic variation.
