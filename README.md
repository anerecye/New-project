# ClinVar Arrhythmia Evaluability Audit

This repository studies how public ClinVar pathogenic and likely pathogenic
(P/LP) assertions behave when they are re-used outside their original submission
context. The current framing centers two questions:

1. For how many public arrhythmia assertions can allele-resolved population
   constraint actually be recovered?
2. When that constraint is available, what disease-model heterogeneity becomes
   visible beneath a shared public P/LP label?

The main manuscript now argues that a public pathogenic label is not, by itself,
a portable disease state. In this project, ancestry-aware population frequency,
allele-resolved evaluability, and mechanism reveal hidden interpretation
regimes beneath the same exported label.

Technical filenames, script names, and data columns still retain the historical
`vital_` prefix for compatibility with cached outputs and prior runs. The
user-facing documentation below follows the newer evaluability and
frequency-tension framing.

## Key Result

Across 1,731 ClinVar P/LP variants in 20 inherited arrhythmia genes, only
357 variants (20.6%) achieved allele-resolved population context after strict
matching, trim-aware reconciliation, decomposition-aware reconciliation, and
reference-based normalization with `bcftools norm`. The remaining
1,326 variants (76.6%) retained only locus- or region-level context without
allele-resolved AF, and 48 (2.8%) remained unevaluable.

Within the Tier 1 subset with usable AF data (`n = 334`), global AF alone
flagged 13 variants above the `1e-5` review threshold, whereas a
global-or-popmax screen flagged 115. Global-only review therefore missed
102 of 115 ancestry-aware frequency alerts (88.7%).

Those 115 frequency-tension variants resolved into three distinct regimes:

- `1` hard incompatibility with an unqualified dominant high-penetrance reading
- `76` boundary or monitored-interpretation cases
- `38` recessive or carrier-compatible cases

## Repository Structure

- `src/`: analysis and data-preparation scripts
- `notebooks/`: quick Colab demo notebook
- `data/processed/`: cached CSV/TSV outputs used by the analyses
- `data/examples/`: tiny demo CSVs for the cached quick runner
- `figures/`: PNG figures generated from processed outputs
- `requirements.txt`: Python dependencies
- `Dockerfile`: minimal reproducible Python environment for cached/demo analyses

Large raw inputs such as full ClinVar downloads, gnomAD VCFs, indexes, and
generated cache files are intentionally excluded. Derived context tables are
committed under `data/processed/`.

## Try the Cached Demo in 2 Minutes

Open the quick demo notebook:

[Cached scoring demo: score ClinVar variants in 2 minutes](https://colab.research.google.com/github/anerecye/New-project/blob/main/notebooks/vital_demo_colab.ipynb)

The notebook has three paths:

- single VCV: edit one line and run `python run_vital.py --vcv VCV000440850`
- CSV batch: run `python run_vital.py --input data/examples/sample_variants.csv`
- demo table: browse `data/processed/vital_top_suspicious.csv`

The quick demo uses cached score tables and makes no ClinVar or gnomAD API
calls.

## Main Analyses

### Component-Consistency Audit

This audit keeps the original expert-specified score structure and asks whether
the component balance behaves sensibly in data without turning the pipeline into
a fitted prediction model.

```bash
python src/run_vital_score_calibration.py
```

Key outputs:

- `data/processed/vital_cross_disease_3000_restricted_calibration_model_metrics.csv`
- `data/processed/vital_cross_disease_3000_restricted_calibration_feature_dominance.csv`
- `data/processed/vital_cross_disease_3000_restricted_calibration_coefficients.csv`
- `figures/vital_cross_disease_3000_restricted_calibration_models.png`

### Biological Contrast Layer

This layer joins cached score tables to LOEUF and Human Protein Atlas
heart-muscle expression and turns the highest-priority cases into compact
biological contrast summaries.

```bash
python src/run_vital_biological_contrast.py
```

### Severe-Annotation Discordance and Expert Comparator

This audit asks whether severe-looking ClinVar P/LP consequence annotations are
sufficient to override ancestry-aware population-frequency contradiction. The
answer here is no: severe annotation does not, by itself, rescue coherence.

```bash
python src/run_vital_external_truth_claim.py
```

Headline outputs include:

- `92/262` AF-observed LOF/splice arrhythmia assertions with popmax/global
  AF `> 1e-5`
- `57/197` severe annotations still discordant after excluding `CASQ2/TRDN`
- `102/115` ancestry-aware alerts hidden by global-AF-only review

### Real Reclassification Audit

This audit separates real downstream classification change from workflow-burden
metrics. It is intentionally narrow and does not claim to be a broad predictor
of future reclassification.

```bash
python src/run_vital_real_reclassification_audit.py
```

### Temporal Robustness

The temporal analysis rescored archived ClinVar snapshots for January 2023,
January 2024, and April 2026 to check whether the urgent-review queue stays
sparse and whether queue composition changes over time.

```bash
python src/run_vital_temporal_robustness.py
```

### Validation Readiness

This package prepares a blinded expert-review pilot, keeps a separate analysis
key, and reports AC reliability as strata rather than as a single magic cutoff.

```bash
python src/run_vital_validation_readiness.py
```

### Provenance Credibility Filter

This sensitivity layer removes the most obviously fragile extreme-frequency
public records and asks what remains true after that conservative filter.

```bash
python src/run_vital_provenance_credibility_sensitivity.py
```

### Tiered Match Reconciliation

This audit asks how much of the ClinVar-to-gnomAD gap is recoverable by
representation-aware cleanup versus how much is a structural infrastructure
limitation.

```bash
python src/run_vital_bcftools_reference_normalization.py
python src/run_vital_match_reconciliation.py
```

Current cached outputs:

- strict exact exome AF-observed space: `334 / 1,731` (`19.3%`)
- strict exact any-dataset space: `350 / 1,731` (`20.2%`)
- exact/equivalent reconciliation space: `357 / 1,731` (`20.6%`)
- locus/regional context without exact/equivalent AF: `1,326 / 1,731` (`76.6%`)
- still unevaluable after both passes: `48 / 1,731` (`2.8%`)

Tier 2 is not random residue:

- same-locus allele discordance: `638 / 1,326` (`48.1%`)
- no same-locus record, regional-only context: `688 / 1,326` (`51.9%`)
- indels/duplications/insertions: `640 / 1,326` Tier 2 variants (`48.3%`)
  versus `106 / 357` Tier 1 variants (`29.7%`)

Reference-normalization audit:

- all `1,731` arrhythmia variants were normalized against a local GRCh38
  primary-assembly FASTA using `bcftools norm`
- coordinate or allele changes introduced by reference-based normalization:
  `0 / 1,731`
- multiallelic decomposition introduced no additional exact/equivalent recoveries

### Clinical Decision-Risk Proxy

ClinVar does not expose patient counts or realized downstream management
changes, so this layer estimates conservative public-exposure upper bounds.

```bash
python src/run_vital_clinical_decision_risk_proxy.py
```

### Clinical Decision Projection

This layer contrasts the clinical pathway implied by an unqualified Mendelian
reading with the pathway implied by a state-aware reading of the same public
label. It identifies interpretive exposure, not measured patient harm.

```bash
python src/run_vital_clinical_decision_projection.py
```

### Clinical Action Contexts and Regime Distribution

This layer converts the 115 frequency-tension variants into two manuscript-ready
summaries: where they fall in real clinical-action contexts and how they split
into label-state regimes.

```bash
python src/run_vital_clinical_action_contexts.py
```

Headline outputs:

- `103 / 115` tension variants fall in an action-linked arrhythmia core
- `60 / 115` lie in drug-restriction contexts
- `59 / 115` lie in device or intensive-surveillance contexts
- regime split: `1` hard dominant incompatibility, `76` boundary/monitoring,
  `38` recessive/carrier-compatible

### Frequency Constraint and Maximum Credible AF

This layer formalizes the “too common for high-penetrance Mendelian disease”
argument using maximum credible allele frequency, disease-prevalence logic,
Bayesian ACMG-style odds, recessive carrier logic, and allele-count reliability.

```bash
python src/run_vital_frequency_constraint_analysis.py
```

## How to Run

Create an environment and install dependencies:

```bash
python -m venv .venv
python -m pip install -r requirements.txt
```

Regenerate the summary tables and the gnomAD AF distribution figure from the
included processed merge:

```bash
python src/create_metric_tables.py
python src/summarize_gnomad_merge.py
python src/make_publication_figure.py
```

To rerun the API-based public-data workflow:

```bash
python src/arrhythmia_variant_pipeline.py --email your@email.example
```

This script:

- fetches ClinVar records for the 20 inherited-arrhythmia genes through NCBI
  Entrez
- matches variants to the gnomAD GraphQL API using exact `chrom-pos-ref-alt`
  variant IDs
- stores exome and genome AF/AC/AN fields when available
- labels unresolved variants as `allele_discordance` or `no_gnomad_record`

NCBI requires an email address for Entrez calls. You can pass it with `--email`
or set `ENTREZ_EMAIL`.

## Reproducibility and API Independence

The current manuscript is tied to an explicit data freeze:

- April 21, 2026 for the current ClinVar/gnomAD analysis
- gnomAD `v4.1.1` for exome and genome frequency queries
- archived January 2023 ClinVar `variant_summary` for historical validation

External APIs are used for acquisition, but downstream reproducibility does not
depend on those APIs remaining unchanged. Intermediate files are cached under
`data/processed/`, and the local-reference normalization audit uses
`data/external/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa`.

## Docker

Build and run the minimal container:

```bash
docker build -t arrhythmia-evaluability .
docker run --rm arrhythmia-evaluability
docker run --rm arrhythmia-evaluability python run_vital.py --vcv VCV000440850
```

The examples above use a neutral image tag, while script names and cached file
paths retain their historical prefixes for compatibility.
