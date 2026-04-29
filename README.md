# VITAL Actionability Routing Audit

VITAL is a post-classification actionability routing layer for public
pathogenic and likely pathogenic (`P/LP`) variant labels.

Public P/LP labels are often reused as actionability proxies, but they are not
portable actionability units. VITAL tests whether a flattened public label can
survive minimal population, mechanism, inheritance, and evaluability
constraints before downstream systems treat it as direct-actionable.

Core chain:

`label -> baseline actionability -> VITAL constraint -> rerouting -> preserved expert truth`

## What VITAL Does

- Treats public `P/LP` as a trigger for routing, not as a direct-actionability certificate.
- Applies allele-resolved evaluability, ancestry-aware frequency, and disease-model constraints.
- Converts a label-driven baseline into reviewable routing categories.
- Separates pathogenicity classification from downstream actionability.
- Preserves expert-curated positives as compatible or review-level when constraints support that route.

## What VITAL Does Not Do

VITAL is not a variant classifier.

VITAL does not claim that nonpass variants are benign, incorrect, or clinically
irrelevant. A nonpass route means that direct actionability cannot be inferred
from the flattened label alone.

VITAL is not a patient-outcome or malpractice detector. The repository audits
variant-level routing logic, not realized clinical harm.

## Routing Categories

| Route | Meaning | Operational output |
| --- | --- | --- |
| `VITAL_OK` | Label remains compatible with the tested actionability model. | Proceed with standard expert interpretation. |
| `CHECK_POPMAX` | Label may remain valid, but direct actionability requires ancestry-aware frequency review. | Review. |
| `CHECK_MODEL` | Label may remain valid, but inheritance, phase, or mechanism is not portable from the flattened label. | Model-specific routing or model repair. |
| `MODEL_CONFLICT` | Flattened P/LP label is incompatible with the tested dominant high-penetrance model. | Do not direct-route as dominant actionability. |
| `EVAL_LIMITED` | No allele-resolved evaluability; direct actionability cannot be justified from population evidence. | Defer direct actionability until representation/callability/evidence is repaired. |

Legacy export states remain available for compatibility:

- `VITAL-1`: allele-resolved, popmax-compatible, model-defined
- `VITAL-2`: locus-level or representation-limited
- `VITAL-X`: recessive or phase-defined route
- `VITAL-ALERT`: allele-resolved but blocked by population/model tension
- `VITAL-0`: representation or normalization unresolved

## Key Result

Across six disease domains, the mean VITAL nonpass rate is `87.7%`: most
label-driven actionability decisions require rerouting once minimal constraints
are restored.

In the inherited-arrhythmia cohort:

- `1,512/1,731` public P/LP labels (`87.3%`) leave the direct-actionable baseline.
- `1,397/1,731` (`80.7%`) are evaluation-limited.
- `77/1,731` (`4.4%`) require population/frequency review.
- `38/1,731` (`2.2%`) require model-specific rerouting.
- In the high-review subset, `309/365` (`84.7%`) still leave the direct-actionable baseline.

This effect is not driven by low-confidence submissions alone. Among 73
evaluable expert-curated P/LP positives, `69/73` (`94.5%`) remain `VITAL_OK` or
review-level and only `4/73` (`5.5%`) enter hard model-conflict routing.

## Repository Structure

- `src/`: analysis and data-preparation scripts
- `notebooks/`: quick Colab demo notebook
- `examples/`: tiny end-to-end VCF and annotation examples
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

## Full Gene-Set CLI

For cohort-level routing summaries, use the cached full mode:

```bash
python run_vital.py --mode full --genes "MYBPC3,MYH7" --pop gnomAD
```

This returns:

- baseline actionable cohort size
- `VITAL-1` / `VITAL_OK` compatible count and rate
- total alert/reroute burden
- alert subtypes (`evaluation-limiting`, `VITAL-ALERT`, `VITAL-X`, `SV/CNV-required`)
- top alerted variants
- an expert-curated benchmark line for signal-to-noise context

Optional outputs:

```bash
python run_vital.py --mode full --genes "MYBPC3,MYH7" --pop gnomAD \
  --output examples/mybpc3_myh7_variants.csv \
  --summary-output examples/mybpc3_myh7_summary.csv
```

For the cross-domain meta-summary and gene-level intersection plot:

```bash
python src/run_vital_standard_meta_analysis.py
```

Key outputs:

- `data/processed/vital_standard_meta_summary.csv`
- `data/processed/vital_standard_meta_overall.csv`
- `data/processed/vital_standard_gene_problem_matrix.csv`
- `data/processed/vital_standard_gene_problem_intersections.csv`
- `figures/vital_standard_gene_problem_upset.png`

## Baseline And Routing Audit

The formal routing audit is the manuscript-facing layer:

```bash
python src/run_vital_routing_validation.py
```

Key outputs:

- `data/processed/vital_formal_baseline_decision_model.csv`
- `data/processed/vital_routing_validation_calls.csv`
- `data/processed/vital_counterfactual_decision_audit_summary.csv`
- `data/processed/vital_simulated_cds_alert_summary.csv`
- `data/processed/vital_actionability_discordance_audit.csv`
- `data/processed/vital_minimal_repair_logic.csv`
- `data/processed/vital_reason_code_definitions.csv`
- `figures/vital_decision_disruption.png`
- `figures/vital_routing_validation.png`

Baseline rules:

| Rule | Meaning |
| --- | --- |
| `B1_LABEL_DRIVEN_BASELINE` | Public P/LP may be reused as direct-actionable unless additional context is encoded. |
| `B2_ACTIONABILITY_CONTEXT_BASELINE` | If P/LP appears in an action-linked gene/context, the baseline route is direct-actionable. |
| `B3_NO_IMPLIED_BENIGNITY` | Absence of `VITAL_OK` does not imply benignity. It means direct actionability is unsupported by the flattened label alone. |

## Simulated CDS Alert Layer

The simulated clinical decision-support layer converts VITAL routes into alert
classes:

| Alert | VITAL routes | Message |
| --- | --- | --- |
| `YELLOW` | `CHECK_POPMAX`, `CHECK_MODEL` | Actionability requires population/mechanism review. |
| `ORANGE` | `EVAL_LIMITED` | Allele-resolved population evidence unavailable; do not infer population compatibility. |
| `RED` | `MODEL_CONFLICT` | Dominant high-penetrance model incompatible under current constraints. |

Headline cached metrics:

- all arrhythmia P/LP alert rate: `1512/1731` (`87.3%`)
- high-review arrhythmia P/LP alert rate: `309/365` (`84.7%`)
- expert-curated evaluable hard-conflict rate: `4/73` (`5.5%`)
- expert-curated evaluable review-level retention: `69/73` (`94.5%`)

## Autopsy x de novo Counterfactual Audit

The autopsy/de novo layer models phenotype-null negative-autopsy settings where
genetic findings can become the main explanatory anchor.

```bash
python src/run_vital_autopsy_denovo_audit.py
```

The audit asks whether `P/LP + de novo` can be safely promoted to cause-of-death
attribution without allele-level evaluability, ancestry-aware frequency, and
disease-model checks. It treats de novo status as evidence within a coherent
model, not as a substitute for the model.

Headline cached outputs:

- baseline false causal attributions: `7,975/9,348` causal attributions (`85.3%`)
- VITAL false causal attributions: `10/48` supported causal attributions (`20.8%`)
- prevented false attributions: `7,965`
- confirmed de novo override errors under the label baseline: `371`
- gold-standard GoF/DN dominant preservation: `42/42 VITAL_OK`, `0 MODEL_CONFLICT`

Key outputs:

- `data/processed/vital_autopsy_denovo_model_inputs.csv`
- `data/processed/vital_autopsy_denovo_variant_design.csv`
- `data/processed/vital_autopsy_denovo_case_calls.csv`
- `data/processed/vital_autopsy_denovo_summary.csv`
- `data/processed/vital_autopsy_denovo_by_denovo_status.csv`
- `data/processed/vital_autopsy_denovo_by_mechanism.csv`
- `data/processed/vital_autopsy_denovo_penetrance_sensitivity.csv`
- `data/processed/vital_autopsy_denovo_denovo_rate_sensitivity.csv`
- `data/processed/vital_autopsy_denovo_mcaf_sensitivity.csv`
- `data/processed/vital_autopsy_denovo_gold_standard_preservation.csv`
- `figures/vital_autopsy_denovo_decision_flow.png`
- `figures/vital_autopsy_denovo_false_attribution.png`
- `figures/vital_autopsy_denovo_override.png`
- `figures/vital_autopsy_denovo_mechanism.png`

## Certification Fields

Recommended reproducibility fields for downstream exports:

- source database
- source version/date
- variant normalization status
- genome build
- evaluability tier
- population frequency source
- popmax source
- actionability domain
- disease model tested
- inheritance model tested
- mechanism class
- review status
- routing output
- reason code

Machine-readable reason codes:

- `EVAL_NO_ALLELE_MATCH`
- `EVAL_LOCI_ONLY`
- `POPMAX_EXCEEDS_MCAF`
- `RECESSIVE_COMPATIBLE_NOT_DOMINANT`
- `EXPERT_REVIEW_RETAINED`
- `MODEL_CONFLICT_DOMINANT_HP`
- `CHECK_LOW_COUNT`
- `CHECK_MECHANISM_AMBIGUOUS`

## VITAL Annotation MVP

This repository now includes a compact annotation export layer for downstream
VCF or variant-table pipelines. It converts the cached VITAL outputs into
machine-readable variant-level columns:

- `VITAL_certification`
- `VITAL_alert`
- `VITAL_public_use`
- `VITAL_sv_required`
- `VITAL_flag`
- `VITAL_evaluability`
- `VITAL_regime`
- `VITAL_action`
- `VITAL_reason`
- `VITAL_threshold`

The MVP is intentionally narrow:

- only popmax/global frequency tension is used
- routing/certification is emitted explicitly instead of collapsing everything into one pass/fail flag
- Tier 2 and still-unevaluable variants are separated via `VITAL_evaluability`
  and `VITAL_certification` instead of being over-called

This is actionability routing, not clinical reclassification.

**`VITAL_OK` does not mean clinically benign.** It means: no
population-frequency conflict detected under the current MVP threshold/model in
usable allele-resolved space.

**`CHECK_POPMAX` does not mean reclassification.** It means:
ancestry-aware AF exceeds the review trigger and requires model-specific
interpretation.

**`MODEL_CONFLICT` means:** observed AF is incompatible with the tested
unqualified dominant high-penetrance model.

**Required precondition:** input variants must be normalized to `GRCh38` with
`bcftools norm` before VITAL lookup. Non-normalized input may produce false
non-matches.

Build the lookup tables once:

```bash
python src/run_vital_annotation_mvp.py
```

This writes:

- `data/processed/arrhythmia_vital_mvp_lookup.tsv`
- `data/processed/arrhythmia_vital_mvp_annovar.tsv`

The normalized lookup table uses `chr / pos / ref / alt`. The ANNOVAR-style
table uses:

`Chr Start End Ref Alt VITAL_evaluability VITAL_flag VITAL_regime VITAL_certification VITAL_alert VITAL_public_use VITAL_sv_required VITAL_popmax_af VITAL_global_af VITAL_threshold VITAL_reason`

Normalize the input VCF before lookup:

```bash
bcftools norm -f GRCh38.fa -m -both input.vcf.gz -Oz -o input.norm.vcf.gz
```

Then annotate the normalized VCF:

```bash
python src/run_vital_annotation_mvp.py \
  --input examples/input.example.vcf \
  --lookup data/processed/arrhythmia_vital_mvp_lookup.tsv \
  --output examples/output.vital.tsv \
  --annovar-export-output examples/annovar_export.hg38_vital.txt
```

Repository examples:

- `examples/input.example.vcf`
- `examples/output.vital.tsv`
- `examples/annovar_export.hg38_vital.txt`

Expected output:

| VITAL_evaluability | VITAL_flag | VITAL_regime | VITAL_reason | Meaning |
| --- | --- | --- | --- | --- |
| `TIER_1` | `OK` | `OK` | `EXPERT_REVIEW_RETAINED` | usable allele-resolved AF was found and did not cross the MVP review trigger |
| `TIER_1` | `CHECK_POPMAX` | `BOUNDARY` | `POPMAX_EXCEEDS_MCAF` | ancestry-aware AF crossed the review trigger and needs model-specific review |
| `TIER_1` | `CHECK_POPMAX` | `BOUNDARY` | `RECESSIVE_COMPATIBLE_NOT_DOMINANT` | dominant actionability is unsupported, but recessive/carrier logic may remain coherent |
| `TIER_1` | `MODEL_CONFLICT` | `DOMINANT_INCOMPATIBLE` | `MODEL_CONFLICT_DOMINANT_HP` | observed AF is incompatible with the tested unqualified dominant high-penetrance model |
| `TIER_2` | `.` | `.` | `EVAL_LOCI_ONLY` | locus or representation context exists, but usable allele-resolved AF does not |
| `UNEVALUABLE` | `.` | `.` | `EVAL_NO_ALLELE_MATCH` | no reliable normalized allele-resolved lookup result was available |

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

### Routing And Clinical Decision Layer

This layer reframes VITAL as a downstream actionability router rather than a
standalone reclassification engine. It compares a ClinVar-only
`ROUTE_PLP_ACTIONABLE` baseline against VITAL constraint routes, quantifies
`actionability at risk`, formalizes the baseline decision model, runs a
counterfactual decision audit, summarizes burden by workflow context, and checks
expert-panel calibration plus representative case vignettes.

```bash
python src/run_vital_routing_validation.py
```

Key outputs:

- `data/processed/vital_formal_baseline_decision_model.csv`
- `data/processed/vital_counterfactual_decision_audit_summary.csv`
- `data/processed/vital_routing_validation_summary.csv`
- `data/processed/vital_routing_validation_context_summary.csv`
- `data/processed/vital_routing_validation_expert_concordance.csv`
- `data/processed/vital_routing_validation_case_vignettes.csv`
- `data/processed/vital_actionability_discordance_audit.csv`
- `data/processed/vital_simulated_cds_alert_summary.csv`
- `data/processed/vital_minimal_repair_logic.csv`
- `data/processed/vital_reason_code_definitions.csv`
- `figures/vital_decision_disruption.png`
- `figures/vital_routing_validation.png`

### Autopsy x de novo Counterfactual Decision Audit

This layer stress-tests negative-autopsy, phenotype-null interpretation. The
label-driven baseline treats public `P/LP` in an arrhythmia gene as sufficient
for probable genetic cause attribution and treats confirmed de novo status as
high confidence. VITAL requires evaluability, popmax compatibility, mechanism
coherence, and de novo compatibility before causal attribution is supported.

```bash
python src/run_vital_autopsy_denovo_audit.py
```

Primary endpoints:

- false causal attribution
- de novo override error
- prevented false attribution
- CHECK/DEFER burden
- MODEL_CONFLICT routing
- gold-standard dominant positive preservation

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

This layer formalizes the "too common for high-penetrance Mendelian disease"
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
