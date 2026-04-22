# Severe Annotation Does Not Override Population-Frequency Tension in ClinVar Arrhythmia Variants

## Abstract

**Purpose:** Severe variant annotations can coexist with ancestry-aware population-frequency contradiction. We tested whether ClinVar pathogenic/likely pathogenic (P/LP) arrhythmia assertions remain frequency-consistent when evaluated with popmax, allele count (AC), variant representation, and review support rather than global AF alone.

**Methods:** We collapsed ClinVar P/LP records across 20 inherited arrhythmia genes to 1,731 unique variants and cross-referenced them with gnomAD v4.1.1 exome data. Exact allele matches, allele-discordant sites, no-record states, global AF, popmax AF, AC, variant type, and ClinVar review strength were retained. We developed VITAL (Variant Interpretation Tension and Allele Load) as an explainable re-review prioritization framework, not an automated reclassification engine.

**Results:** Only 350/1,731 variants (20.2%) had exact gnomAD allele matches, and 334 had usable AF evidence. Global AF alone identified 13 variants above AF >1e-5, whereas popmax/global screening identified 115; global-AF-only review therefore missed 102/115 (88.7%) ancestry-aware frequency alerts. Severe annotation did not protect against this signal: 92/262 AF-observed LOF/splice assertions (35.1%) had popmax/global AF >1e-5, similar to missense assertions (19/66, 28.8%). Mechanism triage showed that this was a tension universe rather than an error-rate estimate: 35/92 severe-discordant assertions were in carrier-compatible CASQ2/TRDN contexts, while only 1/57 non-recessive severe-discordant assertions was AC-supported. VITAL compressed 115 naive AF alerts to 3 red-priority review cases. In an independent 3,000-variant cross-disease historical audit, VITAL-red captured one real P/LP-to-VUS reclassification (CFAP91 VCV000812096), while an expert-panel RYR1 downgrade without frequency tension remained outside scope.

**Conclusion:** Severe annotation alone is not sufficient to override population-frequency contradiction. Global-AF-only ACMG/AMP implementations can systematically miss ancestry-aware frequency tension; popmax-aware, AC-aware, representation-aware review prioritization is needed to make these conflicts clinically inspectable without automating final classification.

## Introduction

Inherited arrhythmia syndromes, including long QT syndrome, Brugada syndrome, catecholaminergic polymorphic ventricular tachycardia, and related disorders, depend on accurate interpretation of rare variants in ion-channel, calcium-handling, and accessory genes. ClinVar is central to this interpretation because it aggregates public clinical assertions used by laboratories, expert panels, and researchers. A P/LP label in ClinVar can influence diagnosis, cascade testing, surveillance, and treatment decisions, so systematic errors in public assertions matter clinically even when individual case-level evidence is unavailable.

Population frequency is one of the strongest checks on high-penetrance Mendelian pathogenicity. A variant that is too common in a general population database should trigger benign evidence under ACMG/AMP BS1/BA1-style reasoning. The operational problem is that global AF can hide ancestry-specific enrichment. A variant may appear globally rare while exceeding a clinically relevant frequency threshold in one population. In that scenario, a global-AF-only workflow can falsely reassure the reviewer and allow frequency contradiction to remain hidden.

Absence from population databases creates a second failure mode. Non-observation is often treated as evidence of rarity, yet exact representation in gnomAD depends on sequence context, variant type, normalization, coverage, alignment, and calling. This is particularly important for indels, duplications, and splice/complex alleles. Frequency-based interpretation therefore fails both when ancestry-aware signals are missed and when absence is overinterpreted for variants that are difficult to represent.

Here, we ask whether ClinVar P/LP arrhythmia assertions remain coherent when evaluated with exact allele matching, popmax, AC, variant type, detectability, and review support. We introduce VITAL as an explainable re-review prioritization framework. VITAL does not classify variants as benign or pathogenic. It identifies where archived clinical assertions and population-frequency evidence are sufficiently discordant to justify expert review.

## Methods

### Variant selection and collapsing

ClinVar P/LP variants were retrieved for 20 canonical inherited arrhythmia genes: KCNQ1, KCNH2, SCN5A, KCNE1, KCNE2, RYR2, CASQ2, TRDN, CALM1, CALM2, CALM3, ANK2, SCN4B, KCNJ2, HCN4, CACNA1C, CACNB2, CACNA2D1, AKAP9, and SNTA1. VCV XML records were retrieved through the NCBI Entrez API. Records were collapsed to unique GRCh38 variant-level entries using chromosome, position, reference allele, and alternate allele. After quality control, including exclusion of one high-AF single-submitter SCN5A record without expert review, 1,731 unique P/LP variants remained. The primary data freeze was April 21, 2026.

### gnomAD matching and evidence states

Variants were queried against gnomAD v4.1.1 exomes through the gnomAD GraphQL API. Exact matching required chromosome, position, reference allele, and alternate allele concordance. Variants were assigned one of five frequency evidence states: `frequency_observed`, `not_observed_in_gnomAD`, `allele_discordance_no_exact_AF`, `exact_match_without_AF`, or `gnomad_query_error_no_frequency_evidence`. Missing AF was never converted to AF=0. Variants without usable exact frequency evidence were retained as gray no-frequency-evidence records.

### Ancestry-aware frequency analysis

For exact AF-observed variants, global AF, global AC, population AF, population AC, and popmax AF were extracted. A naive global screen used global AF >1e-5. An ancestry-aware screen used global or popmax AF >1e-5. AC-supported frequency evidence required a qualifying global or popmax AC >=20 at the relevant frequency signal. These thresholds were used for review prioritization, not automatic classification.

### VITAL score and red-priority gate

VITAL is a 0-100 frequency-assertion tension score. The score is the clipped sum of seven interpretable components:

`VITAL = min(100, AF_pressure + AC_reliability + popmax_enrichment + variant_type_tension + technical_detectability + gene_context + review_fragility)`.

| Component | Role |
|---|---|
| AF pressure | Scales with max(global AF, popmax AF) above 1e-5 |
| AC reliability | Increases with qualifying AC and saturates at AC >=20 |
| Popmax enrichment | Captures ancestry-specific enrichment over global AF |
| Variant-type tension | Adds context for SNV, indel, duplication, or complex representation |
| Technical detectability | Downweights naive absence assumptions by variant-type detectability |
| Gene context | Captures gene-level frequency-tension background |
| Review fragility | Prioritizes weak/single-submitter assertions over expert-reviewed records |

A red-priority call required VITAL >=70, weak review support, and AC-supported frequency evidence. VITAL is not a benignity classifier and does not model inheritance, penetrance, phenotype match, segregation, functional data, or local laboratory evidence. Final interpretation remains expert-driven.

![VITAL review-prioritization workflow](figures/vital_clinical_workflow.png)

### Supporting analyses

We compared global AF, popmax/global AF, AC-supported frequency screening, and VITAL red-priority calls. We assessed variant-type enrichment among non-overlap records, severe-annotation frequency discordance, LOF subtype patterns, and mechanism triage of severe-discordant variants. We performed supporting audits in a 3,000-variant cross-disease ClinVar P/LP sample, including historical January 2023 to April 2026 reclassification, expert-panel reviewed assertions, and real strict P/LP-to-VUS/B/LB events. Threshold sensitivity, AC-gate sensitivity, weight sensitivity, KCNH2 diagnostics, and time-series details are provided in supplementary outputs.

## Results

### Frequency evidence is sparse and structured

ClinVar parsing identified 1,731 unique arrhythmia P/LP variants. Of these, 350 (20.2%) had exact allele-level matches in gnomAD exomes, and 334 had usable AF evidence. Among the remaining variants, 645 were allele-discordant at the same position, 736 had no gnomAD record, and 16 exact matches lacked usable AF blocks. Thus, most ClinVar arrhythmia P/LP assertions cannot be evaluated by a simple exact-AF lookup.

Table 1 summarizes the evidence states used throughout the analysis.

| Evidence state | Count | Interpretation |
|---|---:|---|
| Exact gnomAD match | 350 | Allele represented in gnomAD exomes |
| Usable AF evidence | 334 | Eligible for frequency scoring |
| Allele discordance | 645 | Same position, different alternate allele |
| No gnomAD record | 736 | No record at queried position |
| Exact match without AF | 16 | Exact allele present but no usable AF block |

Among AF-observed variants, most were globally ultra-rare: 321/334 (96.1%) had global AF <=1e-5. This global view was misleading because many globally rare variants were population-enriched by popmax.

### Popmax exposes frequency contradictions missed by global AF

Global AF alone identified only 13 arrhythmia P/LP variants above AF >1e-5. In contrast, popmax/global screening identified 115 variants above the same threshold, including 102 that were globally rare but population-enriched. A global-AF-only ACMG-style screen would therefore miss 102/115 (88.7%) ancestry-aware frequency alerts.

![Population-specific frequency tension missed by global AF](figures/arrhythmia_population_af_outliers.png)

This is an ACMG/AMP implementation failure in a specific scenario: the criteria are not the problem, but global-AF-only implementation hides the population evidence that BS1/BA1-style reasoning is supposed to evaluate.

### Absence is shaped by variant type

Non-overlap was not random. Indels were enriched among variants without exact usable AF evidence: 47.2% of non-overlap variants were indels compared with 28.9% of exact-matched variants (OR=2.20; BH q=1.08e-9). This supports a second operational hazard: absence from gnomAD is not equivalent to rarity for structurally complex or representation-sensitive variants.

![Variant-type detectability bias](figures/arrhythmia_vital_absence_not_rarity.png)

Exome-vs-genome sensitivity did not rescue this problem for duplications. Among 22 exact-matched duplications, none had genome AF >1e-5. The practical conclusion is representation-aware, not mechanistic: absence of an indel or duplication from gnomAD should be interpreted more cautiously than absence of a well-represented SNV.

### Severe annotation does not protect against frequency discordance

The main biological result is that severe annotation does not make a ClinVar P/LP assertion frequency-clean. Among 334 AF-observed arrhythmia variants, 115 (34.4%) had popmax/global AF >1e-5. Among LOF/splice assertions, 92/262 (35.1%) were frequency-discordant, similar to missense assertions (19/66, 28.8%; Fisher OR=1.34, p=0.384). Frameshift, stop-gained, and canonical splice variants all showed naive frequency discordance.

This result does not mean that 35.1% of LOF/splice assertions are wrong. Mechanism triage showed that 35/92 severe-discordant assertions were in CASQ2 or TRDN, where recessive/biallelic disease architecture can make heterozygous carrier frequency biologically compatible. After excluding CASQ2/TRDN, 57/197 severe annotations (28.9%) remained discordant across 8 genes, but only 1/57 non-recessive severe-discordant assertions was AC-supported.

![Severe annotation frequency discordance](figures/vital_external_truth_biological_claim.png)

| Mechanism triage class | N (% of severe-discordant) | AC-supported | Interpretation |
|---|---:|---:|---|
| Carrier-compatible recessive/biallelic architecture | 35 (38.0%) | 6 | Not misclassification by itself |
| Non-recessive AC-supported high-penetrance tension | 1 (1.1%) | 1 | Highest-priority adjudication |
| Non-recessive constrained-gene low-AC surveillance | 24 (26.1%) | 0 | Monitor; insufficient AC for urgent action |
| Non-recessive other low-AC or unresolved surveillance | 32 (34.8%) | 0 | Routine mechanism review |

The claim is therefore narrower and stronger than "LOF variants are often benign": severe annotation is not sufficient to override population-frequency contradiction. It identifies a mechanism hypothesis that must be reconciled with ancestry, AC, inheritance, penetrance, transcript context, and review strength.

### VITAL compresses urgent re-review burden

A naive popmax/global AF >1e-5 screen flagged 115 arrhythmia P/LP variants. Adding AC support reduced the actionable frequency set to 9. VITAL red-priority reduced the urgent review queue to 3 variants by requiring score >=70, weak review, and AC-supported frequency evidence. This is a 97.4% compression relative to naive AF screening.

![VITAL workload compression](figures/arrhythmia_vital_score_model.png)

Table 2 compares the operational layers.

| Review layer | Flagged variants | Practical meaning |
|---|---:|---|
| Global AF >1e-5 | 13 | Misses most ancestry-aware signals |
| Popmax/global AF >1e-5 | 115 | Sensitive but high review burden |
| AC-supported frequency signal | 9 | Reduces one-allele noise |
| VITAL red-priority | 3 | Short urgent expert-review queue |

Workflow-concordance benchmarking supported the selected red gate but was not treated as diagnostic accuracy. The proxy positives were weak-review, AC-supported frequency signals; VITAL red captured the 3 operational positives at cutoffs 65-70 while avoiding proxy false-positive urgent review triggers. Sensitivity details are provided in supplementary tables.

### Red-priority cases show distinct biological tension types

The 3 red-priority arrhythmia variants were not interchangeable AF outliers. They represented three different interpretation problems: a high-frequency haplotype/susceptibility assertion, a recessive-carrier-compatible LOF assertion, and a borderline constrained-gene splice assertion.

| Variant | AF/AC signal | Review support | VITAL | Interpretation question |
|---|---|---|---:|---|
| SCN5A VCV000440850, c.[3919C>T;694G>A] | AFR popmax AF=5.68e-3; AC=190 | No assertion criteria | 96.1 | Haplotype/drug-response susceptibility vs broad Mendelian pathogenic label |
| TRDN VCV001325231, c.1050del | Global AC=40; AMR popmax AF=2.18e-4 | Single submitter | 74.3 | Recessive carrier frequency vs affected-state pathogenicity |
| KCNH2 VCV004535537, c.2398+2T>G | Global AC=24; ASJ popmax AF=1.32e-4 | Single submitter | 70.3 | Borderline splice assertion in constrained dominant LQTS gene |

![Biological contrast of red-priority cases](figures/vital_biological_contrast_cases.png)

Manual review of live ClinVar records confirmed that these cases remained weakly supported or single-submitter assertions with unresolved mechanism-specific interpretation questions. VITAL does not declare them benign. It identifies them as the records a laboratory should inspect first.

### Cross-disease and historical audits support scope, not broad prediction

In an independent 3,000-variant cross-disease current ClinVar P/LP sample, naive popmax/global AF screening flagged 332 variants (11.1%; 95% CI 10.0%-12.2%), while VITAL red-priority contained 3 (0.10%; 95% CI 0.034%-0.294%). Among 204 expert-panel reviewed P/LP assertions in that set, 39 had naive AF tension and 13 were AC-supported, but none were VITAL red-priority. Curated high-frequency exceptions therefore remained visible without being converted into urgent downgrade calls.

Historical audit from January 2023 to April 2026 confirmed high specificity and low recall. In the 3,000-variant set, 23 records underwent real strict P/LP-to-VUS/B/LB reclassification; VITAL-red captured 1/23, CFAP91 VCV000812096, a weak/no-assertion record with popmax AF=3.50e-4, AC=152, and VITAL=76.5. One expert-panel strict downgrade, RYR1 VCV001214001, had no frequency tension (max AF=8.99e-7; VITAL=0.0) and was correctly outside scope. This is the boundary of the framework: VITAL captures frequency-assertion hazards, not all reclassification mechanisms.

## Discussion

The central finding is that severe annotation does not override population evidence. Up to 35.1% of AF-observed LOF/splice ClinVar P/LP arrhythmia assertions showed ancestry-aware frequency discordance. That figure is not an error rate: many cases are compatible with recessive carrier architecture or lack AC support for immediate action. The clinical risk is subtler. A severe-looking ClinVar assertion can be mechanistically plausible and still too frequent for the high-penetrance label attached to it.

Global-AF-only workflows fail operationally because they hide ancestry-specific signals. In this arrhythmia audit, global AF detected 13 frequency alerts, while popmax/global analysis detected 115. If a laboratory applies BS1/BA1-style reasoning with global AF alone, most ancestry-aware contradictions remain invisible. Similarly, absence from gnomAD cannot be treated uniformly as rarity because indels and duplications are underrepresented by exact allele matching.

VITAL's contribution is workflow compression with explainability. It transforms 115 naive AF alerts into 3 urgent review cases while preserving gray no-frequency-evidence states and component-level explanations. This is not a replacement for ACMG/AMP interpretation, ClinGen specification, functional evidence, phenotype review, or laboratory judgment. It is a way to decide what to inspect first when frequency and clinical assertion are visibly in tension.

The three red-priority arrhythmia cases illustrate why a single AF threshold is insufficient. SCN5A VCV000440850 is most consistent with a context-dependent haplotype or susceptibility assertion rather than a straightforward high-penetrance monogenic label. TRDN VCV001325231 shows that a frameshift can be frequency-tolerable under recessive/biallelic architecture. KCNH2 VCV004535537 is the borderline case: biologically uncomfortable in a constrained dominant gene, but still unresolved and weight-sensitive.

Prior ClinVar-gnomAD studies have described rarity patterns and gene constraint. This analysis adds three operational distinctions: exact allele absence is separated from allele discordance and no-record states; popmax plus AC is placed at the center of frequency review; and review burden is treated as a clinical workflow endpoint rather than an incidental byproduct.

## Limitations

First, VITAL weights are expert-specified to preserve interpretability and were not trained as an optimized prediction model. Sensitivity analyses showed stable prioritization around the red gate, but disease-specific deployment should recalibrate thresholds and weights against local disease architecture.

Second, historical validation is sparse. VITAL-red captured 1/23 strict future P/LP-to-VUS/B/LB events in the 3,000-variant audit, confirming that it is not a sensitive predictor of all future reclassification. Its intended endpoint is high-specificity frequency-assertion prioritization.

Third, the score does not directly model penetrance, phase, zygosity, inheritance, phenotype specificity, segregation, MAVE/functional data, transcript rescue, NMD, or private laboratory evidence. These remain expert review tasks.

Fourth, workflow-concordance benchmarking uses proxy labels and is partly circular because weak review and AC-supported frequency evidence help define both the operational positives and the red gate. These metrics should be read as burden diagnostics, not diagnostic accuracy.

Fifth, public ClinVar captures classification-level changes, not patient-level downstream effects. We can identify P/LP-to-VUS/B/LB movement, but not diagnosis reversal, cascade-testing changes, surveillance changes, or therapy changes.

## Clinical and research implications

1. Use popmax, not global AF alone, when auditing ClinVar P/LP assertions for frequency tension.

2. Require AC support before escalating frequency contradictions to urgent review.

3. Do not equate absence from gnomAD with rarity, especially for indels, duplications, and complex alleles.

4. Do not treat severe annotation as self-validating; frameshift, stop-gained, and splice labels do not cancel population-frequency contradiction.

5. Separate re-review prioritization from final classification. VITAL tells laboratories what to inspect first, not what to believe automatically.

VITAL is implemented as a reproducible batch-scoring workflow with cached intermediate files and exportable review queues.

## Data availability

Code, cached intermediate files, figures, and supplementary outputs required to reproduce downstream analyses are available in the project repository. The April 21, 2026 data freeze used ClinVar records retrieved through NCBI Entrez and bulk variant_summary files, gnomAD v4.1.1 exome/genome data queried through the gnomAD GraphQL API, and UCSC annotation resources. Machine-readable supplementary tables include VITAL scores, frequency flags, severe-annotation discordance, mechanism triage, red-priority case summaries, historical audits, expert-panel comparator outputs, and sensitivity analyses.
