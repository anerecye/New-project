# Population-Frequency Tension and Clinical Reclassification Risk Prioritization in ClinVar Pathogenic Arrhythmia Variants

## Abstract

**Purpose:** Clinical variant databases are widely used to support rare disease interpretation, but the relationship between clinically cataloged variants and population allele-frequency spectra remains incompletely characterized. We examined this relationship for inherited arrhythmia-associated variants by cross-referencing ClinVar pathogenic/likely pathogenic (P/LP) variants across 20 canonical arrhythmia genes with gnomAD v4.1.1 exome and ancestry-specific data, and developed a clinical reclassification risk prioritization framework designed to reduce false-positive variant re-evaluation burden in clinical genomics workflows.

**Methods:** After collapsing ClinVar records to unique variant-level entries and applying quality control, 1,731 P/LP variants were evaluated against gnomAD exomes. Exact allele matches, allele discordance at the same genomic position, and variants with no gnomAD record were separated explicitly. We assessed global and ancestry-specific AF, allele count (AC), variant type, gene-level frequency constraint, ClinVar review strength, technical detectability, exome-vs-genome sensitivity, and KCNH2-specific non-overlap diagnostics. We then formalized these features as VITAL (Variant Interpretation Tension and Allele Load), a 0-100 clinical reclassification risk score incorporating global AF, popmax AF, AC reliability, popmax enrichment, variant type, gene-specific constraint, technical detectability, and review fragility. Variants without frequency evidence were not assigned AF=0; they were carried forward under explicit gray no-frequency-evidence labels. VITAL was evaluated as a decision-support workflow and pathogenicity tension continuum for prioritizing human review, not as an automated classifier. Weight-profile sensitivity and historical threshold calibration were performed without outcome refitting.

**Results:** Of 1,731 arrhythmia P/LP variants, 350 (20.2%) had exact allele-level gnomAD matches and 334 (19.3%) had usable exact AF data. Among AF-covered variants, 321/334 (96.1%) had global AF <=1e-5, 12/334 (3.6%) had global AF between 1e-5 and 1e-4, and 1/334 (0.3%) exceeded 1e-4. Popmax analysis identified 115 variants with global or population-specific AF >1e-5, compared with only 13 variants detected by global AF alone. Among 328 variants with complete popmax data, 102 (31.1%) were globally rare but population-enriched. Indels were significantly enriched in the non-overlap set (47.2% vs 28.9% of exact matches; OR=2.20; BH q=1.08e-9). Absence was variant-type biased: exact AF evidence was available for 238/970 SNVs (24.5%) but only 71/529 deletions (13.4%; p=2.36e-7 vs SNVs) and 21/194 duplications (10.8%; p=1.16e-5 vs SNVs), supporting the conclusion that absence from gnomAD is not equivalent to rarity for structurally complex variants. A naive popmax/global AF >1e-5 screen flagged 115 variants, whereas the AC-supported frequency screen flagged 9 and VITAL red prioritized 3 variants, a 97.4% reduction in actionable calls. LOF subtype analysis showed non-equivalent behavior across frameshift, stop_gained, and canonical_splice assertions; high-tension variants included 2 frameshift variants, 1 stop_gained variant, and 1 canonical_splice variant. Weight sensitivity generated no new red variants; SCN5A was anchor, TRDN near-stable, and KCNH2 borderline. In a proxy operational benchmark intended primarily to measure false-positive re-evaluation burden, VITAL red produced 0 false-positive calls among 109 high-confidence retained-P/LP proxy variants, while the naive popmax/global AF screen produced 53; perfect precision/recall values in this benchmark reflect design concordance, not independent clinical accuracy. Historical analyses were retained as secondary audit layers rather than primary validation: red sets remained very small (arrhythmia n=2; pooled multi-domain n=4), strict reclassification was not enriched, and broad/expanded instability endpoints provided only qualitative directional enrichment. External stress tests assessed false-positive burden outside the arrhythmia context; their small size means absence of red calls cannot be interpreted as sensitivity evidence.

**Conclusion:** Clinically cataloged arrhythmia P/LP variants are strongly concentrated at ultra-rare population frequencies, but population-frequency interpretation is shaped by ancestry, AC reliability, variant type, gene context, technical detectability, and review quality. VITAL converts frequency-based contradictions into a pathogenicity tension continuum and a short, auditable clinical reclassification risk queue rather than a black-box reclassification rule. As of the April 21, 2026 data freeze, the framework highlighted three arrhythmia P/LP assertions for expert review while suppressing more than 100 naive AF alerts, preserving explicit no-frequency-evidence states, and leaving final clinical decision-making with qualified reviewers.

## Introduction

Inherited arrhythmia syndromes, including long QT syndrome, Brugada syndrome, catecholaminergic polymorphic ventricular tachycardia, and related conditions, are caused by rare pathogenic variants in a defined set of ion channel, calcium-handling, and accessory genes. Clinical interpretation of rare variants in these genes relies heavily on structured variant databases, of which ClinVar is the largest publicly accessible resource. ClinVar aggregates variant classifications from clinical laboratories, research submitters, and expert panels, providing pathogenic and likely pathogenic designations that inform diagnostic and therapeutic decisions.

Population allele frequency is a cornerstone criterion in variant interpretation frameworks such as ACMG/AMP. A variant observed at excessive frequency in a general population database such as gnomAD is unlikely to cause a highly penetrant Mendelian disorder and may support benign evidence under BS1 or BA1. Conversely, ultra-rare or absent variants can be consistent with pathogenicity. However, absence and rarity are not interchangeable. Population databases are structured by sequencing technology, variant representation, ancestry composition, allele count, and calling pipelines. These factors are especially important for indels, duplications, and other structurally complex alleles.

Several practical challenges complicate the use of population frequency in ClinVar auditing. ClinVar records may not have corresponding exact allele-level entries in gnomAD. A variant may be absent from the population cohort, represented differently, affected by normalization, or located at a position where gnomAD reports a different alternate allele. Even when AF is available, global AF can hide population-specific enrichment, and low allele counts can produce unstable frequency estimates. In addition, ClinVar review strength varies substantially: an assertion from an expert panel is not equivalent to a single-submitter record with minimal criteria.

Here, we report a systematic cross-reference of ClinVar P/LP variants across 20 canonical inherited arrhythmia genes with gnomAD v4.1.1 exome data. We characterize exact matches, allele-discordant sites, no-record variants, ancestry-specific frequency enrichment, variant-type detectability, and gene-level heterogeneity. We then extend the analysis from description to clinical decision support by introducing VITAL (Variant Interpretation Tension and Allele Load), a clinical reclassification risk prioritization framework that integrates AF, popmax, AC, variant type, gene context, technical detectability, and ClinVar review fragility. The goal is not to automate reclassification, but to reduce false-positive variant re-evaluation burden and help laboratories focus limited expert review time on the records most likely to affect clinical decision-making. Finally, we test VITAL against baseline ACMG-style frequency screens, historical ClinVar reclassification from 2023 to 2026, external disease panels, and a negative-control gene set.

## Methods

### Variant selection and ClinVar collapsing

ClinVar P/LP variants were retrieved for 20 canonical inherited arrhythmia-associated genes: KCNQ1, KCNH2, SCN5A, KCNE1, KCNE2, RYR2, CASQ2, TRDN, CALM1, CALM2, CALM3, ANK2, SCN4B, KCNJ2, HCN4, CACNA1C, CACNB2, CACNA2D1, AKAP9, and SNTA1. VCV XML records were retrieved programmatically through the NCBI Entrez API. Clinical significance was determined from the ClinVar germline classification field; Pathogenic and Likely pathogenic records were retained.

Records were collapsed to unique variant-level entries using a composite GRCh38 key of chromosome, genomic position, reference allele, and alternate allele. VCV records lacking GRCh38 VCF coordinates were excluded and logged. After collapsing and quality control, including exclusion of SCN5A VCV000440849 because of a high population AF from a single-submitter assertion without expert review, 1,731 unique arrhythmia P/LP variants were retained.

The primary data freeze was April 21, 2026. Throughout the manuscript, "current" refers to this April 2026 ClinVar and gnomAD analysis snapshot. ClinVar assertions can change after data freeze; therefore, all red-priority variants should be rechecked against the live ClinVar record before publication, clinical reporting, or any patient-facing decision.

### gnomAD matching and frequency evidence states

Each variant was queried against gnomAD v4.1.1 exomes through the gnomAD GraphQL API. Exact matching required concordance of chromosome, position, reference allele, and alternate allele. Exact matches with usable exome AF values were treated as frequency-observed records. Variants without usable exact AF were not imputed as AF=0. Instead, all variants were carried forward under explicit frequency evidence states:

- `frequency_observed`: exact allele match with usable gnomAD AF.
- `not_observed_in_gnomAD`: no gnomAD record at the queried position.
- `allele_discordance_no_exact_AF`: gnomAD reports a different alternate allele at the same position.
- `exact_match_without_AF`: exact allele match but no usable exome AF block.
- `gnomad_query_error_no_frequency_evidence`: transient query failure, kept outside frequency scoring.

This design prevents silent conversion of missing frequency evidence into false rarity. Variants without observed AF retain missing AF fields and are assigned a gray no-frequency-evidence VITAL band.

### Ancestry-specific frequency analysis

For exact AF-covered variants, global AF, global AC, population-specific AF, population-specific AC, and popmax AF were extracted where available. Populations included African/African American, Amish, Ashkenazi Jewish, East Asian, Finnish, Latino/Admixed American, Middle Eastern, Non-Finnish European, South Asian, and remaining ancestry groups. A variant was considered population-enriched when global AF <=1e-5 but popmax AF >1e-5.

### Reclassification-risk baseline screens

We evaluated several simple frequency-based baselines: global AF >1e-5, global AF >1e-4, popmax or global AF >1e-5, popmax or global AF >1e-4, and AF >1e-5 with qualifying AC >=20 in either the global or popmax frequency source. These baselines approximate increasingly strict ACMG-style frequency screens but do not incorporate review quality, variant type, gene context, or technical detectability.

### VITAL clinical reclassification risk prioritization framework

VITAL is a 0-100 clinical reclassification risk score. It is not intended to declare variants benign, and it is not a replacement for clinical interpretation. It prioritizes ClinVar P/LP assertions that are under enough population-frequency, allele-count, technical, gene-context, and review-quality pressure to warrant expert re-review.

The VITAL score combines seven interpretable components:

- AF pressure from global AF and popmax AF, scaled across the range from 1e-5 to 1e-3.
- AC reliability, saturated at AC >=20.
- Popmax enrichment over global AF.
- Variant-type tension, reflecting whether the observed population signal is surprising for SNV, indel, duplication, or complex variant categories.
- Gene-specific frequency constraint, estimated from the gene's observed frequency-positive burden.
- Technical detectability, derived from empirical exact-match/non-overlap behavior by variant type and a prior for short-read detectability.
- ClinVar review fragility, based on review status and submitter count.

The score is the clipped sum of component scores:

`VITAL = min(100, max(0, AF_pressure + AC_reliability + popmax_enrichment + variant_type_tension + technical_detectability + gene_constraint + review_fragility))`.

The component maxima were expert-specified rather than empirically fitted. They were chosen to make high-confidence frequency contradiction the dominant signal, AC reliability a required support term, and review fragility/technical detectability visible but not sufficient alone. VITAL should therefore be interpreted as an expert-weighted clinical workflow framework, not as a trained risk model. These weights should not be assumed to transfer unchanged to other disease areas; disease-specific implementations require local calibration or, at minimum, explicit weight-sensitivity reporting. For exact AF-observed variants, let `F = max(global_AF, popmax_AF)`, `ACq = max(global_AC if global_AF > 1e-5 else 0, popmax_AC if popmax_AF > 1e-5 else 0)`, and `R = (popmax_AF + 1e-12) / (global_AF + 1e-12)`. Components were calculated as follows:

| Component | Maximum | Calculation |
|---|---:|---|
| AF pressure | 45 | `45 * clip(log10((F + 1e-12) / 1e-5) / 2, 0, 1)` |
| AC reliability | 20 | `20 * clip(log1p(ACq) / log1p(20), 0, 1)`; set to 0 when `F <= 1e-5` |
| Popmax enrichment | 10 | `10 * clip(log10(R) / 2, 0, 1)`; set to 0 when `F <= 1e-5` |
| Variant-type tension | 6 | SNV=0, MNV/substitution=1, deletion=3, insertion=3, duplication=6, other/unresolved=2 |
| Technical detectability | 8 | `8 * technical_detectability_index`, derived from empirical exact-match behavior and variant-type priors; set to 0 when `F <= 1e-5` |
| Gene constraint | 10 | `10 * gene_frequency_constraint_proxy`; set to 0 when `F <= 1e-5` |
| Review fragility | 10 | practice guideline/expert panel=0, multiple submitters no conflicts=3, other review=5, single submitter=8, weak/no assertion or missing review=10; set to 0 when `F <= 1e-5` |

The red actionability label was not defined by the score alone. A red call required `VITAL >=70`, weak review support, and AC-supported frequency evidence at the operational AC gate. Variants without exact usable frequency evidence were not scored as zero-frequency variants; they were assigned the gray no-frequency-evidence band.

The final score is reported with component-level breakdowns for explainability. Bands were defined as green frequency-consistent, blue low-support signal, yellow watchlist, orange high-tension, red reclassification-priority, and gray no-frequency-evidence. A VITAL red call requires a high VITAL score, AC-supported frequency evidence, and weak review support. Thus, red is a review-priority label, not a clinical benign classification.

Conceptually, VITAL defines a variant pathogenicity tension continuum rather than a binary filter. Each variant is positioned in a multidimensional discordance space spanning population frequency, AC reliability, ancestry enrichment, review fragility, gene context, variant representation, and predicted functional severity. In the present implementation, functional context is represented lightly using HGVS-derived functional class and the internal gene-frequency constraint proxy. Future implementations could add orthogonal functional axes such as MAVE scores where available, gene constraint metrics such as LOEUF or missense Z, and protein-domain or transcript-specific annotations. We did not use these external functional resources in the current analysis because coverage across the arrhythmia variant set is incomplete; they are proposed as natural extensions of the same continuum framework.

For descriptive subtype analysis, loss-of-function (LOF) variants were further annotated from ClinVar title/HGVS strings as frameshift, stop_gained, canonical_splice, or other_LOF. Frameshift calls required protein-level `fs` or "frameshift" text; stop_gained calls required `Ter`, `*`, "stop gained", or "nonsense"; canonical_splice calls required a +/-1 or +/-2 donor/acceptor pattern or explicit splice donor/acceptor text. This subtype annotation was not used in the VITAL score or red actionability gate. It was added only to test whether LOF subclasses contributed equally to the frequency-function discordance signal.

In the intended clinical workflow, VITAL sits upstream of expert interpretation. It converts population-frequency alerts into a shorter review queue, while final classification remains with laboratory directors, variant scientists, disease experts, and ACMG/AMP evidence review. The framework is therefore a burden-reduction and decision-support tool, not an autonomous diagnostic system.

![VITAL clinical workflow](figures/vital_clinical_workflow.png)

The AC threshold is applied as an actionability filter, not as part of the scoring function. This design separates signal detection (the VITAL score) from operational decision thresholds (the AC gate). To evaluate threshold sensitivity, we repeated current and historical analyses with AC gates of >=5, >=10, >=20, and >=50, recording both the number of flags and the composition of the VITAL-red set at each threshold. Variants entering or leaving the red set were explicitly tracked relative to the prespecified AC>=20 gate.

To evaluate sensitivity to expert-specified component weights, we recalculated VITAL scores under five alternative weighting profiles without refitting to outcomes: balanced equal components, frequency-dominant, review-fragility-dominant, technical-detectability-dominant, and reduced-AF-pressure. For each profile, component scores were normalized to their original maxima and rescaled to the alternative maxima, then red calls were recomputed using the same `score >=70`, AC-supported frequency, and weak-review gate. We recorded red-set overlap with the primary model, gained/lost red variants, compression versus naive AF screening, proxy false-positive burden, and Spearman rank correlation with the primary score.

For historical threshold calibration, we evaluated score cutoffs from 40 to 95 in two modes: score-only and red-gate-compatible. The red-gate-compatible mode applied the same operational actionability gate as the validator: score threshold, AC-supported frequency evidence, and weak review. Calibration QA retained both the field-mapping-inconsistent run and the corrected run. The initial calibration run produced zero red-gate-compatible flags because `weak_review_signal` was absent from the historical prediction table and defaulted to false. The corrected calibration derived weak review from `review_score <=1` or `submitter_count <=1`; it also falls back to `max_frequency_signal >1e-5` plus `qualifying_frequency_ac >=20` when `frequency_signal_ac_ge_20` is absent. This independent field-mapping check was treated as a reproducibility audit, not as a model change.

### Operational proxy benchmark and threshold sweep

Because prospective truth labels are unavailable for current ClinVar assertions, we defined an operational proxy benchmark to test whether VITAL reduces false-positive re-evaluation burden without simply flagging every frequency-positive variant. Positives were weak-review, AC-supported frequency signals. High-confidence current P/LP proxy negatives were variants with stronger review status and no comparable frequency contradiction. This benchmark is not an independent validation of clinical pathogenicity or future reclassification, and it contains an unavoidable degree of circularity because weak review and AC-supported frequency evidence are also part of VITAL's red actionability gate. Its purpose is narrower: to compare how much review burden different frequency screens would generate in the same operational universe. We therefore interpret precision/recall and ROC/PR curves as workflow diagnostics, not as definitive estimates of clinical accuracy. In particular, precision=1.00 and recall=1.00 in this benchmark should be understood as design-concordance metrics under a proxy label definition, not as evidence that VITAL has perfect clinical performance. We compared VITAL red against baseline frequency screens using TP, FP, FN, TN, precision, recall, specificity, false-positive rate, and false-negative rate. A threshold sweep evaluated score cutoffs from 40 to 95.

### Historical reclassification enrichment analysis

To explore whether VITAL enriches for future ClinVar instability rather than only describing current data, we downloaded the January 2023 ClinVar variant_summary archive and scored baseline P/LP variants using the same VITAL framework. The primary historical analysis used the arrhythmia panel. We then expanded the historical dataset to cardiomyopathy, epilepsy, hearing-loss genes, and a random ClinVar P/LP sample, matching the external-domain stress-test design. Baseline 2023 records were compared with the April 21, 2026 ClinVar snapshot. This analysis was treated as preliminary because the red set was small by design. Given the intentionally small red set, historical validation is underpowered and serves only as a directional consistency check rather than a definitive predictive evaluation. Three endpoints were evaluated:

- strict endpoint: P/LP to B/LB or VUS.
- broad endpoint: P/LP to non-P/LP, conflicting, other, or absent from follow-up.
- expanded endpoint: broad endpoint plus aggregate clinical-significance text change or review-status change.

We reported ROC AUC, average precision, method-comparison tables, threshold sweeps, precision/recall, and enrichment of future events among VITAL-red variants. Given the small number of red variants, the historical analysis was interpreted qualitatively as a sparse enrichment test rather than as a stable effect-size estimate or definitive proof of predictive performance. To make this instability explicit, we performed a one-event perturbation analysis showing how the apparent enrichment changes if one red event is added, removed, or replaced by a non-event.

We also performed an exploratory historical calibration analysis to test whether the prespecified VITAL red cutoff of 70 coincided with the threshold that maximized enrichment in the 2023-to-2026 data. This analysis was explicitly not used to refit the score threshold. Because the red set is intentionally small, maxima across thresholds were interpreted as descriptive stress tests rather than evidence for an optimized predictive model.

For the expanded multi-domain analysis, metrics were calculated in two ways: within each domain and after pooled deduplication. This separation was prespecified because disease architecture, baseline reclassification rate, and likely error mechanisms differ across arrhythmia, cardiomyopathy, epilepsy, hearing loss, and random ClinVar records. Pooled deduplication was performed only after within-domain sanity checks. Within each domain, duplicate variation IDs or duplicate variant keys would be collapsed by retaining the row with the highest VITAL score. For pooled analysis, cross-domain duplicates were resolved using a fixed domain priority (arrhythmia > cardiomyopathy > epilepsy > hearing loss > random ClinVar P/LP), with maximum VITAL score used only as a tie-breaker within the same priority level. Dedupe decision tables recorded variation_id, domains_present, kept_domain, dropped_domains, and reason_kept. Stratified bootstrap resampled rows within each domain after deduplication; because red variants remained sparse, bootstrap summaries were interpreted as binary-noise sensitivity checks rather than smooth effect-size estimators.

### External-domain stress tests and descriptive comparator

To evaluate whether VITAL overfits to arrhythmia genes, we ran small external stress tests on three disease domains and one random ClinVar P/LP sample. Panels were intentionally limited to approximately 300 variants per domain for speed and interpretability:

- cardiomyopathy genes: MYH7, MYBPC3, TNNT2, TNNI3, TPM1, ACTC1, MYL2, MYL3, LMNA, FLNC, PLN, DSP, PKP2, DSG2, DSC2.
- epilepsy genes: SCN1A, SCN2A, SCN8A, KCNQ2, KCNQ3, STXBP1, CDKL5, PCDH19, SLC2A1, GABRA1, GABRG2, DEPDC5.
- hearing-loss genes: GJB2, GJB6, SLC26A4, MYO7A, OTOF, TECTA, TMC1, CDH23, USH2A, MYO15A.
- random ClinVar P/LP sample: 500 variants sampled outside the fixed arrhythmia panel.

The existing BRCA1, BRCA2, MLH1, and APC cache was retained as a non-arrhythmia descriptive comparator, not as a definitive negative control. These cancer-predisposition genes have incomplete penetrance, founder effects, and disease-mechanism differences that make them methodologically imperfect as a clean external validation set for arrhythmia interpretation. We therefore use this panel only to ask whether VITAL generates an obvious high-score tail in a familiar non-arrhythmia Mendelian context. These external analyses were designed as false-positive-burden and portability stress tests, not sensitivity analyses: without independent truth labels in these domains, absence of VITAL-red calls cannot establish that no true reclassification candidates were missed.

### Exome-vs-genome sensitivity analysis

For all 350 exact allele-level arrhythmia matches, we queried gnomAD v4.1 genomes using the same GraphQL API. We compared exome and genome AF/AC values and flagged genome-only recoverability. Special attention was given to duplications, where genome data might theoretically improve recovery over exome data.

### Variant type, KCNH2 diagnostics, repeats, GC content, and submission dates

Variants were classified as SNV, deletion, insertion, duplication, MNV/substitution, or other complex events using VCF allele lengths and HGVS descriptions. Duplications were subclassified by size: short (1-10 bp), medium (11-50 bp), and long (>50 bp).

KCNH2 was selected for diagnostic dissection because it has a large number of ClinVar P/LP assertions, high non-overlap, and a prominent duplication burden. We evaluated whether KCNH2 non-overlap was explained by variant type, duplication size, repeat overlap, local GC content, comparison with SCN5A, or submission timing. Local GC content was calculated in a 100-bp window around each variant using UCSC hg38 sequence. RepeatMasker and simpleRepeat annotations were queried around variant breakpoints. ClinVar DateCreated, DateLastUpdated, and MostRecentSubmission fields were parsed from VCV XML.

### Statistical analysis

Frequency classes were defined as ultra-rare (AF <=1e-5), rare (1e-5 < AF <=1e-4), and higher-frequency (AF >1e-4). Fisher's exact tests were used for categorical enrichment analyses, including variant-type non-overlap, gene-level frequency outlier burden, and KCNH2 diagnostic comparisons. Mann-Whitney U tests were used for GC content comparisons. P-values were corrected using the Benjamini-Hochberg false discovery rate procedure where appropriate. Analyses were implemented in Python 3 using requests, pandas, scipy, statsmodels, biopython, seaborn, matplotlib, and tqdm.

Confidence intervals were reported for core binary performance metrics. Wilson 95% confidence intervals were used for proportions, including precision, recall, specificity, false-positive rate, and event rates. For sparse historical enrichment ratios, Haldane-corrected log risk-ratio intervals were used to keep intervals finite in zero-cell settings. These intervals are intentionally shown even when wide, because they make the uncertainty from small red-set counts explicit. One-event perturbation and stratified bootstrap were used as descriptive sensitivity analyses only; neither was treated as formal prospective validation.

## Results

### Variant retrieval and gnomAD match structure

ClinVar parsing identified 1,731 unique P/LP variants across the 20 arrhythmia genes after collapsing and quality control. Of these, 350 (20.2%) had exact allele-level matches in gnomAD v4.1.1 exomes. Among the exact matches, 334 (19.3% of all variants) had usable AF evidence and were eligible for continuous frequency scoring. Another 16 exact matches lacked usable AF blocks and were retained as exact-match-without-AF records.

The remaining 1,381 variants lacked exact usable AF evidence. Among these, 645 (37.3% of total) were allele-discordant, meaning gnomAD reported a different alternate allele at the same position, and 736 (42.5%) had no gnomAD record at the queried position. No current arrhythmia variants were lost to transient query-error status. This split is clinically important: allele discordance is evidence about representation at a locus, not evidence that the ClinVar allele itself is present or absent.

The BRCA/MMR/APC control cache contained 2,000 sampled P/LP variants. It had 187 exact allele-level matches (9.3%), of which 173 had usable exact AF evidence. The VITAL audit table carried all 2,000 control variants forward, including 894 allele-discordant variants and 915 no-record variants.

### Global allele-frequency distribution

Among 334 AF-covered arrhythmia variants, 321 (96.1%) had global AF <=1e-5, 12 (3.6%) had global AF between 1e-5 and 1e-4, and 1 (0.3%) exceeded 1e-4. The median nonzero global AF was 1.37e-6, and 108 AF-covered variants had global AF equal to zero despite exact gnomAD representation. If all exact allele-level matches are used as the denominator, the ultra-rare count is 321/350 (91.7%), consistent with strong concentration of ClinVar P/LP assertions at extremely low population frequency.

The control gene set showed a similar broad pattern among frequency-observed variants: 168/173 control variants (97.1%) had global AF <=1e-5, 5/173 (2.9%) had global AF between 1e-5 and 1e-4, and none exceeded 1e-4. This supports the general expectation that ClinVar P/LP variants in high-penetrance Mendelian disease genes are heavily skewed toward ultra-rare frequencies, while preserving the caveat that cancer-predisposition genes differ in penetrance and founder architecture.

### Popmax reveals frequency contradictions missed by global AF

Global AF alone identified only 13 arrhythmia P/LP variants with AF >1e-5 and 1 with AF >1e-4. In contrast, popmax/global screening identified 115 variants with AF >1e-5 and 9 with AF >1e-4. Among 328 variants with complete popmax data, 102 (31.1%) were globally rare (global AF <=1e-5) but population-enriched (popmax AF >1e-5). This shows that global AF compresses ancestry-specific signals and that popmax is essential for detecting population-specific frequency tension.

### Variant type drives absence and detectability bias

Indels were significantly enriched in the non-overlap set. Insertions, deletions, and duplications together accounted for 652/1,381 non-overlap variants (47.2%) but only 101/350 exact allele-level matches (28.9%; OR=2.20; BH q=1.08e-9). This indicates that exact absence from gnomAD is partly a function of variant representation and detectability rather than biological rarity alone.

The VITAL absence/detectability analysis sharpened this result. Exact AF evidence was available for 238/970 SNVs (24.5%), but only 71/529 deletions (13.4%), 21/194 duplications (10.8%), and 4/30 insertions (13.3%). Compared with SNVs, deletions and duplications were significantly more likely to lack exact AF evidence (deletions: OR=2.10, p=2.36e-7; duplications: OR=2.68, p=1.16e-5). Thus, frequency-based criteria implicitly assume detectability, and that assumption does not hold equally across variant classes.

### Baseline frequency screens generate excess re-evaluation burden

A naive ACMG-style popmax/global AF >1e-5 screen flagged 115 arrhythmia P/LP variants (6.6% of the full analytical set). Requiring a high-frequency threshold of AF >1e-4 reduced this to 9 variants. Requiring AC >=20 in the relevant frequency source also produced 9 AC-supported frequency flags. These AC-supported signals are the frequency contradictions most suitable for risk-prioritized review because they are not driven by one or two observations.

However, AC alone was not enough. Among the 9 AC-supported frequency signals, 6 were stronger-review multiple-submitter/no-conflict assertions. VITAL therefore combined frequency pressure with review fragility, variant type, gene context, and technical detectability. VITAL red prioritized 3 variants, representing a 97.4% reduction in action-priority calls relative to the 115 naive AF flags.

### Operational proxy benchmark against baseline frequency screens

The operational proxy benchmark tested whether VITAL reduces re-evaluation burden relative to baseline frequency screens. It should not be read as independent clinical validation: the positive class was defined using weak review and AC-supported frequency evidence, which also contribute to VITAL's red actionability label. This circularity inflates the apparent precision/recall of VITAL red relative to what would be expected in a fully independent truth set. The resulting precision=1.00 and recall=1.00 are therefore artifacts of design concordance under this proxy task, not estimates of clinical diagnostic accuracy. The benchmark remains informative because it quantifies the practical workflow question under a shared proxy definition: how many high-confidence retained-P/LP proxy records would each method send into unnecessary urgent review?

Under this proxy benchmark, VITAL red had 3 true positives, 0 false positives, 0 false negatives, and 109 true negatives, with precision 1.00 (95% CI 0.44-1.00), recall 1.00 (95% CI 0.44-1.00), specificity 1.00 (95% CI 0.966-1.00), and false-positive rate 0.00 (95% CI 0.00-0.034). The intervals are wide for precision and recall because only three operational positives exist. The popmax/global AF >1e-5 baseline had the same recall but produced 53 false positives, with precision 0.054 (95% CI 0.018-0.146) and false-positive rate 0.486 (95% CI 0.394-0.579). The AF >1e-5 plus AC >=20 baseline still produced 6 false positives.

The threshold sweep showed that this result was not dependent on a single brittle cutoff. At VITAL score >=60, recall remained 1.00 with 2 false positives. At cutoffs of 65 and 70, precision and recall were both 1.00 with 0 false positives. At cutoffs >=75, precision remained 1.00 but recall fell to 0.33.

### Clinical workflow impact and false-positive burden reduction

VITAL's practical value is workload compression for clinical genomics teams. A naive popmax/global AF screen would send 115 current arrhythmia P/LP assertions into manual re-evaluation. VITAL red sends 3. This suppresses 112 potential false-positive review triggers while preserving an auditable explanation for why each non-red variant was withheld from urgent review, such as AC below 20, stronger review status, lower score, or no usable frequency evidence.

In the operational benchmark, the naive popmax/global AF >1e-5 screen produced 53 false positives among 109 high-confidence retained-P/LP proxy negatives, whereas VITAL red produced 0. Interpreted as a laboratory workflow, this is not just a statistical improvement: it prevents nearly half of benchmark-negative records from entering an unnecessary urgent re-evaluation queue. If a first-pass manual review requires approximately 30-60 minutes per variant, suppressing 112 naive alerts corresponds to roughly 56-112 reviewer-hours avoided in this 1,731-variant audit. This time estimate is illustrative, but it makes the clinical impact explicit: VITAL reduces review burden while preserving human decision-making.

### AC threshold sensitivity and red-set stability

The AC sensitivity analysis showed that the current red set is not an artifact of one arbitrary AC threshold, but it is appropriately sensitive at the low and high extremes. At AC>=5, 17 variants had AC-supported frequency signals and 4 were VITAL red. At AC>=10, 15 variants had AC-supported signals and the same 4 were VITAL red. At the prespecified AC>=20 threshold, 9 variants had AC-supported signals and 3 were VITAL red. At AC>=50, only 3 variants had AC-supported signals and 1 remained VITAL red.

The composition analysis is the key result. The AC>=20 red-core consisted of SCN5A VCV000440850, TRDN VCV001325231, and KCNH2 VCV004535537. All three remained red at AC>=5 and AC>=10. CACNB2 VCV003774534 appeared only at lower AC gates (AC>=5 and AC>=10) and was excluded at AC>=20 because its qualifying AC was below 20. At AC>=50, only SCN5A VCV000440850 remained red, while TRDN and KCNH2 were lost because their AC support did not reach the more stringent operational gate. Thus, the clinically relevant core is stable across AC>=5 to AC>=20, while the framework transparently shows which variants depend on lower-count evidence.

### Weight-profile sensitivity and red-set stability classes

Expert-weight sensitivity did not introduce any new red variants. Across five alternative weighting profiles, the red queue ranged from 1 to 3 variants, proxy false-positive count remained 0, and compression versus the naive AF screen remained 97.4%-99.1%. Rank correlation with the primary score was high across profiles (Spearman rho 0.990-1.000), showing that the score ordering is broadly stable even though marginal actionability calls can change.

The primary red variants were therefore assigned stability classes. SCN5A VCV000440850 was an anchor variant, retained as red under all five alternative profiles. TRDN VCV001325231 was near-stable, retained under 4/5 profiles and lost only under the balanced-equal profile. KCNH2 VCV004535537 was borderline, retained under 2/5 profiles and lost under balanced-equal, frequency-dominant, and reduced-AF-pressure profiles. This distinction is clinically important: the three variants are not equivalent in weight robustness, and KCNH2 should be presented as a threshold-adjacent priority signal rather than as an anchor call.

### Variant pathogenicity tension continuum

VITAL reorganizes population-frequency signals along a pathogenicity tension continuum rather than simply applying another binary threshold. Across exact AF-observed variants, increasing 20-point VITAL bands concentrated frequency-function discordance signals. The low-tension band (0-20; n=219) had no popmax AF >1e-4 signals and no AC-supported frequency contradictions. The 40-60 band had a median maximum frequency signal of 3.01e-5, 9.8% AC-supported frequency signals, and 39.3% indel/duplication representation. The 60-80 band (n=5) had 100% popmax AF >1e-4, 40.0% AC-supported frequency signals, 40.0% indel/duplication representation, and 100% canonical-or-atypical functional annotations (LOF, splice/intronic, or unresolved/composite). The 80-100 band contained the SCN5A composite assertion, with popmax AF=5.68e-3, qualifying AC=214, and weak review. The trend is not monotonic for every individual component because the upper bands are intentionally small, but the continuum shows that high VITAL scores represent convergence of frequency pressure, functional severity or atypical annotation, and review/AC context.

| VITAL band | N | Popmax AF >1e-4 | AC-supported AF | Weak/single review | Indel/duplication | LOF/splice/unresolved | Median max AF | Median AC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0-20 | 219 | 0.0% | 0.0% | 74.4% | 27.4% | 78.5% | 8.99e-7 | 0 |
| 20-40 | 48 | 0.0% | 0.0% | 39.6% | 20.8% | 79.2% | 1.85e-5 | 1 |
| 40-60 | 61 | 4.9% | 9.8% | 63.9% | 39.3% | 85.2% | 3.01e-5 | 1 |
| 60-80 | 5 | 100.0% | 40.0% | 60.0% | 40.0% | 100.0% | 1.50e-4 | 19 |
| 80-100 | 1 | 100.0% | 100.0% | 100.0% | 0.0% | 100.0% | 5.68e-3 | 214 |

This framing defines a clinically useful category of frequency-function discordant variants: ClinVar P/LP assertions whose population frequency is difficult to reconcile with a canonical high-penetrance interpretation. These variants may be clinically overcalled, incompletely penetrant, hypomorphic, transcript/domain-specific, ancestry-enriched, or true disease alleles with lower-than-assumed effect size. VITAL does not decide among these explanations. Its contribution is to identify where the tension is strongest and make the basis for that tension auditable.

The same pattern is clearer when comparing signal layers rather than thresholds. A naive AF >1e-5 screen captures 115 variants with median VITAL 41.1 and only 7.8% AC-supported frequency signals. Adding AC support reduces this to 9 variants with median VITAL 52.6. The VITAL 60-80 band contains only 5 variants but has 100% popmax AF >1e-4 and 100% LOF/splice/unresolved annotations; the VITAL-red layer contains 3 variants, all with AC-supported frequency, weak/single review, and canonical-or-atypical functional annotation. Thus, VITAL does not merely cut down a list. It reorganizes the frequency signal into progressively smaller strata in which frequency, function, and review evidence are increasingly discordant.

| Signal layer | N | Median VITAL | Popmax AF >1e-4 | AC-supported AF | Weak/single review | LOF/splice/unresolved |
|---|---:|---:|---:|---:|---:|---:|
| All AF-observed | 334 | 3.0 | 2.7% | 2.7% | 67.4% | 80.2% |
| Naive AF >1e-5 | 115 | 41.1 | 7.8% | 7.8% | 53.9% | 83.5% |
| AC-supported AF >1e-5 | 9 | 52.6 | 33.3% | 100.0% | 33.3% | 88.9% |
| VITAL 60-80 | 5 | 70.3 | 100.0% | 40.0% | 60.0% | 100.0% |
| VITAL 80-100 | 1 | 96.1 | 100.0% | 100.0% | 100.0% | 100.0% |
| VITAL red | 3 | 74.3 | 100.0% | 100.0% | 100.0% | 100.0% |

### LOF subtype discordance

Subtyping LOF variants showed that frameshift, stop_gained, and canonical_splice assertions do not contribute equally to VITAL signal. Among AF-observed LOF variants, 82 were frameshift, 108 were stop_gained, and 60 were canonical_splice. Naive AF >1e-5 flags were common in all three subclasses (29/82 frameshift, 33/108 stop_gained, and 23/60 canonical_splice), but AC-supported frequency contradictions were much rarer (1/82, 3/108, and 2/60, respectively). High-tension VITAL scores (60-100) were observed in 2 frameshift variants, 1 stop_gained variant, and 1 canonical_splice variant.

| LOF subtype | AF-observed LOF N | Naive AF flags | AC-supported AF | VITAL 60-100 | VITAL red | Max VITAL | Example high-tension variants |
|---|---:|---:|---:|---:|---:|---:|---|
| Frameshift | 82 | 29 (35.4%) | 1 (1.2%) | 2 (2.4%) | 1 | 74.3 | TRDN VCV001325231; TRDN VCV001074440 |
| Stop_gained | 108 | 33 (30.6%) | 3 (2.8%) | 1 (0.9%) | 0 | 61.0 | CASQ2 VCV002565185 |
| Canonical_splice | 60 | 23 (38.3%) | 2 (3.3%) | 1 (1.7%) | 1 | 70.3 | KCNH2 VCV004535537 |

This subtype split directly addresses the concern that all LOF variants might be driving VITAL in the same way. They are not. The 60-80 high-tension band contained two frameshift variants, one stop_gained variant, one canonical_splice variant, and one non-LOF/composite assertion. The TRDN frameshift variant VCV001325231 was VITAL-red despite a classically severe predicted consequence, with AMR popmax AF=2.18e-4 and global AC=40. A second TRDN frameshift, VCV001074440, reached the high-tension band but did not pass the red actionability gate. These observations do not prove reduced function, penetrance, or mechanism. They show that even supposedly severe variant classes can occupy frequency-function discordant zones and should not be treated as uniformly textbook high-penetrance evidence without AC, ancestry, review, phenotype, and functional context.

### VITAL red variants and explainable case-level prioritization

The three VITAL-red variants were the highest-priority re-review candidates as of the April 21, 2026 data freeze:

| Gene | ClinVar ID | Variant summary | Key frequency signal | Review support | VITAL | Weight stability |
|---|---|---|---|---|---:|---|
| SCN5A | VCV000440850 | c.[3919C>T;694G>A] | global AF=1.46e-4, global AC=214, AFR popmax AF=5.68e-3, popmax AC=190 | no assertion criteria provided | 96.1 | anchor |
| TRDN | VCV001325231 | c.1050del (p.Glu351fs) | global AF=2.77e-5, global AC=40, AMR popmax AF=2.18e-4 | single submitter | 74.3 | near-stable |
| KCNH2 | VCV004535537 | c.2398+2T>G | global AF=4.37e-5, global AC=24, ASJ popmax AF=1.32e-4 | single submitter | 70.3 | borderline |

These examples illustrate why VITAL is not simply an AF threshold. The naive AF screen flags 115 variants; VITAL prioritizes only 3 because the strongest actionable signals require AC support and fragile review. For example, CACNB2 VCV003774534 has a high VITAL score (71.3) and popmax AF >1e-4 but is not red because qualifying AC is below 20. Conversely, TRDN VCV001325231 remains red despite being a deletion because global AC=40 supports the frequency signal and the assertion is single-submitter.

### Review fragility is an explicit result

Review quality was not hidden inside a black-box score. Of the 9 AC-supported frequency signals, 6 came from multiple-submitter/no-conflict records and were not automatically marked red. All 3 VITAL-red variants came from review-fragile records: 2/3 single-submitter and 1/3 weak/no-assertion. This makes review fragility auditable and clinically interpretable.

### Historical ClinVar analysis: preliminary enrichment in a sparse red set

This analysis is included in the main text as a transparency and reproducibility audit, not as a primary validation claim. Its evidentiary weight is intentionally lower than the current workflow-compression, explainability, and detectability-bias results. The red set is too small for stable inference, and all historical estimates below should be read as secondary, exploratory stress tests.

The January 2023 ClinVar snapshot contained 1,669 baseline arrhythmia P/LP variants suitable for historical validation. By the current 2026 snapshot, 17 met the strict endpoint of P/LP to B/LB or VUS, and 113 met the broader endpoint of P/LP to non-P/LP, conflicting, other, or missing follow-up classification.

VITAL red identified only 2 baseline variants. This small red-set size is a design feature of the framework: VITAL is intended to create a short extreme-priority review queue, not a broad classifier. For the strict endpoint, neither red variant became B/LB or VUS (0/2; 95% CI 0.00-0.66). For the broader endpoint, 1/2 VITAL-red variants later destabilized, compared with 113/1,669 variants overall. Thus, the broad-event rate was 50.0% among VITAL-red variants (95% CI 9.5%-90.5%) versus 6.8% overall (95% CI 5.7%-8.1%), a 7.38-fold enrichment. Compared with non-red variants, the broad-event enrichment was 7.44-fold (95% CI 2.36-23.31). Restricting to the 455 frequency-observed baseline variants gave 1/2 VITAL-red events versus 38/455 overall, a 5.99-fold enrichment and 6.12-fold enrichment versus non-red variants (95% CI 1.87-19.55).

The historical result should therefore be interpreted as exploratory, qualitative enrichment rather than definitive predictive validation. A one-event perturbation analysis illustrates the instability of the apparent fold enrichment:

| Scenario | Red events / red N | Red event rate | Fold enrichment vs overall broad-event rate |
|---|---:|---:|---:|
| Observed historical result | 1/2 | 50.0% | 7.38x |
| If the single red event were absent | 0/2 | 0.0% | 0.00x |
| If one additional red event occurred | 2/2 | 100.0% | 14.77x |
| If one additional red non-event were included | 1/3 | 33.3% | 4.92x |
| If one additional red event and one non-event were included | 2/3 | 66.7% | 9.85x |

This sensitivity is expected for an extreme-priority subset with n=2. The historical snapshot is therefore best read as a directional consistency check: VITAL red may concentrate future instability, but the estimate is not stable enough to support a stand-alone prediction claim. The primary evidence for clinical utility remains the current workflow result: large reduction in false-positive re-evaluation burden with component-level explainability.

Historical AC sensitivity supported the same modular interpretation. For the broad endpoint, AC>=5 produced 40 AC-supported frequency flags and 3 VITAL-red variants (ANK2 VCV000207942, KCNE1 VCV001202620, and TRDN VCV001325231), with 1/3 broad future events. AC>=10 and AC>=20 selected the same two VITAL-red variants (KCNE1 VCV001202620 and TRDN VCV001325231), with 1/2 broad future events at both thresholds. AC>=50 selected only KCNE1 VCV001202620, which was the broad future event. The AC>=10 and AC>=20 red-set identity was therefore stable, while AC>=5 added one lower-AC non-event and AC>=50 became more stringent by losing TRDN. The corresponding AC-supported baseline screens produced substantial false-positive burden across thresholds: 36 false positives at AC>=5, 34 at AC>=10, 25 at AC>=20, and 7 at AC>=50 for the broad endpoint.

Historical threshold calibration identified and documented a field-mapping inconsistency in the first calibration pass. The field-mapping-inconsistent run produced 0 red-gate-compatible flags at all thresholds because the historical prediction table lacked `weak_review_signal` and the calibration script defaulted the field to false. This run is retained as an audit artifact. After applying validator-consistent mapping (`review_score <=1` or `submitter_count <=1` for weak review, with AC support derived from `frequency_signal_ac_ge_20` or the `qualifying_frequency_ac >=20` fallback), the threshold-70 historical red set matched the validator: KCNE1 VCV001202620 and TRDN VCV001325231. KCNH2 VCV004535537 was not red in the 2023 snapshot because it was not present as the same allele-specific ClinVar record: the January 2023 table contained the same splice position as VCV000560686, c.2398+2T>A (A>T), with AC=0 and non-red status, whereas the current red record VCV004535537, c.2398+2T>G (A>C), was created/submitted on December 14, 2025 and has current global AC=24 with single-submitter review. Thus, KCNH2's current red status reflects entry of a later allele-specific assertion with sufficient AC support, not review-score drift within the same VCV record; this supports treating it as a temporally new and weight-borderline signal.

The calibration result did not show that 70 is the enrichment-maximizing historical threshold. For the arrhythmia broad endpoint, the maximum red-gate-compatible enrichment occurred at threshold 80 (1/1 broad events; 14.77x enrichment), whereas the prespecified threshold 70 flagged 2 variants with 1 event (7.38x enrichment). For the expanded endpoint, lower thresholds selected more events, again reflecting the instability expected from very small red sets. We therefore retain 70 as a prespecified operational threshold, not as a fitted optimum. The calibration analysis is useful because it confirms internal consistency of the validator, exposes threshold sensitivity, and documents that enrichment maxima are sparse and unstable.

### Multi-domain historical validation and alignment sanity checks

The expanded historical dataset combined arrhythmia (n=1,669), cardiomyopathy (n=300), epilepsy (n=300), hearing loss (n=300), and random ClinVar P/LP variants (n=500), for 3,069 raw baseline rows. Within-domain duplicate checks found 0 duplicate variation IDs and 0 duplicate variant keys. Cross-domain leakage occurred only through the random ClinVar sample: 6 variants overlapped a fixed disease panel and were removed from pooled analysis using the prespecified domain-priority rule. No VITAL-red rows were removed by within-domain or pooled deduplication (red count remained 4 before and after dedupe), leaving 3,063 pooled-deduplicated variants.

Within-domain results were heterogeneous. Arrhythmia contributed 2 VITAL-red variants; random ClinVar P/LP contributed 2; cardiomyopathy, epilepsy, and hearing loss contributed 0 red variants at the prespecified AC>=20 red gate. For the broad endpoint, arrhythmia had 1/2 red events and random ClinVar P/LP had 1/2 red events. For the expanded endpoint, arrhythmia had 2/2 red events and random ClinVar P/LP had 1/2. The fixed external disease panels therefore increased the historical denominator and event count but did not create additional red calls, consistent with conservative behavior outside the arrhythmia development context.

Historical threshold calibration showed the same conservative pattern at the prespecified threshold. At score threshold 70 with the red-gate-compatible actionability filter, cardiomyopathy, epilepsy, and hearing loss each had 0 red flags. Lower exploratory thresholds did not produce explosive queues: cardiomyopathy and epilepsy remained at 0 red-gate-compatible flags, and hearing loss reached at most 2 flags only at lower thresholds. Thus, adding external historical domains increased denominator size and endpoint heterogeneity without turning VITAL into a broad cross-domain flag generator.

Pooled-deduplicated results were directionally consistent but still sparse. The strict endpoint remained negative: 0/4 VITAL-red variants reached P/LP-to-B/LB-or-VUS despite 24 strict events overall. For the broad endpoint, 2/4 VITAL-red variants reached broad instability compared with 174/3,063 overall (red event rate 50.0%, Wilson 95% CI 15.0%-85.0%; enrichment 8.80x vs overall). For the expanded endpoint, 3/4 VITAL-red variants had broad instability, aggregate clinical-significance change, or review-status change compared with 818/3,063 overall (red event rate 75.0%, Wilson 95% CI 30.1%-95.4%; enrichment 2.81x vs overall). These results support a qualitative enrichment signal under broader instability definitions but still do not justify a stand-alone predictive claim because recall remained very low (1.15% for broad and 0.37% for expanded).

Stratified bootstrap confirmed that the pooled historical signal remains sparse and discrete. For the strict endpoint, 100% of bootstrap samples with at least one red variant had 0 red events. For the broad endpoint, 10.9% of usable bootstrap samples had 0 red events and 30.7% had exactly 1 red event; the median red event rate was 50.0%, but the 95% bootstrap interval spanned 0%-100%. For the expanded endpoint, 2.5% had 0 red events and 15.9% had exactly 1 red event; the median red event rate was 75.0%, with a 14%-100% interval. This is binary-noise sensitivity, not smooth statistical stability.

Alignment sanity checks identified 18 baseline variants missing by VariationID in the 2026 snapshot. Of these, 13 appeared to have disappeared completely or lacked a current GRCh38 key match, while 5 had the same variant key present under a current record and were therefore classified as mapping failures or VariationID changes. The highest-scoring missing case was a random ClinVar MUC4 variant (VITAL=55.5, yellow watchlist), not a red-priority call. Thus, missing follow-up records contributed to broad endpoint counts but did not drive the VITAL-red historical signal.

### External disease panels and ratio compression

Small external-domain stress tests showed that VITAL does not generate large false-positive red queues outside the arrhythmia context:

| Domain | N | Exact AF rows | Naive AF flags | AC-supported flags | VITAL red | Score >=70 | Max score |
|---|---:|---:|---:|---:|---:|---:|---:|
| Cardiomyopathy | 300 | 64 | 23 | 2 | 1 | 1 | 80.8 |
| Epilepsy | 300 | 20 | 5 | 0 | 0 | 0 | 59.9 |
| Hearing loss | 300 | 120 | 76 | 16 | 0 | 3 | 77.4 |
| Random ClinVar P/LP | 500 | 167 | 94 | 20 | 0 | 7 | 85.1 |
| BRCA/MMR/APC descriptive comparator | 2,000 | 173 | 5 | 2 | 0 | 0 | 57.8 |

This is a false-positive-burden and workload-compression result, not a sensitivity result. The panels are deliberately small, and approximately 300 variants per domain is insufficient to determine whether VITAL would recover true reclassification candidates outside arrhythmia. In the random ClinVar P/LP sample, a naive AF screen produced 94/500 actionable flags (18.8%), whereas VITAL red produced 0/500. In hearing-loss genes, naive AF produced 76/300 flags (25.3%) and VITAL red again produced 0/300. This represents a >10-fold reduction in actionable review flags without manual curation. Importantly, VITAL did not flatten external datasets into uniformly green calls: hearing-loss and random ClinVar samples retained high-score yellow/orange tails, but those tails did not become red-priority calls without the full combination of AC support and review fragility. Because independent reclassification truth labels were not available for these external domains, these tests cannot determine whether VITAL missed true reclassification candidates. They show only that the framework does not manufacture large red queues outside the arrhythmia development context.

The BRCA/MMR/APC descriptive comparator had 0 red calls, 0 variants with VITAL score >=60, 0 variants with score >=70, and a maximum score of 57.8. This is reassuring as a portability check, but it should not be overinterpreted as external validation. Because these genes differ from arrhythmia genes in penetrance, founder architecture, and disease mechanism, the result only supports a narrower conclusion: VITAL does not manufacture red-priority calls in this non-arrhythmia comparator set.

### Exome-vs-genome sensitivity

All 350 exact allele-level arrhythmia matches were queried in gnomAD v4.1 genome data. Genome data did not materially change the duplication observation. Among 22 exact-matched duplications, only 2 had any genome AC >0, and none had genome AF >1e-5. Across all exact matches, genome-only AC-positive recovery occurred for 14 variants, but the central frequency contradictions were already visible in exome data. These results indicate that switching from exome to genome gnomAD data does not rescue duplication representation at clinically meaningful frequency thresholds in the current release, but they do not by themselves establish the underlying mechanism.

### KCNH2 diagnostic dissection

KCNH2 had high non-overlap: 406/478 variants (84.9%) lacked exact usable gnomAD representation. However, the updated analysis does not support a simple claim that KCNH2 is uniquely elevated by overall non-overlap alone. SCN5A had a nearly identical overall non-overlap rate (331/394, 84.0%; KCNH2 vs SCN5A p=0.708). The informative signal is instead the variant-type composition of KCNH2. We use "driven by duplications" in a statistical sense: duplications account for most of the observed KCNH2 excess signal. This should not be read as a resolved mechanistic explanation for why those duplications are absent or differently represented in gnomAD.

KCNH2 duplications showed 97/102 non-overlap (95.1%) compared with 75/92 duplications in other arrhythmia genes (81.5%; OR=4.40, p=0.00265). SNVs and deletions did not show comparable KCNH2-specific enrichment. Removing all KCNH2 duplications reduced the excess non-overlap count from +24.6 to +9.0, demonstrating that duplications account for most of the KCNH2 diagnostic signal. The mechanism remains unresolved and could include technical representation artifacts, local sequence-context effects not captured by repeat/GC annotations, true biological depletion from population cohorts, ClinVar ascertainment bias, or some combination of these factors. The analysis therefore identifies a reproducible statistical pattern and a clinical caution, not a resolved molecular or sequencing mechanism.

Duplication size analysis showed that short duplications dominated the KCNH2 duplication set: 90/102 duplications were 1-10 bp, and 85/90 (94.4%) were non-overlap. All medium duplications (9/9) and long duplications (3/3) were non-overlap. Repeat annotations did not explain the signal: only 16/97 non-overlap KCNH2 duplications (16.5%) overlapped a RepeatMasker feature within a 10-bp window, and RepeatMasker overlap was not enriched among non-overlap variants overall.

Local GC content also did not provide a sufficient explanation for KCNH2 absence. Exact-matched KCNH2 variants had median local GC100 of 0.706, whereas non-overlap variants had median local GC100 of 0.664, opposite the expectation if high GC alone drove non-overlap. The SCN5A comparison further argued against a simple GC explanation: KCNH2 and SCN5A had similar non-overlap rates despite different local sequence and variant-type structures.

Submission-date analysis showed only a weak, non-significant trend toward lower non-overlap in more recent KCNH2 variants. Variants created in or before 2022 had non-overlap of 253/291 (86.9%), compared with 153/187 (81.8%) for variants created after 2022 (Fisher p=0.082). This suggests possible improvement over time but does not support a strong temporal conclusion.

## Discussion

This study provides a population-frequency and technical-detectability audit of ClinVar P/LP variants in inherited arrhythmia genes, and extends that audit into an interpretable clinical reclassification risk prioritization framework. The central finding is not simply that most ClinVar P/LP variants are rare. That result is expected. The stronger finding is that population-frequency evidence becomes clinically useful only when interpreted through ancestry, AC support, variant representation, gene context, technical detectability, and review quality.

First, arrhythmia P/LP variants with usable gnomAD AF are strongly concentrated at ultra-rare global frequencies. Among AF-covered variants, 96.1% had AF <=1e-5. This supports the broad validity of frequency constraint in high-penetrance arrhythmia genes. At the same time, global AF alone missed most frequency tension. Popmax identified 115 variants above 1e-5, while global AF identified only 13. Nearly one-third of variants with complete popmax data were globally rare but population-enriched. This reinforces the need for ancestry-aware interpretation.

Second, non-overlap is not equivalent to rarity. More than 80% of arrhythmia P/LP variants lacked usable exact AF evidence, but that set combined allele-discordant sites, no-record variants, and exact matches without AF. Indels and duplications were much more likely than SNVs to lack exact frequency evidence. This is the conceptual reason VITAL includes technical detectability and preserves a gray no-frequency-evidence state instead of treating missing AF as zero.

Third, naive AF screening creates excessive manual review burden. A popmax/global AF >1e-5 screen flagged 115 arrhythmia P/LP variants, but only 9 had AC-supported frequency evidence and only 3 met VITAL red criteria. This 97.4% compression in actionable calls is not cosmetic. It turns an unfocused frequency alert list into a short clinical reclassification risk queue with explicit reasons and component-level explainability. In practical workflow terms, the framework removes 112 urgent review triggers relative to a naive AF screen and eliminates 53 proxy false-positive re-evaluation calls in the operational benchmark.

Fourth, the pathogenicity tension continuum reframes high-frequency P/LP variants as a biological and clinical category rather than only a curation error list. High VITAL variants are frequency-function discordant: they carry population-frequency evidence that is difficult to reconcile with a straightforward high-penetrance interpretation, while also retaining functional labels such as LOF, splice/intronic, or unresolved composite annotation. LOF subtype analysis strengthens this point: frameshift, stop_gained, and canonical_splice variants showed distinct discordance profiles, and frameshift variants were not absent from high-tension space. This pattern is consistent with variants that may not conform to canonical high-penetrance expectations, including clinically overcalled assertions, low-penetrance alleles, hypomorphic effects, ancestry-specific enrichment, transcript-specific effects, or context-dependent pathogenicity. The clinical consequence is concrete: overcalling a low-penetrance or frequency-inconsistent variant as fully penetrant P/LP can trigger unnecessary cascade testing, patient anxiety, reproductive counseling consequences, medication or activity restrictions, or even inappropriate escalation of cardiac surveillance and device discussions. VITAL's role is to bring these cases to expert review, not to replace that review.

Fifth, VITAL's review-fragility term is clinically important. Many frequency-positive P/LP records have stronger review support, and VITAL does not automatically demote them. In the current arrhythmia cache, all red-priority calls came from weak or single-submitter records. This design makes VITAL a decision-support framework rather than a replacement for expert classification: it decides which records deserve urgent review, not what their final classification should be. These results highlight that frequency-based prioritization operates orthogonally to inheritance-aware and phenotype-driven interpretation, complementing rather than replacing established clinical frameworks.

Sixth, the AC sensitivity analysis addresses the concern that the framework depends on an arbitrary threshold. Because the score and the AC gate are separated, the same signal can be viewed under different operational requirements. The current AC>=20 red-core remained stable at AC>=5 and AC>=10, while lower and higher gates transparently added or removed variants. Historical AC sensitivity showed the same pattern: AC>=10 and AC>=20 selected the same two baseline red variants, while AC>=5 added one lower-AC non-event and AC>=50 retained only the strongest case. This modular design is preferable to hiding the threshold inside the score.

Seventh, weight sensitivity shows that VITAL is not creating arbitrary new red calls when component weights are perturbed. Five alternative expert-weight profiles produced no gained red variants, maintained 0 proxy false positives, and preserved high rank correlation with the primary score. However, this analysis also prevents overclaiming: SCN5A is an anchor call, TRDN is near-stable, and KCNH2 is borderline. The KCNH2 result is still a valid re-review signal under the primary framework, but it should be communicated as threshold-adjacent and weight-sensitive.

Eighth, the historical analysis is deliberately secondary. It is retained in the main text to document the audit trail, alignment checks, and negative evidence against overclaiming, not to give the impression of a powered predictive validation. In arrhythmia alone, only 2 baseline 2023 variants were VITAL red, and 1/2 later reached the broad instability endpoint, compared with 113/1,669 overall. In the expanded multi-domain analysis, pooled-deduplicated n increased to 3,063 and the VITAL-red set increased to 4, but the same caution remains: strict events were not enriched (0/4), broad instability was 2/4, and the expanded endpoint was 3/4. These numbers are directionally consistent but still sparse. Historical threshold calibration confirmed that the prespecified threshold 70 is not the enrichment-maximizing threshold in the sparse 2023-to-2026 dataset; it remains an operational review threshold rather than a fitted optimum. The appropriately modest interpretation is that VITAL may identify an extreme-priority subset for monitoring; it does not capture most future ClinVar changes and should not be presented as a stand-alone predictor.

Ninth, external-domain stress tests argue against simple arrhythmia overfitting, but they are not a strong external validation set. Cardiomyopathy produced one VITAL-red call among 300 current sampled variants, while epilepsy, hearing loss, random ClinVar P/LP, and the BRCA/MMR/APC descriptive comparator produced none. In the historical baseline, cardiomyopathy, epilepsy, and hearing loss also produced 0 red calls at AC>=20, while the random ClinVar sample produced 2. The hearing-loss and random samples had many naive AF flags but no red-priority explosion. This shows conservative behavior outside the development context and reframes VITAL as workflow burden reduction: fewer false-positive re-evaluation flags without manual curation. The BRCA/MMR/APC result is especially limited because those genes are affected by incomplete penetrance and founder effects.

Finally, the KCNH2 analysis illustrates why gene-level non-overlap signals require diagnostic dissection. KCNH2 is not uniquely explained by global non-overlap when compared with SCN5A. Instead, much of its statistical excess is accounted for by a high duplication burden and poor duplication representation. GC content, repeat annotation, submission date, and genome data do not provide simple alternative explanations, but they also do not resolve mechanism. The practical implication remains important: absence of a KCNH2 duplication from gnomAD should not be treated as strong evidence of rarity without orthogonal validation.

## Comparison with previous work

Previous work has examined ClinVar-gnomAD concordance and cardiac channelopathy frequency thresholds, but most studies emphasize frequency distributions or gene-level constraint rather than explicit detectability, review-quality modeling, and workflow burden. This analysis extends that literature in four ways. First, it separates exact allele observation, allele discordance, and no-record states rather than collapsing all nonmatches into absence. Second, it places popmax and AC reliability at the center of frequency-based reclassification risk prioritization. Third, it adds technical detectability as an explicit model component, addressing the otherwise hidden assumption that absent variants are equally detectable across variant classes. Fourth, it probes score behavior using a sparse historical ClinVar snapshot and small external disease-domain stress tests, while treating both as supportive but non-definitive validation layers.

## Limitations

This study has several limitations. First, VITAL is a clinical reclassification risk prioritization framework, not a clinical truth model. Red calls should trigger expert review, not automatic reclassification. The final decision must remain with human reviewers applying ACMG/AMP criteria, phenotype, segregation, functional data, disease mechanism, inheritance, and laboratory policy. Second, the component maxima are expert-specified rather than empirically learned. Weight-profile sensitivity showed that no new red variants appeared under five alternative profiles, but KCNH2 was unstable in 3/5 profiles; future disease-specific implementations should recalibrate and report weight sensitivity rather than assuming these weights transfer unchanged. Third, the current operational benchmark uses proxy labels because definitive prospective truth is unavailable, and those proxy labels are partly circular: weak ClinVar review support and AC-supported frequency tension help define the positive class and are also used in the red actionability gate. As a result, ROC/PR curves and precision/recall estimates should not be read as independent clinical performance metrics; the perfect precision and recall values are artifacts of design concordance under the proxy benchmark and are useful only for comparing review burden under different frequency-screening strategies. Fourth, the historical analysis is underpowered for strong predictive claims. The expanded multi-domain analysis increases the denominator but does not eliminate sparsity of red calls; stratified bootstrap behaves as binary-noise sensitivity rather than a smooth stable effect-size estimator. The appropriate interpretation is preliminary qualitative enrichment rather than proof of future reclassification prediction. Fifth, historical threshold calibration is descriptive and was not used to refit the model. The enrichment-maximizing threshold varied by endpoint and domain, confirming that the prespecified threshold 70 should be treated as an operational review gate rather than an optimized predictive cutoff. Sixth, pooled multi-domain analysis mixes different disease architectures, baseline event rates, and error mechanisms. We therefore report within-domain metrics alongside pooled-deduplicated metrics, and the pooled estimate should be read as a sensitivity analysis rather than a primary endpoint. Seventh, cross-domain deduplication requires an operational priority rule. We used a prespecified domain priority followed by maximum VITAL score as a tie-breaker, and no red calls were removed by deduplication; nevertheless, a variant that is red in one domain and non-red in another could be affected by this rule in future datasets. Eighth, the historical endpoint is imperfect: some variants disappear from follow-up, change VariationID, or change terminology without representing true biological reclassification. We classify missing follow-up records as disappeared completely versus mapping failure when possible, but some ambiguity remains. Ninth, external disease-domain panels were intentionally small and should be interpreted as false-positive-burden and portability stress tests, not sensitivity analyses; without independent truth labels, they cannot show whether true reclassification candidates were missed. Tenth, the BRCA/MMR/APC panel is retained only as a descriptive comparator, not a clean negative control, because cancer-predisposition genes have founder effects, incomplete penetrance, and disease-mechanism differences from inherited arrhythmia genes. Eleventh, LOF subtype annotation was derived from ClinVar title/HGVS strings and was used only descriptively; it does not replace transcript-aware consequence annotation, NMD prediction, exon/domain context, MAVE data, or curated functional assays. Twelfth, the KCNH2 duplication result is statistical rather than mechanistic: duplications account for much of the observed non-overlap signal, but the underlying cause remains unresolved. Thirteenth, short-read gnomAD data remain limited for structurally complex variants, and even genome data may not fully resolve duplications or complex indels. Finally, the primary data freeze was April 21, 2026. ClinVar assertions can change after that date, so all red-priority variants should be rechecked against live ClinVar records before publication, clinical reporting, or any patient-facing decision.

## Clinical and research implications

These results support several practical recommendations for clinical genomics workflows:

1. Use ancestry-aware frequency thresholds. Popmax identifies frequency tension that global AF misses.

2. Require AC support before acting on frequency contradictions. AC >=20 is not the only possible threshold, but it prevents overreaction to one- or two-allele observations.

3. Do not treat absence as rarity for all variant types. Indels and duplications require special caution because exact gnomAD representation is systematically lower.

4. Separate review prioritization from clinical classification. VITAL red should initiate expert re-review, not automatically assign benignity.

5. Report explainable component scores. AF pressure, AC reliability, detectability, gene context, and review fragility should be visible for each variant.

6. Use external stress tests and historical snapshots. Frameworks that audit clinical databases should show that they behave conservatively outside the development context and, where possible, show preliminary enrichment for future instability with uncertainty reported explicitly.

7. Track review burden as a clinical endpoint. A framework that removes false-positive re-evaluation calls can improve laboratory throughput, shorten review queues, and focus expert decision-making on variants most likely to change patient-facing interpretation.

## Conclusion

Clinically cataloged P/LP arrhythmia variants are strongly concentrated at ultra-rare population frequencies, but the interpretation of frequency evidence is not reducible to a single AF cutoff. Popmax, AC, variant type, technical detectability, gene context, and ClinVar review strength all materially change which variants should be prioritized for re-review. VITAL formalizes these signals into an explainable 0-100 clinical reclassification risk score, identifies 3 high-priority arrhythmia assertions as of the April 21, 2026 data freeze, suppresses more than 100 naive AF alerts, and preserves explicit no-frequency-evidence states for variants that should not be treated as AF=0. Secondary historical audits across 3,063 pooled-deduplicated variants remain underpowered and should be interpreted only as qualitative consistency checks. The intended clinical value is therefore not prediction in isolation; it is reducing false-positive variant re-evaluation burden without replacing clinician or laboratory decision-making.

## Data availability

Analyses were conducted in April 2026 using gnomAD v4.1.1 exome and genome data queried through the gnomAD GraphQL API. ClinVar data were retrieved through the NCBI Entrez API and the ClinVar bulk variant summary file, including the January 2023 archived snapshot for historical enrichment analysis. UCSC genome sequence and annotation APIs were used for GC and repeat-context analyses. Pipeline source code, cached intermediate files, figures, and machine-readable outputs are available in the project repository.

Key output files include:

- `data/processed/arrhythmia_vital_scores.csv`: per-variant VITAL scores, bands, reasons, and component signals.
- `data/processed/arrhythmia_vital_component_breakdown.csv`: long-format explainability table.
- `data/processed/arrhythmia_vital_method_comparison.csv`: baseline comparison against AF-only and AF+AC screens.
- `data/processed/arrhythmia_vital_threshold_sweep.csv`: sensitivity/specificity sweep by VITAL cutoff.
- `data/processed/arrhythmia_vital_ac_threshold_sensitivity.csv`: AC>=5/10/20/50 actionability-gate sensitivity with red-set composition.
- `data/processed/arrhythmia_vital_weight_sensitivity_summary.csv`: red-set, compression, false-positive, and rank-correlation summary across expert-weight profiles.
- `data/processed/arrhythmia_vital_weight_sensitivity_variant_matrix.csv`: current red variants classified as anchor, near-stable, or borderline across alternative weight profiles.
- `data/processed/arrhythmia_vital_weight_sensitivity_variant_profile_table.csv`: long-format 5 alternative profiles x 3 red variants weight-sensitivity table.
- `data/processed/arrhythmia_vital_pathogenicity_tension_continuum.csv`: per-variant pathogenicity tension zone and frequency-function discordance class.
- `data/processed/arrhythmia_vital_frequency_function_discordance_summary.csv`: 20-point VITAL-band continuum summary.
- `data/processed/arrhythmia_vital_signal_reorganization_summary.csv`: naive AF, AC-supported AF, VITAL-band, and VITAL-red signal reorganization summary.
- `data/processed/arrhythmia_vital_lof_subtype_discordance_summary.csv`: frameshift, stop_gained, and canonical_splice discordance summary across the VITAL continuum.
- `data/processed/arrhythmia_vital_acmg_disagreement.csv`: naive AF flags split by VITAL agreement or non-red status.
- `data/processed/arrhythmia_vital_top_suspicious.csv`: top VITAL-priority variants.
- `data/processed/arrhythmia_vital_absence_detectability_bias.csv`: variant-type evidence for absence/detectability bias.
- `data/processed/arrhythmia_2023_01_to_current_vital_historical_validation.csv`: arrhythmia 2023-to-current historical enrichment analysis.
- `data/processed/arrhythmia_2023_01_to_current_vital_historical_enrichment.csv`: historical enrichment analysis.
- `data/processed/arrhythmia_2023_01_to_current_vital_historical_ac_threshold_sensitivity.csv`: historical AC threshold sensitivity and red-set stability.
- `data/processed/arrhythmia_2023_01_to_current_vital_historical_threshold_calibration.csv`: exploratory historical threshold calibration.
- `data/processed/arrhythmia_2023_01_to_current_vital_historical_threshold_calibration_field_mapping_inconsistent.csv`: retained pre-fix calibration audit artifact showing the zero-flag field-mapping inconsistency.
- `data/processed/arrhythmia_2023_01_to_current_vital_historical_threshold_calibration_field_mapping_audit.csv`: corrected calibration field-mapping audit, including KCNH2 historical red status.
- `data/processed/arrhythmia_2023_01_to_current_vital_historical_threshold_calibration_mapping_fix_summary.csv`: before/after summary of the calibration field-mapping fix.
- `data/processed/combined_2023_01_to_current_vital_historical_endpoint_summary.csv`: within-domain and pooled-deduplicated strict/broad/expanded historical endpoint summary.
- `data/processed/combined_2023_01_to_current_vital_historical_sanity_domain_summary.csv`: baseline n, red counts, duplicates, and missing follow-up counts by domain.
- `data/processed/combined_2023_01_to_current_vital_historical_alignment_breakdown.csv`: missing VariationID breakdown into disappeared versus mapping-failure categories.
- `data/processed/combined_2023_01_to_current_vital_historical_cross_domain_dedupe_decisions.csv`: cross-domain duplicate decisions with kept domain and reason_kept.
- `data/processed/combined_2023_01_to_current_vital_historical_stratified_bootstrap_summary.csv`: stratified bootstrap binary-noise sensitivity summary.
- `data/processed/combined_2023_01_to_current_vital_historical_threshold_calibration.csv`: pooled and within-domain historical threshold calibration.
- `data/processed/combined_2023_01_to_current_vital_historical_threshold_calibration_mapping_fix_summary.csv`: pooled and within-domain before/after field-mapping fix summary.
- `data/processed/vital_external_panel_summary.csv`: external-domain and control stress-test summary.
- `data/processed/vital_external_panel_score_distribution.csv`: external score distribution table.
- `figures/vital_clinical_workflow.png`: clinician-facing workflow schematic linking ClinVar input, VITAL risk prioritization, review queue, and expert decision-making.
- `figures/vital_external_panel_score_distribution.png`: external score distribution figure.

Supplementary tables include:

- Supplementary Table S1: ancestry-specific enrichment of globally rare variants by gene and population.
- Supplementary Table S2: complete list of naive frequency-flagged variants with AF/AC fields and reclassification rationale.
- Supplementary Table S3: KCNH2 diagnostic data, including duplication size, repeats, GC content, and submission dates.
- Supplementary Table S4: exome-vs-genome sensitivity analysis.
- Supplementary Table S5: per-gene frequency outlier and non-overlap summaries.
- Supplementary Table S6: VITAL scores.
- Supplementary Table S7: VITAL method comparison.
- Supplementary Table S8: review fragility summary.
- Supplementary Table S9: VITAL threshold sweep.
- Supplementary Table S10: ACMG/VITAL disagreement table.
- Supplementary Table S11: top suspicious variants.
- Supplementary Table S12: absence/detectability bias.
- Supplementary Table S13: VITAL expert-weight sensitivity across five alternative profiles and three current red variants.
