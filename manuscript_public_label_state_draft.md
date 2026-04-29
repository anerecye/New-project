# Population- and mechanism-aware constraints expose unstable actionability of public pathogenic variant labels

## Abstract

Public pathogenic/likely-pathogenic (P/LP) variant labels are widely reused as portable disease-state claims. After downstream flattening of condition- and submission-level context, however, the disease model implied by a P/LP label may no longer be recoverable. This study tests whether flattened labels reliably preserve dominant disease-model compatibility under population-frequency stress and whether such labels can be safely reused as direct-actionability proxies without an explicit constraint layer.

We examined 1,731 ClinVar P/LP variants across 20 inherited arrhythmia genes, of which 17 contributed P/LP records. Only 357 variants (20.6%) achieved allele-resolved frequency recovery in gnomAD v4.1.1 exomes, and primary frequency conclusions were therefore restricted to the 334 variants with usable allele-resolved AF. Within that evaluable subset, ancestry-aware popmax screening identified frequency alerts missed by 102 of 115 global-only queries (88.7%), including 103 alerts in potentially consequential interpretation contexts such as cascade testing, medication guidance, surveillance intensity, or device-related management.

Applying maximum-credible-allele-frequency (MCAF) stress-testing classified the evaluable exome subset into operationally distinct regimes. Here, decomposition denotes only disease-model incompatibility under population-frequency stress, not a claim that the variant is benign or incorrectly classified. Within the usable-AF Tier 1 subset, incompatibility with an unqualified dominant high-penetrance model was mechanism-dependent: gain-of-function/dominant-negative genes showed zero observed decomposition among evaluable variants (0/34; 95% CI 0-10.2%), haploinsufficient genes showed intermediate decomposition, and recessive genes showed the highest decomposition. This mechanism gradient is therefore a conditional allele-resolved result, not an estimate for non-evaluable variants. Supplementary HGDP and hereditary-cancer analyses were retained only as secondary stress checks rather than as co-equal burden estimates.

We then translated these population and mechanism constraints into VITAL, a post-label actionability routing layer. Under a label-driven baseline, all 1,731 arrhythmia P/LP variants entered a direct-actionable route. After minimal portability constraints, only 219/1,731 labels (12.7%) could be directly routed as actionability-compatible. The remaining 1,512/1,731 decisions (87.3%) were rerouted, but this aggregate intentionally combines distinct reasons: 1,397 (80.7%) were evidence-unavailable EVAL_LIMITED deferrals, 77 (4.4%) required ancestry-aware frequency review, and 38 (2.2%) required disease-model repair. Most failures therefore reflected missing allele-resolved evaluability rather than evidence conflict.

These findings operationalize a measurable downstream reuse failure mode with two separable components: evidence unavailable for most labels, and evidence conflicting or model-incomplete for a smaller allele-resolved subset. Within the evaluable subset, the portability of a P/LP label as a disease-state claim depends on disease model, molecular mechanism, and population-frequency resolution, none of which are currently encoded in the exported label. Population-structured frequency review should therefore be treated as a minimum evidence layer in downstream interpretation workflows, while preserving explicit uncertainty for non-evaluable variants. VITAL does not reclassify variants; it reroutes actionability when flattened labels no longer carry the constraints required for direct reuse.

## Central Claim

Public P/LP is classification evidence, not portable actionability. Direct actionability should be granted only when the label survives explicit allele-resolved evaluability, ancestry-aware frequency, and disease-model constraints.

## Introduction

Public clinical variant assertions now function as shared infrastructure. A pathogenic or likely pathogenic label exported from ClinVar is routinely reimported into laboratory interpretation pipelines, cascade-testing workflows, decision-support systems, review spreadsheets, VCF annotation tools, screening filters, and research analyses as though it were a portable unit of clinical meaning. This pattern of reuse rests on an implicit assumption: that the public label preserves the disease model it was originally intended to describe.

That assumption is structurally unsafe. A public P/LP assertion records a classification outcome but not necessarily the disease-model parameters required for reliable downstream reuse, including inheritance, penetrance, ancestry context, representation stability, molecular mechanism, and whether population evidence was evaluable at all. This study therefore does not adjudicate original ClinVar submissions. It models a downstream reuse failure mode in which exported P/LP alleles are consumed after condition- and submission-level context has been flattened.

This flattening scenario is routine. Laboratory information systems, VCF annotation pipelines, decision-support tools, review spreadsheets, and screening filters often import ClinVar as a categorical flag. Concrete examples include VEP or ANNOVAR consequence tables that append ClinVar significance to a VCF row, laboratory spreadsheets that sort candidate variants by P/LP status, clinical decision-support triggers that watch for actionable genes plus P/LP labels, and research screening filters that retain "pathogenic ClinVar" variants before phenotype-specific review. In each case, the consumer may receive a P/LP label without the disease model that justified it.

The concern is practical rather than semantic. In potentially consequential interpretation contexts such as syndrome diagnosis, family cascade testing, medication guidance, surveillance intensity, or device-related management, the key question is whether an exported label remains interpretable once ancestry-aware frequency, allele-resolved representation, and mechanism are restored downstream.

Inherited arrhythmia genes provide a stringent test domain for this problem. They span dominant, recessive, and mixed architectures; variable penetrance; ancestry-specific enrichment; founder effects; and clinically consequential interpretation contexts. They therefore expose both the biological and operational consequences of label reuse. This is not because arrhythmia genes are unique, but because they make the infrastructure problem easier to detect.

We examine three linked questions. First, how much of the public arrhythmia label space is population-evaluable at all under current representation systems? Second, within the exome Tier 1 usable-AF subset that serves as the primary analysis, how often does a label remain compatible with a dominant high-penetrance disease model when one moves from global AF to ancestry-aware popmax and then to frequency stress-testing? Third, when the same outputs are routed as actionability decisions, how often does a label-driven baseline require review, deferral, or model repair?

We show that most public arrhythmia assertions cannot be connected to allele-resolved population frequency even after exhaustive reconciliation and reference-based normalization. Among the minority for which allele-resolved context is available, ancestry-aware frequency exposes distinct interpretation regimes concealed beneath the same public P/LP label. Within that evaluable subset, the rate of disease-model decomposition is not random: it is structured by molecular mechanism and variant consequence. Outside that subset, the dominant finding is a transportability boundary rather than an estimable decomposition rate. Safe downstream reuse therefore requires more explicit metadata than current public labels provide.

## Methods

### Cohort Definition

We analyzed ClinVar P/LP assertions across 20 inherited arrhythmia genes: KCNQ1, KCNH2, SCN5A, KCNE1, KCNE2, RYR2, CASQ2, TRDN, CALM1, CALM2, CALM3, ANK2, SCN4B, KCNJ2, HCN4, CACNA1C, CACNB2, CACNA2D1, AKAP9, and SNTA1. Variants were collapsed to unique GRCh38 alleles by chromosome, position, reference allele, and alternate allele, yielding 1,731 unique variant records.

Throughout this manuscript, "P/LP variant record" refers operationally to a unique chr-pos-ref-alt allele carrying at least one public P/LP classification in ClinVar after collapsing across condition-specific submissions and submitters. This unit preserves allele identity but does not retain condition-level or submission-level semantics in the April 24, 2026 ClinVar data freeze used for this manuscript.

Of the 20 genes queried, 17 contributed P/LP variants to the final dataset. KCNE2, SCN4B, and SNTA1 had no P/LP variant records in the ClinVar freeze and are therefore absent from downstream analyses. Where this manuscript refers to "17 genes" in mechanism-stratified analyses or cross-domain comparisons, it reflects this data-driven reduction from the 20-gene query set. This deliberate flattening models the reuse scenario under investigation: downstream consumption of exported P/LP alleles after condition and submission context has been lost.

### Exome Reconciliation and Evaluability Tiering

Variants were cross-referenced against gnomAD v4.1.1 exomes using a multi-step reconciliation pipeline. Matching began with strict allele identity and then added trim-aware and decomposition-aware rescue. All variants were additionally normalized against a local GRCh38 reference using bcftools norm. Each variant was then assigned to one of three exome evaluability tiers:

| Tier | Definition | N | % |
| --- | --- | ---: | ---: |
| 1 | Exact or equivalent allele-resolved match with usable allele-level AF | 357 | 20.6 |
| 2 | Locus or regional context only; no allele-resolved AF | 1,326 | 76.6 |
| 3 | Genuinely unevaluable after all reconciliation steps | 48 | 2.8 |

Tier 1 contains the only variants used for allele-level frequency conclusions in the exome analysis. Tier 2 was never interpreted as population-consistent by implication. Because evaluability was not random with respect to variant class or gene architecture, all mechanism-stratified decomposition estimates are conditional on usable allele-resolved AF and are not transported to Tier 2.

### Exome Frequency Screening and Regime Assignment

For Tier 1 variants with usable AF (n = 334), we extracted global AF, popmax AF, allele count, review-strength metadata, and variant-class annotations. We used 1 x 10^-5 as a review trigger rather than as a pathogenicity boundary because it provides a conservative operational screen for dominant high-penetrance arrhythmia claims while remaining below frequencies that would already imply obvious incompatibility. We then compared global-AF-only review with global-or-popmax review.

Regime assignment within the 115 ancestry-aware alerts followed three categories: hard incompatibility with an unqualified dominant high-penetrance reading, boundary or monitoring status, and recessive or carrier-compatible architecture.

Throughout this manuscript, "disease-model decomposition" (or simply "decomposition"; equivalently, disease-model incompatibility under population-frequency stress) refers operationally to the failure of a flattened P/LP label to remain compatible with an unqualified dominant disease-state model under specified population-frequency constraints. A variant "decomposes" when its popmax allele frequency exceeds the maximum credible allele frequency (MCAF) for the tested disease model, indicating that the implied dominant high-penetrance interpretation is not sustainable at that population-frequency resolution. Decomposition does not imply that the variant is benign or that its original classification was incorrect; it indicates that the disease model implied by an unqualified P/LP label requires additional specification, such as recessive inheritance, reduced penetrance, or population-specific context, to remain coherent.

### MCAF Threshold Rationale

MCAF thresholds were treated as disease-model stress ceilings, not as pathogenicity cutoffs. The governing form was:

`MCAF = prevalence x allelic contribution / (penetrance x 2)`

The analysis therefore asks whether an observed population frequency is compatible with a specified disease model, not whether the variant is pathogenic. Three thresholds were used to bracket biologically interpretable dominant arrhythmia scenarios: 2.5 x 10^-5 as a strict high-penetrance dominant ceiling, 1 x 10^-4 as a reduced-penetrance or founder-compatible boundary, and 1 x 10^-3 as a deliberately permissive sensitivity extreme. These are not universal ClinVar thresholds and are not intended to be shared unchanged across all arrhythmia genes. They are stress-test ceilings for the unqualified dominant high-penetrance model; gene- or condition-specific deployment would require disease-specific prevalence, allelic contribution, and penetrance parameters. The 1 x 10^-5 exome screen was used only to nominate variants for review before MCAF bracketing, not as the final model-incompatibility definition.

This threshold strategy is internally checked in three ways. First, absolute alert counts change under threshold sensitivity, but the ordinal mechanism gradient is preserved within the usable-AF subset. Second, the GoF/dominant-negative class behaves as an internal allele-resolved control: at biologically interpretable MCAF thresholds, it shows 0/34 decomposition among evaluable variants. Third, hard conflict is separated from boundary or recessive-compatible routing rather than collapsed into a single "too frequent" label. These safeguards prevent MCAF from becoming an uncalibrated reclassification rule, but they do not remove the Tier 2 transportability boundary.

### Baseline Decision Model and VITAL Routing

To connect the frequency-stress analysis with downstream reuse, we formalized a simple label-driven baseline. This baseline is not a clinical recommendation; it is the null model used to audit flattened reuse.

The baseline is intentionally concrete rather than adversarial. It reflects common downstream workflows in which a ClinVar P/LP field is appended to a variant row and then used for triage before full model reconstruction: VEP/ANNOVAR-style annotation, laboratory candidate-variant spreadsheets, clinical decision-support triggers for actionable genes, and research filters that retain P/LP variants for review.

| Rule | Definition |
| --- | --- |
| B1_LABEL_DRIVEN_BASELINE | If a variant is exported as public ClinVar P/LP, downstream systems may treat it as eligible for direct actionability unless additional context is explicitly encoded. |
| B2_ACTIONABILITY_CONTEXT_BASELINE | If a P/LP variant appears in a gene linked to intervention, surveillance, cascade testing, device consideration, or therapy selection, the baseline route is ROUTE_PLP_ACTIONABLE. |
| B3_NO_IMPLIED_BENIGNITY | Absence of VITAL_OK does not imply benignity. It implies that direct actionability is not supported by the flattened label alone. |

VITAL was implemented as a post-label routing layer:

| VITAL route | Definition | Operational routing |
| --- | --- | --- |
| VITAL_OK | Label remains compatible with the tested actionability model. | Proceed with standard expert interpretation. |
| CHECK_POPMAX | Label may remain valid, but direct actionability requires ancestry-aware frequency review. | Review. |
| CHECK_MODEL | Label may remain valid, but inheritance, phase, or mechanism is not portable from the flattened label. | Model-specific routing or repair. |
| MODEL_CONFLICT | Flattened P/LP label is incompatible with the tested dominant high-penetrance disease model. | Do not direct-route as dominant actionability. |
| EVAL_LIMITED | No allele-resolved evaluability; direct actionability cannot be justified from population evidence. | Defer direct actionability pending representation, callability, or orthogonal evidence. |

Nonpass does not mean wrong. Nonpass means that direct actionability cannot proceed from the flattened label without review, deferral, or model repair.

### MVP Annotation Layer

To make the framework directly reusable in downstream annotation workflows, we implemented a minimal VITAL annotation layer. The MVP annotator accepts normalized VCF, CSV, or TSV inputs and maps variants against a precomputed lookup table derived from the arrhythmia analysis. Outputs include VITAL_evaluability, VITAL_flag, VITAL_regime, VITAL_popmax_af, VITAL_global_af, VITAL_threshold, and VITAL_reason.

For pipeline compatibility, the tool also generates an ANNOVAR-style export using Chr, Start, End, Ref, and Alt coordinates. All lookup-based interpretation requires prior variant normalization against GRCh38 using bcftools norm. The MVP intentionally restricts actionable flags to three categories: OK, CHECK_POPMAX, and MODEL_CONFLICT. VITAL_OK does not imply benign status, CHECK_POPMAX does not imply reclassification, and MODEL_CONFLICT denotes incompatibility only with the tested unqualified dominant high-penetrance model.

### Molecular Mechanism Classification

Genes were classified by molecular mechanism into three categories: gain-of-function/dominant-negative (GoF/DN), haploinsufficient (HI), and recessive (AR). These labels represent the dominant disease-model context being stress-tested for each gene, not a claim that every variant in a given gene acts exclusively through that mechanism. Gene-level assignments, their sources (ClinGen, OMIM, published literature), and caveats for genes with mixed or debated mechanisms, including SCN5A and KCNQ1, are detailed in supplementary tables.

For each category, decomposition rate was defined as the proportion of population-evaluable variants falling outside the dominant-compatible regime. Association between mechanism class and decomposition was assessed within the arrhythmia usable-AF subset by logistic regression with LOEUF score (gnomAD v4.1) and variant consequence as covariates. Because the GoF/DN mechanism class exhibits complete separation at biologically interpretable thresholds, Firth-type penalized logistic regression was used to obtain finite, bias-corrected estimates. Leave-one-gene-out sensitivity analysis was used to assess robustness to gene-level clustering.

### Validation and Calibration Layers

The validation architecture was designed as calibration rather than as a single oversized truth set. Three layers were used in the main text.

First, high-review ClinVar assertions tested whether routing instability persisted after excluding the weakest public assertions. Second, benign/likely benign controls tested whether VITAL behaved as a disease-model compatibility layer rather than as a pathogenicity detector. Third, published case-report reconstructions and repository-derived routing examples tested whether routing categories mapped onto recognizable clinical interpretation failure modes.

This validation suite is intentionally modest. It is not a blinded clinician adjudication study. Its purpose is to check calibration, specificity of hard-conflict behavior, clinical interpretability, and workflow-level decision change.

### Clinical-Action Context and Actionability Discordance Audit

To connect population tension with real downstream use environments, alerted variants were classified into potentially consequential interpretation contexts: cascade testing, drug restriction, intensive surveillance, device-related management, syndrome diagnosis, or carrier/recessive routing. These context labels indicate exposure environments, not measured downstream outcomes.

We also built an Actionability Discordance Audit (ADS) from repository-derived routing calls and representative case vignettes. Inclusion required a public P/LP assertion, action-associated context, normalized variant identity, and a VITAL route of CHECK_POPMAX, CHECK_MODEL, MODEL_CONFLICT, or EVAL_LIMITED. The ADS is not a harm registry; it is a variant-level audit of cases where direct actionability is not fully supported once constraints are restored.

### Analysis Hierarchy

The main text uses one primary denominator chain:

`1,731 public arrhythmia P/LP labels -> 357 allele-resolved exome matches -> 334 usable-AF variants -> 115 ancestry-aware frequency alerts`

All primary disease-model conclusions come from the 334 usable-AF variants and should be read as conditional on allele-level evaluability. All primary actionability-routing conclusions use the 1,731-label arrhythmia baseline. Secondary stress tests are archived separately and do not contribute to the primary burden estimates.

| Layer | Main denominator | Role |
| --- | ---: | --- |
| Public arrhythmia label universe | 1,731 | Primary routing denominator |
| Allele-resolved exome recovery | 357 | Evaluability boundary |
| Usable-AF exome subset | 334 | Primary population/MCAF denominator |
| Popmax alert set | 115 | Primary regime and context analysis |

## Results

### 1. A Structural Evaluability Boundary Limits Allele-Resolved Population Review

Across 1,731 unique arrhythmia P/LP variants, only 357 (20.6%) reached exact or equivalent allele-resolved population context in gnomAD v4.1.1 exomes. Of these, 334 retained usable allele-frequency data for downstream screening. The remaining 1,326 variants occupied Tier 2 locus- or region-context space, and 48 remained genuinely unevaluable after all reconciliation steps.

Strict initial matching recovered 350 variants. Trim-aware and decomposition-aware reconciliation increased recovery by only 7 variants. Full reference-based normalization using bcftools norm changed 0 of 1,731 representations. Together, these negative findings matter more than a cosmetic recovery gain: the dominant barrier is not an avoidable formatting miss, but the structural inability to connect most public arrhythmia assertions to allele-resolved population frequency.

![Allele-level non-observation is shaped by variant class](figures/arrhythmia_vital_absence_not_rarity.png)

### 2. Tier 2 Reflects Structured Representation Limits, Not Random Absence

Tier 2 is not a homogeneous gray zone. Within its 1,326 variants, 638 (48.1%) showed allele discordance at the queried locus, whereas 688 (51.9%) lacked even a same-locus gnomAD record. Tier 2 was also enriched for representation-sensitive classes: indels, duplications, and insertions accounted for 640 of 1,326 Tier 2 variants (48.3%), compared with 106 of 357 Tier 1 variants (29.7%), corresponding to an odds ratio of 2.21. Duplications showed the same pattern (170/1,326 Tier 2 variants, 12.8%, versus 22/357 Tier 1 variants, 6.2%; odds ratio 2.24).

The mechanism denominator was affected by the same evaluability structure. Only 34/248 GoF/dominant-negative labels (13.7%) had usable AF, compared with 92/564 haploinsufficient labels (16.3%), 70/116 recessive labels (60.3%), and 130/779 mixed-mechanism labels (16.7%). The evaluable GoF/dominant-negative subset was entirely SNV-based and came from RYR2, KCNJ2, and CALM2; non-evaluable CALM1, CALM3, and most RYR2 labels therefore cannot be assigned the same 0% decomposition estimate.

| Mechanism class | Total P/LP labels | Usable-AF Tier 1 labels | Not usable for MCAF | Transportability implication |
| --- | ---: | ---: | ---: | --- |
| Gain-of-function / dominant-negative | 248 | 34 (13.7%) | 214 (86.3%) | 0/34 decomposition is an allele-resolved SNV-space result only |
| Haploinsufficient | 564 | 92 (16.3%) | 472 (83.7%) | Frequency-regime estimates do not cover most labels |
| Recessive | 116 | 70 (60.3%) | 46 (39.7%) | Recessive labels are more evaluable and therefore not directly comparable to GoF/DN |
| Mixed or context-dependent | 779 | 130 (16.7%) | 649 (83.3%) | Model-specific inference requires phenotype, inheritance, and allele-level context |

This means that non-observation at the allele level is not a uniform proxy for rarity, and that mechanism-stratified decomposition is vulnerable to verifiability bias if generalized outside the usable-AF subset. For precisely the classes most vulnerable to representation instability, absence of an allele-resolved match often reflects discordance between public assertion and aggregation systems rather than genuine population absence. This argument applies most strongly to indels, duplications, and splice-disrupting variants, where representation instability can cause genuine alleles to be absent from frequency databases for technical reasons. For SNVs in well-covered callable regions, allele absence in gnomAD may constitute meaningful rarity evidence with a computable upper bound. This study conservatively treats all non-recovered variants as not allele-level frequency-evaluable rather than absent from the population, and does not apply coverage-based negative-lookup inference.

### 3. Global AF Alone Suppresses Ancestry-Localized Tension in the Evaluable Exome Subset

Within the 334 Tier 1 variants with usable AF, 321 (96.1%) had global AF <= 1 x 10^-5, so a global-AF-only review would have flagged just 13 variants. Replacing global AF with the maximum of global AF and popmax AF increased the alert set to 115. Global-only review therefore missed 102 of 115 ancestry-aware alerts (88.7%).

This was not a trivial methodological refinement. Of the 115 alerted variants, 103 fell in genes linked to clinically consequential interpretation contexts. Sixty occurred in drug-restriction contexts, 59 in intensive-surveillance or device-related contexts, and 103 in cascade-testing contexts. Global AF alone therefore suppresses precisely the ancestry-localized signal most relevant to downstream disease-model constraint in this domain.

![Population-specific frequency signals missed by global AF](figures/arrhythmia_population_af_outliers.png)

### 4. Exome-Resolved Frequency Tension Spans Three Interpretation Regimes Beneath a Shared Public Label

The 115 ancestry-aware exome alerts did not represent one kind of problem. They resolved into three distinct interpretation regimes:

| Regime | N variants | Representative locus | Defining feature |
| --- | ---: | --- | --- |
| Hard dominant incompatibility | 1 | SCN5A VCV000440850 | AF incompatible with an unqualified dominant high-penetrance reading |
| Boundary / monitoring | 76 | KCNH2 VCV004535537 | Exceeds a strict ceiling but remains compatible with a lower-penetrance dominant interpretation |
| Recessive / carrier-compatible | 38 | TRDN VCV001325231 | Dominant reading implausible; recessive or carrier logic remains coherent |

This MCAF-based regime classification is a structural portability stress-test, not a variant reclassification framework. It does not adjudicate whether any individual variant is pathogenic or benign; it tests whether the dominant high-penetrance disease model implied by an unqualified P/LP label remains compatible with observed population frequencies. A hard-incompatible regime does not reclassify the variant as benign; it identifies a label that requires explicit disease-model specification, such as recessive inheritance, reduced penetrance, or population-specific context, to sustain a coherent interpretation.

SCN5A VCV000440850 remains the clearest incompatibility case in the exome-resolved subset. KCNH2 VCV004535537 illustrates a parameter-sensitive monitoring regime rather than a categorical collapse, and TRDN VCV001325231 shows how a generic P/LP label can flatten a carrier-compatible recessive allele into an apparently dominant claim unless the disease model is made explicit.

### 5. VITAL Converts Frequency and Evaluability Stress Into Actionability Routing

Under the formal label-driven baseline, all 1,731 arrhythmia P/LP variants entered ROUTE_PLP_ACTIONABLE. After VITAL routing, only 219/1,731 (12.7%) remained VITAL_OK. The remaining 1,512/1,731 (87.3%) were rerouted away from direct actionability.

The 87.3% nonpass value should not be read as an 87.3% evidence-conflict rate. Most rerouted decisions were evidence-unavailable deferrals: 1,397/1,731 (80.7%) were EVAL_LIMITED because the asserted allele could not be evaluated at allele level in population data. A smaller subset carried evidence-based tension: 77/1,731 (4.4%) required ancestry-aware frequency review and 38/1,731 (2.2%) required model repair. Thus, VITAL separates absence of evaluability from positive evidence against a flattened model. Across the six-domain meta-analysis, the mean nonpass rate was 87.7%, again interpreted as a portability/routing burden rather than a conflict burden.

![Label-to-actionability routing audit](figures/vital_decision_disruption.png)

The central result is therefore not that 87% of labels are false, nor that 87% conflict with population evidence. It is that 87% of label-driven actionability decisions require rerouting once minimal constraints are restored, with most rerouting caused by evidence unavailability and a smaller component caused by frequency or model tension.

### 6. Routing Instability Persists Among High-Review Assertions

Routing instability persisted after restricting to high-review ClinVar assertions. In the review-score >= 2 subset, 309/365 variants (84.7%) left the direct-actionable baseline. Within that high-review subset, 256/365 (70.1%) were evaluation-limited, 33/365 (9.0%) required population review, and 20/365 (5.5%) required model-specific rerouting.

This pattern shows that the effect is not driven by low-confidence submissions alone. The problem is not simply that public databases contain noisy labels. The problem is that even strong label objects are not automatically portable actionability objects.

### 7. In the Evaluable Subset, Molecular Mechanism Predicts Disease-Model Stability

Decomposition rates varied systematically and predictably across gene classes defined by molecular mechanism within the usable-AF Tier 1 subset. These mechanism labels define the disease-model context being stress-tested for each gene, not the intrinsic mechanism of every variant in that gene.

| Mechanism class | Decomposition rate | Representative genes | Interpretation |
| --- | --- | --- | --- |
| Gain-of-function / dominant-negative | 0% | RYR2, CALM1-CALM3, KCNJ2 | Internal allele-resolved control: evaluable SNVs remain below the tested MCAF ceiling |
| Haploinsufficient | approximately 29% | SCN5A HI context, ANK2 | Moderate selection permits some enrichment in specific ancestry groups |
| Recessive | 37-54% | TRDN, CASQ2, KCNE1 | Weak heterozygous selection allows carrier frequencies inconsistent with dominant models |

The GoF/dominant-negative class at 0% decomposition at MCAF >= 2.5 x 10^-5 constitutes an internal allele-resolved control supporting the interpretation that decomposition at biologically interpretable thresholds reflects genuine biology rather than a uniform pipeline artifact (0/34 evaluable variants; 95% CI 0-10.2%; rule-of-three upper bound 8.8%). It does not establish 0% decomposition for all GoF/dominant-negative labels in ClinVar, because 214/248 GoF/DN labels lacked usable allele-level AF. At sub-biological thresholds below 2.5 x 10^-5, GoF/DN variants are also flagged, confirming that threshold choice determines whether the pipeline detects genuine disease-model incompatibility or normal population variation.

Variant consequence interacts with inheritance mode in the same direction. Missense variants in autosomal dominant genes decomposed at 6.7%, compared with 26.2% for loss-of-function variants in the same genes, a four-fold difference consistent with mechanism-dependent purifying selection. Loss-of-function in autosomal recessive genes reached 34.0%, and splice-disrupting variants reached 54.5%.

Logistic regression confirmed that LOEUF score and recessive or haploinsufficient mechanism class increased decomposition odds, while missense consequence was associated with lower decomposition odds relative to loss-of-function (OR = 0.83). Because the GoF/DN mechanism class exhibits complete separation, Firth-type penalized logistic regression was used to obtain finite estimates.

Within the usable-AF subset, leave-one-gene-out sensitivity analysis preserved the mechanism gradient in every iteration. The overall decomposition rate ranged from 17.2% to 21.3% (overall 19.5%; maximum single-gene shift 2.3 percentage points, driven by TRDN exclusion). The ordering GoF/DN < HI < AR was preserved without exception. Excluding limited-evidence genes (AKAP9, CACNA2D1, CACNB2, SNTA1), mixed-mechanism genes (SCN5A, KCNQ1, CACNA1C, KCNE1), or all eight disputed/mixed genes still preserved the gradient.

ClinVar review status analysis also showed that decomposition is not driven by low-confidence submissions. Among Tier 1 variants at MCAF = 2.5 x 10^-5, no-assertion-criteria variants decomposed at 3/12 (25.0%), single-submitter variants at 28/213 (13.1%), and multiple-submitters-no-conflicts variants at 34/109 (31.2%). The higher rate in the multiple-submitter subset is consistent with enrichment of well-characterized, population-common alleles that attract multiple submissions precisely because they are frequent enough to be independently observed. Restricting to Pathogenic (not Likely pathogenic) classification yielded 21/148 (14.2%), lower than the overall rate but preserving the mechanism gradient. Restricting to multiple-submitter (>=2-star) variants yielded 34/109 (31.2%), higher than overall, with the mechanism gradient steepened.

### 8. Calibration Uses Review Status, Controls, and Case Reconstruction

No single validation layer is large enough to carry the manuscript alone, so validation is presented as a calibration suite.

| Validation layer | Question | Result | Interpretation |
| --- | --- | --- | --- |
| High-review ClinVar assertions | Is routing instability driven by weak submissions? | 309/365 high-review assertions (84.7%) left the direct-actionable baseline. | No; high review status does not restore portability. |
| B/LB controls | Does VITAL behave as a pathogenicity classifier? | Common B/LB controls were incompatible with a forced dominant high-penetrance model. | VITAL detects model incompatibility, not pathogenicity status. |
| Published case reconstructions | Do routes map onto real clinical interpretation modes? | Six portability modes were reconstructed from arrhythmia case reports. | VITAL routes are clinically interpretable, not only abstract states. |

This validation architecture is sufficient for the manuscript's main claim because the claim is routing, not final clinical adjudication. A larger blinded external validation set would be the natural next benchmark, but the current data already show that VITAL separates direct actionability from review, model repair, and evaluability-limited deferral.

### 9. Concrete Variant Examples Reconstruct Label Portability

To connect the routing framework to real clinical interpretation, we reconstructed published arrhythmia case reports and repository-derived routing calls as variant-model-action objects rather than as patient-level outcome claims. Each reconstruction mapped the observed genotype to its current public label, the disease model asserted or implied in the source context, the clinical action context, and the resulting VITAL route.

Each case was processed through the VITAL routing logic: (1) assignment of the label-driven baseline, (2) allele-level evaluability check, (3) ancestry-aware frequency constraint, (4) disease-model compatibility assessment, and (5) final route assignment.

These reconstructions do not evaluate outcomes; they evaluate whether causal interpretation remains coherent once constraints are restored.

| Example | Source type | Label-driven baseline | VITAL route | Portability lesson |
| --- | --- | --- | --- | --- |
| SCN5A p.Asp1790Gly | Published case | Direct AD arrhythmia actionability | VITAL_OK + model specification | Preserved actionability requires explicit LQT3/Brugada overlap and medication-sensitive context |
| KCNH2 p.Ala614Val | Published case | Direct AD LQT2 actionability | VITAL_OK | Dominant actionability is preserved when allele, model, and management context remain aligned |
| KCNQ1 p.Thr322Met | Published case | Direct AD KCNQ1 actionability | CHECK_MODEL | Homozygous JLNS, heterozygous Romano-Ward-like disease, borderline QTc, and unaffected carrier states are not interchangeable |
| SCN5A p.Arg1193Gln | Historical case-association control | No current P/LP baseline; historical disease claim | NO_CURRENT_PLP_BASELINE; historical claim -> MODEL_CONFLICT | Ancestry-aware population frequency breaks an unqualified dominant high-penetrance historical claim |
| KCNH2 p.Arg356Cys | Published case | No P/LP baseline | NO_PLP_BASELINE / CHECK_MODEL | Drug-induced QT prolongation is clinically meaningful but does not make a VUS portable hereditary actionability |
| KCNQ1 large duplication | Published case | Sequencing-negative does not exclude actionability | EVAL_LIMITED / SV_BLOCKED | CNV interpretation requires copy-number-aware representation rather than SNV/indel AF logic |
| RYR2 c.1006-44_1007delinsATTTTG | Published case | Direct RYR2 CPVT actionability | EVAL_LIMITED + functional support | RNA and segregation evidence can repair interpretability when population evaluability is limited |
| SCN5A VCV000440850 | Repository routing call | ROUTE_PLP_ACTIONABLE | MODEL_CONFLICT | Direct dominant actionability is not supported under current population-frequency constraints |
| KCNH2 VCV004535537 | Repository routing call | ROUTE_PLP_ACTIONABLE | CHECK_POPMAX | Actionability requires ancestry-aware frequency review, not global AF alone |
| TRDN VCV001325231 | Repository routing call | ROUTE_PLP_ACTIONABLE | CHECK_MODEL | Recessive or carrier-state logic replaces flattened dominant actionability |
| KCNH2 VCV000405355 | Repository routing call | ROUTE_PLP_ACTIONABLE | EVAL_LIMITED | Direct actionability cannot be justified from allele-resolved population evidence |

The examples separate six portability modes: preservation, model compression, population conflict, event-genotype dissociation, evaluability failure, and functional rescue. The point is not that VITAL assigns a binary truth value to each report. Rather, it reconstructs whether a public label can travel from "variant observed in a case" to "actionable causal object" without losing the disease model, inheritance state, allele-frequency constraint, or evidentiary basis that made the interpretation coherent in the first place.

### 10. Repair Logic Separates Pathogenicity, Evaluability, and Actionability

The repair layer converts a nonpass route into an explicit next step rather than a vague warning.

| VITAL route | Meaning | Recommended next step |
| --- | --- | --- |
| VITAL_OK | Compatible with tested actionability model. | Proceed with standard expert interpretation. |
| CHECK_POPMAX | Ancestry/frequency tension. | Ancestry-aware review and penetrance/model reassessment. |
| CHECK_MODEL | Inheritance or mechanism not portable. | Specify inheritance, phase, and mechanism before action. |
| MODEL_CONFLICT | Tested dominant model incompatible. | Do not use flattened P/LP as a dominant-actionability claim. |
| EVAL_LIMITED | Not allele-resolved. | Defer direct actionability; require representation/callability/orthogonal evidence. |

This repair logic is the operational point of VITAL: it separates pathogenicity, evaluability, and actionability instead of letting a flattened public label stand in for all three.

## Discussion

### A Public P/LP Label Is a Compressed Object, Not a Self-Contained Disease Claim

The main result of this study is a measurable downstream reuse failure mode. Even within the exome Tier 1 usable-AF subset, one exported P/LP label can map to incompatible operational readings once ancestry-aware frequency and mechanism are restored.

This instability arises for three structural reasons. Public labels omit the disease-model parameters needed for reuse; population constraint is ancestry-resolved rather than scalar; and the response to frequency stress is mechanism-dependent. In this framework, decomposition means only incompatibility with the tested dominant high-penetrance model under observed frequency stress, not that the variant is benign, false, or clinically reclassified.

### The Evaluability Boundary Is Structural, Not a Pipeline Artifact

The exome reconciliation result is especially informative because of what it did not show. Reference-based normalization changed 0 of 1,731 variants. Exhaustive reconciliation rescued only 7 additional exact or equivalent exome matches. This localizes the failure clearly: the main barrier is not missed left-alignment or sloppy trimming, but a structural disconnect between public clinical assertions and population aggregation systems.

This distinction matters for governance. If the evaluability gap were a normalization artifact, it could be solved by better preprocessing. Because it is structural, it requires changes at the level of public assertion metadata: variants must carry explicit evaluability tier information so that downstream users know whether population constraint can be applied at the allele level.

### Population Structure Is Not a Refinement Layer; It Is Part of the Claim

Global-only screening missed 88.7% of ancestry-aware exome alerts. In a domain shaped by ancestry-localized enrichment, founder effects, and mixed inheritance architectures, global AF is not a sufficient stand-in for the population signal needed to constrain a disease model. Popmax does not add decorative granularity; it restores biologically relevant structure that global AF suppresses. Among the 102 missed alerts, the dominant signal came from non-European populations, precisely the populations most underrepresented in historical variant databases and most exposed to systematic under-detection of frequency incompatibility.

### In the Evaluable Subset, Decomposition Is Mechanism-Structured

The molecular mechanism findings sharpen the interpretation of the regime structure within the allele-resolved subset. GoF and dominant-negative genes show 0% decomposition at MCAF thresholds >= 2.5 x 10^-5 among 34 usable-AF variants, not as a trivial result but as an internal control that supports the framework at biologically interpretable thresholds. This does not exclude higher or lower decomposition in the non-evaluable GoF/DN majority. Instead, it means that where allele-level population evidence is available, decomposition is absent where selection theory predicts it should be absent and present in proportion to how much constraint is relaxed in haploinsufficient and recessive genes.

### Arrhythmia Genes Act Here as a Stress-Test Domain, Not as a Special Case

This manuscript does not claim that arrhythmia genes are uniquely affected. Rather, they provide a stringent stress-test domain in which ancestry structure, variable penetrance, founder effects, and potentially consequential interpretation contexts coexist. The core claim remains that public P/LP labels lose disease-model portability unless evaluability, inheritance, mechanism, and ancestry-aware AF are made explicit.

### VITAL Is a Routing Layer, Not a Reclassification Engine

VITAL does not claim that nonpass variants are benign, incorrect, or clinically irrelevant. It claims that direct actionability cannot be inferred from the flattened label alone. A public pathogenicity label can remain true in its expert context while still being insufficient as a direct-actionability object in a flattened downstream workflow.

## Limitations

First, all exome disease-model conclusions are restricted to the 334 Tier 1 variants with usable allele-resolved AF. The unevaluable majority cannot be assumed to follow the same regime structure. Because Tier 1 was enriched for SNVs and depleted of indels and duplications relative to Tier 2, direct extrapolation across variant classes would be inappropriate. The 0/34 GoF/DN decomposition result is therefore a statement about evaluable GoF/DN SNVs, not about all GoF/DN assertions.

Second, this is a study of disease-model compatibility and actionability routing, not patient-level outcomes. The clinical-action categories used here identify potentially consequential interpretation contexts, not realized downstream clinical consequences.

Third, frequency-based analysis can exclude a particular disease-model reading, but it cannot on its own establish the correct alternative interpretation. Low penetrance, haplotypic context, allelic phase, transcript-specific rescue, segregation evidence, and functional validation remain outside the scope of this analysis unless orthogonal data are introduced.

Fourth, the 1 x 10^-5 trigger used in the exome screening layer is operational rather than ontological. The qualitative result, especially the large popmax gain over global AF, was stable across a broad threshold range, but absolute alert counts depend on the screening threshold selected. The MCAF ceilings are disease-model stress parameters, not universal pathogenicity cutoffs; clinical deployment should recalibrate them by gene, disease prevalence, penetrance, allelic contribution, and ascertainment context.

Fifth, the current structural-variant layer is a blocking requirement rather than a fully deployed SV/CNV engine. CNVs, repeat-mediated loci, paralogous genes, haplotypes, and phase-dependent mechanisms require specialized representation before direct portability can be granted.

## Conclusion

Flattened public P/LP labels lose disease-model portability in downstream actionability workflows unless allele evaluability, ancestry-aware AF, and disease model are restored. Across 1,731 ClinVar P/LP variants in 20 inherited arrhythmia genes, only 219 labels (12.7%) could be directly routed as actionability-compatible after minimal portability constraints. Most non-direct routes reflected missing allele-resolved evaluability rather than evidence conflict.

Within the 334 usable-AF variants, popmax recovered ancestry-localized frequency tension missed by global AF, and MCAF stress-testing separated compatible, review-level, and model-conflict regimes. The mechanism gradient was conditional on allele-level evaluability: evaluable GoF/DN SNVs showed 0/34 decomposition, whereas haploinsufficient and recessive architectures showed higher incompatibility under flattened dominant models.

The main claim is therefore not that public P/LP labels are wrong. It is that a flattened label is not a self-contained actionability object. Direct downstream actionability requires restoration of allele-level evaluability, ancestry-aware frequency, and disease-model context.

## Tool Availability

The VITAL pipeline (Variant Interpretation Through Ancestry-aware Labeling) is available as an open-source computational tool in the project repository. In addition to generating the analytical outputs reported here, VITAL includes a minimal annotation layer for downstream reuse. The MVP annotator accepts VCF, CSV, or TSV inputs, supports lookup-based annotation, emits VITAL evaluability and disease-model flags, and can export an ANNOVAR-style table with Chr, Start, End, Ref, Alt, VITAL_evaluability, VITAL_flag, VITAL_regime, VITAL_popmax_af, VITAL_global_af, VITAL_threshold, and VITAL_reason.

The tool is designed as a guardrail layer rather than a variant reclassification system: VITAL_OK does not mean benign, CHECK_POPMAX does not mean reclassification, and MODEL_CONFLICT indicates incompatibility with the tested unqualified dominant high-penetrance model. The repository also exposes the routing audit layer (python src/run_vital_routing_validation.py), the cached cohort-level CLI (python run_vital.py --mode full --genes "MYBPC3,MYH7" --pop gnomAD), and machine-readable repair/reason-code tables.

## Data Availability

Code, cached intermediate files, processed tables, and figure assets required to reproduce the analyses summarized here are available in the repository. The April 24, 2026 data freeze used ClinVar public assertion snapshots together with gnomAD v4.1.1 exome outputs cached in the project workspace. Machine-readable supplementary tables include evaluability tier classifications, frequency flags, regime assignments, molecular mechanism analysis outputs, routing validation outputs, and secondary stress-test outputs. Repository: https://github.com/anerecye/New-project
