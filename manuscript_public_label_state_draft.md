# A Public Pathogenic Label Is Not a Disease State: Population Frequency Reveals Hidden Interpretation Regimes in ClinVar Arrhythmia Assertions

## Abstract

**Purpose:** Public ClinVar pathogenic/likely pathogenic (P/LP) assertions are routinely reused as though they encode a single, high-penetrance Mendelian disease state. We tested whether that reading remains coherent for inherited arrhythmia variants once ancestry-aware population frequency, allele-resolved variant representation, allele count support, assertion provenance, and mechanism-specific interpretation are considered together.

**Methods:** We collapsed ClinVar P/LP records across 20 inherited arrhythmia genes to 1,731 unique GRCh38 variants and cross-referenced them against gnomAD v4.1.1 exomes. We first performed strict exact matching and then a second trim- and decomposition-aware reconciliation pass to recover equivalent allele representations and eliminate query failures. Variants were assigned to one of three interpretability tiers: allele-resolved exact/equivalent population context, locus/regional context without allele-resolved AF, or genuinely unevaluable status after both passes. For the canonical exome-resolved subset used for direct ancestry-aware comparison, we retained global allele frequency (AF), popmax AF, gnomAD allele count (`AC_count`), variant class, ClinVar review strength, and submitter structure. We compared global-only versus ancestry-aware screening and applied a frequency plausibility framework based on maximum credible allele frequency logic to test whether observed population frequency remained compatible with high-penetrance dominant, lower-penetrance dominant, or recessive carrier-compatible interpretations. We also performed provenance-aware sensitivity analyses using a pre-specified credibility filter that excluded weak/no-assertion or single-submitter P/LP assertions lacking expert review when they also showed extreme population-frequency inconsistency, defined here as `max_frequency_signal > 1x10^-4` with `qualifying_frequency_ac >= 20`; no separate internal-conflict flag survived into this final non-conflicting P/LP cohort.

**Results:** After reconciliation, allele-resolved exact/equivalent recovery increased only marginally from 350 to 357 variants. Of 1,731 public arrhythmia P/LP assertions, 357 (20.6%) achieved allele-resolved exact/equivalent population context, 1,326 (76.6%) retained only locus/regional context without allele-resolved AF, and 48 (2.8%) remained genuinely unevaluable; all API/query-error states were eliminated. This indicates that the principal limitation is not cosmetic matching failure but the operational unavailability of allele-resolved population frequency for most public P/LP variants under current representation regimes. Within the canonical exome-resolved subset, global AF alone identified 13 alerts whereas global-or-popmax screening identified 115, so global-only review missed 102/115 ancestry-aware contradictions; the 7 recovered equivalent alleles were too few and too rare to change this contrast materially, and only one additional genome-only exact record crossed the global AF threshold. Frequency plausibility analysis exposed at least three interpretation regimes hidden beneath the same public P/LP label. SCN5A `VCV000440850` exceeded a strict dominant high-penetrance AF ceiling by 113.5-fold, a Brugada-like 20% penetrance ceiling by 45.4-fold, and even a deliberately permissive low-penetrance dominant ceiling by 5.7-fold, indicating hard incompatibility with an unqualified dominant high-penetrance Mendelian reading of the exported label. KCNH2 `VCV004535537` exceeded a strict high-penetrance ceiling by 2.6-fold but remained compatible with a more permissive low-penetrance model, placing it in a parameter-sensitive boundary regime. TRDN `VCV001325231` was incompatible with a dominant reading but remained population-plausible under recessive carrier logic, with carrier frequency approximately `4.36x10^-4` and expected homozygote frequency approximately `4.76x10^-8`. Provenance-aware sensitivity analyses excluding the 3 low-credibility extreme-frequency records left the broader structure intact: exact/equivalent reconciliation changed only from `357/1731` to `354/1728`, locus/regional context remained `1,326/1,728`, genuine unevaluable status remained `48/1,728`, and the ancestry-aware alert excess persisted in the remaining AF-observed core (`10` global-only versus `112` global-or-popmax alerts), although the 3 illustrative red-priority cases were removed by construction.

**Conclusion:** The central problem is not merely that some P/LP variants are too frequent. It is that a single public P/LP label can fail to preserve the biological state required to constrain downstream interpretation, while allele-resolved population frequency remains operationally unavailable for most variants even after reconciliation. In inherited arrhythmia genes, public P/LP assertions are therefore both biologically under-specified and operationally under-constrained. Population frequency and mechanism together reveal at least three biologically distinct regimes hidden beneath the same exported label: dominant high-penetrance incompatibility, parameter-sensitive low-penetrance boundary cases, and carrier-compatible recessive architecture.

## Introduction

ClinVar is a foundational public resource for clinical variant interpretation. Its pathogenic and likely pathogenic assertions are reused by diagnostic laboratories, researchers, curation efforts, and downstream reviewers as compact signals of disease relevance. In practice, however, a public P/LP label is often consumed as though it encodes a single, uniform biological claim: that the variant supports a coherent, relatively high-confidence Mendelian disease state.

That assumption is not trivial, and inherited arrhythmia genes are an unusually stringent setting in which to test it. Disorders such as long QT syndrome, Brugada syndrome, catecholaminergic polymorphic ventricular tachycardia, and related channelopathies sit at the intersection of rare disease interpretation, incomplete penetrance, ancestry-specific enrichment, founder effects, susceptibility architecture, recessive disease mechanisms, and complex haplotypic backgrounds. The clinical stakes are equally high. Variant interpretation in this setting can influence diagnostic closure, cascade testing, medication guidance, surveillance, and device-related decisions.

Most discussions of public variant databases focus on whether individual assertions are correct. That is not the primary question here. We ask a narrower and, in some ways, more consequential one: what biological state does a public P/LP label actually preserve once it leaves its original submission context? A label may be reused as shorthand for dominant high-penetrance disease even when the underlying evidence is more consistent with lower penetrance, susceptibility, recessive carrier architecture, haplotypic context, or representation-sensitive uncertainty.

Population frequency is one of the strongest constraints on such over-reading. Variant interpretation frameworks already recognize that alleles too common for a disorder should contribute evidence against a high-penetrance Mendelian interpretation. Yet frequency review fails in at least two distinct ways. First, global AF can conceal ancestry-specific enrichment that becomes visible only under population-specific review. Second, allele-resolved population frequency is often operationally unavailable even when the locus is not empty of information, because exact observation depends on representation, normalization, variant class, calling sensitivity, sequence context, and database aggregation conventions. This matters particularly for indels, duplications, splice-disrupting variants, and complex alleles.

A second challenge concerns assertion provenance. Public P/LP records are not homogeneous objects. They differ in review strength, submitter structure, and internal consistency. Some high-frequency P/LP assertions may reflect lower-credibility single-submitter interpretations rather than a stable cross-laboratory pathogenic claim. Any attempt to use population frequency as a constraint on exported labels must therefore distinguish between biological under-specification and provenance-related fragility.

We therefore tested whether public ClinVar P/LP arrhythmia assertions remain interpretable as a single high-penetrance Mendelian state once allele-resolved matching, trim- and decomposition-aware reconciliation, ancestry-aware frequency, assertion provenance, allele count support, and mechanism-specific disease models are considered together. Our claim is deliberately specific. We do not argue that ClinVar is broadly unreliable, nor do we introduce a new pathogenicity classifier. Instead, we show that, in inherited arrhythmia genes, the same exported public P/LP label can hide biologically distinct and interpretation-relevant states that cannot be safely treated as one thing. We also show that even after aggressive reconciliation, allele-resolved population frequency remains unavailable for most variants. The consequence is structural: the public P/LP label is not only biologically under-specified, but is often attached to an object for which one of the field's key constraining signals cannot be cleanly applied.

## Methods

### Study design and cohort definition

We assembled ClinVar pathogenic and likely pathogenic records across 20 canonical inherited arrhythmia genes: `KCNQ1`, `KCNH2`, `SCN5A`, `KCNE1`, `KCNE2`, `RYR2`, `CASQ2`, `TRDN`, `CALM1`, `CALM2`, `CALM3`, `ANK2`, `SCN4B`, `KCNJ2`, `HCN4`, `CACNA1C`, `CACNB2`, `CACNA2D1`, `AKAP9`, and `SNTA1`. VCV XML records were retrieved through the NCBI Entrez API and collapsed to unique GRCh38 variant-level entries using chromosome, position, reference allele, and alternate allele. After initial de-duplication and quality control, 1,731 unique P/LP variants remained. The primary data freeze was April 21, 2026.

### gnomAD reconciliation and interpretability tiers

Variants were queried against gnomAD v4.1.1 exomes through the gnomAD GraphQL API. We used a two-stage reconciliation strategy. First, we performed strict exact matching requiring concordance of chromosome, position, reference allele, and alternate allele. Second, for unresolved variants, we performed a slower trim- and decomposition-aware reconciliation pass to recover equivalent allele representations and eliminate technical query failures.

After both passes, each queried variant was assigned to one of three interpretability tiers. Tier 1 comprised variants with allele-resolved exact or equivalent population context. Tier 2 comprised variants with locus/regional context but without allele-resolved AF, including same-site allele discordance or other non-equivalent local population context insufficient for direct allele-frequency interpretation. Tier 3 comprised variants that remained genuinely unevaluable after both passes. Missing AF was never converted to zero, and Tier 2 variants were not treated as population-consistent by implication.

Because the current runtime did not provide full `bcftools norm -f GRCh38.fa -m -both` left-alignment against a local reference FASTA, the reconciliation layer was implemented conservatively using trim-normalization of shared prefix/suffix, local `+/-2 bp` indel-payload equivalence, unphased decomposition-aware matching for short multi-base substitutions, and a second slow pass to eliminate transient API/query failures. Tier 2 evidence was retained as context, not as allele-level AF.

### Assertion provenance and credibility filtering

For each ClinVar record, we retained `review_status`, derived `review_strength`, numeric `review_score`, `submitter_count`, and `submitter_count_source`. The arrhythmia XML/API cohort was already restricted to retained P/LP assertions without surviving conflicting-classification labels, so no separate internal-conflict field remained available for the main cohort.

Because public P/LP records differ in provenance quality, we analyzed the full cohort as the primary dataset and then performed pre-specified sensitivity analyses using a credibility-filtered cohort. The credibility filter excluded weak/no-assertion or single-submitter P/LP assertions lacking expert review when they also showed extreme population-frequency inconsistency, defined operationally as `max_frequency_signal > 1x10^-4` with `qualifying_frequency_ac >= 20`. This excluded 3 records (`SCN5A VCV000440850`, `TRDN VCV001325231`, and `KCNH2 VCV004535537`) and was designed to test robustness to fragile public assertions rather than to redefine pathogenicity.

### Ancestry-aware frequency analysis

For the canonical exact exome AF-observed subset (`n=334`), we extracted global AF, global `AC_count`, population-specific AF, population-specific `AC_count`, and popmax AF. We compared two screening strategies. A naive strategy used global AF `> 1x10^-5`. An ancestry-aware strategy used global AF or popmax AF `> 1x10^-5`. Because this threshold was intended as a screening alert rather than a classification boundary, we also performed threshold sensitivity analyses across score cutoffs `40-85` and AC gates `10`, `20`, `50`, `100`, and `200`. `AC_count` was treated as support for the reliability of the observed population signal, not as a disease threshold.

Equivalent recoveries were tracked separately. The 5 local-window equivalent indels were queried directly by recovered gnomAD variant ID and all remained below the alert threshold (`0`, `7.97x10^-6`, `1.37x10^-6`, `6.85x10^-7`, and `6.84x10^-7`). Among the 16 exact genome/no-exome records, only one additional genome-only exact record crossed global AF `> 1x10^-5`. Thus the direct ancestry-aware headline contrast remained the canonical exome-based `13` versus `115`.

### Frequency plausibility framework

We used population frequency not as a benignity classifier but as a constraint on the disease model implied by the public label. For dominant interpretations, maximum credible allele frequency was approximated using standard rare-disease logic:

`AF_max ~= (P x AC_contrib) / (Pen x model_factor)`

where `P` is disease prevalence, `AC_contrib` is the maximum plausible proportion of cases attributable to the variant, `Pen` is penetrance, and `model_factor = 2` for a heterozygous diploid dominant model. We evaluated both stricter high-penetrance and deliberately more permissive low-penetrance dominant scenarios to avoid over-constraint. We then asked whether observed AF remained compatible with an unqualified dominant Mendelian reading across plausible parameter space.

For recessive interpretations, we used standard rare-allele approximations. For a rare recessive allele with frequency `q ~= AF`:

`Carrier frequency ~= 2q`

`Expected homozygote frequency ~= q^2`

This allowed explicit discrimination between dominant-model incompatibility and recessive carrier-compatible plausibility.

### AC_count support

We summarized `AC_count` support as a graded reliability feature rather than imposing a single hard threshold. In the principal illustrative cases, `SCN5A` had `AC_count = 214`, `TRDN` had `AC_count = 40`, and `KCNH2` had `AC_count = 24`. Approximate Poisson relative standard errors were `6.8%`, `15.8%`, and `20.4%`, respectively, with corresponding approximate 95% relative margins of `13.4%`, `31.0%`, and `40.0%`. This was used to distinguish strong contradiction from boundary-level frequency evidence.

### Interpretation routing

Once public label, ancestry-aware frequency, assertion provenance, disease mechanism, and interpretability tier were considered together, variants were conceptually routed into one of several review-relevant states: dominant high-penetrance disease, lower-penetrance or susceptibility signal, carrier-compatible recessive architecture, annotation-sensitive borderline case, representation-limited locus/regional context, or genuinely unevaluable state. This routing was not used for automatic reclassification. It was used to determine what kind of expert review a public label still required.

## Results

### Even after reconciliation, allele-resolved population-frequency evidence remained available for only a minority of public arrhythmia P/LP assertions

ClinVar parsing yielded 1,731 unique arrhythmia P/LP variants. Under strict first-pass matching, 350 variants achieved exact allele-level recovery. After trim- and decomposition-aware reconciliation, this increased only to 357. The marginal gain of seven variants is the key operational result: it shows that the main limitation is not superficial pipeline failure, but the inability to obtain allele-resolved population frequency for most public P/LP variants under current representation regimes.

Using the reconciliation-aware tier system, `357/1,731` variants (`20.6%`) achieved allele-resolved exact/equivalent population context, `1,326/1,731` (`76.6%`) retained only locus/regional context without allele-resolved AF, and `48/1,731` (`2.8%`) remained genuinely unevaluable after both passes. All API/query-error states were eliminated. Frequency-based disease-model analyses in this manuscript are therefore restricted to the canonical exome AF-observed subset and the few additional allele-resolved recoveries tracked separately, whereas Tier 2 variants are interpreted as representation-limited rather than population-consistent.

This was not a minor technical inconvenience. It defined the evaluability boundary of the study. Even after aggressive reconciliation, roughly four of five public arrhythmia P/LP assertions remained outside clean allele-resolved AF review. The problem is therefore structural, not cosmetic.

### Global AF alone missed most ancestry-aware frequency contradictions

Within the canonical exome AF-observed subset (`n=334`), global AF alone identified `13` variants whereas ancestry-aware global-or-popmax screening identified `115`; global-only review therefore missed `102/115` ancestry-aware alerts. The 7 recovered equivalent alleles were too few and too rare to change this contrast materially, and only one additional genome-only exact record crossed the global AF threshold.

This result does not indict population-frequency frameworks themselves. It indicts a common implementation shortcut. Frequency logic cannot constrain interpretation properly when ancestry-aware structure is ignored.

### Representation-limited non-observation was structured by variant class

The absence of allele-resolved AF evidence was not random. Variants outside the exact/equivalent tier were enriched for representation-sensitive classes, particularly indels and duplications. The combined indel/duplication fraction was `47.1%` outside exact/equivalent recovery versus `29.7%` inside it (`OR = 2.11`, `p = 1.49x10^-9`). This supports a separate operational warning: non-observation at the allele-resolved level cannot be treated as generic rarity, especially for technically difficult or representation-unstable variant classes.

The practical lesson is not that another matching pass automatically resolves the problem. It is that representation-aware caution is required. A missing SNV and a missing duplication do not mean the same thing, and a locus/regional signal is not interchangeable with an allele-resolved frequency observation.

### Assertion provenance contributed to frequency tension but did not fully explain it

We next asked whether frequency-discordant P/LP assertions were enriched among lower-credibility public records. Within the canonical AF-observed subset, ancestry-aware frequency discordance was present across provenance classes: `53/109` (`48.6%`) multiple-submitters-no-conflicts records, `57/213` (`26.8%`) single-submitter records, and `5/12` (`41.7%`) weak/no-assertion records were frequency-discordant. Thus, lower-review-strength submissions contributed to the problem, but frequency tension was not confined to that provenance class. All 3 red-priority cases sat inside the fragile provenance layer by design, but the broader ancestry-aware/global mismatch persisted well beyond it.

This point matters because it distinguishes two interpretations of the dataset. One possibility is that frequent P/LP signals are mainly trivial submitter noise. The other is that under-specification persists even among more stable public assertions. Our data support the latter at least in part.

### Severe annotation did not rescue coherence by itself

Among the canonical AF-observed variants, `92/262` loss-of-function/splice assertions (`35.1%`) still showed ancestry-aware frequency tension, similar to missense assertions (`19/66`, `28.8%`). Severe annotation therefore did not protect a public P/LP assertion from population-frequency conflict.

This is a narrower point than "severe variants are often benign," which we do not claim. The actual conclusion is sharper: consequence severity does not specify the disease state required for coherent interpretation. A severe variant can still be compatible with recessive carrier architecture, lower penetrance, founder enrichment, haplotypic context, or annotation-sensitive uncertainty.

### Frequency plausibility exposed three distinct interpretation regimes hidden beneath the same public label

The central result of the study was not the number of frequency-flagged variants. It was the structure of the conflicts that emerged once frequency was treated as a constraint on the implied disease model. Within the allele-resolved subset, the same public P/LP label concealed at least three qualitatively different interpretation regimes.

#### 1. Hard incompatibility with an unqualified dominant high-penetrance reading: SCN5A `VCV000440850`

SCN5A `VCV000440850`, `c.[3919C>T;694G>A]`, showed `AFR popmax AF = 5.68x10^-3` with `AC_count = 214`. Under a strict dominant high-penetrance ceiling, observed AF exceeded the maximum credible allele frequency by `113.5-fold`. Under a Brugada-like scenario with penetrance set to `20%`, the excess remained `45.4-fold`. Even under a deliberately permissive low-penetrance dominant model, observed AF still exceeded the ceiling by `5.7-fold`.

Across a broad and intentionally generous parameter space, no unqualified dominant high-penetrance interpretation remained population-plausible. This is therefore not merely a high-frequency outlier. It is a hard incompatibility between the observed population signal and an unqualified dominant high-penetrance Mendelian reading of the exported public P/LP assertion. This does not exclude alternative biological explanations, including low penetrance, haplotypic context, complex allelic architecture, or disease association detectable only with orthogonal functional or segregation data.

#### 2. Parameter-sensitive boundary conflict: KCNH2 `VCV004535537`

KCNH2 `VCV004535537`, `c.2398+2T>G`, had global `AC_count = 24` and `ASJ popmax AF = 1.32x10^-4`. Under a strict high-penetrance long-QT-like ceiling, observed AF exceeded the plausible maximum by `2.6-fold`. Under a more permissive low-penetrance dominant ceiling, however, the contradiction no longer persisted.

This places KCNH2 in a boundary regime rather than the same contradiction class as SCN5A. The case still matters because a generic public P/LP reading can imply stronger disease certainty than population plausibility supports. But the correct conclusion is not absolute incompatibility. It is that the variant occupies a parameter-sensitive, mechanism-sensitive zone that requires qualified expert review rather than automatic trust in splice severity or public label status.

#### 3. Dominant-model failure with recessive carrier-compatible plausibility: TRDN `VCV001325231`

TRDN `VCV001325231`, `c.1050del`, had global `AC_count = 40` and `AMR popmax AF = 2.18x10^-4`. Under a dominant Mendelian interpretation, the observed population signal was not plausible. Under recessive carrier logic, however, the same observation became coherent: carrier frequency was approximately `4.36x10^-4` and expected homozygote frequency was approximately `4.76x10^-8`.

TRDN therefore illustrates a distinct interpretation failure. The issue is not simply that the public label is too frequent. The issue is that the wrong model class is being imported. A generic P/LP label can flatten a carrier-compatible recessive allele into an apparently dominant disease-causing statement unless inheritance state is made explicit.

### AC_count graded confidence across these conflicts

The three principal cases also differed in signal stability. SCN5A had `AC_count = 214`, supporting a robust population signal. TRDN had `AC_count = 40`, consistent with moderate but still meaningful support. KCNH2 had `AC_count = 24`, reinforcing its treatment as a boundary case rather than a hard contradiction.

This gradient matters. Frequency tension is not monolithic. Some cases support strong incompatibility with an unqualified dominant high-penetrance reading, some support inheritance-routed reinterpretation, and others identify boundary conditions where the disease model itself becomes the critical variable.

### The overall structure of the signal persisted after provenance-aware filtering

We repeated the main analyses after applying the pre-specified credibility filter to weak/no-assertion or single-submitter P/LP assertions lacking expert review and showing extreme AF inconsistency (`max_frequency_signal > 1x10^-4` with `qualifying_frequency_ac >= 20`). This excluded exactly 3 records and removed the red-priority exemplars by construction. Even so, the higher-level architecture of the results persisted: exact/equivalent reconciliation changed only from `357/1,731` to `354/1,728`, locus/regional context remained `1,326/1,728`, genuine unevaluable status remained `48/1,728`, and the ancestry-aware alert gap persisted in the remaining AF-observed core (`10` global-only versus `112` global-or-popmax alerts).

This indicates that the core findings are not driven solely by a small number of obviously fragile public assertions. What disappears under the filter is the urgent exemplar set, which is precisely what the filter was designed to remove. What remains is the broader structure of under-specification and ancestry-aware frequency tension.

## Discussion

The main result of this study is not that ClinVar contains imperfect assertions. That observation is too shallow to be useful. The stronger result is that a public pathogenic label is not, by itself, a disease state. In inherited arrhythmia genes, the same public P/LP label can remain compatible with multiple biologically distinct interpretations once ancestry-aware population frequency, allele-resolved variant representation, mechanism, and assertion provenance are reintroduced.

Those interpretations are not interchangeable. In this study, they fell into at least three clinically and conceptually distinct regimes: hard incompatibility with an unqualified dominant high-penetrance disease reading, parameter-sensitive low-penetrance boundary cases, and recessive carrier-compatible architecture. That distinction matters because these regimes do not support the same downstream actions. They imply different expectations for penetrance, family testing, counseling, and interpretive caution. A public label that fails to preserve this state information is not merely incomplete. It is liable to be over-read.

The reconciliation result sharpens this argument in an important way. Before any biological interpretation begins, one must first be able to obtain allele-resolved population frequency for the variant under review. Here, a second trim- and decomposition-aware reconciliation pass increased exact/equivalent recovery only from 350 to 357 variants and reduced query-error states to zero. That is a strong negative finding. It shows that the principal limitation is not a fragile initial pipeline or an easily fixable formatting defect. Rather, allele-resolved population frequency is operationally unavailable for most public P/LP variants under the current representation and aggregation regime.

This reframes the role of population frequency. Frequency is not used here as a blunt benignity rule, nor as a shortcut for reclassification. It is used as a plausibility constraint on the disease model being imported into downstream interpretation. But that constraint can function only when allele-resolved population evidence is actually available. The same public P/LP label is therefore doubly compressed: it does not preserve a unique biological disease state, and for most variants it is not attached to an object for which a key constraining signal can be cleanly applied.

The ancestry-aware findings reinforce this point from the opposite direction. Among the variants for which allele-resolved AF was available, global AF alone missed most frequency contradictions visible under popmax review. This means that even when analysts invoke population logic in principle, they may still fail to constrain interpretation in practice if they ignore ancestry structure. The system therefore fails in two separable ways. For a minority of variants, allele-resolved AF exists but is easily misused if ancestry-aware structure is ignored. For the majority, allele-resolved AF is not cleanly available in the first place, even after reconciliation.

Assertion provenance adds a third layer. Some frequency tension is attributable to lower-credibility public records, especially single-submitter and weak/no-assertion claims. But provenance alone does not dissolve the problem. Even after pre-specified filtering of the 3 most extreme fragile records, the broader structure of under-specification persisted. This matters because it means the phenomenon is not reducible to trivial database noise.

The representation finding is especially important. Public databases and downstream users often treat non-observation as an argument for rarity. Our data show that this inference is unsafe in precisely the variant classes where technical representation is least stable. A locus/regional signal without allele-resolved equivalence is not the same as rarity, and it is not the same as allele-level absence. What appears to be a frequency signal may instead be a representation artifact.

The severe-annotation result makes the same point from another angle. A severe consequence term does not rescue an under-specified public label. Severe variants can still occupy lower-penetrance, recessive, founder-enriched, haplotypic, or otherwise state-dependent disease models. Annotation severity is not a substitute for explicit biological specification, just as locus-level context is not a substitute for allele-resolved AF.

Taken together, these findings suggest that public P/LP assertions should not be consumed as though they encode a single kind of pathogenic claim. At minimum, downstream review should ask two distinct questions. First, what disease state is actually being implied: dominant high-penetrance disease, lower-penetrance susceptibility, recessive carrier architecture, haplotypic/complex-allele context, or another state-dependent interpretation? Second, what level of population interpretability is actually available: allele-resolved AF, locus/regional context without allele resolution, or genuinely unevaluable status? Without both distinctions, the same public object can support mutually inconsistent interpretations while appearing deceptively simple.

This argument is deliberately narrower than a broad critique of public databases. We do not claim that ClinVar is unreliable overall, nor that public assertions should be discarded. The point is more structural. Public labels are compressed outputs. When that compression erases the disease state needed to constrain interpretation, and when allele-resolved population frequency remains unavailable for most variants even after reconciliation, the label remains useful only if mechanism and population logic are explicitly reintroduced downstream.

## Limitations

First, all frequency-based disease-model analyses are restricted to the allele-resolved subset. Most public arrhythmia P/LP assertions did not achieve clean allele-resolved AF even after trim- and decomposition-aware reconciliation and therefore remain representation-limited rather than population-consistent.

Second, the plausibility framework evaluates incompatibility with an unqualified disease-model reading of the public label, not final causal architecture. It does not distinguish low penetrance, haplotypic context, allelic phase, segregation-supported pathogenicity, transcript-specific rescue, or functionally validated effects in the absence of orthogonal data. Accordingly, the strongest frequency conflicts should be interpreted as model-level incompatibility with an unqualified high-penetrance dominant reading, not as definitive exclusion of all biologically relevant disease association.

Third, frequency ceilings depend on assumptions about prevalence, penetrance, and allelic contribution. To minimize the risk of over-constraint, we intentionally evaluated permissive parameter ranges. The resulting analyses should therefore be interpreted as plausibility tests, not automatic reclassification rules.

Fourth, the global AF `> 1x10^-5` threshold was used as a screening alert rather than a universal classification boundary. Its appropriateness may vary across genes, diseases, ancestries, and penetrance regimes. Threshold sensitivity analyses are therefore important for calibrating the stability of specific alert counts.

Fifth, the principal illustrative cases are intentionally few. This manuscript is designed to establish the existence and structure of under-specified public label states, not to estimate the full prevalence of clinically consequential downstream misinterpretation across all genes or all clinical settings.

Sixth, ClinVar captures classification-level assertions rather than patient-level outcomes. We therefore identify interpretive risk and state ambiguity, not directly measured rates of downstream misdiagnosis or management change.

Seventh, Tier 2 locus/regional context is informative but not interchangeable with allele-resolved AF. It should be interpreted as structured partial context, not as a substitute for allele-level frequency evidence.

## Clinical and Research Implications

For inherited arrhythmia variant review, popmax should be used routinely alongside global AF when evaluating frequency tension. `AC_count` should be treated as graded support for the reliability of the observed population signal, not as a substitute for disease modeling. Non-observation at the allele-resolved level should not be treated as generic rarity, particularly for indels, duplications, splice-disrupting variants, and complex alleles.

More broadly, public P/LP assertions should be interpreted as compressed public objects rather than fully specified disease statements. The practical question is not only whether a variant has been called pathogenic, but what kind of pathogenic claim the public label can still support once population frequency and mechanism are taken seriously, and whether allele-resolved population frequency is even operationally available for that variant under current representation constraints.

The interpretability tier framework proposed here offers a practical route for downstream review. A variant with allele-resolved AF can be subjected to direct ancestry-aware plausibility analysis. A variant with only locus/regional context requires representation-aware caution and should not be over-read as rare. A genuinely unevaluable variant marks the boundary at which public population-frequency logic cannot be cleanly applied.

For Tier 2 variants, non-observation should not be treated as rarity without additional review. At minimum, reviewers should assess:

1. whether the variant class is representation-unstable, particularly for indels, duplications, splice-disrupting variants, or complex alleles;
2. whether same-locus population variation is present despite lack of allele equivalence;
3. whether transcript-aware normalization or decomposition could alter representation;
4. whether haplotypic or complex-allele structure is plausible; and
5. whether assertion provenance is limited to a single submitter without stronger review support.

That distinction is likely relevant beyond arrhythmia genes. Any setting that combines incomplete penetrance, ancestry-specific enrichment, founder effects, recessive architecture, technically unstable representation, and public reuse of compressed pathogenic labels may be vulnerable to the same interpretive failure.

## Conclusion

The central problem exposed here is not simply that some ClinVar P/LP variants are too frequent. It is that a single exported public label can fail to preserve the biological state required to constrain downstream interpretation. In inherited arrhythmia genes, ancestry-aware population frequency and mechanism reveal at least three different interpretation regimes hidden beneath the same public P/LP assertion: hard incompatibility with an unqualified dominant high-penetrance disease reading, parameter-sensitive low-penetrance boundary cases, and carrier-compatible recessive architecture.

But the problem is deeper than disease-model ambiguity alone. Even after trim- and decomposition-aware reconciliation, exact/equivalent recovery increased only from 350 to 357 variants. Thus, the primary limitation is not cosmetic matching failure but the operational unavailability of allele-resolved population frequency for most public P/LP variants under current representation regimes.

This limitation is not fully explained by lower-credibility public assertions. Provenance-aware filtering removes the 3 most extreme fragile records by design, but the broader structure of under-specification persists.

A public pathogenic label is therefore not, on its own, a disease state, and for most variants it is not attached to an object that can be cleanly constrained by allele-resolved population frequency. Until both the relevant disease state and the relevant level of population interpretability are made explicit, public P/LP assertions should be treated as under-specified.
