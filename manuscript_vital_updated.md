# A Public Pathogenic Label Is Not a Disease State: Population Evaluability Reveals Hidden Interpretation Regimes in ClinVar Arrhythmia Assertions

## Abstract

Public pathogenic and likely pathogenic (P/LP) assertions are widely reused across clinical databases, laboratory pipelines, and research workflows as though they carry a stable and portable disease claim. In practice, however, a public P/LP label is a compressed clinical object. Once separated from its original submission context, it does not reliably preserve the inheritance model, penetrance assumptions, or evaluability conditions required for responsible downstream interpretation.

We examined this problem in a high-stakes domain by analyzing 1,731 unique ClinVar P/LP variants across 20 inherited arrhythmia genes. We integrated ancestry-aware population frequency, allele-resolved variant representation, allele count support, and disease mechanism, and classified each variant into one of three evaluability tiers according to the depth of population constraint that could technically be applied.

After strict initial matching, exhaustive trim-aware and decomposition-aware reconciliation, and full reference-based normalization with `bcftools norm`, only 357 of 1,731 variants (20.6%) achieved allele-resolved population context. The remaining 1,326 (76.6%) retained only locus- or region-level context without allele-resolved allele frequency, and 48 (2.8%) remained unevaluable. Reference-based normalization changed 0 of 1,731 variant representations, indicating that the dominant limitation was not residual normalization failure but the structural unavailability of allele-resolved population constraint under current representation and aggregation regimes.

Disease-model analysis was possible only for the Tier 1 subset with usable allele frequency data (`n = 334`, 19.3% of the full cohort). Within this subset, a screen based on global allele frequency (AF) alone identified 13 variants above the `1x10^-5` review threshold. A global-or-popmax screen identified 115, meaning that global-only review missed 102 of 115 ancestry-aware frequency alerts (88.7%). Of these 115 alerts, 103 fell within genes linked to clinically actionable interpretation contexts, including drug restriction, intensive surveillance, and device-related management. These 115 variants, despite sharing the same public P/LP label, resolved into three biologically distinct interpretation regimes: hard incompatibility with an unqualified dominant high-penetrance reading (`n = 1`), boundary or monitoring status (`n = 76`), and recessive or carrier-compatible architecture (`n = 38`).

Two limitations define the scope of inference. First, all disease-model findings are restricted to the minority of variants for which allele-resolved population constraint was technically recoverable. The prevalence of the same regime structure in the unevaluable majority cannot be estimated from these data. Second, the clinical-action classification used here identifies exposure to consequential decision environments, not measured downstream harm.

These findings support a governance conclusion rather than a simple technical one. A public P/LP label does not reliably preserve the disease-state definition required for safe downstream reuse. For most variants in this study, one of the key constraining signals, allele-resolved population frequency, was not operationally available. Where such constraint was available, the same exported label supported mutually inconsistent biological readings. Responsible reuse therefore requires explicit specification of the disease model being invoked and explicit declaration of the level of population evaluability available for the asserted allele.

## Introduction

Public clinical variant assertions now function as shared infrastructure. A pathogenic or likely pathogenic label exported from ClinVar is routinely reimported into laboratory interpretation pipelines, cascade testing workflows, decision-support systems, and research analyses as though it were a portable unit of clinical meaning. This pattern of reuse rests on an implicit assumption: that the public label preserves the disease state it was originally intended to describe.

That assumption is structurally unsafe. A public P/LP assertion is not a self-sufficient disease claim. It is a compressed object that records a classification outcome while omitting key components of the interpretive model that produced it, including the relevant inheritance logic, penetrance assumptions, ancestry context, and the evaluability conditions under which the assertion was judged. Once the label circulates independently of those parameters, downstream users are forced to reconstruct a disease model that the public object does not itself encode.

This problem is not merely theoretical. In domains where variant interpretation influences drug restriction, surveillance intensity, cascade testing, or device management, importing the wrong disease model can change the practical meaning of a public label. The central question is therefore not whether an individual assertion may have been internally coherent in its original submission context. The question is whether the exported label remains interpretable as a coherent disease-state claim once ancestry-aware population structure and mechanism are reintroduced at the point of downstream use.

Inherited arrhythmia genes provide an unusually stringent test domain for this problem. They span dominant, recessive, and mixed architectures; variable penetrance; ancestry-specific enrichment patterns; founder effects; and clinically consequential interpretation contexts. They therefore constitute a stress-test setting in which both the biological and operational consequences of label reuse become especially visible. This is not because arrhythmia genes are unique, but because they make the underlying infrastructure problem easier to detect. If disease-state ambiguity and population evaluability limits are visible here, they are unlikely to be confined to this gene space.

In this study, we examined 1,731 unique ClinVar P/LP variants across 20 inherited arrhythmia genes and asked two related questions. First, for how many public assertions can allele-resolved population constraint actually be recovered under current public representation and aggregation systems? Second, when that constraint is available, what kinds of disease-model heterogeneity become visible beneath a shared public P/LP label?

We show that most public arrhythmia assertions cannot be connected to allele-resolved population frequency even after exhaustive reconciliation and reference-based normalization. Among the minority of variants for which allele-resolved population context can be established, ancestry-aware frequency exposes three biologically distinct interpretation regimes concealed beneath the same exported label. These findings do not estimate the full scale of the problem across all public assertions. Rather, they establish its structure and show that safe downstream reuse requires more explicit metadata than current public labels provide.

## Results

### Evaluability tiering: allele-resolved population constraint is unavailable for most variants

We assembled 1,731 unique ClinVar P/LP variants across 20 canonical inherited arrhythmia genes and classified each into one of three evaluability tiers according to the depth of allele-resolved population constraint that could be established against gnomAD v4.1.1 exomes after exhaustive reconciliation (Table 1).

Table 1. Evaluability tier distribution (`n = 1,731 variants`)

| Tier | Definition | N | % |
|---|---|---:|---:|
| 1 | Exact or equivalent allele-resolved match in gnomAD v4.1.1 exomes | 357 | 20.6 |
| 2 | Locus or regional context only; no allele-resolved AF | 1,326 | 76.6 |
| 3 | Genuinely unevaluable after all reconciliation steps | 48 | 2.8 |

Under strict initial matching, 350 variants achieved exact allele-level recovery. Trim-aware and decomposition-aware reconciliation increased this number to 357. Full reference-based normalization using `bcftools norm` against a local GRCh38 reference changed 0 of 1,731 variant representations. The combination of a marginal reconciliation gain (7 variants, 0.4%) and zero normalization-induced changes argues strongly against residual left-alignment artefacts or incomplete decomposition as primary explanations for missing allele-level matches. The dominant limitation is structural: under current representation and aggregation regimes, most public arrhythmia assertions cannot be connected to allele-resolved population frequency.

This finding defines an operational boundary on downstream interpretation. For most variants in this cohort, population frequency cannot be applied at the asserted allele level, not because it was overlooked analytically, but because the relevant representation is not recoverable from current public infrastructure.

### Tier 2 reflects structured representation limits, not random absence

The 1,326 Tier 2 variants were not uniformly uninformative. Within this group, 638 (48.1%) showed allele discordance at the queried locus, meaning that gnomAD contained a record at the same position but not for the exact ClinVar allele. The remaining 688 (51.9%) lacked even a same-locus gnomAD record and retained only regional context.

Tier 2 was also enriched for representation-sensitive variant classes. Indels, duplications, and insertions collectively accounted for 640 of 1,326 Tier 2 variants (48.3%), compared with 106 of 357 Tier 1 variants (29.7%), corresponding to an odds ratio of 2.21. This enrichment has a direct interpretive consequence: non-observation at the allele level is not a uniform proxy for rarity. For representation-sensitive classes, especially indels and duplications, the absence of an allele-resolved match often reflects discordance between public assertion and population aggregation systems rather than genuine population absence.

Accordingly, Tier 2 should not be treated as weak positive evidence of rarity. Its primary value is diagnostic: it identifies where rarity logic cannot safely be applied because allele-resolved evaluability has not been achieved.

![Allele-level non-observation is shaped by variant class](figures/arrhythmia_vital_absence_not_rarity.png)

### Global AF alone misses most ancestry-aware frequency signals

Disease-model analysis was restricted to Tier 1 variants with usable AF data (`n = 334`). Within this subset, 321 variants (96.1%) had global AF `<= 1x10^-5`, so a global-AF-only screen would have flagged 13 variants for closer review. Replacing global AF with the maximum of global AF and popmax AF increased the number of flagged variants to 115. Global-only review therefore missed 102 of 115 ancestry-aware frequency alerts (88.7%).

The threshold of `1x10^-5` is used here as a screening trigger rather than as a biological boundary separating pathogenic from benign variants. Sensitivity analysis indicated that the qualitative structure of the signal, specifically the large proportional gain from incorporating popmax, was robust to reasonable threshold variation, although absolute counts shifted accordingly.

Of the 115 alerted variants, 103 occurred in genes associated with clinically consequential interpretation contexts, including drug restriction, intensive surveillance, and device-related management. Some variants fell into more than one decision context. These classifications identify exposure to meaningful clinical decision environments, not measured downstream outcomes.

The central result of this analysis is therefore not merely that popmax increases alert counts. It is that global AF alone systematically suppresses ancestry-localized signals in a gene set where ancestry structure, founder effects, and non-uniform inheritance architecture are biologically relevant.

![Population-specific frequency signals missed by global AF](figures/arrhythmia_population_af_outliers.png)

### Three interpretation regimes beneath a shared public label

The principal finding was not the number of frequency alerts, but the biological heterogeneity revealed when ancestry-aware frequency was used to constrain the disease model implicitly carried by the public P/LP label. The 115 frequency-tension variants resolved into three distinct interpretation regimes (Table 2).

Table 2. Interpretation regimes among frequency-tension variants (`n = 115`)

| Regime | N variants | Representative locus | Defining feature |
|---|---:|---|---|
| Hard dominant incompatibility | 1 | SCN5A `VCV000440850` | AF incompatible with an unqualified dominant high-penetrance reading |
| Boundary / monitoring | 76 | KCNH2 `VCV004535537` | Exceeds strict ceiling but remains compatible with lower-penetrance dominant interpretation |
| Recessive / carrier-compatible | 38 | TRDN `VCV001325231` | Dominant reading implausible; recessive or carrier logic remains coherent |

#### Hard dominant incompatibility

One variant, SCN5A `VCV000440850`, `c.[3919C>T;694G>A]`, showed population signals incompatible with any unqualified dominant high-penetrance reading across a broad and deliberately permissive parameter space. Its AFR popmax AF was `5.68x10^-3` (allele count 214), exceeding a strict dominant high-penetrance ceiling by 113.5-fold, a Brugada-appropriate 20% penetrance ceiling by 45.4-fold, and a deliberately permissive low-penetrance dominant ceiling by 5.7-fold.

This does not exclude all possible biological relevance. Low penetrance, haplotypic dependence, or more complex allelic architecture remain conceivable. What it does exclude is the most common downstream reading of a standalone public P/LP label as a portable dominant Mendelian disease claim without qualification.

#### Boundary and monitoring cases

Seventy-six variants occupied an intermediate zone. These variants were incompatible with a strict high-penetrance dominant disease model but not contradicted to the same degree as the SCN5A outlier. KCNH2 `VCV004535537`, `c.2398+2T>G`, illustrates this class. It showed a global allele count of 24 and an ASJ popmax AF of `1.32x10^-4`, exceeding a strict long-QT high-penetrance ceiling by 2.6-fold while remaining compatible with a more permissive low-penetrance dominant model.

The allele count of 24 is part of the interpretation, not a side note. At this level, the signal is meaningful but still carries material statistical uncertainty. These variants therefore warrant qualified review and monitored interpretation rather than automated assignment of full dominant disease weight.

#### Recessive and carrier-compatible architecture

Thirty-eight variants were more coherent under recessive or carrier architecture than under a dominant disease interpretation. TRDN `VCV001325231`, `c.1050del`, illustrates this regime. It had a global allele count of 40 and an AMR popmax AF of `2.18x10^-4`, rendering a dominant Mendelian reading implausible while remaining fully coherent under recessive logic, with an implied carrier frequency of approximately `4.36x10^-4` and an expected homozygote frequency of approximately `4.76x10^-8`.

The core interpretive failure in this regime is not excessive frequency in isolation. It is the importation of the wrong inheritance-class reading from a public label that does not specify which disease model it encodes.

Across all three regimes, allele count acted as a graded modifier of confidence in the population signal. The SCN5A variant, with allele count 214, supported a robust and stable incompatibility conclusion. The TRDN variant, with allele count 40, supported a meaningful but more cautious inference. The KCNH2 variant, with allele count 24, was most appropriately treated as a monitoring case. Allele count is therefore not a binary quality filter but a continuous modifier of evidentiary weight.

## Discussion

### The evaluability boundary is structural, not technical

The most informative result of the reconciliation analysis is what further technical effort did not achieve. Despite strict matching, trim-aware reconciliation, decomposition-aware reconciliation, and full reference-based normalization, only 357 of 1,731 variants (20.6%) reached allele-resolved population context, and normalization contributed zero additional matches. This localizes the problem clearly. The dominant barrier is not an avoidable pipeline failure or a remediable normalization oversight. It is a structural property of the current interface between public clinical assertions and population aggregation systems.

This distinction matters because it changes what downstream users should infer from missing allele-level population context. In this setting, lack of allele-resolved recovery is often an infrastructure limitation rather than an interpretable biological signal. Analytical diligence cannot fully repair a representational interface that does not preserve allele-level evaluability across systems.

### Inherited arrhythmia genes function here as a stress-test domain

This study does not claim that inherited arrhythmia genes are uniquely affected. The point of this gene set is methodological. These genes combine ancestry structure, variable penetrance, dominant and recessive mechanisms, founder effects, and clinically consequential interpretation contexts. They therefore offer a stringent domain in which failures of disease-model portability and population evaluability are easier to detect.

The present findings should therefore be read as proof of structure rather than a direct estimate of prevalence across all disease areas. The study demonstrates that the problem exists in a domain where the stakes of interpretive slippage are unusually visible. Whether the same regime distribution holds elsewhere requires separate empirical study.

### The popmax gap is not a minor methodological detail

The finding that global-AF-only review missed 88.7% of ancestry-aware frequency alerts is not an artifact of an eccentric threshold choice. It reflects a general problem that arises whenever ancestry-heterogeneous frequencies are collapsed into a single global summary statistic. When a variant is rare globally but enriched in a specific ancestry group, global AF dilutes precisely the signal that matters most for disease-model constraint. Popmax restores that signal.

In inherited arrhythmia genes, where ancestry-localized enrichment, founder effects, and recessive carrier architecture are biologically relevant features rather than edge cases, global-only screening is not merely incomplete. It is misaligned with the biology it is meant to constrain.

### The label compression problem has both biological and operational dimensions

The three interpretation regimes identified here illustrate two distinct dimensions of label compression. The first is biological. The same public P/LP label may remain compatible with dominant high-penetrance disease, low-penetrance susceptibility, recessive or carrier architecture, or a biologically implausible unqualified dominant claim. The disease-state model is therefore not contained in the label itself.

The second dimension is operational. For most variants in this cohort, the principal quantitative constraint on that disease model, allele-resolved population frequency, was not technically recoverable. A downstream user facing a Tier 2 assertion therefore lacks both an explicit disease model and a clean population constraint. The public label is most compressed precisely where interpretive support is weakest.

This is why the problem cannot be reduced to a generic statement that "context matters." The more consequential claim is that the shared infrastructure surrounding public labels behaves as though that missing context were optional, even in settings where downstream interpretation depends on it.

## Limitations

First, all disease-model analyses were restricted to Tier 1 variants with usable allele-resolved AF data (`n = 334`, 19.3% of the full cohort). The interpretation regimes described here therefore cannot be assumed to represent the unevaluable majority. Because Tier 1 was enriched for SNVs and depleted of indels and duplications relative to Tier 2, direct extrapolation across variant classes would be inappropriate.

Second, this is a study of disease-model compatibility, not a study of patient-level outcomes. The clinical-action categories used here identify exposure to consequential decision environments but do not quantify realized downstream harm.

Third, frequency-based analysis can exclude a particular disease-model reading, but it cannot on its own establish the correct alternative interpretation. Low penetrance, haplotypic context, allelic phase, transcript-specific effects, segregation evidence, and functional validation remain outside the scope of this analysis unless orthogonal data are introduced.

Fourth, the `1x10^-5` threshold was used as a review trigger rather than as a biological boundary. The qualitative structure of the findings was robust to threshold variation, but the absolute number of flagged variants depends on the screening threshold chosen.

Fifth, Tier 2 does not provide allele-level frequency evidence. Its interpretive value lies in identifying where allele-resolved rarity logic cannot be safely applied, not in resolving the status of individual variants.

## Governance Implications

Three governance-relevant conclusions follow from these findings.

First, popmax should be required alongside global AF in inherited arrhythmia variant review. In this study, global-only screening missed 88.7% of ancestry-aware frequency alerts, and 103 of those alerts occurred in genes linked to clinically consequential interpretation contexts. Omitting popmax in such settings is not a small methodological omission. It is a systematic failure to apply available population constraint.

Second, public variant assertions intended for downstream reuse should carry three additional fields alongside pathogenicity category. These are:

1. the disease-state model being asserted, including inheritance class and penetrance logic, for example dominant high penetrance, dominant reduced penetrance, recessive, or unspecified;
2. the evaluability tier of the asserted allele with respect to population constraint, using a schema analogous to the one developed here; and
3. the allele count supporting any cited population signal, together with explicit acknowledgment of the uncertainty associated with low-count observations.

These additions do not require new biological evidence. They require explicit curation of information that downstream users already need but that the public label does not currently preserve.

Third, downstream users should treat non-observation in population databases as a property of representation rather than biology for indels, duplications, insertions, and splice-disrupting variants. The enrichment of such classes in Tier 2 indicates that absence of an allele-level gnomAD match cannot be naively interpreted as evidence of rarity. Evaluability metadata must therefore remain conceptually distinct from pathogenicity classification.

## Conclusion

Across 1,731 ClinVar P/LP variants in 20 inherited arrhythmia genes, two problems converge. The first is operational. Even after exhaustive reconciliation and reference-based normalization, 76.6% of variants lacked allele-resolved population context, meaning that one of the main quantitative constraints on disease-model evaluation was not technically recoverable for most public assertions in this setting. The second is biological. Among the minority of variants for which allele-resolved frequency could be established, ancestry-aware population data exposed three distinct interpretation regimes beneath the same exported P/LP label: hard incompatibility with an unqualified dominant high-penetrance reading, boundary cases requiring monitored interpretation, and recessive or carrier-compatible architecture.

These findings should be read as establishing a structural problem, not as claiming a complete estimate of its prevalence. The unevaluable majority cannot yet be analyzed in the same way, and that inability to assess scope is itself part of the governance failure documented here.

A public pathogenic label does not reliably preserve the disease-state definition required for safe downstream reuse. Responsible reuse requires, at minimum, explicit specification of the disease model being invoked and explicit declaration of the degree of population evaluability available for the asserted allele. Without both, a compact public object can circulate as though it carried stable clinical meaning while supporting mutually inconsistent biological readings.
