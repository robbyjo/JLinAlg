# Multi-signal colocalization

Colocalization asks whether two association signals in the same region are
best explained by one shared causal variant. It is a follow-up to fine mapping
or cis-MR, not evidence that one trait causes the other.

## Align the inputs

Fit each trait with variants from the same genome build and ancestry-matched
LD. Variant identifiers must be unambiguous. `ColocSusie` intersects IDs in the
first input's order, so allele alignment and duplicate removal must happen
before fine mapping. Sample sizes, phenotype scales, and SuSiE prior choices
should be recorded with the result.

```java
SusieResult exposure = Susie.fitSummary(
    exposureZ, exposureLd, exposureN, variantIds,
    SusieOptions.defaults(), BackendPolicy.PREFERRED);
SusieResult outcome = Susie.fitSummary(
    outcomeZ, outcomeLd, outcomeN, variantIds,
    SusieOptions.defaults(), BackendPolicy.PREFERRED);

ColocSusieResult coloc = ColocSusie.analyze(
    exposure, outcome, ColocOptions.defaults());
```

`ColocSusieInput` can instead consume a signal-by-variant matrix of log Bayes
factors. This is useful when the SuSiE fits were produced elsewhere.

## Priors, overlap, and trimming

The defaults match `coloc::coloc.susie`: trait-specific priors `p1=p2=1e-4`,
shared prior `p12=5e-6`, and removal of signals with less than 0.5 posterior
mass on the common variant set. `ColocOptions` exposes all four values and
optional positive variant-specific prior weights. Sensitivity analysis should
vary plausible priors; a high H4 under only one aggressive prior is fragile.

Each trait's weights are normalized on the common variants. The conditional
shared-variant posterior includes the **product** of the two normalized prior
weights, consistently with H4. Multiplying either weight vector by a positive
constant does not change the result, including very large or small finite
scalings. Overlap trimming and the reported per-trait lead remain based on the
input Bayes factors with uniform variant weights; they are not a new weighted
fine-mapping fit.

Inspect `commonVariants()` and `skippedSignalPairs()`. Heavy trimming usually
indicates mismatched variant coverage or a signal whose posterior lies outside
the overlap, not absence of biological sharing.

## Interpret H0-H4 and the shared variant

```java
for (int index = 0; index < coloc.signalPairs().size(); index++) {
    ColocSignalPair pair = coloc.signalPairs().get(index);
    double[] conditionalShared = coloc.sharedVariantPosterior(index);
    System.out.printf("signals %d/%d H3=%.3f H4=%.3f lead=%s%n",
        pair.trait1EffectIndex(), pair.trait2EffectIndex(),
        pair.posteriorH3(), pair.posteriorH4(),
        coloc.commonVariants().get(maxIndex(conditionalShared)));
}
```

- H0: neither trait is associated in the region.
- H1/H2: only trait 1/trait 2 is associated.
- H3: both are associated, with distinct causal variants.
- H4: both are associated and share a causal variant.

The five values are posterior probabilities for each signal pair. The
variant-level posterior is conditional on H4 and must be interpreted together
with H4, not as an unconditional causal probability.

When the traits have disjoint finite Bayes-factor support, H4 is zero and the
conditional shared-variant vector is all zeros (an undefined conditional
distribution), but H0–H3 are still evaluated. In particular, strong disjoint
signals favor H3, not H0. Distinct configurations are accumulated directly in
log space in O(P) time per signal pair; no subtraction of nearly equal shared
and total evidence is used. This preserves tiny nonzero H3 when a shared SNP
dominates. Log evidence outside double precision's representable range is
rejected instead of returning NaN posteriors.

## xWAS and MR follow-up

For a retained cis-MR/xWAS hit, reconstruct or retain the harmonized regional
summary statistics, fine-map exposure and outcome separately, then run the
analysis above. Report the MR estimate, H3/H4, common-variant count, skipped
signals, priors, and the leading H4-conditional variants. Strong MR with H3
favored over H4 is a warning that LD between distinct signals may explain the
association.

The implementation is regression-tested against `coloc::coloc.susie` 5.2.3
for single- and multi-signal fixtures. See the
[verification guide](../performance-benchmarks.md) and the compact
[SuSiE example](susie-and-sem.md).

The post-`dd45e80` boundary audit adds independent base-R enumeration, weighted
posterior and scaling checks, disjoint-support tests, and warm checksum-gated
timings. Reproduce with
`src/test/resources/r-reference/run-genetic-audit.ps1`; raw evidence and
workload limitations are in
[the genetic audit report](../../src/benchmark/resources/genetic-audit/AUDIT.md).
