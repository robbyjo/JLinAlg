# CLI-only multi-signal colocalization

Colocalization asks whether two association signals in the same region are
best explained by a shared causal variant. It does not estimate a causal
effect and does not replace Mendelian randomization.

> **Version availability:** the `coloc` CLI command was added after v0.3.4.
> Use current `main` (`jlinalg-0.3.5-SNAPSHOT.jar`) or a later release; it
> is not present in the published v0.3.4 JAR.

## 1. Prepare two comparable regional analyses

Before fitting either trait, verify:

- the same genome build and unambiguous variant identifiers;
- allele-aligned summary statistics;
- an appropriate LD matrix for each trait's ancestry and sample;
- a defensible common locus boundary;
- recorded sample sizes, phenotype scales, and SuSiE options.

Fine-map each trait separately:

```powershell
java -jar jlinalg-<version>.jar susie `
  --summary expression-locus.tsv --ld expression-ld.tsv `
  --sample-size 1200 --effects 10 --out expression-susie.tsv

java -jar jlinalg-<version>.jar susie `
  --summary disease-locus.tsv --ld disease-ld.tsv `
  --sample-size 85000 --effects 10 --out disease-susie.tsv
```

Check that both logs report convergence. The two files needed below are
`expression-susie.tsv.effects.tsv` and
`disease-susie.tsv.effects.tsv`.

## 2. Understand the colocalization input

The `susie` effects sidecar is long format. The required columns are:

```text
effect_index,variant_id,log_bayes_factor
1,rs1001,15.22
1,rs1002,4.01
1,rs1003,-0.30
```

The sidecar contains all variants for every credible effect, not only members
of that effect's credible set. A custom CSV/TSV from another SuSiE
implementation can be used if it has `effect_index`, `variant_id`, and
`log_bayes_factor` (aliases `effect`, `SNP`/`rsid`, and
`lbf`/`logbf` are recognized). Each effect within one file must contain
the same unique variant set.

Trait 1 and trait 2 may contain different regional variant sets. JLinAlg
intersects them by exact identifier in trait-1 order and reports the common
count. It does not infer equivalent variants from coordinates or alleles.

## 3. Run colocalization

```powershell
java -jar jlinalg-<version>.jar coloc `
  --trait1 expression-susie.tsv.effects.tsv `
  --trait2 disease-susie.tsv.effects.tsv `
  --out expression-disease-coloc.csv
```

Defaults match the usual `coloc.susie` priors:

- trait-1 association prior `p1 = 1e-4`;
- trait-2 association prior `p2 = 1e-4`;
- shared association prior `p12 = 5e-6`;
- minimum posterior overlap 0.5, with low-overlap signal pairs trimmed.

Specify a sensitivity run explicitly:

```powershell
java -jar jlinalg-<version>.jar coloc `
  --trait1 expression-susie.tsv.effects.tsv `
  --trait2 disease-susie.tsv.effects.tsv `
  --trait1-prior 1e-4 --trait2-prior 1e-4 --shared-prior 1e-6 `
  --minimum-overlap 0.5 --out coloc-sensitive.tsv
```

`--no-trim` evaluates all signal pairs, including pairs with little
posterior mass on the common variants. This is a diagnostic option, not a
remedy for badly mismatched variant coverage.

## 4. Read the outputs

With `--out expression-disease-coloc.csv`, JLinAlg writes:

1. `expression-disease-coloc.csv`, one row per retained signal pair, with
   source effect indices, common-variant count, lead variants, and posterior
   probabilities `posterior_h0` through `posterior_h4`.
2. `expression-disease-coloc.csv.variants.tsv`, one row per common variant
   per retained signal pair, with `posterior_shared_given_h4`.
3. `expression-disease-coloc.csv.log`, with input paths, signal counts,
   common variants, retained and skipped pair counts, priors, overlap rule,
   and output paths.

Use `--variant-output FILE` and `--log FILE` for explicit sidecar names.
The main and variant outputs use commas only when their own filename ends in
`.csv`; otherwise they use tabs. Existing files are protected unless
`--overwrite` is supplied.

## 5. Interpret H0 through H4

- H0: neither trait has an association in the region.
- H1: only trait 1 is associated.
- H2: only trait 2 is associated.
- H3: both are associated, but with distinct causal variants.
- H4: both are associated and share a causal variant.

These five posterior probabilities apply to one pair of SuSiE signals and sum
to one. A high H4 supports a shared regional signal under the supplied priors,
LD, variant set, and fine-mapping models. A high H3 supports distinct signals.

`posterior_shared_given_h4` is conditional on H4. For example, a variant
posterior of 0.90 is not an unconditional 90% probability of a shared causal
variant if H4 itself is only 0.05. Report both values.

## 6. Diagnostics and sensitivity

Inspect `common_variants` and `skipped_signal_pairs` in the log. Heavy
trimming commonly means that one trait's posterior lies mostly outside the
shared variant set. Fix coverage or harmonization rather than immediately
disabling trimming.

Repeat the analysis over prespecified plausible priors. Confirm that the
conclusion is not driven by one unusually permissive shared prior. Also
compare locus boundaries and appropriate LD references. Do not combine
results across genome builds without explicit liftover and allele checks.

For MR or xWAS follow-up, report the causal estimate and its diagnostics
separately from H3/H4, common-variant count, skipped pairs, priors, and leading
H4-conditional variants. Strong MR evidence with H3 favored over H4 warns
that LD between distinct signals may explain the association.

The [multi-signal colocalization vignette](colocalization.md) documents the
algorithm and R validation. The [CLI-only SuSiE tutorial](cli-susie-tutorial.md)
explains how to create the input sidecars.
