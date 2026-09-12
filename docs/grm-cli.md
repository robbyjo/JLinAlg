# Build a genomic relationship matrix from the CLI

The `grm` command constructs an additive genomic relationship matrix. The existing
`--grm` analysis option reads that matrix for Gaussian/non-Gaussian mixed models
or a Cox kinship frailty. Construction is new in source builds after v0.3.5;
build with `./gradlew executableJar` (Windows: `.\gradlew.bat executableJar`).

## Construct and use a GRM

Run this small example from the repository root, with fresh output paths:

```shell
java -jar build/cli/jlinalg-0.3.5.jar grm --genotypes examples/rare-meta/cohort-a.vcf --maf 0.01 --call-rate 0.95 --block-size 256 --backend cpu --out build/grm-demo/cohort.tsv
java -jar build/cli/jlinalg-0.3.5.jar --pheno examples/rare-meta/cohort-a.tsv --id sample --formula "trait ~ 1" --grm build/grm-demo/cohort.tsv --backend cpu --out build/grm-demo/model.tsv
```

The six-person, two-variant fixture demonstrates file interoperability; it is
too small for meaningful relatedness estimation. For an analysis, supply a
prespecified set of quality-controlled, LD-pruned diploid variants. LD pruning,
ancestry adjustment, sample selection, and leave-one-chromosome-out (LOCO) GRMs
are not performed automatically. Construct each desired subset separately.

## Inputs and options

`--genotypes` accepts VCF, VCF.gz/VCF.bgz, BCF, BGEN layout 2, or a plain CSV/TSV
matrix with variants in rows and samples in columns. Format follows the suffix.
Delimited input can have `id,sample1,sample2,...` or
`id,chromosome,position,ref,alt,sample1,...` headers, with commas replaced by tabs
for TSV. Values are alternate-allele dosages in `[0,2]`; `NA` denotes missing.
All source samples are included, in source order. For BGEN without embedded IDs,
add `--bgen-samples cohort.sample` (Oxford sample file or one ID per line).
VCF/BCF uses the existing additive dosage reader. Supply biallelic diploid data
for this diploid relationship estimator; the command does not infer ploidy corrections.

| Option | Default | Meaning |
|---|---|---|
| `--genotypes FILE` | required | Source genotype file |
| `--out FILE` | required | Labeled `.tsv` or `.csv` GRM |
| `--maf` | `0.01` | Minimum minor allele frequency, in `[0,0.5)` |
| `--call-rate` | `0.95` | Minimum called-sample fraction, in `(0,1]` |
| `--block-size` | `256` | Positive number of source variants per block |
| `--backend` | `cpu` | Compute backend policy; explicit accelerators must be available |
| `--bgen-samples FILE` | none | External sample IDs for BGEN |

The command rejects existing result or log paths. `OUT.log` records UTC start/end,
readable elapsed time, source and format, filters, sample count, considered/used/
excluded variant counts, and selected backend. Errors retain a failed timing log
but do not publish a partial matrix. `grm --help` prints the command syntax.

## Matrix definition and limits

For each retained variant, estimate `p` from called dosages in the source sample.
Set `z_i=(g_i-2p)/sqrt(2p(1-p))`; missing dosages are mean imputed, so their
standardized value is zero. All-missing and monomorphic variants are excluded,
as are variants below the MAF or call-rate cutoff. The matrix is
`K = sum(z z') / number_of_retained_variants`.

This is an additive relationship matrix; kinship coefficients are `K/2`.
The denominator is common to all sample pairs. It is not the pairwise-observed
denominator used by some external GRM implementations. Do not expect equality
with an external tool using different filters, missingness, scaling, or weights.

Genotypes are streamed in blocks and filtered-block sums are divided by the
total retained variant count, including when blocks retain different numbers.
The matrix itself remains dense: one double matrix needs `8*N*N` bytes, with
additional construction copies and block workspace. Output is also quadratic;
the log's matrix-byte count is not a peak-memory estimate. Large sample counts
can exceed heap or Java array limits. No sparse approximation is silently used.

The Java equivalent is
`GenomicRelationshipMatrix.fromSource(source, options, backendPolicy, blockSize)`.
It shares normalization and filters with `fromVariantDosages`. Regression tests
use analytically derived matrix entries, missing calls, excluded/unequal blocks,
CSV/TSV round trips, and a VCF-to-mixed-model CLI run. See the
[mixed-model CLI guide](vignettes/command-line.md) for phenotype/sample alignment,
repeated observations, and existing labeled/GCTA GRM input.
