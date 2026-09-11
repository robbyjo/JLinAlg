# CLI-only SuSiE fine mapping

SuSiE fine mapping separates association evidence in one locus into a small
number of single effects and reports variant posterior inclusion
probabilities (PIPs) and credible sets. This command uses summary statistics,
an LD correlation matrix, and the analysis sample size. No Java code is
required.

> **Release availability:** the `susie` CLI command is included in v0.3.5
> and later.

## 1. Define and harmonize one locus

Choose a locus window and one genome build. Use the same effect-allele
orientation for summary statistics and LD. Remove duplicate variant IDs.
External LD should be ancestry-matched to the association study; in-sample LD
is preferable when available.

The summary file may be CSV or TSV. It needs one row per variant and either a
z score:

```text
variant_id,z
rs1001,8.42
rs1002,2.13
rs1003,-6.90
```

or beta and positive SE:

```text
SNP,beta,se,effect_allele,other_allele
rs1001,0.168,0.020,A,G
rs1002,0.041,0.019,C,T
rs1003,-0.138,0.020,T,C
```

The allele columns are useful provenance but the SuSiE command reads the
variant ID and z, or computes z as beta/SE. It does not flip alleles or
perform liftover. Standard names are detected; unusual names can be specified
with `--variant-column`, `--z-column`, `--beta-column`, and
`--se-column`.

## 2. Prepare the LD matrix

An unlabeled LD CSV/TSV is accepted only when its rows and columns are in the
exact summary-row order:

```text
1.0,0.31,0.02
0.31,1.0,0.08
0.02,0.08,1.0
```

A labeled matrix is safer because JLinAlg checks identifiers and reorders it
to the summary file:

```text
variant_id,rs1003,rs1001,rs1002
rs1003,1.0,0.02,0.08
rs1001,0.02,1.0,0.31
rs1002,0.08,0.31,1.0
```

The first row and first column must contain the same unique variant set as the
summary file. The matrix must be finite, symmetric, positive semidefinite, and
have a valid correlation diagonal. Singular LD is allowed; inconsistent or
indefinite LD is rejected.

If LD is generated with PLINK or another program, preserve the variant order
file beside the matrix. Either convert it to the labeled layout above or
reorder the summary file exactly before using an unlabeled matrix.

## 3. Run SuSiE

```powershell
java -jar jlinalg-<version>.jar susie `
  --summary locus-summary.csv --ld locus-ld.csv --sample-size 125000 `
  --effects 10 --coverage 0.95 --min-purity 0.5 `
  --max-iterations 200 --tolerance 1e-6 `
  --out locus-susie.csv
```

`--effects` is the maximum number of single effects and is capped by the
number of variants. The defaults are 10 effects, 200 iterations, tolerance
1e-6, prior variance 0.2, 95% credible-set coverage, minimum purity 0.5, and
estimated residual variance. Use `--fixed-residual-variance` when the
analysis plan requires the summary-statistic residual variance to remain
fixed. `--backend cpu` gives the portable path; `preferred` uses the
configured backend policy.

## 4. Read the three output files

With `--out locus-susie.csv`, JLinAlg writes:

1. `locus-susie.csv` in CSV because the name ends in `.csv`. It has
   `variant_id`, `pip`, and `posterior_mean`.
2. `locus-susie.csv.effects.tsv`. For every retained credible effect, it has
   every locus variant's `alpha`, effect posterior mean, and
   `log_bayes_factor`, plus credible-set membership, coverage, and purity.
   `effect_index` is one-based. This file is accepted directly by
   `jlinalg coloc`.
3. `locus-susie.csv.log`, containing inputs, sample size, variant/effect
   counts, credible-set count, iterations, convergence, objective, residual
   variance, backend, and output paths.

Choose explicit sidecar paths if preferred:

```powershell
java -jar jlinalg-<version>.jar susie `
  --summary locus-summary.tsv --ld locus-ld.tsv --sample-size 125000 `
  --out locus-pip.tsv --effects-output locus-effects.tsv `
  --log locus-susie.log
```

Output delimiters follow the output suffix. Existing files are protected
unless `--overwrite` is supplied.

## 5. Check the fit before using it

Require `converged=true` in the log before treating the result as final.
Inspect:

- PIP for evidence that a variant participates in at least one effect.
- Credible-set coverage and membership.
- Purity, the minimum absolute within-set LD; diffuse low-purity sets should
  not be presented as well-resolved signals.
- Sensitivity to the maximum effect count, prior variance, and LD source.

In summary mode, posterior means are on the standardized summary-statistic
scale. They are not automatically original-unit GWAS effects. PIPs and
credible sets also depend on the chosen locus and on the supplied LD.

If no effect passes the purity rule, the effects sidecar contains only its
header. That is a valid indication that no coloc-ready credible signal was
retained, not permission to lower purity after seeing the desired result.

## 6. Common failures

- **LD and summary identifier sets differ:** use a labeled matrix or explicitly
  intersect both inputs before fitting.
- **LD is not symmetric or positive semidefinite:** verify allele order,
  variant order, missing-value handling, and LD calculation; do not silently
  replace the matrix with identity.
- **Duplicate IDs:** create unambiguous build-aware IDs such as
  `chr:position:ref:alt`.
- **No convergence:** increase iterations only after checking LD and locus
  definition; compare results across sensible priors.
- **Very different ancestry:** obtain a better matched LD reference rather
  than interpreting unstable credible sets.

Continue with [CLI-only colocalization](cli-colocalization-tutorial.md).
The [SuSiE API vignette](susie-and-sem.md) explains individual-level and
sufficient-statistic modes and the numerical validation in more depth.
