# Meta-analysis from separate cohort files

These commands are available in source builds after v0.3.5. Build the executable
with `./gradlew executableJar` (`.\gradlew.bat executableJar` on Windows).
The examples below use the resulting `build/cli/jlinalg-0.3.5.jar`, rather than
the previously published v0.3.5 release asset.

## Inputs and cohort order

Supply one TSV or CSV per independent cohort. Gzip-compressed `.tsv.gz` and
`.csv.gz` files are supported. Each file has unique feature IDs, effect
estimates, and positive standard errors:

```text
feature_id  beta  se
gene1       0.2   0.1
single      0.4   0.1
```

The displayed columns should be tab-separated in a TSV. IDs can identify
genes, probes, SNPs, or other features. Use `--id-column`, `--effect-column`,
and `--se-column` when the shared column names differ from `feature_id`,
`beta`, and `se`. Extra columns are ignored. Files may have different row
orders and feature sets. The command sorts inputs on disk and joins their
union by exact, case-sensitive feature ID, after trimming surrounding spaces.
Duplicate IDs within a cohort are errors, including duplicates with missing
values. A header-only cohort file is allowed.

Effect units and effect-allele orientations must already agree. Ratios such
as odds ratios should be converted to log effects and their corresponding
SEs before pooling. This CLI does not align alleles or model overlapping
cohort samples.

Repeated `--cohort NAME=FILE` arguments define the direction-string order.
Names are unique and use letters, digits, dot, dash, or underscore. An absent
feature or a missing beta/SE (`NA`, `NaN`, `.`, or empty) contributes `?`.
Positive, negative, and exactly zero effects contribute `+`, `-`, and `0`.
Missing pairs are excluded from fitting, while malformed numbers and
nonpositive SEs fail with a file/line diagnostic.

## Fixed- and random-effects pooling

From the repository root, run the supplied six-cohort example:

```powershell
java -jar build/cli/jlinalg-0.3.5.jar meta-analysis --cohort c1=examples/meta-analysis/c1.tsv --cohort c2=examples/meta-analysis/c2.tsv --cohort c3=examples/meta-analysis/c3.tsv --cohort c4=examples/meta-analysis/c4.tsv --cohort c5=examples/meta-analysis/c5.tsv --cohort c6=examples/meta-analysis/c6.tsv --model fixed --out fixed-meta.tsv
```

`gene1` has five complete cohorts and direction `++--?-`. `single` has
direction `+?????` and status `single_cohort`. Repeat with `--model random`
and a fresh output path for random-effects pooling. Random effects, REML,
normal coefficient inference, and 95% intervals are the defaults.

| Option | Choices / default |
| --- | --- |
| `--model` | `fixed`, `random` (default) |
| `--tau-estimator` | `reml` (default), `dl`, `pm`; ignored by fixed effects |
| `--inference` | `normal` (default), `t`, `hk`, `modified-hk` |
| `--confidence` | Confidence level in `(0,1)`; default `0.95` |
| `--min-cohorts` | Positive complete-cohort count; default `1` |
| `--max-iterations` | Heterogeneity optimizer iterations; default `200` |
| `--tolerance` | Relative heterogeneity tolerance; default `1e-10` |

`--min-cohorts 2` leaves the singleton row in the output with status
`below_min_cohorts` and `NA` estimates. Keeping these rows makes the excluded
features auditable. A row with no complete cohorts is also excluded.

With the default minimum of one, a single cohort passes through with its
supplied beta and SE and normal inference, even if a random-effects model
or Hartung–Knapp inference was requested. No pooling has occurred: tau², Q,
I², and prediction intervals are `NA`. The status and log identify this case.

## Meta-regression

Add a small moderator table containing a `cohort` column and numeric
moderator columns. Rows are matched by cohort name, so their order may differ
from the command line. Every configured cohort must have a moderator row;
a missing moderator value excludes that cohort from all regression rows.
The regression direction string and `n_cohorts` describe these complete cases.
Encode categorical moderators as numeric contrasts before running.

```powershell
java -jar build/cli/jlinalg-0.3.5.jar meta-regression --cohort c1=examples/meta-analysis/c1.tsv --cohort c2=examples/meta-analysis/c2.tsv --cohort c3=examples/meta-analysis/c3.tsv --cohort c4=examples/meta-analysis/c4.tsv --cohort c5=examples/meta-analysis/c5.tsv --cohort c6=examples/meta-analysis/c6.tsv --moderator-file examples/meta-analysis/moderators.tsv --moderators age --model random --out meta-regression.tsv
```

Use `--moderators age,proportion_female` for multiple moderators.
An intercept is included unless `--no-intercept` is supplied. Fixed and random
models, tau estimators, and inference options are shared with ordinary pooling.
Each fitted feature produces one output row per coefficient. Regression needs
more complete cohorts than coefficients, irrespective of `--min-cohorts`.
Unestimable rows receive `insufficient_cohorts` or `rank_deficient` with `NA`
estimates. Numerical fitting failures stop the command instead of publishing
an incomplete result.

## Output and throughput

The output includes feature ID, coefficient term, model, complete-cohort count,
direction, status, beta, SE, statistic, p-value, negative log10 p-value,
confidence limits, tau², Q/df/p-value, I², and prediction limits. Regression
also reports the omnibus moderator Q and its p-value. Its residual Q is the
zero-heterogeneity diagnostic, matching the Java API; moderator Q uses the
unscaled model covariance and remains an asymptotic chi-square test when
coefficient inference uses Hartung–Knapp.

Results are sorted lexicographically by feature ID. Output is CSV for `.csv`
paths and TSV otherwise, with an `OUT.log` sidecar recording settings and
cohort order. Existing result/log files are rejected. Choose a fresh path for
another run. The complete result is moved into place only after all rows
succeed; temporary sort files are then removed.

The log is opened and flushed at the start of execution. It records UTC
`started`/`finished`, numeric `elapsed_ms`, a readable `elapsed` duration,
and completion status. A failed run retains its timing log even when no result
is published; use a fresh output path when retrying.

`--sort-chunk-rows 100000` limits the number of cohort rows sorted in memory
at once. Merges open at most 32 runs per pass. `--block-rows 4096` bounds the
feature batch; use smaller values under memory pressure. Scratch disk space
is required beside the output, and temporary sorted data can exceed the input
size. Input parsing and joining are serial. `--threads N` parallelizes ordinary
batch fits in chunks of 2,048 features; meta-regression currently uses serial
CPU fitting and requires `--threads 1`.

The ordinary CLI passes primitive arrays into `PreparedMetaAnalysisBatch`
and does not construct a `List<MetaStudy>` per feature. Meta-regression uses
the new primitive-array overload. See the [array benchmark](../meta-cli-validation.md)
for measured allocation/time results and reproducible validation commands.
