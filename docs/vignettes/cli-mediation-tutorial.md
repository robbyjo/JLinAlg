# CLI-only mediation analysis

JLinAlg's `mediation` command fits Gaussian linear mediation with optional
numeric covariates, grouped random intercepts, or pedigree relatedness. It
reports the treatment-to-mediator path, mediator-to-outcome path, indirect
effect, direct effect, and total effect. No Java code is required.

> **Version availability:** the `mediation` CLI command was added after
> v0.3.4. Use current `main` (`jlinalg-0.3.5-SNAPSHOT.jar`) or a later
> release; it is not present in the published v0.3.4 JAR.

## 1. State the three variables

For treatment or exposure X, mediator M, outcome Y, and covariates C, the
command fits:

```text
M = intercept + a X + alpha C + error
Y = intercept + direct X + b M + gamma C + error
Y = intercept + total X + delta C + error
indirect = a * b
```

The input is observation-by-variable CSV or TSV. Y, X, M, and covariates must
be numeric:

```text
sample_id,treatment,mediator,outcome,age,sex01,subject
O001,0.2,1.1,2.8,51,0,P001
O002,0.8,1.7,3.9,48,1,P002
O003,1.1,2.2,4.7,61,0,P003
```

Categorical covariates must currently be encoded into numeric indicator
columns before using this specialized command. Do not include an intercept;
JLinAlg supplies it.

## 2. Ordinary Gaussian mediation

```powershell
java -jar jlinalg-<version>.jar mediation --input mediation.csv `
  --outcome outcome --treatment treatment --mediator mediator `
  --covariates age,sex01 --confidence 0.95 `
  --out mediation-effects.csv
```

The same complete-case sample is used for all three component models. Rows
missing the outcome, treatment, mediator, or any selected covariate are
omitted and counted on screen and in the log.

## 3. Understand the five output rows

- `a`: treatment-to-mediator coefficient.
- `b`: mediator-to-outcome coefficient adjusted for treatment.
- `indirect`: product `a*b`.
- `direct`: treatment coefficient in the outcome model containing M.
- `total`: treatment coefficient in the outcome model without M.

Every row reports estimate, standard error, statistic, p-value, confidence
limits, and degrees of freedom. The indirect SE is the first-order Sobel
delta-method value:

```text
sqrt(b^2 Var(a) + a^2 Var(b))
```

Its interval and p-value are asymptotic normal. The command does not bootstrap,
simulate, or sample a posterior distribution.

## 4. Repeated observations or clusters

Add a random-intercept column with `--group`:

```powershell
java -jar jlinalg-<version>.jar mediation --input repeated.tsv `
  --outcome outcome --treatment treatment --mediator mediator `
  --covariates age,sex01 --group subject `
  --out mediation-subject-reml.tsv
```

This fits the same grouped random-intercept structure in the mediator,
outcome, and total-effect components using sparse REML. Repeat `--group` or
pass comma-separated columns for multiple independent random intercepts:

```text
--group subject --group clinic
```

Check `converged=true` in the log. The three component models estimate their
own covariance parameters and are warm-started in sequence.

## 5. Pedigree mediation

The mediation table needs a subject column, which may repeat:

```text
observation,subject,treatment,mediator,outcome,age
O101,1001,0.2,1.1,2.8,51
O102,1002,0.8,1.7,3.9,48
O103,1003,1.1,2.2,4.7,29
O104,1003,1.3,2.5,5.0,31
O105,9009,0.5,1.4,3.3,44
```

The pedigree file names individuals and their two parents:

```text
family_id,member_id,parent1_id,parent2_id
F10,1001,0,0
F10,1002,0,0
F10,1003,1001,1002
```

Invoke the columns explicitly:

```powershell
java -jar jlinalg-<version>.jar mediation --input mediation-pedigree.csv `
  --outcome outcome --treatment treatment --mediator mediator `
  --covariates age --individual-id subject `
  --pedigree pedigree.csv --pedigree-id member_id `
  --sire-id parent1_id --dam-id parent2_id `
  --out mediation-pedigree.tsv
```

`--individual-id subject` selects the subject key in the mediation table.
`--pedigree-id`, `--sire-id`, and `--dam-id` select the corresponding
pedigree columns. Blank, `NA`, `0`, `.`, and `-9` parent values mean
unknown.

If `--pedigree-family-id` is present, known subject values must be qualified
as `family_id:member_id`, for example `F10:1003`. Omit this option when
member IDs are already globally unique. A phenotype subject absent from the
pedigree, such as 9009, is added as an unrelated singleton and reported.
Repeated rows with the same absent subject share one singleton founder.

Ordinary grouped effects can accompany pedigree relatedness, for example
`--group clinic`. Pedigree ancestry uses Henderson's sparse additive
relationship precision in all three REML component models.

## 6. Output files and checks

`--out mediation-effects.csv` writes:

- `mediation-effects.csv`, using its extension to choose CSV or TSV;
- `mediation-effects.csv.log`, with model type, variable names, selected
  covariates/groups, input and retained rows, singleton counts, confidence,
  backend, convergence, and output path.

Use `--log FILE` for another log path and `--overwrite` to replace existing
outputs deliberately.

Before causal interpretation, assess treatment-outcome, treatment-mediator,
and mediator-outcome confounding assumptions, temporal order, measurement
error, and model form. The Sobel calculation assumes zero cross-model
covariance between a and b. This is not a general multilevel causal mediation
model with correlated equation errors or random slopes.

The [mediation methods vignette](mediation.md) documents the Java component
fits and numerical validation. Return to the
[progressive CLI association tutorial](cli-association-tutorial.md) for the
broader model workflow.
