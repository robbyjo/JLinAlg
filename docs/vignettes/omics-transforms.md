# Omics transforms from the command line

The main JLinAlg command accepts a deterministic transform pipeline for
feature-by-sample CSV or TSV omics matrices. Use it for expression,
methylation, protein, or other numeric molecular features before an OLS or GLM
association scan.

```powershell
java -jar build/cli/jlinalg-0.3.2.jar `
  --omics methylation.tsv `
  --pheno phenotype.tsv `
  --id IID `
  --formula "trait ~ age + sex + <omics>" `
  --transform "<omics> = mvalue(epsilon=1e-6) | zscore()" `
  --out ewas-results.tsv
```

`--transform` applies only to the `<omics>` term. It is not available for
VCF, BCF, or BGEN genotype inputs.

## Input and pipeline syntax

The omics file has one feature per row and one sample per column:

```text
feature_id  S1    S2    S3    S4
cg000001    0.12  0.20  NA    0.85
cg000002    0.44  0.51  0.49  0.55
```

Use the literal target `<omics>`, an equals sign, and one or more stages
separated by `|`:

```text
--transform "<omics> = stage(parameters) | next_stage()"
```

Quote the complete specification so that the shell does not interpret `<`,
`>`, or `|`. Whitespace around the target-assignment `=` and stage separator
`|` is optional; the examples include it to make those boundaries clear.
Stages run from left to right. The transform is applied
independently to every feature row after sample IDs have been aligned to the
phenotype file; it never pools values across features. Repeating
`--transform` appends more stages in command-line order, although one quoted
pipeline is usually easier to audit. Omitting the option is the same as
`identity()`.

## Built-in stages

| Stage | What it does | Requirements and defaults |
| --- | --- | --- |
| `identity()` | Copies values without changing them. | `identity` is also accepted. |
| `winsor(p=0.01)` | Clamps finite values to symmetric empirical quantiles. | `p` sets lower `p` and upper `1-p`. Quantiles use linear interpolation. |
| `winsor(lower=0.01,upper=0.99)` | Clamps finite values to separate empirical quantiles. | Defaults are `lower=0` and `upper=1`; require `0 <= lower <= upper <= 1`. |
| `winsor_mad(k=4)` | Clamps finite values to `median +/- k * MAD`, with `MAD = median(abs(x-median)) / qnorm(0.75)`. | `k` defaults to 4 and must be finite and nonnegative. Missing values do not enter either median; zero MAD leaves the row unchanged. `winsor_mad` and `winsor_mad()` are also accepted. |
| `log1p()` | Computes `ln(1+x)`. | Every finite value must be greater than `-1`. `log1p` is also accepted. |
| `log(offset=1)` | Computes `ln(x+offset)`. | The offset must be finite and every shifted finite value must be positive. |
| `zscore()` | Centers and scales using the row's sample standard deviation. | Requires at least two varying finite values. `zscore`, `zscore()`, and `zscore(ddof=1)` are accepted. |
| `int()` | Applies a tie-aware Blom rank inverse-normal transform. | Average ranks are used for ties; `int`, `int()`, and `int(method=blom)` are accepted. |
| `mvalue(epsilon=1e-6)` | Converts methylation beta values to base-2 log odds after clipping to `[epsilon, 1-epsilon]`. | Finite inputs must be in `[0,1]`; epsilon defaults to `1e-6` and must be in `(0,0.5)`. |

Names are case-insensitive. Use the parameter names shown above. Custom
expression syntax such as `expr(...)` is deliberately not enabled.

## Choose a pipeline

### Robust median/MAD clipping

Clamp outliers relative to a robust row center and scale, then standardize:

```text
--transform "<omics> = winsor_mad(k=4) | zscore()"
```

This matches R's `mad(..., constant = 1/qnorm(0.75), na.rm = TRUE)`
scaling. When the raw MAD is zero, the winsorization stage returns the row
unchanged rather than clipping values to the median.

### Methylation beta values

Convert beta values to M values, then standardize the tested effect to one
post-transform standard deviation:

```text
--transform "<omics> = mvalue(epsilon=1e-6) | zscore()"
```

M values are often preferable for modeling methylation intensity, while the
final z-score makes the reported coefficient a change in outcome per one
standard deviation of the transformed feature.

### Nonnegative expression measurements

Compress a long right tail and standardize:

```text
--transform "<omics> = log1p() | zscore()"
```

Use `log(offset=...)` instead when the required pseudocount is not one. For
example, values whose minimum is `-0.25` can use
`log(offset=0.250001)`, provided every shifted value is strictly positive.

### Protein abundance with outliers

Limit extreme observations before converting ranks to normal scores:

```text
--transform "<omics> = winsor(p=0.01) | int(method=blom)"
```

Rank inverse-normalization changes the estimand to the rank-normalized scale.
Choose it as an analysis decision, not as an automatic repair for poor data.
Applying `zscore()` after `int()` is generally redundant.

### Asymmetric clipping

Use separate tails when the scientific preprocessing rule is asymmetric:

```text
--transform "<omics> = winsor(lower=0.005,upper=0.975) | zscore()"
```

## Missing values and failures

Blank cells plus `.`, `NA`, `N/A`, `null`, and `NaN` are read as missing.
Built-in transforms leave them missing and estimate row statistics from finite
values only. After the complete pipeline, the CLI mean-imputes remaining
missing values within that transformed feature row.

Domain checks are strict. For example, `log1p()` rejects `-1`, `log()` rejects
a nonpositive shifted value, `mvalue()` rejects a beta value outside `[0,1]`,
and `zscore()` rejects a constant row. A transform-domain error or a row with
no finite values stops the run with an error. A feature that passes the
transform but is non-estimable after model adjustment uses the normal
structured failure output and accounting.

For reproducibility, record the exact transform string with the input version
and interpret coefficients on the final transformed scale. JLinAlg's output
manifest and log accompany the result file, but the transform string should
also be part of the analysis protocol.

## Trusted custom transforms

A trusted Java JAR can add a named stage by implementing
`OmicsTransformProvider`:

```java
package example;

import java.util.Map;
import org.jlinalg.pipeline.OmicsTransform;
import org.jlinalg.pipeline.OmicsTransformProvider;

public final class ScaleTransformProvider
        implements OmicsTransformProvider {
    @Override public String name() { return "scale"; }

    @Override
    public OmicsTransform create(Map<String, String> parameters) {
        double factor = Double.parseDouble(
            parameters.getOrDefault("factor", "1"));
        return values -> {
            double[] result = values.clone();
            for (int i = 0; i < result.length; i++)
                if (Double.isFinite(result[i])) result[i] *= factor;
            return result;
        };
    }
}
```

Register the provider in the JAR at
`META-INF/services/org.jlinalg.pipeline.OmicsTransformProvider` with this
single line:

```text
example.ScaleTransformProvider
```

Then load and use it like a built-in stage:

```text
--transform-plugin transforms.jar --transform "<omics> = scale(factor=2) | zscore()"
```

Providers must return a row with the same sample count. Plugin JARs execute
arbitrary JVM code in the analysis process, so load only artifacts you trust.
