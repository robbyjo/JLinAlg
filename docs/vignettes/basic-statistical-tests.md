# Basic statistical tests

`StatisticalTests` provides one camel-case entry point and one immutable
`StatisticalTestResult` shape for the common tests in base R's `stats`
package. Tests already present in JDistlib are delegated to the pinned
JDistlib implementation; JLinAlg adds the missing table, blocked-design,
correlation, ANOVA, proportion, and power calculations.

```java
import org.jlinalg.stats.Alternative;
import org.jlinalg.stats.StatisticalTestResult;
import org.jlinalg.stats.StatisticalTests;

StatisticalTestResult welch = StatisticalTests.studentT(
    new double[] {4.1, 5.2, 5.8, 6.0},
    new double[] {3.7, 4.0, 4.9, 5.1},
    0.0, false, Alternative.TWO_SIDED, 0.95);

System.out.println(welch.statistic());
System.out.println(welch.parameters().get("df"));
System.out.println(welch.pValue());
System.out.println(welch.confidenceInterval().orElseThrow());
```

Every result exposes:

- `method()` and `statisticName()`;
- `statistic()` and named `parameters()` such as degrees of freedom;
- `pValue()` and `pValueMethod()`;
- named `estimates()` and `nullValues()`;
- `alternative()` and an optional `confidenceInterval()`.

## Available families

The facade covers Ansari-Bradley, Bartlett, binomial, Pearson chi-square,
correlation, Fisher exact, Fligner-Killeen, Friedman, Kolmogorov-Smirnov,
Kruskal-Wallis, Mantel-Haenszel, McNemar, Mood, one-way ANOVA, Poisson, ANOVA
power, proportion, proportion trend, Quade, Shapiro-Wilk, Student t, variance,
and Wilcoxon signed-rank/rank-sum tests.

Correlation supports Pearson, Kendall, and Spearman methods. The default
Kendall overload uses the exact permutation distribution below 50 untied
pairs and the normal approximation otherwise. An explicit `exact=true`
request uses a normalized, tail-symmetric recurrence; calculations exceeding
50 million recurrence cells fail with an instruction to use `exact=false`.
No silent approximation is substituted for that explicit request.
Untied Spearman uses exact enumeration through nine pairs, AS 89 through
1,289 pairs, and Student t from 1,290 pairs, matching base R's cutoffs.
Tied rank tests use the corresponding asymptotic test, and `pValueMethod()`
records that choice. Positive and negative zero count as ties.

Version 0.3.0 centers and rescales Pearson/ANOVA calculations and rescales
Student/Welch/variance-test inputs to avoid overflow or underflow when units
change. It preserves small exact Kendall
upper-tail probabilities, and replaces quadratic Kendall pair counting with
an O(n log n) tie-aware algorithm. Equal-variance ANOVA permits an individual
constant group when the pooled within-group variance is positive; Welch's
test still requires positive variance in every group.

Reproduce the independent base-R accuracy fixtures and paired timing with
`Rscript src/test/R/core-statistics-audit.R` and
`./gradlew benchmarkCoreStatisticsAudit`. Details and caveats are in the
[v0.3.0 audit report](../release-0.3.0-audit.md).

```java
StatisticalTestResult fisher = StatisticalTests.fisherExact(new long[][] {
    {1, 9},
    {11, 3}
});

StatisticalTestResult friedman = StatisticalTests.friedman(new double[][] {
    {5.4, 5.5, 5.55},
    {5.85, 5.7, 5.75},
    {5.2, 5.6, 5.5}
});
```

Two-by-two Fisher tests return the conditional maximum-likelihood odds ratio
and exact confidence interval. Larger rectangular tables use bounded
Fisher-Freeman-Halton enumeration; the overload accepting `maximumTables`
lets callers choose the resource bound. Mantel-Haenszel accepts
`tables[stratum][row][column]`, returning common-odds-ratio inference for
2-by-2 strata and the generalized Cochran-Mantel-Haenszel statistic otherwise.

`powerAnova` mirrors `power.anova.test`: exactly one boxed argument is `null`
and is solved from the remaining values.

```java
AnovaPowerResult required = StatisticalTests.powerAnova(
    4.0, null, 1.0, 3.0, 0.05, 0.80);
System.out.println(required.observationsPerGroup());
```

All sample APIs reject `NaN` and infinity. Compact complete cases before the
call when omission is intended. Chi-square approximations and continuity
corrections are explicit in `pValueMethod()`; inspect that field rather than
assuming every returned p-value is exact.
