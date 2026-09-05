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

Correlation supports Pearson, Kendall, and Spearman methods. Untied Kendall
tests use the exact permutation distribution when requested. Untied Spearman
tests use exact enumeration through nine pairs and base R's AS 89 tail
approximation for larger samples. Tied rank tests use the corresponding
asymptotic test, and `pValueMethod()` records that choice.

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
