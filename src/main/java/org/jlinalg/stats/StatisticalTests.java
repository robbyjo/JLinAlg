/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.stats;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import jdistlib.Beta;
import jdistlib.ChiSquare;
import jdistlib.F;
import jdistlib.NonCentralF;
import jdistlib.Normal;
import jdistlib.T;
import jdistlib.disttest.DistributionTest;
import jdistlib.disttest.NormalityTest;
import jdistlib.disttest.TestKind;
import jdistlib.generic.GenericDistribution;
import jdistlib.math.MathFunctions;

/**
 * Basic hypothesis tests with a uniform Java API modeled after base R's
 * {@code stats} tests. Existing JDistlib implementations are delegated to;
 * this class supplies the missing tests and consistent result metadata.
 * Inputs must be finite: callers that want complete-case omission should
 * compact their data explicitly before invoking a test.
 */
public final class StatisticalTests {
    private static final double DEFAULT_CONFIDENCE_LEVEL = 0.95;
    private static final double FISHER_RELATIVE_ERROR = 1e-7;

    private StatisticalTests() { }

    /** Ansari-Bradley two-sample scale test, using JDistlib. */
    public static StatisticalTestResult ansariBradley(
            double[] first, double[] second) {
        return ansariBradley(first, second, Alternative.TWO_SIDED, false);
    }

    /** Ansari-Bradley two-sample scale test, using JDistlib. */
    public static StatisticalTestResult ansariBradley(
            double[] first, double[] second, Alternative alternative,
            boolean forceExact) {
        double[] x = StatisticsSupport.sample(first, 1, "first sample");
        double[] y = StatisticsSupport.sample(second, 1, "second sample");
        double[] raw = DistributionTest.ansari_bradley_test(
            x, y, forceExact, kind(alternative));
        boolean exact = forceExact || (x.length < 50 && y.length < 50
            && !hasTies(concatenate(x, y)));
        return result("Ansari-Bradley test", "AB", raw[0], Map.of(), raw[1],
            Map.of(), Map.of("ratio of scales", 1.0), Optional.empty(),
            alternative, exact ? "exact" : "normal approximation");
    }

    /** Bartlett test of homogeneity of variances, using JDistlib. */
    public static StatisticalTestResult bartlett(double[] values, int[] groups) {
        double[] x = StatisticsSupport.sample(values, 2, "values");
        int[] labels = StatisticsSupport.groups(groups, x.length);
        int count = groupCount(labels, 2);
        double[] raw = DistributionTest.bartlett_test(x, labels);
        return result("Bartlett test of homogeneity of variances", "Bartlett's K-squared",
            raw[0], Map.of("df", count - 1.0), raw[1], Map.of(), Map.of(),
            Optional.empty(), Alternative.TWO_SIDED, "chi-square approximation");
    }

    /** Exact binomial test, using JDistlib, with a Clopper-Pearson interval. */
    public static StatisticalTestResult binomial(
            int successes, int trials, double probability) {
        return binomial(successes, trials, probability,
            Alternative.TWO_SIDED, DEFAULT_CONFIDENCE_LEVEL);
    }

    /** Exact binomial test, using JDistlib, with a Clopper-Pearson interval. */
    public static StatisticalTestResult binomial(
            int successes, int trials, double probability,
            Alternative alternative, double confidenceLevel) {
        if (trials <= 0 || successes < 0 || successes > trials
                || !Double.isFinite(probability)
                || probability < 0.0 || probability > 1.0) {
            throw new IllegalArgumentException("invalid binomial counts or null probability");
        }
        StatisticsSupport.probability(confidenceLevel, "confidence level");
        double[] raw = DistributionTest.binomial_test(
            successes, trials, probability, kind(alternative));
        double alpha = alternative == Alternative.TWO_SIDED
            ? (1.0 - confidenceLevel) / 2.0 : 1.0 - confidenceLevel;
        double lower = successes == 0 ? 0.0
            : Beta.quantile(alpha, successes, trials - successes + 1.0, true, false);
        double upper = successes == trials ? 1.0
            : Beta.quantile(1.0 - alpha, successes + 1.0,
                trials - successes, true, false);
        if (alternative == Alternative.LESS) lower = 0.0;
        if (alternative == Alternative.GREATER) upper = 1.0;
        return result("Exact binomial test", "number of successes", raw[0],
            Map.of("trials", (double) trials), raw[1],
            Map.of("probability of success", successes / (double) trials),
            Map.of("probability of success", probability),
            Optional.of(new ConfidenceInterval(lower, upper, confidenceLevel)),
            alternative, "exact binomial");
    }

    /** Pearson chi-square goodness-of-fit test, using JDistlib. */
    public static StatisticalTestResult chiSquareGoodnessOfFit(
            long[] observed, double[] probabilities) {
        return chiSquareGoodnessOfFit(observed, probabilities, 0);
    }

    /** Pearson chi-square goodness-of-fit test, using JDistlib. */
    public static StatisticalTestResult chiSquareGoodnessOfFit(
            long[] observed, double[] probabilities, int estimatedParameters) {
        double[] raw = DistributionTest.chi_square_goodness_of_fit_test(
            observed == null ? null : observed.clone(),
            probabilities == null ? null : probabilities.clone(), estimatedParameters);
        return result("Chi-squared test for given probabilities", "X-squared", raw[0],
            Map.of("df", raw[2]), raw[1], Map.of(), Map.of(), Optional.empty(),
            Alternative.TWO_SIDED, "chi-square approximation");
    }

    /** Pearson chi-square independence test without continuity correction. */
    public static StatisticalTestResult chiSquareIndependence(long[][] counts) {
        long[][] copy = copy(counts);
        double[] raw = DistributionTest.chi_square_independence_test(copy);
        return result("Pearson's Chi-squared test", "X-squared", raw[0],
            Map.of("df", raw[2]), raw[1], Map.of(), Map.of(), Optional.empty(),
            Alternative.TWO_SIDED, "chi-square approximation without correction");
    }

    /** Pearson product-moment correlation test. */
    public static StatisticalTestResult correlation(double[] first, double[] second) {
        return correlation(first, second, CorrelationMethod.PEARSON,
            Alternative.TWO_SIDED, DEFAULT_CONFIDENCE_LEVEL, true, false);
    }

    /** Pearson, Kendall, or Spearman correlation test. */
    public static StatisticalTestResult correlation(
            double[] first, double[] second, CorrelationMethod method,
            Alternative alternative) {
        return correlation(first, second, method, alternative,
            DEFAULT_CONFIDENCE_LEVEL,
            method != CorrelationMethod.KENDALL || (first != null && first.length < 50), false);
    }

    /**
     * Pearson, Kendall, or Spearman correlation test.
     * Rank tests use exact/AS 89 probabilities when requested and there are no
     * ties (Spearman switches to Student t at n >= 1290, as in R).
     * Explicit Kendall enumeration is bounded to 50 million recurrence cells;
     * requests above that bound fail rather than silently approximate.
     * Continuity correction applies only to asymptotic rank tests.
     */
    public static StatisticalTestResult correlation(
            double[] first, double[] second, CorrelationMethod method,
            Alternative alternative, double confidenceLevel,
            boolean exact, boolean continuityCorrection) {
        StatisticsSupport.paired(first, second, 2);
        double[] x = first.clone();
        double[] y = second.clone();
        Objects.requireNonNull(method, "method");
        Objects.requireNonNull(alternative, "alternative");
        StatisticsSupport.probability(confidenceLevel, "confidence level");
        return switch (method) {
            case PEARSON -> pearson(x, y, alternative, confidenceLevel);
            case KENDALL -> kendall(x, y, alternative, exact, continuityCorrection);
            case SPEARMAN -> spearman(x, y, alternative, exact, continuityCorrection);
        };
    }

    private static StatisticalTestResult pearson(
            double[] first, double[] second, Alternative alternative,
            double confidenceLevel) {
        if (first.length < 3) {
            throw new IllegalArgumentException("Pearson correlation needs at least 3 pairs");
        }
        double correlation = StatisticsSupport.correlation(first, second);
        double degrees = first.length - 2.0;
        double statistic = Math.sqrt(degrees) * correlation
            / Math.sqrt(1.0 - correlation * correlation);
        double pValue = tPValue(statistic, degrees, alternative);
        Optional<ConfidenceInterval> interval = Optional.empty();
        if (first.length > 3 && Double.isFinite(correlation)) {
            double z = atanh(correlation);
            double sigma = 1.0 / Math.sqrt(first.length - 3.0);
            double critical = Normal.quantile(
                alternative == Alternative.TWO_SIDED
                    ? (1.0 + confidenceLevel) / 2.0 : confidenceLevel,
                0.0, 1.0, true, false);
            double lower = alternative == Alternative.LESS
                ? Double.NEGATIVE_INFINITY : z - sigma * critical;
            double upper = alternative == Alternative.GREATER
                ? Double.POSITIVE_INFINITY : z + sigma * critical;
            interval = Optional.of(new ConfidenceInterval(
                Math.tanh(lower), Math.tanh(upper), confidenceLevel));
        }
        return result("Pearson's product-moment correlation", "t", statistic,
            Map.of("df", degrees), pValue, Map.of("correlation", correlation),
            Map.of("correlation", 0.0), interval, alternative, "Student t");
    }

    private static StatisticalTestResult kendall(
            double[] first, double[] second, Alternative alternative,
            boolean exact, boolean correction) {
        int size = first.length;
        long[] pairs = StatisticsSupport.kendallPairs(first, second);
        long concordant = pairs[0];
        long discordant = pairs[1];
        Map<Double, Integer> xCounts = StatisticsSupport.tieCounts(first);
        Map<Double, Integer> yCounts = StatisticsSupport.tieCounts(second);
        double pairCount = size * (size - 1.0) / 2.0;
        double firstTies = tiePairs(xCounts);
        double secondTies = tiePairs(yCounts);
        double denominator = Math.sqrt(
            (pairCount - firstTies) * (pairCount - secondTies));
        double tau = (concordant - discordant) / denominator;
        boolean ties = !xCounts.isEmpty() || !yCounts.isEmpty();
        if (!Double.isFinite(tau)) {
            return result("Kendall's rank correlation tau", "z", Double.NaN,
                Map.of(), Double.NaN, Map.of("tau", Double.NaN),
                Map.of("tau", 0.0), Optional.empty(), alternative,
                "undefined for a constant sample");
        }
        if (exact && !ties) {
            long total = size * (size - 1L) / 2;
            long observed = concordant;
            double lower = kendallCumulative(observed, size);
            // Symmetry avoids subtracting a tiny upper tail from one.
            double upper = kendallCumulative(total - observed, size);
            double pValue = switch (alternative) {
                case LESS -> lower;
                case GREATER -> upper;
                case TWO_SIDED -> Math.min(1.0, 2.0
                    * (observed > total / 2.0 ? upper : lower));
            };
            return result("Kendall's rank correlation tau", "T", observed,
                Map.of(), pValue, Map.of("tau", tau), Map.of("tau", 0.0),
                Optional.empty(), alternative, "exact permutation");
        }
        double score = tau * Math.sqrt(
            (pairCount - firstTies) * (pairCount - secondTies));
        double v0 = size * (size - 1.0) * (2.0 * size + 5.0);
        double vt = tieVarianceTerm(xCounts);
        double vu = tieVarianceTerm(yCounts);
        double v1 = tieSecondTerm(xCounts) * tieSecondTerm(yCounts);
        double v2 = tieThirdTerm(xCounts) * tieThirdTerm(yCounts);
        double variance = (v0 - vt - vu) / 18.0
            + v1 / (2.0 * size * (size - 1.0))
            + (size > 2 ? v2 / (9.0 * size * (size - 1.0) * (size - 2.0)) : 0.0);
        if (correction && score != 0.0) score -= Math.copySign(1.0, score);
        double statistic = score / Math.sqrt(variance);
        return result("Kendall's rank correlation tau", "z", statistic,
            Map.of(), StatisticsSupport.normalPValue(statistic, alternative),
            Map.of("tau", tau), Map.of("tau", 0.0), Optional.empty(),
            alternative, correction ? "normal approximation with continuity correction"
                : "normal approximation");
    }

    private static StatisticalTestResult spearman(
            double[] first, double[] second, Alternative alternative,
            boolean exact, boolean correction) {
        double[] firstRanks = StatisticsSupport.ranks(first);
        double[] secondRanks = StatisticsSupport.ranks(second);
        double rho = StatisticsSupport.correlation(firstRanks, secondRanks);
        if (!Double.isFinite(rho)) {
            return result("Spearman's rank correlation rho", "S", Double.NaN,
                Map.of(), Double.NaN, Map.of("rho", Double.NaN), Map.of("rho", 0.0),
                Optional.empty(), alternative, "undefined for a constant sample");
        }
        int size = first.length;
        double q = (Math.pow(size, 3.0) - size) * (1.0 - rho) / 6.0;
        boolean ties = hasTies(first) || hasTies(second);
        if (exact && !ties && size < 1290) {
            double lower = spearmanTail(q, size, true);
            double upper = spearmanTail(q, size, false);
            double pValue = switch (alternative) {
                case GREATER -> lower;
                case LESS -> upper;
                case TWO_SIDED -> Math.min(1.0, 2.0 * (q
                    > (Math.pow(size, 3.0) - size) / 6.0 ? upper : lower));
            };
            return result("Spearman's rank correlation rho", "S", q, Map.of(),
                pValue, Map.of("rho", rho), Map.of("rho", 0.0), Optional.empty(),
                alternative, size <= 9 ? "exact permutation" : "AS 89 approximation");
        }
        double denominator = size * (size * (double) size - 1.0) / 6.0;
        if (correction) denominator += 1.0;
        double adjustedRho = 1.0 - q / denominator;
        double statistic = adjustedRho / Math.sqrt(
            (1.0 - adjustedRho * adjustedRho) / (size - 2.0));
        double pValue = tPValue(statistic, size - 2.0, alternative);
        return result("Spearman's rank correlation rho", "S", q,
            Map.of("df", size - 2.0), pValue, Map.of("rho", rho),
            Map.of("rho", 0.0), Optional.empty(), alternative,
            correction ? "asymptotic Student t with continuity correction"
                : "asymptotic Student t");
    }

    /** Fisher's exact conditional test for a rectangular contingency table. */
    public static StatisticalTestResult fisherExact(long[][] table) {
        if (StatisticsSupport.total(table, 2, 2) == 0L) {
            throw new IllegalArgumentException("contingency table total must be positive");
        }
        if (table.length != 2 || table[0].length != 2) {
            return fisherExact(table, 1_000_000L);
        }
        return fisherExact(table, 1.0, Alternative.TWO_SIDED,
            DEFAULT_CONFIDENCE_LEVEL);
    }

    /**
     * Fisher-Freeman-Halton exact test for a general rectangular table.
     * Enumeration stops at {@code maximumTables} to keep the calculation
     * explicitly bounded; use a chi-square test for tables beyond that bound.
     */
    public static StatisticalTestResult fisherExact(
            long[][] table, long maximumTables) {
        if (StatisticsSupport.total(table, 2, 2) == 0L) {
            throw new IllegalArgumentException("contingency table total must be positive");
        }
        if (maximumTables < 1L) {
            throw new IllegalArgumentException("maximum table count must be positive");
        }
        if (table.length == 2 && table[0].length == 2) {
            return fisherExact(table);
        }
        int rows = table.length;
        int columns = table[0].length;
        int[] rowSums = new int[rows];
        int[] columnSums = new int[columns];
        double observedLogWeight = 0.0;
        for (int row = 0; row < rows; row++) {
            for (int column = 0; column < columns; column++) {
                int count = Math.toIntExact(table[row][column]);
                rowSums[row] = Math.addExact(rowSums[row], count);
                columnSums[column] = Math.addExact(columnSums[column], count);
                observedLogWeight -= MathFunctions.lgammafn(count + 1.0);
            }
        }
        ExactTableEnumerator enumerator = new ExactTableEnumerator(
            rowSums, columnSums, observedLogWeight, maximumTables);
        double pValue = enumerator.enumerate();
        return result("Fisher's Exact Test for Count Data",
            "observed conditional probability", enumerator.observedProbability(),
            Map.of("enumerated tables", (double) enumerator.tableCount()), pValue,
            Map.of(), Map.of(), Optional.empty(), Alternative.TWO_SIDED,
            "exact conditional minimum-likelihood enumeration");
    }

    /** Fisher's exact conditional test for a 2-by-2 table. */
    public static StatisticalTestResult fisherExact(
            long[][] table, double nullOddsRatio, Alternative alternative,
            double confidenceLevel) {
        requireTwoByTwo(table);
        if (Double.isNaN(nullOddsRatio) || nullOddsRatio < 0.0) {
            throw new IllegalArgumentException("null odds ratio must be in [0, infinity]");
        }
        StatisticsSupport.probability(confidenceLevel, "confidence level");
        Objects.requireNonNull(alternative, "alternative");
        ConditionalTableDistribution distribution =
            ConditionalTableDistribution.from(table);
        int observed = Math.toIntExact(table[0][0]);
        double[] probabilities = distribution.probabilities(nullOddsRatio);
        double lowerTail = distribution.tail(probabilities, observed, false);
        double upperTail = distribution.tail(probabilities, observed, true);
        double pValue = switch (alternative) {
            case LESS -> lowerTail;
            case GREATER -> upperTail;
            case TWO_SIDED -> distribution.minimumLikelihoodTail(
                probabilities, observed, FISHER_RELATIVE_ERROR);
        };
        double estimate = distribution.conditionalMle(observed);
        double alpha = alternative == Alternative.TWO_SIDED
            ? (1.0 - confidenceLevel) / 2.0 : 1.0 - confidenceLevel;
        double lower = alternative == Alternative.LESS ? 0.0
            : distribution.lowerConfidenceLimit(observed, alpha);
        double upper = alternative == Alternative.GREATER ? Double.POSITIVE_INFINITY
            : distribution.upperConfidenceLimit(observed, alpha);
        return result("Fisher's Exact Test for Count Data", "count", observed,
            Map.of(), pValue, Map.of("odds ratio", estimate),
            Map.of("odds ratio", nullOddsRatio),
            Optional.of(new ConfidenceInterval(lower, upper, confidenceLevel)),
            alternative, "exact conditional minimum-likelihood");
    }

    /** Fligner-Killeen test of homogeneity of variances, using JDistlib. */
    public static StatisticalTestResult flignerKilleen(
            double[] values, int[] groups) {
        double[] x = StatisticsSupport.sample(values, 2, "values");
        int[] labels = StatisticsSupport.groups(groups, x.length);
        int count = groupCount(labels, 2);
        double[] raw = DistributionTest.fligner_test(x, labels);
        return result("Fligner-Killeen test of homogeneity of variances",
            "Fligner-Killeen:med chi-squared", raw[0], Map.of("df", count - 1.0),
            raw[1], Map.of(), Map.of(), Optional.empty(), Alternative.TWO_SIDED,
            "chi-square approximation");
    }

    /** Friedman rank-sum test for a complete blocks-by-groups matrix. */
    public static StatisticalTestResult friedman(double[][] blocksByGroups) {
        StatisticsSupport.validateMatrix(blocksByGroups, 2, 2, "blocked data");
        int blocks = blocksByGroups.length;
        int groups = blocksByGroups[0].length;
        double[] rankSums = new double[groups];
        double tieAdjustment = 0.0;
        for (double[] block : blocksByGroups) {
            double[] ranks = StatisticsSupport.ranks(block);
            for (int group = 0; group < groups; group++) rankSums[group] += ranks[group];
            for (int count : StatisticsSupport.tieCounts(block).values()) {
                tieAdjustment += count * (double) count * count - count;
            }
        }
        double center = blocks * (groups + 1.0) / 2.0;
        double numerator = 0.0;
        for (double sum : rankSums) {
            double difference = sum - center;
            numerator += difference * difference;
        }
        double denominator = blocks * groups * (groups + 1.0)
            - tieAdjustment / (groups - 1.0);
        double statistic = 12.0 * numerator / denominator;
        double degrees = groups - 1.0;
        double pValue = ChiSquare.cumulative(statistic, degrees, false, false);
        return result("Friedman rank sum test", "Friedman chi-squared", statistic,
            Map.of("df", degrees), pValue, Map.of(), Map.of(), Optional.empty(),
            Alternative.TWO_SIDED, "chi-square approximation");
    }

    /** Two-sample Kolmogorov-Smirnov test, using JDistlib. */
    public static StatisticalTestResult kolmogorovSmirnov(
            double[] first, double[] second) {
        return kolmogorovSmirnov(first, second, Alternative.TWO_SIDED, true);
    }

    /** Two-sample Kolmogorov-Smirnov test, using JDistlib. */
    public static StatisticalTestResult kolmogorovSmirnov(
            double[] first, double[] second, Alternative alternative,
            boolean exact) {
        double[] x = StatisticsSupport.sample(first, 1, "first sample");
        double[] y = StatisticsSupport.sample(second, 1, "second sample");
        double[] raw = DistributionTest.kolmogorov_smirnov_test(
            x, y, kind(alternative), exact);
        boolean effectiveExact = exact && !hasTies(concatenate(x, y));
        return result("Two-sample Kolmogorov-Smirnov test", "D", raw[0], Map.of(),
            raw[1], Map.of(), Map.of(), Optional.empty(), alternative,
            effectiveExact ? "exact" : "asymptotic");
    }

    /** One-sample Kolmogorov-Smirnov test against a fixed distribution. */
    public static StatisticalTestResult kolmogorovSmirnov(
            double[] sample, GenericDistribution distribution,
            Alternative alternative, boolean exact) {
        double[] x = StatisticsSupport.sample(sample, 1, "sample");
        Objects.requireNonNull(distribution, "distribution");
        double[] raw = DistributionTest.kolmogorov_smirnov_test(
            x, distribution, kind(alternative), exact);
        boolean effectiveExact = exact && !hasTies(x);
        return result("One-sample Kolmogorov-Smirnov test", "D", raw[0], Map.of(),
            raw[1], Map.of(), Map.of(), Optional.empty(), alternative,
            effectiveExact ? "exact" : "asymptotic");
    }

    /** Kruskal-Wallis rank-sum test, using JDistlib. */
    public static StatisticalTestResult kruskalWallis(
            double[] values, int[] groups) {
        double[] x = StatisticsSupport.sample(values, 2, "values");
        int[] labels = StatisticsSupport.groups(groups, x.length);
        int count = groupCount(labels, 2);
        double[] raw = DistributionTest.kruskal_wallis_test(x, labels);
        return result("Kruskal-Wallis rank sum test", "Kruskal-Wallis chi-squared",
            raw[0], Map.of("df", count - 1.0), raw[1], Map.of(), Map.of(),
            Optional.empty(), Alternative.TWO_SIDED, "chi-square approximation");
    }

    /**
     * Classical Mantel-Haenszel test for 2-by-2 tables across strata.
     * The array layout is {@code tables[stratum][row][column]}.
     */
    public static StatisticalTestResult mantelHaenszel(long[][][] tables) {
        return mantelHaenszel(tables, Alternative.TWO_SIDED, true,
            DEFAULT_CONFIDENCE_LEVEL);
    }

    /** Classical Mantel-Haenszel test for 2-by-2 tables across strata. */
    public static StatisticalTestResult mantelHaenszel(
            long[][][] tables, Alternative alternative, boolean correction,
            double confidenceLevel) {
        validateStrata(tables);
        Objects.requireNonNull(alternative, "alternative");
        StatisticsSupport.probability(confidenceLevel, "confidence level");
        if (tables[0].length != 2 || tables[0][0].length != 2) {
            if (alternative != Alternative.TWO_SIDED) {
                throw new IllegalArgumentException(
                    "general Cochran-Mantel-Haenszel tests are omnibus and two-sided");
            }
            return generalMantelHaenszel(tables);
        }
        double delta = 0.0;
        double variance = 0.0;
        double diagonal = 0.0;
        double offDiagonal = 0.0;
        double varianceOne = 0.0;
        double varianceTwo = 0.0;
        double varianceThree = 0.0;
        for (long[][] table : tables) {
            double a = table[0][0];
            double b = table[0][1];
            double c = table[1][0];
            double d = table[1][1];
            double total = a + b + c + d;
            double rowOne = a + b;
            double rowTwo = c + d;
            double columnOne = a + c;
            double columnTwo = b + d;
            delta += a - rowOne * columnOne / total;
            variance += rowOne * rowTwo * columnOne * columnTwo
                / (total * total * (total - 1.0));
            diagonal += a * d / total;
            offDiagonal += b * c / total;
            varianceOne += (a + d) * a * d / (total * total);
            varianceTwo += ((a + d) * b * c + (b + c) * a * d)
                / (total * total);
            varianceThree += (b + c) * b * c / (total * total);
        }
        double yates = correction && Math.abs(delta) >= 0.5 ? 0.5 : 0.0;
        double statistic = Math.pow(Math.abs(delta) - yates, 2.0) / variance;
        double signedZ = Math.copySign(Math.sqrt(statistic), delta);
        double pValue = alternative == Alternative.TWO_SIDED
            ? ChiSquare.cumulative(statistic, 1.0, false, false)
            : StatisticsSupport.normalPValue(signedZ, alternative);
        double estimate = diagonal / offDiagonal;
        double standardError = Math.sqrt(
            varianceOne / (2.0 * diagonal * diagonal)
            + varianceTwo / (2.0 * diagonal * offDiagonal)
            + varianceThree / (2.0 * offDiagonal * offDiagonal));
        Optional<ConfidenceInterval> interval = Optional.empty();
        if (estimate > 0.0 && Double.isFinite(estimate)
                && Double.isFinite(standardError)) {
            double critical = Normal.quantile(
                alternative == Alternative.TWO_SIDED
                    ? (1.0 + confidenceLevel) / 2.0 : confidenceLevel,
                0.0, 1.0, true, false);
            double lower = alternative == Alternative.LESS ? 0.0
                : estimate * Math.exp(-critical * standardError);
            double upper = alternative == Alternative.GREATER ? Double.POSITIVE_INFINITY
                : estimate * Math.exp(critical * standardError);
            interval = Optional.of(new ConfidenceInterval(lower, upper, confidenceLevel));
        }
        String pMethod = correction
            ? "chi-square approximation with continuity correction"
            : "chi-square approximation without continuity correction";
        return result("Mantel-Haenszel chi-squared test", "Mantel-Haenszel X-squared",
            statistic, Map.of("df", 1.0), pValue,
            Map.of("common odds ratio", estimate),
            Map.of("common odds ratio", 1.0),
            interval, alternative, alternative == Alternative.TWO_SIDED ? pMethod
                : pMethod.replace("chi-square", "normal"));
    }

    /** McNemar chi-square test for a square paired-count table. */
    public static StatisticalTestResult mcnemar(long[][] table) {
        return mcnemar(table, true);
    }

    /** McNemar chi-square test for a square paired-count table. */
    public static StatisticalTestResult mcnemar(long[][] table, boolean correction) {
        StatisticsSupport.total(table, 2, 2);
        int dimension = table.length;
        if (table[0].length != dimension) {
            throw new IllegalArgumentException("McNemar table must be square");
        }
        double statistic = 0.0;
        boolean undefined = false;
        for (int row = 0; row < dimension; row++) {
            for (int column = row + 1; column < dimension; column++) {
                double difference = Math.abs((double) table[row][column]
                    - table[column][row]);
                if (correction && dimension == 2 && difference != 0.0) difference -= 1.0;
                double total = (double) table[row][column] + table[column][row];
                if (total == 0.0) undefined = true;
                else statistic += difference * difference / total;
            }
        }
        if (undefined) statistic = Double.NaN;
        double degrees = dimension * (dimension - 1.0) / 2.0;
        double pValue = Double.isNaN(statistic) ? Double.NaN
            : ChiSquare.cumulative(statistic, degrees, false, false);
        return result("McNemar's Chi-squared test", "McNemar's chi-squared", statistic,
            Map.of("df", degrees), pValue, Map.of(), Map.of(), Optional.empty(),
            Alternative.TWO_SIDED, correction && dimension == 2
                ? "chi-square approximation with continuity correction"
                : "chi-square approximation without continuity correction");
    }

    /** Mood two-sample scale test, using JDistlib. */
    public static StatisticalTestResult mood(double[] first, double[] second) {
        return mood(first, second, Alternative.TWO_SIDED);
    }

    /** Mood two-sample scale test, using JDistlib. */
    public static StatisticalTestResult mood(
            double[] first, double[] second, Alternative alternative) {
        double[] x = StatisticsSupport.sample(first, 1, "first sample");
        double[] y = StatisticsSupport.sample(second, 1, "second sample");
        double[] raw = DistributionTest.mood_test(x, y, kind(alternative));
        return result("Mood two-sample test of scale", "Z", raw[0], Map.of(), raw[1],
            Map.of(), Map.of("ratio of scales", 1.0), Optional.empty(), alternative,
            "normal approximation");
    }

    /** Welch one-way analysis of means. Each row is one independent group. */
    public static StatisticalTestResult oneWayAnova(double[][] groups) {
        return oneWayAnova(groups, false);
    }

    /** One-way ANOVA, optionally assuming equal group variances. */
    public static StatisticalTestResult oneWayAnova(
            double[][] groups, boolean equalVariances) {
        validateGroups(groups);
        // The F statistic is invariant to a common location and scale. Work in
        // bounded coordinates so squares/weights cannot overflow or underflow.
        double minimum = groups[0][0], maximum = minimum;
        for (double[] group : groups) for (double value : group) {
            minimum = Math.min(minimum, value);
            maximum = Math.max(maximum, value);
        }
        double range = maximum - minimum;
        double center = Double.isFinite(range) ? minimum + range / 2.0
            : minimum / 2.0 + maximum / 2.0;
        double scale = Math.max(Math.abs(minimum - center), Math.abs(maximum - center));
        if (!(scale > 0.0)) throw new IllegalArgumentException("data have no variation");
        double[][] scaled = new double[groups.length][];
        for (int g = 0; g < groups.length; g++) {
            scaled[g] = new double[groups[g].length];
            for (int i = 0; i < scaled[g].length; i++) scaled[g][i] = (groups[g][i] - center) / scale;
        }
        groups = scaled;
        int groupCount = groups.length;
        int totalCount = 0;
        double[] means = new double[groupCount];
        double[] variances = new double[groupCount];
        double[] weights = new double[groupCount];
        double grandSum = 0.0;
        for (int index = 0; index < groupCount; index++) {
            means[index] = StatisticsSupport.mean(groups[index]);
            variances[index] = StatisticsSupport.variance(groups[index]);
            if ((!equalVariances && !(variances[index] > 0.0))
                    || !Double.isFinite(variances[index])) {
                throw new IllegalArgumentException("each group must have positive finite variance");
            }
            weights[index] = groups[index].length / variances[index];
            totalCount += groups[index].length;
            grandSum += means[index] * groups[index].length;
        }
        double statistic;
        double denominatorDegrees;
        String method = "One-way analysis of means";
        if (equalVariances) {
            double grandMean = grandSum / totalCount;
            double between = 0.0;
            double within = 0.0;
            for (int index = 0; index < groupCount; index++) {
                double difference = means[index] - grandMean;
                between += groups[index].length * difference * difference;
                within += (groups[index].length - 1.0) * variances[index];
            }
            denominatorDegrees = totalCount - groupCount;
            if (!(within > 0.0)) throw new IllegalArgumentException("pooled within-group variance must be positive");
            statistic = (between / (groupCount - 1.0))
                / (within / denominatorDegrees);
        } else {
            double sumWeights = 0.0;
            double weightedMeanNumerator = 0.0;
            for (int index = 0; index < groupCount; index++) {
                sumWeights += weights[index];
                weightedMeanNumerator += weights[index] * means[index];
            }
            double weightedMean = weightedMeanNumerator / sumWeights;
            double weightedBetween = 0.0;
            double tmp = 0.0;
            for (int index = 0; index < groupCount; index++) {
                double difference = means[index] - weightedMean;
                weightedBetween += weights[index] * difference * difference;
                double complement = 1.0 - weights[index] / sumWeights;
                tmp += complement * complement / (groups[index].length - 1.0);
            }
            tmp /= groupCount * (double) groupCount - 1.0;
            statistic = weightedBetween
                / ((groupCount - 1.0) * (1.0 + 2.0 * (groupCount - 2.0) * tmp));
            denominatorDegrees = 1.0 / (3.0 * tmp);
            method += " (not assuming equal variances)";
        }
        double numeratorDegrees = groupCount - 1.0;
        double pValue = F.cumulative(
            statistic, numeratorDegrees, denominatorDegrees, false, false);
        return result(method, "F", statistic,
            orderedMap("num df", numeratorDegrees, "denom df", denominatorDegrees),
            pValue, Map.of(), Map.of(), Optional.empty(), Alternative.TWO_SIDED,
            "F distribution");
    }

    /** Exact one-sample Poisson rate test, using JDistlib. */
    public static StatisticalTestResult poisson(
            int events, double exposure, double rate) {
        return poisson(events, exposure, rate, Alternative.TWO_SIDED);
    }

    /** Exact one-sample Poisson rate test, using JDistlib. */
    public static StatisticalTestResult poisson(
            int events, double exposure, double rate, Alternative alternative) {
        if (events < 0 || !(exposure > 0.0) || rate < 0.0
                || !Double.isFinite(exposure) || !Double.isFinite(rate)) {
            throw new IllegalArgumentException("invalid Poisson count, exposure, or rate");
        }
        double[] raw = DistributionTest.poisson_test(
            events, exposure, rate, kind(alternative));
        return result("Exact Poisson test", "number of events", raw[0],
            Map.of("exposure", exposure), raw[1],
            Map.of("event rate", events / exposure), Map.of("event rate", rate),
            Optional.empty(), alternative, "exact Poisson");
    }

    /** Exact comparison of two Poisson rates, using JDistlib. */
    public static StatisticalTestResult poisson(
            int firstEvents, int secondEvents, double firstExposure,
            double secondExposure, double rateRatio, Alternative alternative) {
        if (firstEvents < 0 || secondEvents < 0 || !(firstExposure > 0.0)
                || !(secondExposure > 0.0) || rateRatio < 0.0
                || !Double.isFinite(firstExposure) || !Double.isFinite(secondExposure)
                || !Double.isFinite(rateRatio)) {
            throw new IllegalArgumentException("invalid Poisson comparison input");
        }
        double[] raw = DistributionTest.poisson_test(firstEvents, secondEvents,
            firstExposure, secondExposure, rateRatio, kind(alternative));
        double estimate = (firstEvents / firstExposure) / (secondEvents / secondExposure);
        return result("Comparison of Poisson rates", "first event count", raw[0],
            Map.of("total events", (double) firstEvents + secondEvents), raw[1],
            Map.of("rate ratio", estimate), Map.of("rate ratio", rateRatio),
            Optional.empty(), alternative, "exact conditional binomial");
    }

    /** Balanced one-way ANOVA power with all inputs supplied except power. */
    public static AnovaPowerResult powerAnova(
            int groups, double observationsPerGroup, double betweenVariance,
            double withinVariance, double significanceLevel) {
        return powerAnova((double) groups, observationsPerGroup, betweenVariance,
            withinVariance, significanceLevel, null);
    }

    /**
     * Base-R-style balanced one-way ANOVA power solver. Exactly one argument
     * must be {@code null}; that quantity is solved from the other five.
     */
    public static AnovaPowerResult powerAnova(
            Double groups, Double observationsPerGroup, Double betweenVariance,
            Double withinVariance, Double significanceLevel, Double power) {
        Object[] values = {groups, observationsPerGroup, betweenVariance,
            withinVariance, significanceLevel, power};
        int missing = 0;
        for (Object value : values) if (value == null) missing++;
        if (missing != 1) {
            throw new IllegalArgumentException("exactly one ANOVA power argument must be null");
        }
        validatePowerInputs(groups, observationsPerGroup, betweenVariance,
            withinVariance, significanceLevel, power);
        double solvedGroups = groups == null
            ? solvePower(value -> anovaPower(value, observationsPerGroup,
                betweenVariance, withinVariance, significanceLevel) - power,
                2.0, 100.0) : groups;
        double solvedN = observationsPerGroup == null
            ? solvePower(value -> anovaPower(solvedGroups, value,
                betweenVariance, withinVariance, significanceLevel) - power,
                2.0, 100_000.0) : observationsPerGroup;
        double solvedBetween = betweenVariance == null
            ? solvePower(value -> anovaPower(solvedGroups, solvedN, value,
                withinVariance, significanceLevel) - power,
                withinVariance * 1e-7, withinVariance * 1e7) : betweenVariance;
        double solvedWithin = withinVariance == null
            ? solvePower(value -> anovaPower(solvedGroups, solvedN, betweenVariance,
                value, significanceLevel) - power,
                betweenVariance * 1e-7, betweenVariance * 1e7) : withinVariance;
        double solvedSignificance = significanceLevel == null
            ? solvePower(value -> anovaPower(solvedGroups, solvedN, solvedBetween,
                solvedWithin, value) - power, 1e-10, 1.0 - 1e-10)
            : significanceLevel;
        double solvedPower = power == null
            ? anovaPower(solvedGroups, solvedN, solvedBetween, solvedWithin,
                solvedSignificance) : power;
        return new AnovaPowerResult(solvedGroups, solvedN, solvedBetween,
            solvedWithin, solvedSignificance, solvedPower);
    }

    /** Score test for one or more proportions with continuity correction. */
    public static StatisticalTestResult proportion(
            int[] successes, int[] trials) {
        return proportion(successes, trials, null, Alternative.TWO_SIDED,
            DEFAULT_CONFIDENCE_LEVEL, true);
    }

    /** Score test for one or more proportions. */
    public static StatisticalTestResult proportion(
            int[] successes, int[] trials, double[] probabilities,
            Alternative alternative, double confidenceLevel, boolean correction) {
        validateProportions(successes, trials, probabilities);
        Objects.requireNonNull(alternative, "alternative");
        StatisticsSupport.probability(confidenceLevel, "confidence level");
        int groups = successes.length;
        double[] estimates = new double[groups];
        for (int index = 0; index < groups; index++) {
            estimates[index] = successes[index] / (double) trials[index];
        }
        double[] nullProbabilities;
        if (probabilities == null) {
            if (groups == 1) nullProbabilities = new double[] {0.5};
            else {
                long totalSuccess = 0L;
                long totalTrials = 0L;
                for (int index = 0; index < groups; index++) {
                    totalSuccess += successes[index];
                    totalTrials += trials[index];
                }
                nullProbabilities = new double[groups];
                Arrays.fill(nullProbabilities, totalSuccess / (double) totalTrials);
            }
        } else nullProbabilities = probabilities.clone();
        if (groups > 2 || (groups == 2 && probabilities != null)) {
            alternative = Alternative.TWO_SIDED;
        }
        double yates = correction && groups <= 2 ? 0.5 : 0.0;
        Optional<ConfidenceInterval> interval = Optional.empty();
        if (groups == 1) {
            yates = Math.min(yates,
                Math.abs(successes[0] - trials[0] * nullProbabilities[0]));
            interval = Optional.of(proportionInterval(trials[0],
                estimates[0], yates, alternative, confidenceLevel));
        } else if (groups == 2 && probabilities == null) {
            double difference = estimates[0] - estimates[1];
            yates = Math.min(yates,
                Math.abs(difference) / (1.0 / trials[0] + 1.0 / trials[1]));
            double critical = Normal.quantile(
                alternative == Alternative.TWO_SIDED
                    ? (1.0 + confidenceLevel) / 2.0 : confidenceLevel,
                0.0, 1.0, true, false);
            double standardError = Math.sqrt(
                estimates[0] * (1.0 - estimates[0]) / trials[0]
                + estimates[1] * (1.0 - estimates[1]) / trials[1]);
            double width = critical * standardError
                + yates * (1.0 / trials[0] + 1.0 / trials[1]);
            double lower = alternative == Alternative.LESS ? -1.0
                : Math.max(-1.0, difference - width);
            double upper = alternative == Alternative.GREATER ? 1.0
                : Math.min(1.0, difference + width);
            interval = Optional.of(new ConfidenceInterval(lower, upper, confidenceLevel));
        }
        double statistic = 0.0;
        for (int index = 0; index < groups; index++) {
            double expectedSuccess = trials[index] * nullProbabilities[index];
            double expectedFailure = trials[index] * (1.0 - nullProbabilities[index]);
            statistic += Math.pow(Math.abs(successes[index] - expectedSuccess) - yates, 2.0)
                / expectedSuccess;
            statistic += Math.pow(Math.abs(
                trials[index] - successes[index] - expectedFailure) - yates, 2.0)
                / expectedFailure;
        }
        double degrees = probabilities == null && groups > 1
            ? groups - 1.0 : groups;
        double pValue;
        if (alternative == Alternative.TWO_SIDED) {
            pValue = ChiSquare.cumulative(statistic, degrees, false, false);
        } else {
            double direction = groups == 1
                ? Math.signum(estimates[0] - nullProbabilities[0])
                : Math.signum(estimates[0] - estimates[1]);
            pValue = StatisticsSupport.normalPValue(
                direction * Math.sqrt(statistic), alternative);
        }
        Map<String, Double> estimateMap = indexedMap("proportion", estimates);
        Map<String, Double> nullMap = probabilities == null && groups > 1
            ? Map.of() : indexedMap("proportion", nullProbabilities);
        String method = groups == 1 ? "1-sample proportions test"
            : groups + "-sample test for "
                + (probabilities == null ? "equality of proportions" : "given proportions");
        return result(method, "X-squared", statistic, Map.of("df", degrees), pValue,
            estimateMap, nullMap, interval, alternative, correction && yates != 0.0
                ? "chi-square approximation with continuity correction"
                : "chi-square approximation without continuity correction");
    }

    /** Cochran-Armitage chi-square trend test using scores 1 through k. */
    public static StatisticalTestResult proportionTrend(
            int[] successes, int[] trials) {
        double[] scores = new double[successes == null ? 0 : successes.length];
        for (int index = 0; index < scores.length; index++) scores[index] = index + 1.0;
        return proportionTrend(successes, trials, scores);
    }

    /** Cochran-Armitage chi-square trend test with caller-supplied scores. */
    public static StatisticalTestResult proportionTrend(
            int[] successes, int[] trials, double[] scores) {
        validateProportions(successes, trials, null);
        if (scores == null || scores.length != successes.length) {
            throw new IllegalArgumentException("scores must match the number of groups");
        }
        for (double score : scores) if (!Double.isFinite(score)) {
            throw new IllegalArgumentException("scores must be finite");
        }
        scores = StatisticsSupport.centeredUnitScale(scores);
        long totalSuccess = 0L;
        long totalTrials = 0L;
        double weightedScore = 0.0;
        for (int index = 0; index < successes.length; index++) {
            totalSuccess += successes[index];
            totalTrials += trials[index];
            weightedScore += trials[index] * scores[index];
        }
        double pooled = totalSuccess / (double) totalTrials;
        double meanScore = weightedScore / totalTrials;
        double numerator = 0.0;
        double scoreVariance = 0.0;
        for (int index = 0; index < successes.length; index++) {
            double centered = scores[index] - meanScore;
            numerator += centered * successes[index];
            scoreVariance += trials[index] * centered * centered;
        }
        double statistic = numerator * numerator
            / (pooled * (1.0 - pooled) * scoreVariance);
        double pValue = ChiSquare.cumulative(statistic, 1.0, false, false);
        return result("Chi-squared Test for Trend in Proportions", "X-squared",
            statistic, Map.of("df", 1.0), pValue, Map.of(), Map.of(),
            Optional.empty(), Alternative.TWO_SIDED, "chi-square approximation");
    }

    /** Quade test for a complete blocks-by-groups matrix. */
    public static StatisticalTestResult quade(double[][] blocksByGroups) {
        StatisticsSupport.validateMatrix(blocksByGroups, 2, 2, "blocked data");
        int blocks = blocksByGroups.length;
        int groups = blocksByGroups[0].length;
        double[] ranges = new double[blocks];
        double[][] withinRanks = new double[blocks][groups];
        for (int block = 0; block < blocks; block++) {
            withinRanks[block] = StatisticsSupport.ranks(blocksByGroups[block]);
            double minimum = Arrays.stream(blocksByGroups[block]).min().orElseThrow();
            double maximum = Arrays.stream(blocksByGroups[block]).max().orElseThrow();
            ranges[block] = maximum - minimum;
        }
        double[] rangeRanks = StatisticsSupport.ranks(ranges);
        double center = (groups + 1.0) / 2.0;
        double totalSquares = 0.0;
        double[] columnSums = new double[groups];
        for (int block = 0; block < blocks; block++) {
            for (int group = 0; group < groups; group++) {
                double score = rangeRanks[block] * (withinRanks[block][group] - center);
                totalSquares += score * score;
                columnSums[group] += score;
            }
        }
        double treatmentSquares = 0.0;
        for (double sum : columnSums) treatmentSquares += sum * sum;
        treatmentSquares /= blocks;
        double statistic;
        double numeratorDegrees;
        double denominatorDegrees;
        double pValue;
        String pMethod;
        if (totalSquares == treatmentSquares) {
            statistic = Double.NaN;
            numeratorDegrees = Double.NaN;
            denominatorDegrees = Double.NaN;
            pValue = Math.pow(MathFunctions.gammafn(groups + 1.0), 1.0 - blocks);
            pMethod = "Conover zero-denominator convention";
        } else {
            statistic = (blocks - 1.0) * treatmentSquares
                / (totalSquares - treatmentSquares);
            numeratorDegrees = groups - 1.0;
            denominatorDegrees = (blocks - 1.0) * (groups - 1.0);
            pValue = F.cumulative(
                statistic, numeratorDegrees, denominatorDegrees, false, false);
            pMethod = "F distribution";
        }
        return result("Quade test", "Quade F", statistic,
            orderedMap("num df", numeratorDegrees, "denom df", denominatorDegrees),
            pValue, Map.of(), Map.of(), Optional.empty(), Alternative.TWO_SIDED,
            pMethod);
    }

    /** Shapiro-Wilk normality test, using JDistlib. */
    public static StatisticalTestResult shapiroWilk(double[] sample) {
        double[] x = StatisticsSupport.sample(sample, 3, "sample");
        if (x.length > 5000) {
            throw new IllegalArgumentException("Shapiro-Wilk supports 3 to 5000 observations");
        }
        double statistic = NormalityTest.shapiro_wilk_statistic(x);
        double pValue = NormalityTest.shapiro_wilk_pvalue(statistic, x.length);
        return result("Shapiro-Wilk normality test", "W", statistic, Map.of(),
            pValue, Map.of(), Map.of(), Optional.empty(), Alternative.TWO_SIDED,
            "Royston approximation");
    }

    /** One-sample Student t test, using JDistlib for tail probabilities. */
    public static StatisticalTestResult studentT(double[] sample, double mean) {
        return studentT(sample, mean, Alternative.TWO_SIDED,
            DEFAULT_CONFIDENCE_LEVEL);
    }

    /** One-sample Student t test. */
    public static StatisticalTestResult studentT(
            double[] sample, double mean, Alternative alternative,
            double confidenceLevel) {
        double[] x = StatisticsSupport.sample(sample, 2, "sample");
        if (!Double.isFinite(mean)) throw new IllegalArgumentException("null mean must be finite");
        StatisticsSupport.probability(confidenceLevel, "confidence level");
        Objects.requireNonNull(alternative, "alternative");
        double scale = StatisticsSupport.rescale(0.0, x);
        double estimate = StatisticsSupport.mean(x);
        double standardError = Math.sqrt(StatisticsSupport.variance(x) / x.length);
        if (!(standardError > 0.0)) throw new IllegalArgumentException("sample must have positive variance");
        double degrees = x.length - 1.0;
        double statistic = (estimate - mean / scale) / standardError;
        ConfidenceInterval interval = tInterval(
            estimate, standardError, degrees, alternative, confidenceLevel);
        return result("One Sample t-test", "t", statistic, Map.of("df", degrees),
            tPValue(statistic, degrees, alternative), Map.of("mean", estimate * scale), Map.of("mean", mean),
            Optional.of(scaleInterval(interval, scale)), alternative, "Student t");
    }

    /** Welch or pooled two-sample t test. */
    public static StatisticalTestResult studentT(
            double[] first, double[] second, double meanDifference,
            boolean equalVariances, Alternative alternative,
            double confidenceLevel) {
        double[] x = StatisticsSupport.sample(first, 2, "first sample");
        double[] y = StatisticsSupport.sample(second, 2, "second sample");
        if (!Double.isFinite(meanDifference)) {
            throw new IllegalArgumentException("null mean difference must be finite");
        }
        StatisticsSupport.probability(confidenceLevel, "confidence level");
        Objects.requireNonNull(alternative, "alternative");
        double scale = StatisticsSupport.rescale(0.0, x, y);
        double firstVariance = StatisticsSupport.variance(x);
        double secondVariance = StatisticsSupport.variance(y);
        double standardError;
        double degrees;
        if (equalVariances) {
            degrees = x.length + y.length - 2.0;
            double pooled = ((x.length - 1.0) * firstVariance
                + (y.length - 1.0) * secondVariance) / degrees;
            standardError = Math.sqrt(pooled * (1.0 / x.length + 1.0 / y.length));
        } else {
            double firstTerm = firstVariance / x.length;
            double secondTerm = secondVariance / y.length;
            standardError = Math.sqrt(firstTerm + secondTerm);
            degrees = Math.pow(firstTerm + secondTerm, 2.0)
                / (firstTerm * firstTerm / (x.length - 1.0)
                    + secondTerm * secondTerm / (y.length - 1.0));
        }
        if (!(standardError > 0.0)) throw new IllegalArgumentException("pooled standard error must be positive");
        double estimate = StatisticsSupport.mean(x) - StatisticsSupport.mean(y);
        double statistic = (estimate - meanDifference / scale) / standardError;
        return result(equalVariances ? "Two Sample t-test" : "Welch Two Sample t-test",
            "t", statistic, Map.of("df", degrees), tPValue(statistic, degrees, alternative),
            Map.of("difference in means", estimate * scale),
            Map.of("difference in means", meanDifference),
            Optional.of(scaleInterval(tInterval(estimate, standardError, degrees,
                alternative, confidenceLevel), scale)), alternative, "Student t");
    }

    /** Paired t test, implemented as a one-sample test of differences. */
    public static StatisticalTestResult pairedT(
            double[] first, double[] second, double meanDifference,
            Alternative alternative, double confidenceLevel) {
        double[] differences = StatisticsSupport.differences(first, second);
        StatisticalTestResult base = studentT(
            differences, meanDifference, alternative, confidenceLevel);
        return result("Paired t-test", base.statisticName(), base.statistic(),
            base.parameters(), base.pValue(),
            Map.of("mean difference", base.estimates().get("mean")),
            Map.of("mean difference", meanDifference), base.confidenceInterval(),
            alternative, base.pValueMethod());
    }

    /** F test for the ratio of two variances, using JDistlib. */
    public static StatisticalTestResult variance(
            double[] first, double[] second) {
        return variance(first, second, 1.0, Alternative.TWO_SIDED,
            DEFAULT_CONFIDENCE_LEVEL);
    }

    /** F test for the ratio of two variances, using JDistlib. */
    public static StatisticalTestResult variance(
            double[] first, double[] second, double ratio,
            Alternative alternative, double confidenceLevel) {
        double[] x = StatisticsSupport.sample(first, 2, "first sample");
        double[] y = StatisticsSupport.sample(second, 2, "second sample");
        if (!(ratio > 0.0) || !Double.isFinite(ratio)) {
            throw new IllegalArgumentException("null variance ratio must be positive");
        }
        StatisticsSupport.probability(confidenceLevel, "confidence level");
        StatisticsSupport.rescale(0.0, x, y);
        if (!(StatisticsSupport.variance(x) > 0.0)
                || !(StatisticsSupport.variance(y) > 0.0)) {
            throw new IllegalArgumentException("both samples must have positive variance");
        }
        double[] raw = DistributionTest.var_test(x, y, ratio, kind(alternative));
        double estimate = StatisticsSupport.variance(x) / StatisticsSupport.variance(y);
        double firstDegrees = x.length - 1.0;
        double secondDegrees = y.length - 1.0;
        double alpha = alternative == Alternative.TWO_SIDED
            ? (1.0 - confidenceLevel) / 2.0 : 1.0 - confidenceLevel;
        double lower = alternative == Alternative.LESS ? 0.0
            : estimate / F.quantile(1.0 - alpha,
                firstDegrees, secondDegrees, true, false);
        double upper = alternative == Alternative.GREATER ? Double.POSITIVE_INFINITY
            : estimate / F.quantile(alpha,
                firstDegrees, secondDegrees, true, false);
        return result("F test to compare two variances", "F", raw[0],
            orderedMap("num df", firstDegrees, "denom df", secondDegrees), raw[1],
            Map.of("ratio of variances", estimate), Map.of("ratio of variances", ratio),
            Optional.of(new ConfidenceInterval(lower, upper, confidenceLevel)),
            alternative, "F distribution");
    }

    /** One-sample Wilcoxon signed-rank test, using JDistlib. */
    public static StatisticalTestResult wilcoxonSignedRank(
            double[] sample, double location, Alternative alternative,
            boolean continuityCorrection) {
        double[] x = StatisticsSupport.sample(sample, 1, "sample");
        if (!Double.isFinite(location)) {
            throw new IllegalArgumentException("null location must be finite");
        }
        double[] raw = DistributionTest.wilcoxon_test(
            x, location, continuityCorrection, kind(alternative));
        boolean exact = !hasTies(absoluteDifferences(x, location))
            && Arrays.stream(x).noneMatch(value -> value == location);
        return result("Wilcoxon signed rank test", "V", raw[0], Map.of(), raw[1],
            Map.of(), Map.of("location shift", location), Optional.empty(), alternative,
            exact ? "exact" : continuityCorrection
                ? "normal approximation with continuity correction"
                : "normal approximation");
    }

    /** Two-sample Wilcoxon rank-sum (Mann-Whitney) test, using JDistlib. */
    public static StatisticalTestResult wilcoxonRankSum(
            double[] first, double[] second, double locationShift,
            Alternative alternative, boolean continuityCorrection) {
        double[] x = StatisticsSupport.sample(first, 1, "first sample");
        double[] y = StatisticsSupport.sample(second, 1, "second sample");
        if (!Double.isFinite(locationShift)) {
            throw new IllegalArgumentException("null location shift must be finite");
        }
        double[] raw = DistributionTest.mann_whitney_u_test(x, y, locationShift,
            continuityCorrection, false, kind(alternative));
        double[] shifted = x.clone();
        for (int index = 0; index < shifted.length; index++) {
            shifted[index] -= locationShift;
        }
        boolean exact = !hasTies(concatenate(shifted, y));
        return result("Wilcoxon rank sum test", "W", raw[0], Map.of(), raw[1],
            Map.of(), Map.of("location shift", locationShift), Optional.empty(),
            alternative, exact ? "exact" : continuityCorrection
                ? "normal approximation with continuity correction"
                : "normal approximation");
    }

    private static StatisticalTestResult result(
            String method, String statisticName, double statistic,
            Map<String, Double> parameters, double pValue,
            Map<String, Double> estimates, Map<String, Double> nullValues,
            Optional<ConfidenceInterval> confidenceInterval,
            Alternative alternative, String pValueMethod) {
        double bounded = Double.isNaN(pValue) ? Double.NaN
            : Math.max(0.0, Math.min(1.0, pValue));
        return new StatisticalTestResult(method, statisticName, statistic,
            parameters, bounded, estimates, nullValues, confidenceInterval,
            alternative, pValueMethod);
    }

    private static TestKind kind(Alternative alternative) {
        Objects.requireNonNull(alternative, "alternative");
        return switch (alternative) {
            case TWO_SIDED -> TestKind.TWO_SIDED;
            case LESS -> TestKind.LOWER;
            case GREATER -> TestKind.GREATER;
        };
    }

    private static double tPValue(
            double statistic, double degrees, Alternative alternative) {
        return switch (alternative) {
            case LESS -> T.cumulative(statistic, degrees, true, false);
            case GREATER -> T.cumulative(statistic, degrees, false, false);
            case TWO_SIDED -> Math.min(1.0, 2.0 * T.cumulative(
                Math.abs(statistic), degrees, false, false));
        };
    }

    private static ConfidenceInterval tInterval(
            double estimate, double standardError, double degrees,
            Alternative alternative, double confidenceLevel) {
        double critical = T.quantile(
            alternative == Alternative.TWO_SIDED
                ? (1.0 + confidenceLevel) / 2.0 : confidenceLevel,
            degrees, true, false);
        double lower = alternative == Alternative.LESS ? Double.NEGATIVE_INFINITY
            : estimate - critical * standardError;
        double upper = alternative == Alternative.GREATER ? Double.POSITIVE_INFINITY
            : estimate + critical * standardError;
        return new ConfidenceInterval(lower, upper, confidenceLevel);
    }

    private static double atanh(double value) {
        return 0.5 * Math.log((1.0 + value) / (1.0 - value));
    }

    private static ConfidenceInterval scaleInterval(ConfidenceInterval interval, double scale) {
        return new ConfidenceInterval(interval.lower() * scale, interval.upper() * scale, interval.level());
    }

    private static boolean hasTies(double[] values) {
        double[] sorted = values.clone();
        Arrays.sort(sorted);
        for (int index = 1; index < sorted.length; index++) {
            if (sorted[index - 1] == sorted[index]) return true;
        }
        return false;
    }

    private static double[] concatenate(double[] first, double[] second) {
        double[] result = Arrays.copyOf(first, first.length + second.length);
        System.arraycopy(second, 0, result, first.length, second.length);
        return result;
    }

    private static double[] absoluteDifferences(double[] values, double location) {
        double[] result = values.clone();
        for (int index = 0; index < result.length; index++) {
            result[index] = Math.abs(result[index] - location);
        }
        return result;
    }

    private static int groupCount(int[] labels, int minimumPerGroup) {
        Map<Integer, Integer> counts = new LinkedHashMap<>();
        for (int label : labels) counts.merge(label, 1, Integer::sum);
        if (counts.size() < 2) {
            throw new IllegalArgumentException("at least two groups are required");
        }
        for (int count : counts.values()) {
            if (count < minimumPerGroup) {
                throw new IllegalArgumentException(
                    "each group must contain at least " + minimumPerGroup + " observations");
            }
        }
        return counts.size();
    }

    private static double tiePairs(Map<Double, Integer> counts) {
        double result = 0.0;
        for (int count : counts.values()) result += count * (count - 1.0) / 2.0;
        return result;
    }

    private static double tieVarianceTerm(Map<Double, Integer> counts) {
        double result = 0.0;
        for (int count : counts.values()) {
            result += count * (count - 1.0) * (2.0 * count + 5.0);
        }
        return result;
    }

    private static double tieSecondTerm(Map<Double, Integer> counts) {
        double result = 0.0;
        for (int count : counts.values()) result += count * (count - 1.0);
        return result;
    }

    private static double tieThirdTerm(Map<Double, Integer> counts) {
        double result = 0.0;
        for (int count : counts.values()) {
            result += count * (count - 1.0) * (count - 2.0);
        }
        return result;
    }

    private static double kendallCumulative(long concordant, int size) {
        long maximum = size * (size - 1L) / 2;
        if (concordant < 0) return 0.0;
        if (concordant >= maximum) return 1.0;
        if (concordant > maximum / 2) return 1.0 - kendallCumulative(maximum - concordant - 1, size);
        if ((concordant + 1) * (double) size > 50_000_000.0) {
            throw new IllegalArgumentException("exact Kendall exceeds 50 million recurrence cells; use exact=false");
        }
        int limit = Math.toIntExact(concordant);
        double[] probabilities = {1.0};
        for (int n = 2; n <= size; n++) {
            double[] next = new double[(int) Math.min(limit, n * (n - 1L) / 2) + 1];
            double window = 0.0;
            for (int k = 0; k < next.length; k++) {
                if (k < probabilities.length) window += probabilities[k];
                if (k >= n && k - n < probabilities.length) window -= probabilities[k - n];
                next[k] = Math.max(0.0, window / n);
            }
            probabilities = next;
        }
        double cumulative = 0.0;
        for (double probability : probabilities) cumulative += probability;
        return Math.min(1.0, cumulative);
    }

    private static double spearmanTail(double q, int size, boolean lowerTail) {
        double threshold = Math.rint(q) + (lowerTail ? 2.0 : 0.0);
        double maximum = size * (size * (double) size - 1.0) / 3.0;
        if (threshold <= 0.0) return lowerTail ? 0.0 : 1.0;
        if (threshold > maximum) return lowerTail ? 1.0 : 0.0;
        if (size <= 9) {
            int[] permutation = new int[size];
            for (int index = 0; index < size; index++) permutation[index] = index + 1;
            long[] counts = new long[2];
            enumerateSpearman(permutation, size, threshold, lowerTail, counts);
            return counts[0] / (double) counts[1];
        }
        double n = size;
        double inverseN = 1.0 / n;
        double x = (6.0 * (threshold - 1.0) * inverseN / (n * n - 1.0) - 1.0)
            * Math.sqrt(n - 1.0);
        double square = x * x;
        double correction = x * inverseN * (0.2274
            + inverseN * (0.2531 + 0.1745 * inverseN)
            + square * (-0.0758 + inverseN * (0.1033 + 0.3932 * inverseN)
                - square * inverseN * (0.0879 + 0.0151 * inverseN
                    - square * (0.0072 - 0.0831 * inverseN
                        + square * inverseN * (0.0131 - 0.00046 * square)))));
        double edgeworth = correction / Math.exp(square / 2.0);
        double probability = Normal.cumulative(
            x, 0.0, 1.0, lowerTail, false) + (lowerTail ? -edgeworth : edgeworth);
        return Math.max(0.0, Math.min(1.0, probability));
    }

    private static void enumerateSpearman(
            int[] permutation, int remaining, double threshold,
            boolean lowerTail, long[] counts) {
        if (remaining == 1) {
            int statistic = 0;
            for (int index = 0; index < permutation.length; index++) {
                int difference = index + 1 - permutation[index];
                statistic += difference * difference;
            }
            if (lowerTail ? statistic < threshold : statistic >= threshold) counts[0]++;
            counts[1]++;
            return;
        }
        enumerateSpearman(permutation, remaining - 1, threshold, lowerTail, counts);
        for (int index = 0; index < remaining - 1; index++) {
            int swap = remaining % 2 == 0 ? index : 0;
            int temporary = permutation[swap];
            permutation[swap] = permutation[remaining - 1];
            permutation[remaining - 1] = temporary;
            enumerateSpearman(permutation, remaining - 1, threshold, lowerTail, counts);
        }
    }

    private static void requireTwoByTwo(long[][] table) {
        if (StatisticsSupport.total(table, 2, 2) == 0L) {
            throw new IllegalArgumentException("contingency table total must be positive");
        }
        if (table.length != 2 || table[0].length != 2) {
            throw new IllegalArgumentException("Fisher exact test currently requires a 2-by-2 table");
        }
        for (long[] row : table) for (long count : row) {
            if (count > Integer.MAX_VALUE) {
                throw new IllegalArgumentException("Fisher table entries exceed supported range");
            }
        }
    }

    private static long[][] copy(long[][] table) {
        if (table == null) return null;
        long[][] result = new long[table.length][];
        for (int row = 0; row < table.length; row++) {
            result[row] = table[row] == null ? null : table[row].clone();
        }
        return result;
    }

    private static StatisticalTestResult generalMantelHaenszel(long[][][] tables) {
        int rows = tables[0].length;
        int columns = tables[0][0].length;
        int rowDegrees = rows - 1;
        int columnDegrees = columns - 1;
        int dimension = rowDegrees * columnDegrees;
        double[] difference = new double[dimension];
        double[][] covariance = new double[dimension][dimension];
        for (long[][] table : tables) {
            double total = StatisticsSupport.total(table, rows, columns);
            double[] rowSums = new double[rows];
            double[] columnSums = new double[columns];
            for (int row = 0; row < rows; row++) {
                for (int column = 0; column < columns; column++) {
                    rowSums[row] += table[row][column];
                    columnSums[column] += table[row][column];
                }
            }
            for (int row = 0; row < rowDegrees; row++) {
                for (int column = 0; column < columnDegrees; column++) {
                    int index = row * columnDegrees + column;
                    difference[index] += table[row][column]
                        - rowSums[row] * columnSums[column] / total;
                    for (int otherRow = 0; otherRow < rowDegrees; otherRow++) {
                        double rowFactor = (row == otherRow ? total * rowSums[row] : 0.0)
                            - rowSums[row] * rowSums[otherRow];
                        for (int otherColumn = 0;
                                otherColumn < columnDegrees; otherColumn++) {
                            int other = otherRow * columnDegrees + otherColumn;
                            double columnFactor = column == otherColumn
                                ? total * columnSums[column]
                                : 0.0;
                            columnFactor -= columnSums[column] * columnSums[otherColumn];
                            covariance[index][other] += rowFactor * columnFactor
                                / (total * total * (total - 1.0));
                        }
                    }
                }
            }
        }
        double[] solved = solve(covariance, difference);
        double statistic = 0.0;
        for (int index = 0; index < dimension; index++) {
            statistic += difference[index] * solved[index];
        }
        double pValue = ChiSquare.cumulative(statistic, dimension, false, false);
        return result("Cochran-Mantel-Haenszel test",
            "Cochran-Mantel-Haenszel M^2", statistic,
            Map.of("df", (double) dimension), pValue, Map.of(), Map.of(),
            Optional.empty(), Alternative.TWO_SIDED, "chi-square approximation");
    }

    private static double[] solve(double[][] matrix, double[] rightHandSide) {
        int size = rightHandSide.length;
        double[][] augmented = new double[size][size + 1];
        for (int row = 0; row < size; row++) {
            System.arraycopy(matrix[row], 0, augmented[row], 0, size);
            augmented[row][size] = rightHandSide[row];
        }
        for (int pivot = 0; pivot < size; pivot++) {
            int selected = pivot;
            for (int row = pivot + 1; row < size; row++) {
                if (Math.abs(augmented[row][pivot])
                        > Math.abs(augmented[selected][pivot])) selected = row;
            }
            if (!(Math.abs(augmented[selected][pivot]) > 1e-14)) {
                throw new IllegalArgumentException(
                    "Cochran-Mantel-Haenszel covariance matrix is singular");
            }
            double[] temporary = augmented[pivot];
            augmented[pivot] = augmented[selected];
            augmented[selected] = temporary;
            for (int row = pivot + 1; row < size; row++) {
                double factor = augmented[row][pivot] / augmented[pivot][pivot];
                for (int column = pivot; column <= size; column++) {
                    augmented[row][column] -= factor * augmented[pivot][column];
                }
            }
        }
        double[] solution = new double[size];
        for (int row = size - 1; row >= 0; row--) {
            double residual = augmented[row][size];
            for (int column = row + 1; column < size; column++) {
                residual -= augmented[row][column] * solution[column];
            }
            solution[row] = residual / augmented[row][row];
        }
        return solution;
    }

    private static void validateStrata(long[][][] tables) {
        if (tables == null || tables.length == 0) {
            throw new IllegalArgumentException("at least one stratum is required");
        }
        int rows = tables[0] == null ? 0 : tables[0].length;
        int columns = rows == 0 || tables[0][0] == null
            ? 0 : tables[0][0].length;
        if (rows < 2 || columns < 2) {
            throw new IllegalArgumentException("strata must have at least two rows and columns");
        }
        for (long[][] table : tables) {
            StatisticsSupport.total(table, rows, columns);
            if (table.length != rows || table[0].length != columns) {
                throw new IllegalArgumentException("all strata must have the same dimensions");
            }
            if (StatisticsSupport.total(table, 2, 2) < 2L) {
                throw new IllegalArgumentException("each stratum must contain at least two observations");
            }
        }
    }

    private static void validateGroups(double[][] groups) {
        if (groups == null || groups.length < 2) {
            throw new IllegalArgumentException("at least two groups are required");
        }
        for (int index = 0; index < groups.length; index++) {
            StatisticsSupport.sample(groups[index], 2, "group " + (index + 1));
        }
    }

    private static void validateProportions(
            int[] successes, int[] trials, double[] probabilities) {
        if (successes == null || trials == null || successes.length == 0
                || successes.length != trials.length
                || (probabilities != null && probabilities.length != successes.length)) {
            throw new IllegalArgumentException("proportion inputs must have equal positive lengths");
        }
        for (int index = 0; index < successes.length; index++) {
            if (trials[index] <= 0 || successes[index] < 0
                    || successes[index] > trials[index]) {
                throw new IllegalArgumentException("invalid successes or trial count");
            }
            if (probabilities != null) {
                StatisticsSupport.probability(
                    probabilities[index], "null probability");
            }
        }
    }

    private static ConfidenceInterval proportionInterval(
            int trials, double estimate,
            double correction, Alternative alternative, double confidenceLevel) {
        double critical = Normal.quantile(
            alternative == Alternative.TWO_SIDED
                ? (1.0 + confidenceLevel) / 2.0 : confidenceLevel,
            0.0, 1.0, true, false);
        double zTerm = critical * critical / (2.0 * trials);
        double centeredUpper = estimate + correction / trials;
        double upper = centeredUpper >= 1.0 ? 1.0
            : (centeredUpper + zTerm + critical * Math.sqrt(
                centeredUpper * (1.0 - centeredUpper) / trials
                    + zTerm / (2.0 * trials))) / (1.0 + 2.0 * zTerm);
        double centeredLower = estimate - correction / trials;
        double lower = centeredLower <= 0.0 ? 0.0
            : (centeredLower + zTerm - critical * Math.sqrt(
                centeredLower * (1.0 - centeredLower) / trials
                    + zTerm / (2.0 * trials))) / (1.0 + 2.0 * zTerm);
        if (alternative == Alternative.LESS) lower = 0.0;
        if (alternative == Alternative.GREATER) upper = 1.0;
        return new ConfidenceInterval(
            Math.max(0.0, lower), Math.min(1.0, upper), confidenceLevel);
    }

    private static Map<String, Double> indexedMap(String name, double[] values) {
        Map<String, Double> result = new LinkedHashMap<>();
        if (values.length == 1) result.put(name, values[0]);
        else for (int index = 0; index < values.length; index++) {
            result.put(name + " " + (index + 1), values[index]);
        }
        return result;
    }

    private static Map<String, Double> orderedMap(
            String firstName, double firstValue,
            String secondName, double secondValue) {
        Map<String, Double> result = new LinkedHashMap<>();
        result.put(firstName, firstValue);
        result.put(secondName, secondValue);
        return result;
    }

    private static void validatePowerInputs(
            Double groups, Double observations, Double betweenVariance,
            Double withinVariance, Double significance, Double power) {
        if (groups != null && (!(groups >= 2.0) || !Double.isFinite(groups))) {
            throw new IllegalArgumentException("number of groups must be at least 2");
        }
        if (observations != null
                && (!(observations >= 2.0) || !Double.isFinite(observations))) {
            throw new IllegalArgumentException("observations per group must be at least 2");
        }
        if (betweenVariance != null && (!(betweenVariance >= 0.0)
                || !Double.isFinite(betweenVariance))) {
            throw new IllegalArgumentException("between-group variance must be nonnegative");
        }
        if (withinVariance != null && (!(withinVariance > 0.0)
                || !Double.isFinite(withinVariance))) {
            throw new IllegalArgumentException("within-group variance must be positive");
        }
        if (significance != null) StatisticsSupport.probability(
            significance, "significance level");
        if (power != null) StatisticsSupport.probability(power, "power");
    }

    private static double anovaPower(
            double groups, double observations, double betweenVariance,
            double withinVariance, double significance) {
        double numeratorDegrees = groups - 1.0;
        double denominatorDegrees = (observations - 1.0) * groups;
        double critical = F.quantile(
            significance, numeratorDegrees, denominatorDegrees, false, false);
        double noncentrality = numeratorDegrees * observations
            * betweenVariance / withinVariance;
        return NonCentralF.cumulative(critical, numeratorDegrees,
            denominatorDegrees, noncentrality, false, false);
    }

    private static double solvePower(
            StatisticsSupport.DoubleFunction function, double lower, double upper) {
        return StatisticsSupport.bisect(function, lower, upper);
    }

    private static final class ExactTableEnumerator {
        private final int[] rowSums;
        private final int[] columnSums;
        private final double observedLogWeight;
        private final double extremeThreshold;
        private final long maximumTables;
        private long tableCount;
        private double allLogSum = Double.NEGATIVE_INFINITY;
        private double extremeLogSum = Double.NEGATIVE_INFINITY;

        private ExactTableEnumerator(
                int[] rowSums, int[] columnSums, double observedLogWeight,
                long maximumTables) {
            this.rowSums = rowSums.clone();
            this.columnSums = columnSums.clone();
            this.observedLogWeight = observedLogWeight;
            this.extremeThreshold = observedLogWeight
                + Math.log1p(FISHER_RELATIVE_ERROR);
            this.maximumTables = maximumTables;
        }

        double enumerate() {
            enumerateRow(0, columnSums.clone(), 0.0);
            return Math.exp(extremeLogSum - allLogSum);
        }

        long tableCount() { return tableCount; }

        double observedProbability() {
            return Math.exp(observedLogWeight - allLogSum);
        }

        private void enumerateRow(int row, int[] remainingColumns, double logWeight) {
            if (row == rowSums.length - 1) {
                int total = 0;
                double completedWeight = logWeight;
                for (int count : remainingColumns) {
                    total = Math.addExact(total, count);
                    completedWeight -= MathFunctions.lgammafn(count + 1.0);
                }
                if (total == rowSums[row]) record(completedWeight);
                return;
            }
            enumerateCell(row, 0, rowSums[row], remainingColumns, logWeight);
        }

        private void enumerateCell(
                int row, int column, int remainingRow, int[] remainingColumns,
                double logWeight) {
            if (column == remainingColumns.length - 1) {
                if (remainingRow <= remainingColumns[column]) {
                    remainingColumns[column] -= remainingRow;
                    enumerateRow(row + 1, remainingColumns,
                        logWeight - MathFunctions.lgammafn(remainingRow + 1.0));
                    remainingColumns[column] += remainingRow;
                }
                return;
            }
            int maximum = Math.min(remainingRow, remainingColumns[column]);
            for (int count = 0; count <= maximum; count++) {
                remainingColumns[column] -= count;
                enumerateCell(row, column + 1, remainingRow - count,
                    remainingColumns,
                    logWeight - MathFunctions.lgammafn(count + 1.0));
                remainingColumns[column] += count;
            }
        }

        private void record(double logWeight) {
            tableCount++;
            if (tableCount > maximumTables) {
                throw new IllegalArgumentException(
                    "Fisher exact enumeration exceeded maximumTables=" + maximumTables);
            }
            allLogSum = logAdd(allLogSum, logWeight);
            if (logWeight <= extremeThreshold) {
                extremeLogSum = logAdd(extremeLogSum, logWeight);
            }
        }

        private static double logAdd(double first, double second) {
            if (first == Double.NEGATIVE_INFINITY) return second;
            if (second > first) {
                double temporary = first;
                first = second;
                second = temporary;
            }
            return first + Math.log1p(Math.exp(second - first));
        }
    }

    private static final class ConditionalTableDistribution {
        private final int lower;
        private final int upper;
        private final double[] logCentralProbabilities;

        private ConditionalTableDistribution(
                int lower, int upper, double[] logCentralProbabilities) {
            this.lower = lower;
            this.upper = upper;
            this.logCentralProbabilities = logCentralProbabilities;
        }

        static ConditionalTableDistribution from(long[][] table) {
            int firstColumn = Math.toIntExact(table[0][0] + table[1][0]);
            int secondColumn = Math.toIntExact(table[0][1] + table[1][1]);
            int firstRow = Math.toIntExact(table[0][0] + table[0][1]);
            int lower = Math.max(0, firstRow - secondColumn);
            int upper = Math.min(firstRow, firstColumn);
            double normalizer = MathFunctions.lchoose(
                firstColumn + (double) secondColumn, firstRow);
            double[] logProbabilities = new double[upper - lower + 1];
            for (int value = lower; value <= upper; value++) {
                logProbabilities[value - lower] = MathFunctions.lchoose(firstColumn, value)
                    + MathFunctions.lchoose(secondColumn, firstRow - value) - normalizer;
            }
            return new ConditionalTableDistribution(lower, upper, logProbabilities);
        }

        double[] probabilities(double oddsRatio) {
            double[] probabilities = new double[logCentralProbabilities.length];
            if (oddsRatio == 0.0) {
                probabilities[0] = 1.0;
                return probabilities;
            }
            if (Double.isInfinite(oddsRatio)) {
                probabilities[probabilities.length - 1] = 1.0;
                return probabilities;
            }
            return probabilitiesForLogOdds(Math.log(oddsRatio));
        }

        private double[] probabilitiesForLogOdds(double logOdds) {
            double maximum = Double.NEGATIVE_INFINITY;
            for (int index = 0; index < logCentralProbabilities.length; index++) {
                maximum = Math.max(maximum,
                    logCentralProbabilities[index] + (lower + index) * logOdds);
            }
            double sum = 0.0;
            double[] probabilities = new double[logCentralProbabilities.length];
            for (int index = 0; index < probabilities.length; index++) {
                probabilities[index] = Math.exp(logCentralProbabilities[index]
                    + (lower + index) * logOdds - maximum);
                sum += probabilities[index];
            }
            for (int index = 0; index < probabilities.length; index++) {
                probabilities[index] /= sum;
            }
            return probabilities;
        }

        double tail(double[] probabilities, int observed, boolean upperTail) {
            double result = 0.0;
            for (int value = lower; value <= upper; value++) {
                if (upperTail ? value >= observed : value <= observed) {
                    result += probabilities[value - lower];
                }
            }
            return result;
        }

        double minimumLikelihoodTail(
                double[] probabilities, int observed, double relativeError) {
            double observedProbability = probabilities[observed - lower];
            double result = 0.0;
            for (double probability : probabilities) {
                if (probability <= observedProbability * (1.0 + relativeError)) {
                    result += probability;
                }
            }
            return result;
        }

        double conditionalMle(int observed) {
            if (observed == lower) return 0.0;
            if (observed == upper) return Double.POSITIVE_INFINITY;
            double logOdds = StatisticsSupport.bisect(
                value -> expectation(probabilitiesForLogOdds(value)) - observed,
                -700.0, 700.0);
            return Math.exp(logOdds);
        }

        double lowerConfidenceLimit(int observed, double alpha) {
            if (observed == lower) return 0.0;
            double logOdds = StatisticsSupport.bisect(value -> tail(
                probabilitiesForLogOdds(value), observed, true) - alpha,
                -700.0, 700.0);
            return Math.exp(logOdds);
        }

        double upperConfidenceLimit(int observed, double alpha) {
            if (observed == upper) return Double.POSITIVE_INFINITY;
            double logOdds = StatisticsSupport.bisect(value -> tail(
                probabilitiesForLogOdds(value), observed, false) - alpha,
                -700.0, 700.0);
            return Math.exp(logOdds);
        }

        private double expectation(double[] probabilities) {
            double result = 0.0;
            for (int index = 0; index < probabilities.length; index++) {
                result += (lower + index) * probabilities[index];
            }
            return result;
        }
    }
}
