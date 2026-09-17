/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.imputation;

import java.util.ArrayList;
import java.util.List;
import jdistlib.Normal;
import jdistlib.T;

/** Rubin variance pooling with optional Barnard-Rubin finite-sample degrees of freedom. */
public final class RubinPooling {
    private RubinPooling() { }

    /** One pooled scalar estimate. */
    public record Estimate(double estimate, double withinVariance,
            double betweenVariance, double totalVariance,
            double standardError, double degreesOfFreedom,
            double relativeIncreaseInVariance,
            double fractionMissingInformation,
            double statistic, double pValue) { }

    /** Pools one parameter across imputations. Infinite completeDf uses Rubin's large-sample df. */
    public static Estimate pool(double[] estimates, double[] variances,
            double completeDf) {
        if (estimates == null || variances == null || estimates.length < 2
                || estimates.length != variances.length)
            throw new IllegalArgumentException(
                "aligned estimates and variances from at least two imputations are required");
        if (!(completeDf > 0.0) || Double.isNaN(completeDf))
            throw new IllegalArgumentException("completeDf must be positive or infinite");
        int imputations = estimates.length;
        double estimate = 0.0;
        double within = 0.0;
        for (int index = 0; index < imputations; index++) {
            if (!Double.isFinite(estimates[index]) || !Double.isFinite(variances[index])
                    || variances[index] < 0.0)
                throw new IllegalArgumentException(
                    "estimates must be finite and variances finite and nonnegative");
            estimate += estimates[index] / imputations;
            within += variances[index] / imputations;
        }
        double between = 0.0;
        for (double value : estimates) between += (value - estimate) * (value - estimate);
        between /= imputations - 1.0;
        double inflation = (1.0 + 1.0 / imputations) * between;
        double total = within + inflation;
        double relative = within == 0.0
            ? (inflation == 0.0 ? 0.0 : Double.POSITIVE_INFINITY)
            : inflation / within;
        double oldDf = inflation == 0.0 ? Double.POSITIVE_INFINITY
            : (imputations - 1.0) * Math.pow(1.0 + within / inflation, 2.0);
        double lambda = total == 0.0 ? 0.0 : inflation / total;
        double degrees = oldDf;
        if (Double.isFinite(completeDf)) {
            double observedDf = ((completeDf + 1.0) / (completeDf + 3.0))
                * completeDf * (1.0 - lambda);
            degrees = Double.isInfinite(oldDf) ? observedDf
                : 1.0 / (1.0 / oldDf + 1.0 / observedDf);
        }
        double fraction = Double.isInfinite(relative) ? 1.0
            : (relative + 2.0 / (degrees + 3.0)) / (relative + 1.0);
        double standardError = Math.sqrt(total);
        double statistic = standardError == 0.0
            ? (estimate == 0.0 ? Double.NaN
                : Math.copySign(Double.POSITIVE_INFINITY, estimate))
            : estimate / standardError;
        double pValue = Double.isNaN(statistic) ? Double.NaN
            : Double.isInfinite(degrees)
            ? Math.min(1.0, 2.0 * Normal.cumulative(
                Math.abs(statistic), 0.0, 1.0, false, false))
            : Math.min(1.0, 2.0 * T.cumulative(
                Math.abs(statistic), degrees, false, false));
        return new Estimate(estimate, within, between, total, standardError,
            degrees, relative, Math.min(1.0, fraction), statistic, pValue);
    }

    /** Pools aligned parameter columns. */
    public static List<Estimate> pool(double[][] estimates,
            double[][] variances, double completeDf) {
        if (estimates == null || variances == null || estimates.length < 2
                || estimates.length != variances.length || estimates[0].length == 0)
            throw new IllegalArgumentException("aligned imputation-by-parameter matrices are required");
        int parameters = estimates[0].length;
        for (int row = 0; row < estimates.length; row++)
            if (estimates[row].length != parameters || variances[row].length != parameters)
                throw new IllegalArgumentException("pooling matrices must be rectangular and aligned");
        List<Estimate> result = new ArrayList<>();
        for (int parameter = 0; parameter < parameters; parameter++) {
            double[] beta = new double[estimates.length];
            double[] variance = new double[estimates.length];
            for (int imputation = 0; imputation < estimates.length; imputation++) {
                beta[imputation] = estimates[imputation][parameter];
                variance[imputation] = variances[imputation][parameter];
            }
            result.add(pool(beta, variance, completeDf));
        }
        return List.copyOf(result);
    }
}
