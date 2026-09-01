/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.inference;

import java.util.Arrays;
import java.util.Objects;
import jdistlib.Normal;
import jdistlib.T;

/**
 * Immutable coefficient-level association estimates and two-sided Wald tests.
 * Arrays follow the fixed-effect design-column order.
 */
public final class AssociationStatistics {
    private static final double LOG_TWO = Math.log(2.0);
    private static final double LOG_TEN = Math.log(10.0);
    private final double[] beta;
    private final double[] standardErrors;
    private final double[] statistics;
    private final double[] pValues;
    private final double[] negativeLog10PValues;
    private final double[] degreesOfFreedom;
    private final DegreesOfFreedomMethod degreesOfFreedomMethod;
    private final StatisticDistribution statisticDistribution;

    private AssociationStatistics(
            double[] beta,
            double[] standardErrors,
            double[] degreesOfFreedom,
            DegreesOfFreedomMethod degreesOfFreedomMethod,
            StatisticDistribution statisticDistribution) {
        if (beta == null || standardErrors == null || degreesOfFreedom == null
                || beta.length != standardErrors.length
                || beta.length != degreesOfFreedom.length) {
            throw new IllegalArgumentException(
                "beta, standard errors, and degrees of freedom must have equal lengths");
        }
        this.beta = beta.clone();
        this.standardErrors = standardErrors.clone();
        this.degreesOfFreedom = degreesOfFreedom.clone();
        this.degreesOfFreedomMethod = Objects.requireNonNull(
            degreesOfFreedomMethod, "degreesOfFreedomMethod");
        this.statisticDistribution = Objects.requireNonNull(
            statisticDistribution, "statisticDistribution");
        this.statistics = new double[beta.length];
        this.pValues = new double[beta.length];
        this.negativeLog10PValues = new double[beta.length];
        calculateTests();
    }

    /** Creates Student's t tests with scalar denominator degrees of freedom. */
    public static AssociationStatistics studentT(
            double[] beta,
            double[] standardErrors,
            double degreesOfFreedom,
            DegreesOfFreedomMethod method) {
        double[] values = new double[beta.length];
        Arrays.fill(values, degreesOfFreedom);
        return studentT(beta, standardErrors, values, method);
    }

    /** Creates Student's t tests with coefficient-specific denominator DF. */
    public static AssociationStatistics studentT(
            double[] beta,
            double[] standardErrors,
            double[] degreesOfFreedom,
            DegreesOfFreedomMethod method) {
        if (method == DegreesOfFreedomMethod.ASYMPTOTIC) {
            throw new IllegalArgumentException(
                "Student's t inference requires finite-DF method metadata");
        }
        return new AssociationStatistics(beta, standardErrors,
            degreesOfFreedom, method, StatisticDistribution.STUDENT_T);
    }

    /** Creates asymptotic standard-normal Wald tests. */
    public static AssociationStatistics normal(
            double[] beta, double[] standardErrors) {
        double[] degreesOfFreedom = new double[beta.length];
        Arrays.fill(degreesOfFreedom, Double.POSITIVE_INFINITY);
        return new AssociationStatistics(beta, standardErrors,
            degreesOfFreedom, DegreesOfFreedomMethod.ASYMPTOTIC,
            StatisticDistribution.STANDARD_NORMAL);
    }

    private void calculateTests() {
        for (int index = 0; index < beta.length; index++) {
            double estimate = beta[index];
            double standardError = standardErrors[index];
            if (!Double.isFinite(standardError) || standardError < 0.0) {
                statistics[index] = Double.NaN;
                pValues[index] = Double.NaN;
                negativeLog10PValues[index] = Double.NaN;
                continue;
            }
            if (standardError == 0.0) {
                statistics[index] = estimate == 0.0
                    ? Double.NaN
                    : Math.copySign(Double.POSITIVE_INFINITY, estimate);
            } else {
                statistics[index] = estimate / standardError;
            }
            if (Double.isNaN(statistics[index])) {
                pValues[index] = Double.NaN;
                negativeLog10PValues[index] = Double.NaN;
            } else if (statisticDistribution == StatisticDistribution.STANDARD_NORMAL) {
                double magnitude = Math.abs(statistics[index]);
                pValues[index] = Math.min(1.0, 2.0 * Normal.cumulative(
                    magnitude, 0.0, 1.0, false, false));
                negativeLog10PValues[index] = negativeLog10(
                    Normal.cumulative(magnitude, 0.0, 1.0, false, true));
            } else {
                double degrees = degreesOfFreedom[index];
                if (Double.isFinite(degrees) && degrees > 0.0) {
                    double magnitude = Math.abs(statistics[index]);
                    pValues[index] = Math.min(1.0, 2.0 * T.cumulative(
                        magnitude, degrees, false, false));
                    negativeLog10PValues[index] = negativeLog10(
                        T.cumulative(magnitude, degrees, false, true));
                } else {
                    pValues[index] = Double.NaN;
                    negativeLog10PValues[index] = Double.NaN;
                }
            }
        }
    }

    private static double negativeLog10(double logUpperTail) {
        return -Math.min(0.0, LOG_TWO + logUpperTail) / LOG_TEN;
    }

    /** Fixed-effect estimates (beta). */
    public double[] beta() { return beta.clone(); }

    /** Alias emphasizing that beta contains the requested effect sizes. */
    public double[] effectSizes() { return beta.clone(); }

    /** Standard errors corresponding to {@link #beta()}. */
    public double[] standardErrors() { return standardErrors.clone(); }

    /** Wald t or z statistics as identified by {@link #statisticDistribution()}. */
    public double[] statistics() { return statistics.clone(); }

    /** Two-sided p-values. */
    public double[] pValues() { return pValues.clone(); }

    /** Two-sided p-values represented on a caller-selected scale. */
    public double[] pValues(PValueScale scale) {
        Objects.requireNonNull(scale, "scale");
        return switch (scale) {
            case REGULAR -> pValues();
            case NEGATIVE_LOG10 -> negativeLog10PValues();
            case LOG10 -> log10PValues();
        };
    }

    /** Base-10 logarithm of each two-sided p-value. */
    public double[] log10PValues() {
        double[] result = negativeLog10PValues.clone();
        for (int index = 0; index < result.length; index++)
            result[index] = -result[index];
        return result;
    }

    /** GWAS-standard negative base-10 logarithm of each two-sided p-value. */
    public double[] negativeLog10PValues() {
        return negativeLog10PValues.clone();
    }

    /** Denominator DF; positive infinity denotes asymptotic normal inference. */
    public double[] degreesOfFreedom() { return degreesOfFreedom.clone(); }

    public DegreesOfFreedomMethod degreesOfFreedomMethod() {
        return degreesOfFreedomMethod;
    }

    public StatisticDistribution statisticDistribution() {
        return statisticDistribution;
    }
}
