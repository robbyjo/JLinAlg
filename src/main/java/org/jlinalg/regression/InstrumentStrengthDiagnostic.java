/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.regression;

import org.jlinalg.inference.StatisticDistribution;

/**
 * Excluded-instrument diagnostics for one endogenous regressor.
 *
 * <p>The partial R-squared and classical F statistic compare the unrestricted
 * first stage to a regression on the included exogenous regressors alone. The
 * reported Wald statistic uses the covariance estimator selected for the 2SLS
 * fit. For HC covariance its reference law is chi-square and
 * {@link #effectiveFStatistic()} is the chi-square statistic divided by its
 * numerator degrees of freedom. For homoskedastic or clustered covariance the
 * reference law is F and {@code effectiveFStatistic} is that F statistic.</p>
 */
public record InstrumentStrengthDiagnostic(
        int endogenousColumn,
        double partialRSquared,
        double classicalFStatistic,
        double waldStatistic,
        double effectiveFStatistic,
        int numeratorDegreesOfFreedom,
        double denominatorDegreesOfFreedom,
        double pValue,
        StatisticDistribution referenceDistribution) {
}
