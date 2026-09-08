/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import java.util.List;
import org.jlinalg.inference.AssociationStatistics;

/** Direct effects from multivariable summary-data MR. */
public final class MultivariableMrResult {
    private final List<String> exposureNames;
    private final AssociationStatistics statistics;
    private final double[] covariance;
    private final double[] marginalFStatistics;
    private final double intercept;
    private final double interceptStandardError;
    private final double q;
    private final int qDegreesOfFreedom;

    MultivariableMrResult(List<String> names, AssociationStatistics statistics,
            double[] covariance, double[] marginalFStatistics,
            double intercept, double interceptStandardError,
            double q, int qDegreesOfFreedom) {
        this.exposureNames = List.copyOf(names);
        this.statistics = statistics;
        this.covariance = covariance.clone();
        this.marginalFStatistics = marginalFStatistics.clone();
        this.intercept = intercept;
        this.interceptStandardError = interceptStandardError;
        this.q = q;
        this.qDegreesOfFreedom = qDegreesOfFreedom;
    }
    public List<String> exposureNames() { return exposureNames; }
    public AssociationStatistics statistics() { return statistics; }
    public double[] beta() { return statistics.beta(); }
    public double[] standardErrors() { return statistics.standardErrors(); }
    public double[] pValues() { return statistics.pValues(); }
    public double[] covariance() { return covariance.clone(); }
    /** Mean squared marginal exposure z scores; not conditional instrument strength. */
    public double[] marginalFStatistics() { return marginalFStatistics.clone(); }
    /**
     * @deprecated This accessor previously mislabeled marginal statistics as
     * conditional F statistics. Conditional strength requires a specified
     * cross-exposure sampling covariance model, which this fit does not have.
     * Use {@link #marginalFStatistics()} only for marginal strength.
     * @throws UnsupportedOperationException always
     */
    @Deprecated
    public double[] conditionalFStatistics() {
        throw new UnsupportedOperationException("conditional F statistics require cross-exposure sampling covariance; marginalFStatistics() reports marginal strength only");
    }
    public double intercept() { return intercept; }
    public double interceptStandardError() { return interceptStandardError; }
    public double cochranQ() { return q; }
    public int heterogeneityDegreesOfFreedom() { return qDegreesOfFreedom; }
}
