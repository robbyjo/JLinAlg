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
    private final double[] conditionalFStatistics;
    private final double intercept;
    private final double interceptStandardError;
    private final double q;
    private final int qDegreesOfFreedom;

    MultivariableMrResult(List<String> names, AssociationStatistics statistics,
            double[] covariance, double[] marginalFStatistics,
            double[] conditionalFStatistics,
            double intercept, double interceptStandardError,
            double q, int qDegreesOfFreedom) {
        this.exposureNames = List.copyOf(names);
        this.statistics = statistics;
        this.covariance = covariance.clone();
        this.marginalFStatistics = marginalFStatistics.clone();
        this.conditionalFStatistics = conditionalFStatistics == null
            ? null : conditionalFStatistics.clone();
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
     * Sanderson-Windmeijer-style conditional F statistics when this result was
     * created by a covariance-aware generalized fit.
     *
     * @deprecated Prefer the explicitly typed {@link ConditionalStrengthResult}
     * returned by {@link MultivariableMendelianRandomization#conditionalStrength}
     * because it also retains conditional Q and degrees of freedom.
     */
    @Deprecated
    public double[] conditionalFStatistics() {
        if (conditionalFStatistics == null)
            throw new UnsupportedOperationException("conditional F statistics require cross-exposure sampling covariance; use the covariance-aware fit overload");
        return conditionalFStatistics.clone();
    }
    public boolean conditionalStrengthAvailable() {
        return conditionalFStatistics != null;
    }
    public double intercept() { return intercept; }
    public double interceptStandardError() { return interceptStandardError; }
    public double cochranQ() { return q; }
    public int heterogeneityDegreesOfFreedom() { return qDegreesOfFreedom; }
}
