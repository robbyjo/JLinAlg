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
    private final double[] conditionalFStatistics;
    private final double intercept;
    private final double interceptStandardError;
    private final double q;
    private final int qDegreesOfFreedom;

    MultivariableMrResult(List<String> names, AssociationStatistics statistics,
            double[] covariance, double[] conditionalFStatistics,
            double intercept, double interceptStandardError,
            double q, int qDegreesOfFreedom) {
        this.exposureNames = List.copyOf(names);
        this.statistics = statistics;
        this.covariance = covariance.clone();
        this.conditionalFStatistics = conditionalFStatistics.clone();
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
    public double[] conditionalFStatistics() { return conditionalFStatistics.clone(); }
    public double intercept() { return intercept; }
    public double interceptStandardError() { return interceptStandardError; }
    public double cochranQ() { return q; }
    public int heterogeneityDegreesOfFreedom() { return qDegreesOfFreedom; }
}
