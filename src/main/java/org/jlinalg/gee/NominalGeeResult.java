/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

import org.jlinalg.compute.BackendProvenance;

/** Marginal multinomial-logit estimates with a cluster sandwich covariance. */
public final class NominalGeeResult {
    private final int categories;
    private final int predictors;
    private final double[] coefficients;
    private final double[] covariance;
    private final double[] standardErrors;
    private final double[] statistics;
    private final double[] pValues;
    private final double[] fittedProbabilities;
    private final int clusters;
    private final int iterations;
    private final boolean converged;
    private final double degreesOfFreedom;
    private final BackendProvenance backend;

    NominalGeeResult(
            int categories, int predictors,
            double[] coefficients, double[] covariance,
            double[] standardErrors, double[] statistics, double[] pValues,
            double[] fittedProbabilities, int clusters, int iterations,
            boolean converged, double degreesOfFreedom,
            BackendProvenance backend) {
        this.categories = categories;
        this.predictors = predictors;
        this.coefficients = coefficients.clone();
        this.covariance = covariance.clone();
        this.standardErrors = standardErrors.clone();
        this.statistics = statistics.clone();
        this.pValues = pValues.clone();
        this.fittedProbabilities = fittedProbabilities.clone();
        this.clusters = clusters;
        this.iterations = iterations;
        this.converged = converged;
        this.degreesOfFreedom = degreesOfFreedom;
        this.backend = backend;
    }

    public int categories() { return categories; }
    public int predictors() { return predictors; }
    /** Row-major non-reference-category by predictor coefficient matrix. */
    public double[] coefficients() { return coefficients.clone(); }
    public double[] covariance() { return covariance.clone(); }
    public double[] standardErrors() { return standardErrors.clone(); }
    public double[] statistics() { return statistics.clone(); }
    public double[] pValues() { return pValues.clone(); }
    /** Row-major observation-by-category fitted probability matrix. */
    public double[] fittedProbabilities() { return fittedProbabilities.clone(); }
    public int clusters() { return clusters; }
    public int iterations() { return iterations; }
    public boolean converged() { return converged; }
    public double degreesOfFreedom() { return degreesOfFreedom; }
    public BackendProvenance backend() { return backend; }
}
