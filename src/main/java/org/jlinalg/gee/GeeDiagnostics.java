/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

/** Cluster score, leverage, influence, and deletion diagnostics. */
public final class GeeDiagnostics {
    private final int[] clusterIds;
    private final int parameters;
    private final double[] clusterScores;
    private final double[] leverage;
    private final double[] cookDistances;
    private final double[] oneStepDeletedCoefficients;
    private final double[] exactDeletedCoefficients;

    GeeDiagnostics(
            int[] clusterIds,
            int parameters,
            double[] clusterScores,
            double[] leverage,
            double[] cookDistances,
            double[] oneStepDeletedCoefficients,
            double[] exactDeletedCoefficients) {
        this.clusterIds = clusterIds.clone();
        this.parameters = parameters;
        this.clusterScores = clusterScores.clone();
        this.leverage = leverage.clone();
        this.cookDistances = cookDistances.clone();
        this.oneStepDeletedCoefficients = oneStepDeletedCoefficients.clone();
        this.exactDeletedCoefficients = exactDeletedCoefficients.clone();
    }

    public int clusters() { return clusterIds.length; }
    public int parameters() { return parameters; }
    public int[] clusterIds() { return clusterIds.clone(); }
    /** Row-major {@code clusters x parameters} ordinary cluster scores. */
    public double[] clusterScores() { return clusterScores.clone(); }
    /** Trace of the observation-level cluster leverage block. */
    public double[] leverage() { return leverage.clone(); }
    public double[] cookDistances() { return cookDistances.clone(); }
    /** Row-major one-step delete-cluster coefficient approximations. */
    public double[] oneStepDeletedCoefficients() {
        return oneStepDeletedCoefficients.clone();
    }
    /** Exact delete-cluster coefficients, or an empty array when not requested. */
    public double[] exactDeletedCoefficients() {
        return exactDeletedCoefficients.clone();
    }
    public boolean hasExactDeletionFits() {
        return exactDeletedCoefficients.length != 0;
    }
}
