/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.penalized;

/** Deterministic K-fold squared-error selection for a regularization path. */
public final class PenalizedCrossValidationResult {
    private final PenalizedRegressionPath fullDataPath;
    private final double[] meanSquaredErrors;
    private final double[] standardErrors;
    private final int minimumIndex;
    private final int oneStandardErrorIndex;
    private final int folds;
    private final long randomSeed;

    PenalizedCrossValidationResult(
            PenalizedRegressionPath fullDataPath,
            double[] meanSquaredErrors,
            double[] standardErrors,
            int minimumIndex,
            int oneStandardErrorIndex,
            int folds,
            long randomSeed) {
        this.fullDataPath = fullDataPath;
        this.meanSquaredErrors = meanSquaredErrors.clone();
        this.standardErrors = standardErrors.clone();
        this.minimumIndex = minimumIndex;
        this.oneStandardErrorIndex = oneStandardErrorIndex;
        this.folds = folds;
        this.randomSeed = randomSeed;
    }

    public PenalizedRegressionPath fullDataPath() { return fullDataPath; }
    public double[] lambdas() { return fullDataPath.lambdas(); }
    public double[] meanSquaredErrors() { return meanSquaredErrors.clone(); }
    public double[] standardErrors() { return standardErrors.clone(); }
    public int folds() { return folds; }
    public long randomSeed() { return randomSeed; }
    public double lambdaMinimum() { return lambdas()[minimumIndex]; }
    public double lambdaOneStandardError() {
        return lambdas()[oneStandardErrorIndex];
    }
    public PenalizedRegressionResult minimumErrorFit() {
        return fullDataPath.fit(minimumIndex);
    }
    public PenalizedRegressionResult oneStandardErrorFit() {
        return fullDataPath.fit(oneStandardErrorIndex);
    }
}
