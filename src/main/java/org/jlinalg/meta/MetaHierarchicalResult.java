/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.Arrays;
import java.util.List;
import org.jlinalg.inference.AssociationStatistics;

/** Jointly estimated hierarchical random-coefficient fit. Covariances are row-major. */
public record MetaHierarchicalResult(MetaMultilevelRegressionResult model,
        List<String> randomEffectNames, List<double[]> randomCovariances,
        double[][] totalCovariance, AssociationStatistics associationStatistics,
        double[] confidenceLower, double[] confidenceUpper,
        MetaMultilevelRegression.Estimation estimation, int iterations) {
    public MetaHierarchicalResult {
        randomEffectNames = List.copyOf(randomEffectNames);
        randomCovariances = randomCovariances.stream().map(double[]::clone).toList();
        totalCovariance = copy(totalCovariance);
        confidenceLower = confidenceLower.clone(); confidenceUpper = confidenceUpper.clone();
    }
    @Override public List<double[]> randomCovariances() { return randomCovariances.stream().map(double[]::clone).toList(); }
    @Override public double[][] totalCovariance() { return copy(totalCovariance); }
    @Override public double[] confidenceLower() { return confidenceLower.clone(); }
    @Override public double[] confidenceUpper() { return confidenceUpper.clone(); }
    public boolean converged() { return model.converged(); }
    public double[] beta() { return model.beta(); }
    public double logLikelihood() { return model.logLikelihood(); }
    private static double[][] copy(double[][] a) { return Arrays.stream(a).map(double[]::clone).toArray(double[][]::new); }
}
