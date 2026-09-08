/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.List;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;

/** Cluster-robust meta-regression coefficients and sandwich covariance. */
public record MetaClusterRobustResult(
        List<String> coefficientNames,
        AssociationStatistics associationStatistics,
        double[] covariance,
        double tauSquared,
        int clusterCount,
        BackendProvenance backend,
        MetaRobustCorrection correction,
        double[] confidenceLower, double[] confidenceUpper) {
    public MetaClusterRobustResult {
        coefficientNames = List.copyOf(coefficientNames);
        covariance = covariance.clone();
        confidenceLower = confidenceLower.clone(); confidenceUpper = confidenceUpper.clone();
    }
    public double[] beta() { return associationStatistics.beta(); }
    @Override public double[] covariance() { return covariance.clone(); }
    @Override public double[] confidenceLower() { return confidenceLower.clone(); }
    @Override public double[] confidenceUpper() { return confidenceUpper.clone(); }
}
