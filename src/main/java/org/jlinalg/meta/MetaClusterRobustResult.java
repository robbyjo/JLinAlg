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
        BackendProvenance backend) {
    public MetaClusterRobustResult {
        coefficientNames = List.copyOf(coefficientNames);
        covariance = covariance.clone();
    }
    public double[] beta() { return associationStatistics.beta(); }
    @Override public double[] covariance() { return covariance.clone(); }
}
