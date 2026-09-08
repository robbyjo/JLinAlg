/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.List;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;

/** Generalized inverse-covariance meta-analysis result. */
public record MetaCorrelatedResult(
        List<String> studyNames,
        AssociationStatistics associationStatistics,
        double[] covariance,
        double tauSquared,
        double generalizedQ,
        double qDegreesOfFreedom,
        double[] normalizedWeights,
        BackendProvenance backend,
        double confidenceLower, double confidenceUpper) {
    public MetaCorrelatedResult {
        studyNames = List.copyOf(studyNames);
        covariance = covariance.clone();
        normalizedWeights = normalizedWeights.clone();
    }
    public double pooledEffectSize() { return associationStatistics.beta()[0]; }
    public double standardError() { return associationStatistics.standardErrors()[0]; }
    @Override public double[] covariance() { return covariance.clone(); }
    @Override public double[] normalizedWeights() { return normalizedWeights.clone(); }
}
