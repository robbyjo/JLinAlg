/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

/** Optional exact derivatives for a GlmFamily's conditional log likelihood.
 * Both derivatives include the prior weight and hold nuisance parameters fixed.
 * Information means minus the OBSERVED second derivative, not Fisher information.
 * Families without this interface use their working-score contract and a central
 * difference of that score for noncanonical observed curvature. */
public interface LaplaceFamilyDerivatives {
    double linearPredictorScore(double response,double linearPredictor,double priorWeight);
    double linearPredictorInformation(double response,double linearPredictor,double priorWeight);
}
