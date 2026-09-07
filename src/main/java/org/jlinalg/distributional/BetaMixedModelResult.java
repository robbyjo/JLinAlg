/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.distributional;

import java.util.List;
import java.util.Map;
import org.jlinalg.glmm.GlmmLaplaceResult;
import org.jlinalg.inference.AssociationStatistics;

/** Fixed effects, conditional modes, and variance components from a beta GLMM. */
public final class BetaMixedModelResult {
    private final GlmmLaplaceResult fit;
    private final double precision;

    BetaMixedModelResult(GlmmLaplaceResult fit, double precision) {
        this.fit = fit;
        this.precision = precision;
    }

    public double[] meanCoefficients() { return fit.beta(); }
    public double[] standardErrors() { return fit.standardErrors(); }
    public AssociationStatistics associationStatistics() {
        return fit.associationStatistics();
    }
    public double[] fixedEffectCovariance() {
        return fit.fixedEffectCovariance();
    }
    /** Constant beta precision phi, where Var(Y)=mu(1-mu)/(1+phi). */
    public double precision() { return precision; }
    public List<String> componentNames() { return fit.componentNames(); }
    public double[] varianceComponents() { return fit.varianceComponents(); }
    /** Conditional coefficient modes; pedigree terms include unobserved ancestors. */
    public double[] randomEffects(String name) {
        return fit.componentCoefficients(name);
    }
    public Map<String, double[]> randomEffects() {
        return fit.componentCoefficients();
    }
    public double[] componentPredictor(String name) {
        return fit.componentPredictor(name);
    }
    public double[] linearPredictor() { return fit.linearPredictor(); }
    public double[] fittedMeans() { return fit.fittedMeans(); }
    public double marginalLogLikelihood() { return fit.marginalLogLikelihood(); }
    public int outerIterations() { return fit.outerIterations(); }
    public int modeIterations() { return fit.modeIterations(); }
    public boolean converged() { return fit.converged(); }
    public String convergenceMessage() { return fit.convergenceMessage(); }
}
