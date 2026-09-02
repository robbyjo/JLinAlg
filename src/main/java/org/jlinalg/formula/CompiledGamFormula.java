/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.formula;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gam.Gam;
import org.jlinalg.gam.GamResult;
import org.jlinalg.gam.GeneralizedGam;
import org.jlinalg.gam.GeneralizedGamResult;
import org.jlinalg.gam.PSplineTerm;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.glmm.GlmmPqlOptions;
import org.jlinalg.reml.RemlOptions;

/** A fixed and smooth formula compiled once into reusable numerical arrays. */
public final class CompiledGamFormula {
    private final CompiledFormula fixed;
    private final List<PSplineTerm> smoothTerms;

    CompiledGamFormula(
            CompiledFormula fixed, List<PSplineTerm> smoothTerms) {
        this.fixed = fixed;
        this.smoothTerms = List.copyOf(smoothTerms);
    }

    public int rows() { return fixed.rows(); }
    public int parametricColumns() { return fixed.columns(); }
    public List<String> parametricCoefficientNames() {
        return fixed.coefficientNames();
    }
    public List<PSplineTerm> smoothTerms() { return smoothTerms; }

    /** Fits an exact Gaussian REML GAM without rebuilding bases. */
    public GamResult fitGaussian(
            RemlOptions options, BackendPolicy backendPolicy) {
        if (fixed.weightsView() != null) {
            throw new IllegalArgumentException(
                "weighted Gaussian GAM is not implemented yet");
        }
        double[] response = fixed.responseView().clone();
        double[] offset = fixed.offsetView();
        if (offset != null) {
            for (int row = 0; row < response.length; row++) {
                response[row] -= offset[row];
            }
        }
        return Gam.fitGaussian(response, fixed.designView(), fixed.rows(),
            fixed.columns(), smoothTerms, options, backendPolicy);
    }

    /** Fits a generalized GAM without rebuilding bases. */
    public GeneralizedGamResult fit(
            GlmFamily family,
            GlmmPqlOptions options,
            BackendPolicy backendPolicy) {
        return GeneralizedGam.fit(fixed.responseView(), fixed.designView(),
            fixed.rows(), fixed.columns(), smoothTerms, family,
            fixed.weightsView(), fixed.offsetView(), options, backendPolicy);
    }
}
