/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.formula;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gam.GaussianSmoothSelectionResult;
import org.jlinalg.gam.GaussianSmoothSelector;
import org.jlinalg.gam.PenalizedPredictor;
import org.jlinalg.gam.QuadraticPenalizedPredictor;
import org.jlinalg.gam.QuadraticSmoothTerm;
import org.jlinalg.gam.SmoothingSelectionOptions;

/** Reusable formula compilation for multi-penalty {@code s()}, {@code te()}, and {@code ti()}. */
public final class CompiledQuadraticGamFormula {
    private final CompiledFormula fixed;
    private final List<QuadraticSmoothTerm> smoothTerms;

    CompiledQuadraticGamFormula(
            CompiledFormula fixed, List<QuadraticSmoothTerm> smoothTerms) {
        this.fixed = fixed;
        this.smoothTerms = List.copyOf(smoothTerms);
    }

    public int rows() { return fixed.rows(); }
    public int parametricColumns() { return fixed.columns(); }
    public List<String> parametricCoefficientNames() {
        return fixed.coefficientNames();
    }
    public double[] response() { return fixed.response(); }
    public double[] parametricDesign() { return fixed.design(); }
    public List<QuadraticSmoothTerm> smoothTerms() { return smoothTerms; }

    /** Selects all smoothing parameters using GCV, UBRE, or AIC. */
    public GaussianSmoothSelectionResult fitGaussian(
            SmoothingSelectionOptions options,
            BackendPolicy backendPolicy) {
        if (fixed.weightsView() != null || fixed.offsetView() != null) {
            throw new IllegalArgumentException(
                "direct tensor GAM selection currently requires no weights or offset");
        }
        return GaussianSmoothSelector.fit(fixed.responseView(),
            matrix(fixed.designView(), fixed.rows(), fixed.columns()),
            smoothTerms, null, options, backendPolicy);
    }

    /** Creates a fixed-smoothing predictor for distributional models. */
    public PenalizedPredictor predictor(
            List<double[]> smoothingParameters,
            BackendPolicy backendPolicy) {
        return QuadraticPenalizedPredictor.compile(
            matrix(fixed.designView(), fixed.rows(), fixed.columns()),
            smoothTerms, smoothingParameters, backendPolicy);
    }

    private static double[][] matrix(double[] values, int rows, int columns) {
        double[][] result = new double[rows][columns];
        for (int row = 0; row < rows; row++) {
            System.arraycopy(values, row * columns, result[row], 0, columns);
        }
        return result;
    }
}
