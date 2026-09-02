/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.gam.PenalizedPredictor;

/** Distributional fit with AIC-selected parameter-specific smooth penalties. */
public final class DistributionalSmoothingResult {
    private final DistributionalResult fit;
    private final List<PenalizedPredictor> predictors;
    private final List<List<double[]>> smoothingParameters;
    private final double aic;
    private final int evaluations;

    DistributionalSmoothingResult(
            DistributionalResult fit,
            List<PenalizedPredictor> predictors,
            List<List<double[]>> smoothingParameters,
            double aic,
            int evaluations) {
        this.fit = fit;
        this.predictors = List.copyOf(predictors);
        this.smoothingParameters = copy(smoothingParameters);
        this.aic = aic;
        this.evaluations = evaluations;
    }
    public DistributionalResult fit() { return fit; }
    public List<PenalizedPredictor> predictors() { return predictors; }
    public List<List<double[]>> smoothingParameters() {
        return copy(smoothingParameters);
    }
    public double aic() { return aic; }
    public int evaluations() { return evaluations; }

    private static List<List<double[]>> copy(List<List<double[]>> source) {
        List<List<double[]>> result = new ArrayList<>(source.size());
        for (List<double[]> parameter : source) {
            List<double[]> terms = new ArrayList<>(parameter.size());
            for (double[] term : parameter) terms.add(term.clone());
            result.add(List.copyOf(terms));
        }
        return List.copyOf(result);
    }
}
