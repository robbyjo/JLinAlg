/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import java.util.List;
import org.jlinalg.gam.PenalizedPredictor;

/** Prediction from parameter-specific compiled design matrices. */
public final class DistributionalPrediction {
    private DistributionalPrediction() { }

    /** Returns parameter-by-observation fitted distribution parameters. */
    public static double[][] parameters(
            DistributionalResult fit,
            List<PenalizedPredictor> predictors,
            DistributionalFamily family) {
        if (fit == null || predictors == null || family == null
                || predictors.size() != family.parameterCount()) {
            throw new IllegalArgumentException("fit, family, and one predictor per parameter are required");
        }
        int observations = predictors.get(0).observations();
        double[][] result = new double[family.parameterCount()][observations];
        for (int parameter = 0; parameter < family.parameterCount(); parameter++) {
            PenalizedPredictor predictor = predictors.get(parameter);
            double[] coefficients = fit.parameters().get(parameter).coefficients();
            if (predictor.observations() != observations
                    || predictor.columns() != coefficients.length) {
                throw new IllegalArgumentException("prediction design dimensions do not match fit");
            }
            double[] design = predictor.design();
            for (int row = 0; row < observations; row++) {
                double linear = 0.0;
                for (int column = 0; column < coefficients.length; column++) {
                    linear += design[row * coefficients.length + column]
                        * coefficients[column];
                }
                result[parameter][row] = family.inverseLink(parameter, linear);
            }
        }
        return result;
    }

    /**
     * Returns an observation-by-category probability matrix for multinomial
     * and ordinal vector-response families.
     */
    public static double[][] categoryProbabilities(
            DistributionalResult fit,
            List<PenalizedPredictor> predictors,
            DistributionalFamily family) {
        double[][] linked = parameters(fit, predictors, family);
        int observations = linked[0].length;
        double[][] result = null;
        double[] rowParameters = new double[linked.length];
        for (int row = 0; row < observations; row++) {
            for (int parameter = 0; parameter < linked.length; parameter++) {
                rowParameters[parameter] = linked[parameter][row];
            }
            double[] probabilities = family.categoryProbabilities(rowParameters);
            if (probabilities == null || probabilities.length < 2) {
                throw new IllegalArgumentException(
                    "family returned an invalid category-probability vector");
            }
            if (result == null) result = new double[observations][probabilities.length];
            if (probabilities.length != result[row].length) {
                throw new IllegalArgumentException(
                    "family changed its category count between observations");
            }
            double sum = 0.0;
            for (int category = 0; category < probabilities.length; category++) {
                double probability = probabilities[category];
                if (!(probability >= 0.0) || !Double.isFinite(probability)) {
                    throw new IllegalArgumentException(
                        "family returned an invalid category probability");
                }
                result[row][category] = probability;
                sum += probability;
            }
            if (Math.abs(sum - 1.0) > 1e-8) {
                throw new IllegalArgumentException(
                    "category probabilities must sum to one");
            }
        }
        return result;
    }
}
