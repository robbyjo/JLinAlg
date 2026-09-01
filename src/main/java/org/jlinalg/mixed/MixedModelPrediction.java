/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.jlinalg.internal.MatrixOps;

/** Marginal and conditional prediction, including policy for unseen groups. */
public final class MixedModelPrediction {
    private MixedModelPrediction() { }

    /** Predicts from fixed effects only, corresponding to lme4 {@code re.form=NA}. */
    public static double[] marginal(
            LinearMixedModelResult fitted, double[][] fixedEffects) {
        if (fitted == null || fixedEffects == null || fixedEffects.length == 0)
            throw new IllegalArgumentException("fit and fixed-effect design are required");
        double[] beta = fitted.beta();
        double[] fixed = MatrixOps.rowMajor(fixedEffects, fixedEffects.length);
        if (fixedEffects[0].length != beta.length)
            throw new IllegalArgumentException(
                "prediction fixed columns must equal fitted beta length");
        double[] result = new double[fixedEffects.length];
        for (int row = 0; row < result.length; row++)
            for (int column = 0; column < beta.length; column++)
                result[row] += fixed[row * beta.length + column] * beta[column];
        return result;
    }

    /**
     * Adds requested random terms to marginal predictions. Unseen coefficient
     * labels contribute zero only when {@code allowNewLevels} is true.
     */
    public static double[] conditional(
            LinearMixedModelResult fitted,
            double[][] fixedEffects,
            List<RandomEffectTerm> randomEffects,
            boolean allowNewLevels) {
        double[] result = marginal(fitted, fixedEffects);
        if (randomEffects == null)
            throw new IllegalArgumentException("prediction random effects are required");
        for (RandomEffectTerm term : randomEffects) {
            if (term == null || term.observations() != result.length)
                throw new IllegalArgumentException(
                    "prediction random-effect rows must equal prediction rows");
            RandomEffectEstimates estimates = fitted.randomEffects(term.name());
            Map<String, Double> byCoefficient = new HashMap<>();
            List<String> fittedNames = estimates.coefficientNames();
            double[] fittedValues = estimates.estimates();
            for (int index = 0; index < fittedNames.size(); index++)
                byCoefficient.put(fittedNames.get(index), fittedValues[index]);
            double[] termDesign = term.design();
            List<String> names = term.coefficientNames();
            for (int column = 0; column < names.size(); column++) {
                Double value = byCoefficient.get(names.get(column));
                if (value == null) {
                    if (!allowNewLevels)
                        throw new IllegalArgumentException(
                            "new random-effect level in " + term.name()
                                + ": " + names.get(column));
                    continue;
                }
                for (int row = 0; row < result.length; row++)
                    result[row] += termDesign[row * names.size() + column] * value;
            }
        }
        return result;
    }
}
