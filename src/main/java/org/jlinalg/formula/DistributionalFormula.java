/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.formula;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.distributional.DistributionalFamily;
import org.jlinalg.distributional.DistributionalModel;
import org.jlinalg.distributional.DistributionalOptions;
import org.jlinalg.distributional.DistributionalResult;
import org.jlinalg.gam.PenalizedPredictor;

/** User-facing parameter-specific formulas for GAMLSS/VGAM-style models. */
public final class DistributionalFormula {
    private DistributionalFormula() { }

    /**
     * Fits one formula per family parameter. Smooth formulas accept s/te/ti;
     * each smooth needs a fixed positive smoothing vector in formula order.
     */
    public static DistributionalResult fit(
            List<String> formulas,
            List<List<double[]>> smoothingParameters,
            ModelTable table,
            DistributionalFamily family,
            FormulaOptions formulaOptions,
            DistributionalOptions fitOptions,
            BackendPolicy backendPolicy) {
        if (formulas == null || formulas.size() != family.parameterCount()
                || smoothingParameters == null
                || smoothingParameters.size() != formulas.size()) {
            throw new IllegalArgumentException(
                "one formula and smoothing collection are required per parameter");
        }
        List<PenalizedPredictor> predictors = new ArrayList<>(formulas.size());
        double[] response = null;
        for (int parameter = 0; parameter < formulas.size(); parameter++) {
            String formula = formulas.get(parameter);
            boolean smooth = formula.contains("s(") || formula.contains("te(")
                || formula.contains("ti(");
            if (smooth) {
                CompiledQuadraticGamFormula compiled = AdvancedGamFormula.compile(
                    formula, table, formulaOptions);
                response = consistentResponse(response, compiled.response());
                predictors.add(compiled.predictor(
                    smoothingParameters.get(parameter), backendPolicy));
            } else {
                CompiledFormula compiled = Formula.compile(
                    formula, table, formulaOptions);
                response = consistentResponse(response, compiled.response());
                if (compiled.weightsView() != null || compiled.offsetView() != null) {
                    throw new IllegalArgumentException(
                        "distributional formulas do not yet accept weights or offset");
                }
                predictors.add(PenalizedPredictor.linear(
                    matrix(compiled.designView(), compiled.rows(), compiled.columns())));
            }
        }
        return DistributionalModel.fit(response, predictors, family,
            fitOptions, backendPolicy);
    }

    private static double[] consistentResponse(double[] previous, double[] current) {
        if (previous == null) return current;
        if (!java.util.Arrays.equals(previous, current)) {
            throw new IllegalArgumentException(
                "all parameter formulas must use the same response");
        }
        return previous;
    }

    private static double[][] matrix(double[] values, int rows, int columns) {
        double[][] result = new double[rows][columns];
        for (int row = 0; row < rows; row++) {
            System.arraycopy(values, row * columns, result[row], 0, columns);
        }
        return result;
    }
}
