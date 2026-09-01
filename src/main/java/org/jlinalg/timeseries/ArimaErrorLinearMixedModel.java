/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.timeseries;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import jdistlib.math.MultivariableFunction;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.mixed.LinearMixedModel;
import org.jlinalg.mixed.LinearMixedModelResult;
import org.jlinalg.mixed.RandomEffectTerm;

/** Profile REML for Gaussian LMMs with ARMA/SARMA errors. */
public final class ArimaErrorLinearMixedModel {
    private static final double TRANSFORM_BOUND = 3.8;

    private ArimaErrorLinearMixedModel() { }

    public static ArimaErrorLmmResult fit(
            double[] response,
            double[][] fixedEffects,
            List<RandomEffectTerm> randomEffects,
            ArimaOrder errorOrder) {
        return fit(response, fixedEffects, randomEffects, errorOrder,
            ArimaErrorLmmOptions.defaults(), BackendPolicy.PREFERRED);
    }

    /**
     * Profiles REML variance components at every stationary ARMA candidate.
     * Integrated orders difference response, fixed design, and random design.
     */
    public static ArimaErrorLmmResult fit(
            double[] response,
            double[][] fixedEffects,
            List<RandomEffectTerm> randomEffects,
            ArimaOrder errorOrder,
            ArimaErrorLmmOptions options,
            BackendPolicy backendPolicy) {
        if (response == null) {
            throw new IllegalArgumentException("response is required");
        }
        double[] originalResponse = MatrixOps.finiteCopy(response, "response");
        double[] originalFixed = MatrixOps.rowMajor(fixedEffects, response.length);
        if (errorOrder == null || options == null || backendPolicy == null) {
            throw new IllegalArgumentException(
                "errorOrder, options, and backendPolicy are required");
        }
        int fixedColumns = fixedEffects[0].length;
        ModelData model = differenceModel(originalResponse, originalFixed,
            fixedColumns, randomEffects, errorOrder, options.seasonalOrder());
        int parameterCount = errorOrder.autoregressive()
            + errorOrder.movingAverage()
            + options.seasonalOrder().autoregressive()
            + options.seasonalOrder().movingAverage();
        double[] initial = new double[parameterCount];
        double[] lower = new double[parameterCount];
        double[] upper = new double[parameterCount];
        Arrays.fill(lower, -TRANSFORM_BOUND);
        Arrays.fill(upper, TRANSFORM_BOUND);
        MultivariableFunction objective = parameters -> {
            try {
                ArimaMath.Coefficients coefficients = ArimaMath.decode(
                    parameters, errorOrder, options.seasonalOrder(), false, false);
                double[] correlation = ArimaMath.correlationMatrix(
                    model.response().length, coefficients.effectiveAr(),
                    coefficients.effectiveMa());
                LinearMixedModelResult fit =
                    LinearMixedModel.fitWithResidualCorrelation(
                        model.response(), model.fixed(), model.random(),
                        correlation, options.remlOptions(), backendPolicy);
                double likelihood = fit.reml().restrictedLogLikelihood();
                return Double.isFinite(likelihood)
                    ? -likelihood : Double.MAX_VALUE / 4.0;
            } catch (IllegalArgumentException | IllegalStateException exception) {
                return Double.MAX_VALUE / 4.0;
            }
        };
        BoundedOptimizer.Result optimized = Arima.multiStart(
            initial, lower, upper, objective, parameterCount,
            options.maximumFunctionEvaluations(),
            options.optimizationTolerance());
        ArimaMath.Coefficients coefficients = ArimaMath.decode(
            optimized.parameters(), errorOrder, options.seasonalOrder(),
            false, false);
        double[] correlation = ArimaMath.correlationMatrix(
            model.response().length, coefficients.effectiveAr(),
            coefficients.effectiveMa());
        LinearMixedModelResult finalFit =
            LinearMixedModel.fitWithResidualCorrelation(
                model.response(), model.fixed(), model.random(), correlation,
                options.remlOptions(), backendPolicy);
        boolean converged = optimized.converged() && finalFit.reml().converged();
        int loss = response.length - model.response().length;
        return new ArimaErrorLmmResult(finalFit, errorOrder,
            options.seasonalOrder(), coefficients, loss,
            optimized.evaluations(), converged,
            converged ? "ARMA profile and REML tolerances reached"
                : "ARMA profile or inner REML stopped before tolerance");
    }

    private static ModelData differenceModel(
            double[] response,
            double[] fixed,
            int fixedColumns,
            List<RandomEffectTerm> random,
            ArimaOrder order,
            SeasonalArimaOrder seasonal) {
        if (random == null || random.isEmpty()) {
            throw new IllegalArgumentException(
                "at least one random-effect term is required");
        }
        double[] transformedResponse = ArimaMath.difference(
            response, order, seasonal);
        double[] transformedFixed = differenceMatrix(
            fixed, response.length, fixedColumns, order, seasonal);
        List<RandomEffectTerm> transformedRandom = new ArrayList<>(random.size());
        for (RandomEffectTerm term : random) {
            if (term == null || term.observations() != response.length) {
                throw new IllegalArgumentException(
                    "random-effect terms must match response length");
            }
            double[] design = term.design();
            double[] transformed = differenceMatrix(design, response.length,
                term.coefficients(), order, seasonal);
            transformedRandom.add(RandomEffectTerm.of(term.name(),
                twoDimensional(transformed, transformedResponse.length,
                    term.coefficients()), term.coefficientNames()));
        }
        ensureNonzeroColumns(transformedFixed,
            transformedResponse.length, fixedColumns);
        return new ModelData(transformedResponse,
            twoDimensional(transformedFixed, transformedResponse.length,
                fixedColumns), transformedRandom);
    }

    private static double[] differenceMatrix(
            double[] matrix,
            int rows,
            int columns,
            ArimaOrder order,
            SeasonalArimaOrder seasonal) {
        int resultRows = ArimaMath.difference(
            new double[rows], order, seasonal).length;
        double[] result = new double[resultRows * columns];
        for (int column = 0; column < columns; column++) {
            double[] values = new double[rows];
            for (int row = 0; row < rows; row++) {
                values[row] = matrix[row * columns + column];
            }
            double[] differenced = ArimaMath.difference(values, order, seasonal);
            for (int row = 0; row < resultRows; row++) {
                result[row * columns + column] = differenced[row];
            }
        }
        return result;
    }

    private static double[][] twoDimensional(
            double[] matrix, int rows, int columns) {
        double[][] result = new double[rows][columns];
        for (int row = 0; row < rows; row++) {
            System.arraycopy(matrix, row * columns,
                result[row], 0, columns);
        }
        return result;
    }

    private static void ensureNonzeroColumns(
            double[] matrix, int rows, int columns) {
        for (int column = 0; column < columns; column++) {
            boolean nonzero = false;
            for (int row = 0; row < rows; row++) {
                nonzero |= matrix[row * columns + column] != 0.0;
            }
            if (!nonzero) {
                throw new IllegalArgumentException(
                    "differencing annihilates fixed-effect column " + column
                        + "; remove level intercepts or use a drift/trend column");
            }
        }
    }

    private record ModelData(
            double[] response,
            double[][] fixed,
            List<RandomEffectTerm> random) { }
}
