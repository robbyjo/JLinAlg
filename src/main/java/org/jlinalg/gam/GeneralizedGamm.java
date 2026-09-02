/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.glmm.GlmmPql;
import org.jlinalg.glmm.GlmmPqlOptions;
import org.jlinalg.glmm.GlmmPqlResult;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.VarianceComponent;

/** PQL generalized additive mixed models with arbitrary covariance components. */
public final class GeneralizedGamm {
    private GeneralizedGamm() { }

    /** Fits with unit prior weights, zero offset, and default controls. */
    public static GeneralizedGammResult fit(
            double[] response,
            double[][] parametricDesign,
            List<PSplineTerm> smoothTerms,
            GlmFamily family,
            List<VarianceComponent> randomComponents) {
        return fit(response, parametricDesign, smoothTerms, family,
            randomComponents, null, null, GlmmPqlOptions.defaults(),
            BackendPolicy.PREFERRED);
    }

    /** Fits a generalized GAMM with explicit weights, offset, and controls. */
    public static GeneralizedGammResult fit(
            double[] response,
            double[][] parametricDesign,
            List<PSplineTerm> smoothTerms,
            GlmFamily family,
            List<VarianceComponent> randomComponents,
            double[] priorWeights,
            double[] offset,
            GlmmPqlOptions options,
            BackendPolicy backendPolicy) {
        if (response == null) {
            throw new IllegalArgumentException("response is required");
        }
        double[] fixed = MatrixOps.rowMajor(parametricDesign, response.length);
        return fit(response, fixed, response.length, parametricDesign[0].length,
            smoothTerms, family, randomComponents, priorWeights, offset,
            options, backendPolicy);
    }

    /** Contiguous row-major overload for allocation-sensitive callers. */
    public static GeneralizedGammResult fit(
            double[] response,
            double[] parametricDesign,
            int rows,
            int parametricColumns,
            List<PSplineTerm> smoothTerms,
            GlmFamily family,
            List<VarianceComponent> randomComponents,
            double[] priorWeights,
            double[] offset,
            GlmmPqlOptions options,
            BackendPolicy backendPolicy) {
        MatrixOps.validateModelData(
            response, parametricDesign, rows, parametricColumns);
        if (family == null || options == null || backendPolicy == null) {
            throw new IllegalArgumentException(
                "family, options, and backendPolicy are required");
        }
        List<VarianceComponent> additional = randomComponents == null
            ? List.of() : List.copyOf(randomComponents);
        validateComponents(additional, rows);
        double[] weights = weights(priorWeights, rows);
        double[] offsets = offset == null ? new double[rows]
            : MatrixOps.finiteCopy(offset, "offset");
        if (offsets.length != rows) {
            throw new IllegalArgumentException("offset length must equal rows");
        }

        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            PSplineMixedModelCompiler.Compiled compiled =
                PSplineMixedModelCompiler.compile(parametricDesign, rows,
                    parametricColumns, smoothTerms, backend);
            List<VarianceComponent> components = new ArrayList<>(
                compiled.terms().size() + additional.size());
            for (PSplineMixedModelCompiler.Term term : compiled.terms()) {
                components.add(new VarianceComponent(term.term().name(), rows,
                    PSplineMixedModelCompiler.covarianceBasis(term)));
            }
            components.addAll(additional);
            GlmmPqlResult fitted = GlmmPql.fit(response,
                compiled.fixedDesign(), rows, compiled.fixedColumns(), family,
                components, weights, offsets, options, backendPolicy);
            PqlWorkingState.State state = PqlWorkingState.reconstruct(
                response, compiled.fixedDesign(), rows, compiled.fixedColumns(),
                family, weights, offsets, components, fitted, backend);
            return result(parametricColumns, compiled, additional,
                components, fitted, state, backend);
        }
    }

    private static GeneralizedGammResult result(
            int parametricColumns,
            PSplineMixedModelCompiler.Compiled compiled,
            List<VarianceComponent> additional,
            List<VarianceComponent> components,
            GlmmPqlResult fitted,
            PqlWorkingState.State state,
            ComputeBackend backend) {
        double[] variances = fitted.varianceComponents();
        List<SmoothTermEstimate> smoothEstimates = new ArrayList<>();
        int fixedStart = parametricColumns;
        int componentIndex = 0;
        double totalEdf = parametricColumns;
        for (PSplineMixedModelCompiler.Term term : compiled.terms()) {
            double[] fixedCoefficients = Arrays.copyOfRange(fitted.beta(),
                fixedStart, fixedStart + term.fixedColumns());
            fixedStart += term.fixedColumns();
            double variance = variances[componentIndex++];
            double[] randomCoefficients =
                PSplineMixedModelCompiler.randomCoefficients(
                    term, state.projectedResidual(), variance);
            double[] contribution = PSplineMixedModelCompiler.contribution(
                term, fixedCoefficients, randomCoefficients);
            double edf = term.fixedColumns()
                + PSplineMixedModelCompiler.randomEdf(
                    term, state.projection(), variance);
            edf = Math.max(term.fixedColumns(),
                Math.min(term.fixedColumns() + term.randomColumns(), edf));
            totalEdf += edf;
            smoothEstimates.add(new SmoothTermEstimate(term.term(),
                term.fixedTransform(), term.fixedMeans(), fixedCoefficients,
                term.randomTransform(), term.randomMeans(), randomCoefficients,
                contribution, 1.0 / variance, edf));
        }
        Map<String, double[]> randomPredictors = new LinkedHashMap<>();
        for (VarianceComponent component : additional) {
            double variance = variances[componentIndex++];
            double[] contribution = MatrixOps.multiply(backend,
                component.covariance(), component.dimension(),
                component.dimension(), state.projectedResidual());
            for (int row = 0; row < contribution.length; row++) {
                contribution[row] *= variance;
            }
            randomPredictors.put(component.name(), contribution);
            totalEdf += variance * traceProduct(component.covariance(),
                state.projection(), component.dimension());
        }
        return new GeneralizedGammResult(fitted, parametricColumns,
            smoothEstimates, randomPredictors, totalEdf);
    }

    private static double traceProduct(
            double[] left, double[] right, int dimension) {
        double result = 0.0;
        for (int row = 0; row < dimension; row++) {
            for (int column = 0; column < dimension; column++) {
                result += left[row * dimension + column]
                    * right[column * dimension + row];
            }
        }
        return result;
    }

    private static double[] weights(double[] supplied, int rows) {
        if (supplied == null) {
            double[] result = new double[rows];
            Arrays.fill(result, 1.0);
            return result;
        }
        if (supplied.length != rows) {
            throw new IllegalArgumentException("prior weight length must equal rows");
        }
        double[] result = MatrixOps.finiteCopy(supplied, "priorWeights");
        for (double value : result) {
            if (!(value > 0.0)) {
                throw new IllegalArgumentException(
                    "prior weights must be strictly positive");
            }
        }
        return result;
    }

    private static void validateComponents(
            List<VarianceComponent> components, int rows) {
        java.util.HashSet<String> names = new java.util.HashSet<>();
        for (VarianceComponent component : components) {
            if (component == null || component.dimension() != rows
                    || !names.add(component.name())) {
                throw new IllegalArgumentException(
                    "covariance components must have unique names and matching rows");
            }
        }
    }
}
