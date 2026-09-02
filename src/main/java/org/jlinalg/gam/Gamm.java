/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTranspose;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.Reml;
import org.jlinalg.reml.RemlOptions;
import org.jlinalg.reml.RemlResult;
import org.jlinalg.reml.VarianceComponent;

/** Gaussian additive mixed models with arbitrary covariance components. */
public final class Gamm {
    private static final String RESIDUAL = "residual";

    private Gamm() { }

    /** Fits with independent residuals and default REML controls. */
    public static GammResult fitGaussian(
            double[] response,
            double[][] parametricDesign,
            List<PSplineTerm> smoothTerms,
            List<VarianceComponent> randomComponents) {
        return fitGaussian(response, parametricDesign, smoothTerms,
            randomComponents, null, RemlOptions.defaults(),
            BackendPolicy.PREFERRED);
    }

    /**
     * Fits smooths together with grouped, pedigree, GRM, or other covariance
     * components and an optional residual correlation matrix.
     */
    public static GammResult fitGaussian(
            double[] response,
            double[][] parametricDesign,
            List<PSplineTerm> smoothTerms,
            List<VarianceComponent> randomComponents,
            double[] residualCorrelation,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        if (response == null) {
            throw new IllegalArgumentException("response is required");
        }
        double[] fixed = MatrixOps.rowMajor(parametricDesign, response.length);
        return fitGaussian(response, fixed, response.length,
            parametricDesign[0].length, smoothTerms, randomComponents,
            residualCorrelation, options, backendPolicy);
    }

    /** Contiguous row-major overload for allocation-sensitive callers. */
    public static GammResult fitGaussian(
            double[] response,
            double[] parametricDesign,
            int rows,
            int parametricColumns,
            List<PSplineTerm> smoothTerms,
            List<VarianceComponent> randomComponents,
            double[] residualCorrelation,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        MatrixOps.validateModelData(
            response, parametricDesign, rows, parametricColumns);
        if (options == null || backendPolicy == null) {
            throw new IllegalArgumentException(
                "options and backendPolicy are required");
        }
        List<VarianceComponent> additional = randomComponents == null
            ? List.of() : List.copyOf(randomComponents);
        validateComponents(additional, rows);
        double[] residual = residualCorrelation == null
            ? MatrixOps.identity(rows)
            : new VarianceComponent(RESIDUAL, rows, residualCorrelation)
                .covariance();

        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            PSplineMixedModelCompiler.Compiled compiled =
                PSplineMixedModelCompiler.compile(parametricDesign, rows,
                    parametricColumns, smoothTerms, backend);
            List<VarianceComponent> components = new ArrayList<>(
                compiled.terms().size() + additional.size() + 1);
            for (PSplineMixedModelCompiler.Term term : compiled.terms()) {
                components.add(new VarianceComponent(term.term().name(), rows,
                    PSplineMixedModelCompiler.covarianceBasis(term)));
            }
            components.addAll(additional);
            components.add(new VarianceComponent(RESIDUAL, rows, residual));
            RemlResult fitted = Reml.fit(response, compiled.fixedDesign(), rows,
                compiled.fixedColumns(), components, options, backendPolicy);
            return result(response, parametricColumns, compiled,
                additional, components, fitted, backend);
        }
    }

    private static GammResult result(
            double[] response,
            int parametricColumns,
            PSplineMixedModelCompiler.Compiled compiled,
            List<VarianceComponent> additional,
            List<VarianceComponent> components,
            RemlResult fitted,
            ComputeBackend backend) {
        int rows = response.length;
        double[] variances = fitted.varianceComponents();
        double[] covariance = new double[rows * rows];
        for (int component = 0; component < components.size(); component++) {
            double[] basis = components.get(component).covariance();
            for (int index = 0; index < covariance.length; index++) {
                covariance[index] += variances[component] * basis[index];
            }
        }
        CholeskyFactor factor = backend.dpotrf(covariance, rows);
        double[] projectedResidual = factor.solve(fitted.residuals());
        double[] projection = projection(
            factor, compiled.fixedDesign(), rows, compiled.fixedColumns(),
            fitted.fixedEffectCovariance(), backend);

        double[] conditionalFitted = MatrixOps.multiply(backend,
            compiled.fixedDesign(), rows, compiled.fixedColumns(), fitted.beta());
        double residualVariance = variances[variances.length - 1];
        List<SmoothTermEstimate> smoothEstimates = new ArrayList<>();
        int fixedStart = parametricColumns;
        double totalEdf = parametricColumns;
        int componentIndex = 0;
        for (PSplineMixedModelCompiler.Term term : compiled.terms()) {
            double[] fixedCoefficients = java.util.Arrays.copyOfRange(
                fitted.beta(), fixedStart, fixedStart + term.fixedColumns());
            fixedStart += term.fixedColumns();
            double variance = variances[componentIndex++];
            double[] randomCoefficients =
                PSplineMixedModelCompiler.randomCoefficients(
                    term, projectedResidual, variance);
            double[] contribution = PSplineMixedModelCompiler.contribution(
                term, fixedCoefficients, randomCoefficients);
            addRandomContribution(conditionalFitted, term, randomCoefficients);
            double edf = term.fixedColumns()
                + PSplineMixedModelCompiler.randomEdf(
                    term, projection, variance);
            edf = Math.max(term.fixedColumns(),
                Math.min(term.fixedColumns() + term.randomColumns(), edf));
            totalEdf += edf;
            smoothEstimates.add(new SmoothTermEstimate(term.term(),
                term.fixedTransform(), term.fixedMeans(), fixedCoefficients,
                term.randomTransform(), term.randomMeans(), randomCoefficients,
                contribution, residualVariance / variance, edf));
        }

        Map<String, double[]> randomContributions = new LinkedHashMap<>();
        for (VarianceComponent component : additional) {
            double variance = variances[componentIndex++];
            double[] contribution = MatrixOps.multiply(backend,
                component.covariance(), rows, rows, projectedResidual);
            for (int row = 0; row < rows; row++) contribution[row] *= variance;
            addInPlace(conditionalFitted, contribution);
            randomContributions.put(component.name(), contribution);
            totalEdf += variance * traceProduct(
                component.covariance(), projection, rows);
        }
        return new GammResult(fitted, parametricColumns, smoothEstimates,
            randomContributions, conditionalFitted,
            MatrixOps.subtract(response, conditionalFitted), totalEdf);
    }

    private static double[] projection(
            CholeskyFactor factor,
            double[] fixed,
            int rows,
            int columns,
            double[] fixedCovariance,
            ComputeBackend backend) {
        double[] inverse = factor.solve(MatrixOps.identity(rows), rows);
        double[] inverseFixed = factor.solve(fixed, columns);
        double[] temporary = MatrixOps.multiply(backend,
            inverseFixed, rows, columns, fixedCovariance, columns);
        double[] correction = new double[rows * rows];
        backend.dgemm(MatrixTranspose.NONE, MatrixTranspose.TRANSPOSE,
            rows, rows, columns, 1.0,
            temporary, inverseFixed, 0.0, correction);
        return MatrixOps.subtract(inverse, correction);
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

    private static void addInPlace(double[] destination, double[] source) {
        for (int index = 0; index < destination.length; index++) {
            destination[index] += source[index];
        }
    }

    private static void addRandomContribution(
            double[] destination,
            PSplineMixedModelCompiler.Term term,
            double[] coefficients) {
        for (int row = 0; row < destination.length; row++) {
            for (int column = 0; column < term.randomColumns(); column++) {
                destination[row] += term.randomDesign()[
                    row * term.randomColumns() + column]
                    * coefficients[column];
            }
        }
    }

    private static void validateComponents(
            List<VarianceComponent> components, int rows) {
        java.util.HashSet<String> names = new java.util.HashSet<>();
        for (VarianceComponent component : components) {
            if (component == null || component.dimension() != rows
                    || RESIDUAL.equals(component.name())
                    || !names.add(component.name())) {
                throw new IllegalArgumentException(
                    "additional covariance components must have unique names "
                        + "and matching dimensions");
            }
        }
    }
}
