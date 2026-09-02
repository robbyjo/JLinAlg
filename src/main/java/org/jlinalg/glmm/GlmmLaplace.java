/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.glmm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.SymmetricEigenDecomposition;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.VarianceComponent;

/** Dense first-order Laplace GLMM with arbitrary positive-semidefinite covariance bases. */
public final class GlmmLaplace {
    private GlmmLaplace() { }

    public static GlmmLaplaceResult fit(
            double[] response,
            double[][] fixedEffects,
            GlmFamily family,
            List<VarianceComponent> components) {
        return fit(response, MatrixOps.rowMajor(fixedEffects, response.length),
            response.length, fixedEffects[0].length, family, components,
            null, null, GlmmLaplaceOptions.defaults(), BackendPolicy.PREFERRED);
    }

    /** Fits a Laplace marginal likelihood and profiles every covariance variance. */
    public static GlmmLaplaceResult fit(
            double[] response,
            double[] fixedEffects,
            int observations,
            int fixedColumns,
            GlmFamily family,
            List<VarianceComponent> components,
            double[] priorWeights,
            double[] offset,
            GlmmLaplaceOptions options,
            BackendPolicy backendPolicy) {
        MatrixOps.validateModelData(response, fixedEffects, observations, fixedColumns);
        if (family == null || components == null || components.isEmpty()
                || options == null || backendPolicy == null) {
            throw new IllegalArgumentException(
                "family, covariance components, controls, and backend are required");
        }
        double[] weights = weights(priorWeights, observations);
        double[] offsets = offset == null ? new double[observations]
            : MatrixOps.finiteCopy(offset, "offset");
        if (offsets.length != observations) {
            throw new IllegalArgumentException("offset length must match observations");
        }
        for (int row = 0; row < observations; row++) {
            family.validateResponse(response[row], weights[row]);
        }
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            List<ComponentFactor> factors = factors(components, observations, backend);
            double[] logVariances = initial(options, components.size());
            Mode best = mode(response, fixedEffects, observations, fixedColumns,
                family, factors, logVariances, weights, offsets, options, backend);
            int evaluations = 1;
            int outerIterations = 0;
            double step = options.initialLogVarianceStep();
            boolean outerConverged = false;
            for (int sweep = 1; sweep <= options.maximumOuterIterations(); sweep++) {
                outerIterations = sweep;
                boolean improved = false;
                for (int component = 0; component < logVariances.length; component++) {
                    double original = logVariances[component];
                    double selected = original;
                    Mode coordinateBest = best;
                    for (double direction : new double[] {-1.0, 1.0}) {
                        double trial = Math.max(Math.log(options.minimumVariance()),
                            Math.min(Math.log(options.maximumVariance()),
                                original + direction * step));
                        if (trial == original) continue;
                        logVariances[component] = trial;
                        Mode candidate = mode(response, fixedEffects,
                            observations, fixedColumns, family, factors,
                            logVariances, weights, offsets, options, backend);
                        evaluations++;
                        if (candidate.laplaceLogLikelihood()
                                > coordinateBest.laplaceLogLikelihood()) {
                            coordinateBest = candidate;
                            selected = trial;
                        }
                    }
                    logVariances[component] = selected;
                    if (coordinateBest.laplaceLogLikelihood()
                            > best.laplaceLogLikelihood()
                                + options.relativeTolerance()
                                    * (1.0 + Math.abs(best.laplaceLogLikelihood()))) {
                        best = coordinateBest;
                        improved = true;
                    }
                }
                if (!improved) step *= 0.5;
                if (step <= options.relativeTolerance() * 10.0) {
                    outerConverged = true;
                    break;
                }
            }
            best = mode(response, fixedEffects, observations, fixedColumns,
                family, factors, logVariances, weights, offsets, options, backend);
            evaluations++;
            return result(family, components, factors, logVariances,
                best, fixedColumns, outerIterations,
                outerConverged && best.converged());
        }
    }

    private static Mode mode(
            double[] response,
            double[] fixed,
            int observations,
            int fixedColumns,
            GlmFamily family,
            List<ComponentFactor> factors,
            double[] logVariances,
            double[] priorWeights,
            double[] offsets,
            GlmmLaplaceOptions options,
            ComputeBackend backend) {
        int randomColumns = factors.stream().mapToInt(ComponentFactor::columns).sum();
        int columns = fixedColumns + randomColumns;
        double[] design = design(fixed, observations, fixedColumns,
            factors, logVariances, columns);
        double[] coefficients = new double[columns];
        initializeIntercept(response, fixed, observations, fixedColumns,
            family, offsets, coefficients);
        double[] means = new double[observations];
        double[] linear = new double[observations];
        double[] information = null;
        boolean converged = false;
        int iterations = 0;
        for (int iteration = 1; iteration <= options.maximumModeIterations(); iteration++) {
            iterations = iteration;
            double[] workingWeights = new double[observations];
            double[] workingResponse = new double[observations];
            for (int row = 0; row < observations; row++) {
                double eta = offsets[row];
                for (int column = 0; column < columns; column++) {
                    eta += design[row * columns + column] * coefficients[column];
                }
                linear[row] = eta;
                double mean = family.inverseLink(eta);
                means[row] = mean;
                double derivative = family.meanDerivative(eta);
                double variance = family.variance(mean);
                double weight = priorWeights[row] * derivative * derivative / variance;
                workingWeights[row] = Math.max(1e-12, weight);
                workingResponse[row] = eta - offsets[row]
                    + (response[row] - mean) / derivative;
            }
            double[] weightedDesign = design.clone();
            double[] weightedResponse = workingResponse.clone();
            for (int row = 0; row < observations; row++) {
                double squareRoot = Math.sqrt(workingWeights[row]);
                weightedResponse[row] *= squareRoot;
                for (int column = 0; column < columns; column++) {
                    weightedDesign[row * columns + column] *= squareRoot;
                }
            }
            information = MatrixOps.transposeMultiply(backend,
                weightedDesign, observations, columns, weightedDesign, columns);
            for (int column = fixedColumns; column < columns; column++) {
                information[column * columns + column] += 1.0;
            }
            double[] right = MatrixOps.transposeMultiply(backend,
                weightedDesign, observations, columns, weightedResponse, 1);
            double[] candidate = backend.dpotrf(information, columns).solve(right);
            double change = relativeChange(coefficients, candidate);
            coefficients = candidate;
            if (change <= options.relativeTolerance()) {
                converged = true;
                break;
            }
        }
        double conditional = 0.0;
        for (int row = 0; row < observations; row++) {
            double eta = offsets[row];
            for (int column = 0; column < columns; column++) {
                eta += design[row * columns + column] * coefficients[column];
            }
            linear[row] = eta;
            means[row] = family.inverseLink(eta);
            conditional += family.logLikelihood(
                response[row], means[row], priorWeights[row], 1.0);
        }
        if (!Double.isFinite(conditional)) {
            throw new IllegalArgumentException(
                "Laplace fitting requires a family with a finite likelihood");
        }
        double randomQuadratic = 0.0;
        for (int column = fixedColumns; column < columns; column++) {
            randomQuadratic += coefficients[column] * coefficients[column];
        }
        conditional -= 0.5 * randomQuadratic;
        int random = columns - fixedColumns;
        double[] randomInformation = new double[random * random];
        for (int row = 0; row < random; row++) {
            System.arraycopy(information, (fixedColumns + row) * columns + fixedColumns,
                randomInformation, row * random, random);
        }
        double logDeterminant = backend.dpotrf(randomInformation, random)
            .logDeterminant();
        double laplace = conditional - 0.5 * logDeterminant;
        double[] covariance = backend.dpotrf(information, columns)
            .solve(MatrixOps.identity(columns), columns);
        return new Mode(coefficients, covariance, design, linear, means,
            laplace, iterations, converged);
    }

    private static GlmmLaplaceResult result(
            GlmFamily family,
            List<VarianceComponent> components,
            List<ComponentFactor> factors,
            double[] logVariances,
            Mode mode,
            int fixedColumns,
            int outerIterations,
            boolean converged) {
        int totalColumns = mode.coefficients().length;
        double[] beta = Arrays.copyOf(mode.coefficients(), fixedColumns);
        double[] fixedCovariance = new double[fixedColumns * fixedColumns];
        double[] standardErrors = new double[fixedColumns];
        for (int row = 0; row < fixedColumns; row++) {
            for (int column = 0; column < fixedColumns; column++) {
                fixedCovariance[row * fixedColumns + column] =
                    mode.covariance()[row * totalColumns + column];
            }
            standardErrors[row] = Math.sqrt(Math.max(0.0,
                fixedCovariance[row * fixedColumns + row]));
        }
        AssociationStatistics association = AssociationStatistics.normal(
            beta, standardErrors);
        double[] variances = new double[logVariances.length];
        List<String> names = new ArrayList<>(components.size());
        Map<String, double[]> predictors = new LinkedHashMap<>();
        int coefficientStart = fixedColumns;
        for (int component = 0; component < components.size(); component++) {
            variances[component] = Math.exp(logVariances[component]);
            names.add(components.get(component).name());
            ComponentFactor factor = factors.get(component);
            double[] contribution = new double[mode.linear().length];
            for (int row = 0; row < contribution.length; row++) {
                for (int column = 0; column < factor.columns(); column++) {
                    contribution[row] += mode.design()[row * totalColumns
                        + coefficientStart + column]
                        * mode.coefficients()[coefficientStart + column];
                }
            }
            predictors.put(components.get(component).name(), contribution);
            coefficientStart += factor.columns();
        }
        return new GlmmLaplaceResult(family.name(), names, variances,
            association, fixedCovariance, predictors, mode.linear(), mode.means(),
            mode.laplaceLogLikelihood(), outerIterations,
            mode.iterations(), converged);
    }

    private static List<ComponentFactor> factors(
            List<VarianceComponent> components,
            int observations,
            ComputeBackend backend) {
        List<ComponentFactor> result = new ArrayList<>(components.size());
        java.util.HashSet<String> names = new java.util.HashSet<>();
        for (VarianceComponent component : components) {
            if (component == null || component.dimension() != observations
                    || !names.add(component.name())) {
                throw new IllegalArgumentException(
                    "components must have unique names and matching observations");
            }
            SymmetricEigenDecomposition eigen = backend.dsyev(
                component.covariance(), observations);
            double[] values = eigen.eigenvalues();
            double[] vectors = eigen.eigenvectors();
            double maximum = Arrays.stream(values).map(Math::abs).max().orElse(0.0);
            double tolerance = 1e-10 * Math.max(1.0, maximum);
            int rank = 0;
            for (double value : values) {
                if (value < -tolerance) {
                    throw new IllegalArgumentException(
                        "covariance component is not positive semidefinite: "
                            + component.name());
                }
                if (value > tolerance) rank++;
            }
            if (rank == 0) {
                throw new IllegalArgumentException(
                    "covariance component has zero rank: " + component.name());
            }
            double[] factor = new double[observations * rank];
            int destination = 0;
            for (int source = 0; source < observations; source++) {
                if (values[source] <= tolerance) continue;
                double scale = Math.sqrt(values[source]);
                for (int row = 0; row < observations; row++) {
                    factor[row * rank + destination] =
                        vectors[row * observations + source] * scale;
                }
                destination++;
            }
            result.add(new ComponentFactor(factor, observations, rank));
        }
        return List.copyOf(result);
    }

    private static double[] design(
            double[] fixed,
            int observations,
            int fixedColumns,
            List<ComponentFactor> factors,
            double[] logVariances,
            int columns) {
        double[] result = new double[observations * columns];
        for (int row = 0; row < observations; row++) {
            System.arraycopy(fixed, row * fixedColumns,
                result, row * columns, fixedColumns);
        }
        int destination = fixedColumns;
        for (int component = 0; component < factors.size(); component++) {
            ComponentFactor factor = factors.get(component);
            double scale = Math.exp(0.5 * logVariances[component]);
            for (int row = 0; row < observations; row++) {
                for (int column = 0; column < factor.columns(); column++) {
                    result[row * columns + destination + column] = scale
                        * factor.values()[row * factor.columns() + column];
                }
            }
            destination += factor.columns();
        }
        return result;
    }

    private static void initializeIntercept(
            double[] response,
            double[] fixed,
            int observations,
            int fixedColumns,
            GlmFamily family,
            double[] offset,
            double[] coefficients) {
        double constant = fixed[0];
        boolean isConstant = constant != 0.0;
        double mean = 0.0;
        for (int row = 0; row < observations; row++) {
            isConstant &= fixed[row * fixedColumns] == constant;
            mean += family.initialMean(response[row]);
        }
        if (isConstant) {
            double averageOffset = Arrays.stream(offset).average().orElse(0.0);
            coefficients[0] = (family.link(mean / observations)
                - averageOffset) / constant;
        }
    }

    private static double relativeChange(double[] current, double[] candidate) {
        double result = 0.0;
        for (int index = 0; index < current.length; index++) {
            result = Math.max(result, Math.abs(candidate[index] - current[index])
                / (1.0 + Math.abs(current[index])));
        }
        return result;
    }

    private static double[] initial(GlmmLaplaceOptions options, int count) {
        double[] supplied = options.initialVariances();
        if (supplied != null && supplied.length != count) {
            throw new IllegalArgumentException("one initial variance is required per component");
        }
        double[] result = new double[count];
        if (supplied == null) return result;
        for (int index = 0; index < count; index++) result[index] = Math.log(supplied[index]);
        return result;
    }

    private static double[] weights(double[] supplied, int observations) {
        double[] result = supplied == null ? new double[observations]
            : MatrixOps.finiteCopy(supplied, "priorWeights");
        if (supplied == null) Arrays.fill(result, 1.0);
        if (result.length != observations) {
            throw new IllegalArgumentException("prior weight length must match observations");
        }
        for (double value : result) {
            if (!(value > 0.0)) {
                throw new IllegalArgumentException("prior weights must be positive");
            }
        }
        return result;
    }

    private record ComponentFactor(double[] values, int observations, int columns) { }
    private record Mode(
            double[] coefficients,
            double[] covariance,
            double[] design,
            double[] linear,
            double[] means,
            double laplaceLogLikelihood,
            int iterations,
            boolean converged) { }
}
