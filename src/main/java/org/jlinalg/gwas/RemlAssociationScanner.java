/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gwas;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTranspose;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.Reml;
import org.jlinalg.reml.RemlOptions;
import org.jlinalg.reml.RemlResult;
import org.jlinalg.reml.VarianceComponent;

/**
 * P3D/EMMAX-style GWAS or TWAS scanner. Variance components and the null-model
 * projection are computed once; marker blocks then use accelerated matrix
 * multiplication without refitting REML per marker.
 */
public final class RemlAssociationScanner {
    private final double[] projection;
    private final double[] projectedResponse;
    private final int observations;
    private final double degreesOfFreedom;
    private final RemlResult nullModel;
    private final BackendPolicy backendPolicy;

    private RemlAssociationScanner(
            double[] projection,
            double[] projectedResponse,
            int observations,
            double degreesOfFreedom,
            RemlResult nullModel,
            BackendPolicy backendPolicy) {
        this.projection = projection;
        this.projectedResponse = projectedResponse;
        this.observations = observations;
        this.degreesOfFreedom = degreesOfFreedom;
        this.nullModel = nullModel;
        this.backendPolicy = backendPolicy;
    }

    /** Fits and factorizes a reusable mixed-model null model. */
    public static RemlAssociationScanner prepare(
            double[] response,
            double[][] covariates,
            List<VarianceComponent> components,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        if (response == null) throw new IllegalArgumentException("response is required");
        double[] fixed = MatrixOps.rowMajor(covariates, response.length);
        return prepare(response, fixed, response.length, covariates[0].length,
            components, options, backendPolicy);
    }

    /** Contiguous row-major null-model preparation. */
    public static RemlAssociationScanner prepare(
            double[] response,
            double[] covariates,
            int rows,
            int columns,
            List<VarianceComponent> components,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        MatrixOps.validateModelData(response, covariates, rows, columns);
        RemlResult nullModel = Reml.fit(response, covariates, rows, columns,
            components, options, backendPolicy);
        double[] variances = nullModel.varianceComponents();
        double[] covariance = new double[rows * rows];
        for (int component = 0; component < components.size(); component++) {
            double[] basis = components.get(component).covariance();
            for (int index = 0; index < covariance.length; index++) {
                covariance[index] += variances[component] * basis[index];
            }
        }
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            CholeskyFactor factor = backend.dpotrf(covariance, rows);
            double[] inverse = factor.solve(MatrixOps.identity(rows), rows);
            double[] inverseFixed = factor.solve(covariates, columns);
            double[] temporary = MatrixOps.multiply(backend,
                inverseFixed, rows, columns,
                nullModel.fixedEffectCovariance(), columns);
            double[] correction = new double[rows * rows];
            backend.dgemm(MatrixTranspose.NONE, MatrixTranspose.TRANSPOSE,
                rows, rows, columns, 1.0,
                temporary, inverseFixed, 0.0, correction);
            double[] projection = MatrixOps.subtract(inverse, correction);
            double[] projectedResponse = factor.solve(nullModel.residuals());
            double degrees = rows - columns - 1.0;
            if (!(degrees > 0.0)) {
                throw new IllegalArgumentException(
                    "association scan requires positive denominator DF");
            }
            return new RemlAssociationScanner(projection, projectedResponse,
                rows, degrees, nullModel, backendPolicy);
        }
    }

    /** Scans a conventional observation-by-marker dosage matrix. */
    public AssociationScanResult scan(double[][] markers, List<String> names) {
        return scan(markers, names, AssociationScanOptions.defaults());
    }

    /** Scans a dosage matrix with explicit batching and parallelism controls. */
    public AssociationScanResult scan(
            double[][] markers,
            List<String> names,
            AssociationScanOptions options) {
        if (markers == null || markers.length != observations) {
            throw new IllegalArgumentException("marker row count must equal observations");
        }
        double[] rowMajor = MatrixOps.rowMajorUnchecked(markers, observations);
        return scan(rowMajor, markers[0].length, names,
            options);
    }

    /** Scans contiguous row-major dosages in accelerated blocks. */
    public AssociationScanResult scan(
            double[] markers,
            int markerCount,
            List<String> names,
            AssociationScanOptions options) {
        if (markers == null || markerCount < 1
                || markers.length != observations * markerCount
                || options == null) {
            throw new IllegalArgumentException("marker dimensions or options are invalid");
        }
        List<String> markerNames = names == null
            ? defaultNames(markerCount) : List.copyOf(names);
        if (markerNames.size() != markerCount) {
            throw new IllegalArgumentException("one name is required per marker");
        }
        double[] prepared = prepareMarkers(
            markers, markerCount, options.missingPolicy());
        double[] beta = new double[markerCount];
        double[] standardErrors = new double[markerCount];
        int batches = (markerCount + options.batchSize() - 1)
            / options.batchSize();
        execute(batches, options.parallelism(), batch -> {
            int first = batch * options.batchSize();
            try (BackendContext context = BackendContext.select(backendPolicy)) {
                ComputeBackend backend = context.backend();
                int count = Math.min(options.batchSize(), markerCount - first);
                double[] block = columns(prepared, markerCount, first, count);
                double[] projected = MatrixOps.multiply(backend,
                    projection, observations, observations, block, count);
                for (int marker = 0; marker < count; marker++) {
                    double numerator = 0.0;
                    double information = 0.0;
                    for (int row = 0; row < observations; row++) {
                        double dosage = block[row * count + marker];
                        numerator += dosage * projectedResponse[row];
                        information += dosage * projected[row * count + marker];
                    }
                    int destination = first + marker;
                    if (!(information > 1e-14) || !Double.isFinite(information)) {
                        beta[destination] = Double.NaN;
                        standardErrors[destination] = Double.NaN;
                    } else {
                        beta[destination] = numerator / information;
                        standardErrors[destination] = Math.sqrt(1.0 / information);
                    }
                }
            }
        });
        AssociationStatistics statistics = AssociationStatistics.studentT(
            beta, standardErrors, degreesOfFreedom,
            DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION);
        return new AssociationScanResult(markerNames, statistics, nullModel);
    }

    public RemlResult nullModel() { return nullModel; }

    /**
     * Applies the retained mixed-model null projection P to a row-major
     * observation-by-variable matrix without refitting variance components.
     */
    public double[] project(double[] variables, int variableCount) {
        if (variables == null || variableCount < 1
                || variables.length != observations * variableCount)
            throw new IllegalArgumentException(
                "projection matrix dimensions are invalid");
        MatrixOps.requireFinite(variables, "projected variables");
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            return MatrixOps.multiply(context.backend(), projection,
                observations, observations, variables, variableCount);
        }
    }

    /** Returns P y from the retained mixed-model null fit. */
    public double[] projectedResponse() { return projectedResponse.clone(); }
    public int observations() { return observations; }
    public double associationDegreesOfFreedom() { return degreesOfFreedom; }
    public BackendPolicy backendPolicy() { return backendPolicy; }

    private static void execute(
            int batches, int parallelism, BatchOperation operation) {
        ForkJoinPool pool = new ForkJoinPool(Math.min(batches, parallelism));
        try {
            pool.submit(() -> IntStream.range(0, batches).parallel()
                .forEach(operation::run)).get();
        } catch (InterruptedException exception) {
            Thread.currentThread().interrupt();
            throw new IllegalStateException(
                "REML association scan was interrupted", exception);
        } catch (ExecutionException exception) {
            Throwable cause = exception.getCause();
            if (cause instanceof RuntimeException runtime) throw runtime;
            throw new IllegalStateException("REML association scan failed", cause);
        } finally {
            pool.shutdownNow();
        }
    }

    private double[] prepareMarkers(
            double[] markers, int markerCount, GenotypeMissingPolicy policy) {
        double[] result = markers;
        for (int marker = 0; marker < markerCount; marker++) {
            double sum = 0.0;
            int finite = 0;
            for (int row = 0; row < observations; row++) {
                double value = markers[row * markerCount + marker];
                if (Double.isFinite(value)) {
                    sum += value;
                    finite++;
                } else if (policy == GenotypeMissingPolicy.ERROR) {
                    throw new IllegalArgumentException(
                        "non-finite dosage for marker " + marker + ", row " + row);
                }
            }
            if (finite != observations) {
                if (finite == 0) {
                    throw new IllegalArgumentException(
                        "marker has no finite dosages: " + marker);
                }
                if (result == markers) result = markers.clone();
                double mean = sum / finite;
                for (int row = 0; row < observations; row++) {
                    int index = row * markerCount + marker;
                    if (!Double.isFinite(result[index])) result[index] = mean;
                }
            }
        }
        return result;
    }

    private double[] columns(
            double[] matrix, int totalColumns, int first, int count) {
        if (first == 0 && count == totalColumns) return matrix;
        double[] result = new double[observations * count];
        for (int row = 0; row < observations; row++) {
            System.arraycopy(matrix, row * totalColumns + first,
                result, row * count, count);
        }
        return result;
    }

    private static List<String> defaultNames(int count) {
        List<String> result = new ArrayList<>(count);
        for (int index = 0; index < count; index++) {
            result.add("marker" + (index + 1));
        }
        return result;
    }

    @FunctionalInterface
    private interface BatchOperation { void run(int batch); }
}
