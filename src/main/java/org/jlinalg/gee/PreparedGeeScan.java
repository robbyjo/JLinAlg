/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamily;

/**
 * Reuses cluster ordering, the null-model start, and one backend context while
 * scanning many observation-aligned predictors through the same marginal GEE.
 */
public final class PreparedGeeScan implements AutoCloseable {
    private final PreparedGeeData base;
    private final GlmFamily family;
    private final GeeOptions options;
    private final BackendContext context;
    private final double[] initialCoefficients;
    private boolean closed;

    private PreparedGeeScan(
            double[] response,
            double[] baseDesign,
            int rows,
            int columns,
            int[] cluster,
            int[] repeated,
            double[] weights,
            double[] offset,
            GlmFamily family,
            GeeOptions options,
            BackendPolicy backendPolicy) {
        if (family == null || options == null || backendPolicy == null) {
            throw new IllegalArgumentException("family and scan options are required");
        }
        if (options.covariance() != GeeCovariance.NAIVE
                && options.covariance() != GeeCovariance.ROBUST
                && options.covariance() != GeeCovariance.DF_ADJUSTED) {
            throw new IllegalArgumentException(
                "prepared scans support naive, robust, or DF-adjusted covariance");
        }
        if (options.exactClusterDeletion()) {
            throw new IllegalArgumentException(
                "prepared scans do not materialize delete-cluster fits");
        }
        if (options.initialCoefficients() != null) {
            throw new IllegalArgumentException(
                "prepared scans compute their own reusable null-model start");
        }
        this.family = family;
        this.options = options.toBuilder().parallelism(1).build();
        base = Gee.prepare(response, baseDesign, rows, columns, cluster,
            repeated, weights, offset, this.options, family);
        BackendContext selected = BackendContext.select(backendPolicy);
        try {
            double[] start = Gee.startingCoefficients(base, family,
                this.options, backendPolicy);
            GeeResult nullFit = Gee.fitPrepared(base, family, this.options,
                start, selected.backend(), selected.provenance(), false, true);
            initialCoefficients = new double[columns + 1];
            System.arraycopy(nullFit.coefficients(), 0, initialCoefficients,
                0, columns);
            context = selected;
        } catch (RuntimeException | Error exception) {
            selected.close();
            throw exception;
        }
    }

    /** Prepares fixed covariates and clustered observations for predictor scans. */
    public static PreparedGeeScan prepare(
            double[] response,
            double[] baseDesign,
            int rows,
            int columns,
            int[] cluster,
            int[] repeated,
            double[] weights,
            double[] offset,
            GlmFamily family,
            GeeOptions options,
            BackendPolicy backendPolicy) {
        return new PreparedGeeScan(response, baseDesign, rows, columns,
            cluster, repeated, weights, offset, family, options, backendPolicy);
    }

    /** Fits one marginal GEE per predictor; predictors are row-major by observation. */
    public GeeScanResult scan(
            double[] predictors,
            int predictorCount,
            List<String> predictorNames,
            int parallelism) {
        if (closed) throw new IllegalStateException("prepared GEE scan is closed");
        if (predictorCount < 1 || parallelism < 1) {
            throw new IllegalArgumentException(
                "predictorCount and parallelism must be positive");
        }
        List<String> names = names(predictorNames, predictorCount);
        double[] coefficients = new double[predictorCount];
        double[] standardErrors = new double[predictorCount];
        double[] statistics = new double[predictorCount];
        double[] pValues = new double[predictorCount];
        double[] association = new double[predictorCount];
        int[] iterations = new int[predictorCount];
        boolean[] converged = new boolean[predictorCount];
        GeeOptions fitOptions = options.toBuilder()
            .initialCoefficients(initialCoefficients).build();
        long started = System.nanoTime();
        java.util.function.IntConsumer fit = predictor -> {
            PreparedGeeData data = base.withAppendedPredictor(
                predictors, predictorCount, predictor);
            GeeResult result = Gee.fitPrepared(data, family, fitOptions,
                initialCoefficients.clone(), context.backend(),
                context.provenance(), false, true);
            int coefficient = data.columns() - 1;
            coefficients[predictor] = result.coefficients()[coefficient];
            standardErrors[predictor] = result.standardErrors()[coefficient];
            statistics[predictor] = result.statistics()[coefficient];
            pValues[predictor] = result.pValues()[coefficient];
            double[] parameters = result.associationParameters();
            association[predictor] = parameters.length == 0
                ? Double.NaN : parameters[0];
            iterations[predictor] = result.iterations();
            converged[predictor] = result.converged();
        };
        int workers = Math.min(parallelism, predictorCount);
        if (workers == 1) {
            for (int predictor = 0; predictor < predictorCount; predictor++) {
                fit.accept(predictor);
            }
        } else {
            ForkJoinPool pool = new ForkJoinPool(workers);
            try {
                pool.submit(() -> IntStream.range(0, predictorCount)
                    .parallel().forEach(fit)).get();
            } catch (InterruptedException exception) {
                Thread.currentThread().interrupt();
                throw new IllegalStateException("parallel GEE scan interrupted", exception);
            } catch (ExecutionException exception) {
                throw new IllegalStateException(
                    "parallel GEE scan failed", exception.getCause());
            } finally {
                pool.shutdown();
            }
        }
        return new GeeScanResult(names, coefficients, standardErrors,
            statistics, pValues, association, iterations, converged, workers,
            System.nanoTime() - started, context.provenance());
    }

    public int observations() { return base.observations(); }
    public int clusters() { return base.clusters(); }

    @Override
    public void close() {
        if (closed) return;
        closed = true;
        context.close();
    }

    private static List<String> names(List<String> supplied, int count) {
        if (supplied != null) {
            if (supplied.size() != count) {
                throw new IllegalArgumentException(
                    "predictor name count must equal predictorCount");
            }
            return List.copyOf(supplied);
        }
        List<String> result = new ArrayList<>(count);
        for (int index = 0; index < count; index++) {
            result.add("predictor_" + index);
        }
        return List.copyOf(result);
    }
}
