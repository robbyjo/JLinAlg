/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pedigree;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;
import org.jlinalg.inference.BootstrapFailure;
import org.jlinalg.inference.BootstrapOptions;
import org.jlinalg.inference.GaussianBootstrapResult;
import org.jlinalg.mixed.MixedModelSimulationMode;

/** Parametric bootstrap for a retained Gaussian pedigree animal model. */
public final class PedigreeBootstrap {
    private PedigreeBootstrap() { }

    public static GaussianBootstrapResult bootstrap(
            PreparedPedigreeReml model,
            PedigreeRemlResult fitted,
            BootstrapOptions options) {
        if (model == null || fitted == null || options == null)
            throw new IllegalArgumentException(
                "prepared model, fitted result, and bootstrap options are required");
        int simulations = options.simulations();
        double[][] fixed = new double[simulations][];
        double[][] variances = new double[simulations][];
        BootstrapFailure[] failures = new BootstrapFailure[simulations];
        double[] relationshipLower = PedigreeSimulation.relationshipLower(model);
        execute(simulations, options.parallelism(), simulation -> {
            try {
                Random random = new Random(replicateSeed(
                    options.randomSeed(), simulation));
                double[] response = PedigreeSimulation.simulateOne(
                    model, fitted, random, MixedModelSimulationMode.MARGINAL,
                    relationshipLower);
                PedigreeRemlResult refit = model.refit(fitted, response);
                if (!refit.reml().converged()) {
                    failures[simulation] = new BootstrapFailure(simulation,
                        "ConvergenceFailure",
                        refit.reml().convergenceMessage());
                    return;
                }
                fixed[simulation] = refit.beta();
                variances[simulation] = refit.reml().varianceComponents();
            } catch (RuntimeException exception) {
                failures[simulation] = new BootstrapFailure(simulation,
                    exception.getClass().getSimpleName(), message(exception));
            }
        });
        List<double[]> retainedFixed = new ArrayList<>();
        List<double[]> retainedVariances = new ArrayList<>();
        List<BootstrapFailure> retainedFailures = new ArrayList<>();
        for (int simulation = 0; simulation < simulations; simulation++) {
            if (fixed[simulation] == null)
                retainedFailures.add(failures[simulation]);
            else {
                retainedFixed.add(fixed[simulation]);
                retainedVariances.add(variances[simulation]);
            }
        }
        return new GaussianBootstrapResult(simulations, options.randomSeed(),
            options.confidenceLevel(), fitted.beta(),
            fitted.reml().componentNames(),
            fitted.reml().varianceComponents(),
            retainedFixed.toArray(double[][]::new),
            retainedVariances.toArray(double[][]::new), retainedFailures);
    }

    private static long replicateSeed(long seed, int simulation) {
        long value = seed + 0x9E3779B97F4A7C15L * (simulation + 1L);
        value = (value ^ (value >>> 30)) * 0xBF58476D1CE4E5B9L;
        value = (value ^ (value >>> 27)) * 0x94D049BB133111EBL;
        return value ^ (value >>> 31);
    }

    private static String message(RuntimeException exception) {
        return exception.getMessage() == null
            ? exception.getClass().getSimpleName() : exception.getMessage();
    }

    private static void execute(
            int simulations, int parallelism,
            java.util.function.IntConsumer operation) {
        if (parallelism == 1) {
            for (int simulation = 0; simulation < simulations; simulation++)
                operation.accept(simulation);
            return;
        }
        ForkJoinPool pool = new ForkJoinPool(
            Math.min(simulations, parallelism));
        try {
            pool.submit(() -> IntStream.range(0, simulations).parallel()
                .forEach(operation)).join();
        } finally {
            pool.shutdown();
        }
    }
}
