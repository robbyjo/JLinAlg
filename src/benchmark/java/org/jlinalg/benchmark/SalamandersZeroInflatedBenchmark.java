/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Supplier;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.distributional.SparseZeroInflatedMixedModel;
import org.jlinalg.distributional.ZeroInflatedMixedOptions;
import org.jlinalg.distributional.ZeroInflatedMixedResult;
import org.jlinalg.distributional.ZeroInflatedOuterOptimizer;
import org.jlinalg.mixed.RandomEffectTerm;

/** Warmed fit-only timings with separately reported inference for R-exported data. */
public final class SalamandersZeroInflatedBenchmark {
    private SalamandersZeroInflatedBenchmark() { }

    public static void main(String[] args) throws Exception {
        Path directory = Path.of(System.getProperty("jlinalg.benchmark.salamanders.data",
            "build/tmp/salamander-comparison"));
        Data zip = read(directory.resolve("zip.tsv"));
        Data zinb = read(directory.resolve("zinb.tsv"));
        System.out.println("java_version=" + System.getProperty("java.version"));
        System.out.println("gradient_threads=" + Integer.getInteger(
            "jlinalg.benchmark.salamanders.threads", 4));
        System.out.println("rows=" + zip.rows());
        if (args.length == 0 || args[0].equals("zip")) run("zip", zip, false);
        if (args.length == 0 || args[0].equals("zinb")) run("zinb", zinb, true);
    }

    private static void run(String label, Data data, boolean negativeBinomial) {
        RandomEffectTerm site = RandomEffectTerm.randomIntercept("site", data.site());
        ZeroInflatedMixedOptions options = new ZeroInflatedMixedOptions(
            2000, 150, 1e-8, 0.3, 1e-9, 1e5, 1e-6, 1e6, 30.0,
            new double[] {0.15},
            ZeroInflatedOuterOptimizer.AUTO,
            Integer.getInteger("jlinalg.benchmark.salamanders.threads", 4));
        try (SparseZeroInflatedMixedModel.Prepared prepared = negativeBinomial
                ? SparseZeroInflatedMixedModel.prepareNegativeBinomial(
                    data.rows(), List.of(site), null, List.of(), null,
                    List.of(), options, BackendPolicy.CPU)
                : SparseZeroInflatedMixedModel.preparePoisson(
                    data.rows(), List.of(site), null, List.of(), null,
                    List.of(), options, BackendPolicy.CPU)) {
            ZeroInflatedMixedResult inferred = negativeBinomial
                ? prepared.fitWithInference(data.y(), data.x(), data.columns(),
                    data.x(), data.columns(), data.intercept(), 1, null)
                : prepared.fitWithInference(data.y(), data.x(), data.columns(),
                    data.x(), data.columns(), null, 0, null);
            System.out.println(label + "_count=" + Arrays.toString(inferred.countCoefficients()));
            System.out.println(label + "_zero=" + Arrays.toString(inferred.zeroCoefficients()));
            System.out.println(label + "_count_sd=" + Math.sqrt(inferred.varianceComponents()[0]));
            if (negativeBinomial) System.out.println(label + "_size=" + inferred.sizes()[0]);
            System.out.println(label + "_loglik=" + inferred.marginalLogLikelihood());
            System.out.println(label + "_se=" + Arrays.toString(inferred.standardErrors()));
            System.out.println(label + "_converged=" + inferred.converged());
            System.out.println(label + "_evaluations=" + inferred.objectiveEvaluations());
            System.out.println(label + "_mode_iterations=" + inferred.modeIterations());
            requireConverged(inferred);

            Supplier<ZeroInflatedMixedResult> preparedFit = () -> {
                if (negativeBinomial) return prepared.fit(data.y(), data.x(), data.columns(),
                    data.x(), data.columns(), data.intercept(), 1, null);
                else return prepared.fit(data.y(), data.x(), data.columns(), data.x(),
                    data.columns(), null, 0, null);
            };
            System.out.println(label + "_prepared_time=" + timing(preparedFit));
        }

        Supplier<ZeroInflatedMixedResult> fullFit = () -> {
            if (negativeBinomial) return SparseZeroInflatedMixedModel.fitNegativeBinomial(
                data.y(), data.x(), data.columns(), data.x(), data.columns(),
                data.intercept(), 1, List.of(site), null, null, options,
                BackendPolicy.CPU);
            else return SparseZeroInflatedMixedModel.fitPoisson(data.y(), data.x(),
                data.columns(), data.x(), data.columns(), List.of(site), null,
                null, options, BackendPolicy.CPU);
        };
        System.out.println(label + "_full_time=" + timing(fullFit));
    }

    private static void requireConverged(ZeroInflatedMixedResult fit) {
        if (!fit.converged() || !Double.isFinite(fit.marginalLogLikelihood()))
            throw new IllegalStateException(fit.convergenceMessage());
    }

    private static String timing(
            Supplier<ZeroInflatedMixedResult> operation) {
        int warmups = Integer.getInteger("jlinalg.benchmark.salamanders.warmups", 2);
        int runs = Integer.getInteger("jlinalg.benchmark.salamanders.runs", 7);
        if (warmups < 0 || runs < 1)
            throw new IllegalArgumentException("invalid timing controls");
        for (int index = 0; index < warmups; index++) requireConverged(operation.get());
        double[] seconds = new double[runs];
        for (int index = 0; index < runs; index++) {
            long started = System.nanoTime();
            ZeroInflatedMixedResult fit = operation.get();
            seconds[index] = (System.nanoTime() - started) / 1e9;
            requireConverged(fit);
        }
        String samples = Arrays.toString(seconds);
        Arrays.sort(seconds);
        return "median:" + (runs % 2 == 0
            ? (seconds[runs / 2 - 1] + seconds[runs / 2]) / 2.0 : seconds[runs / 2]) + ",min:" + seconds[0]
            + ",max:" + seconds[runs - 1] + ",samples:" + samples;
    }

    private static Data read(Path path) throws Exception {
        List<String> lines = Files.readAllLines(path);
        int columns = lines.get(0).split("\t", -1).length - 2;
        int rows = lines.size() - 1;
        double[] y = new double[rows];
        double[] x = new double[rows * columns];
        double[] intercept = new double[rows];
        List<String> site = new ArrayList<>(rows);
        for (int row = 0; row < rows; row++) {
            String[] fields = lines.get(row + 1).split("\t", -1);
            y[row] = Double.parseDouble(fields[0]);
            site.add(fields[1]);
            intercept[row] = 1.0;
            for (int column = 0; column < columns; column++)
                x[row * columns + column] = Double.parseDouble(fields[column + 2]);
        }
        return new Data(rows, columns, y, x, intercept, List.copyOf(site));
    }

    private record Data(int rows, int columns, double[] y, double[] x,
                        double[] intercept, List<String> site) { }
}
