/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.distributional.BetaMixedModel;
import org.jlinalg.distributional.BetaMixedModelOptions;
import org.jlinalg.distributional.BetaMixedModelResult;
import org.jlinalg.glmm.GlmmLaplaceOptions;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.pedigree.PedigreeIndividual;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;

/** Reproducible grouped and pedigree beta mixed-model macrobenchmark. */
public final class BetaMixedModelBenchmark {
    private BetaMixedModelBenchmark() { }

    public static void main(String[] arguments) {
        int rows = Integer.getInteger(
            "jlinalg.benchmark.beta.mixed.rows", 30_000);
        int groups = Integer.getInteger(
            "jlinalg.benchmark.beta.mixed.groups", 1_000);
        int repetitions = Integer.getInteger(
            "jlinalg.benchmark.beta.mixed.repetitions", 3);
        Data data = data(rows, groups);
        BetaMixedModelOptions options = new BetaMixedModelOptions(
            new GlmmLaplaceOptions(24, 60, 1e-5, 0.5,
                1e-7, 100.0, null), 12.0, 1e-5, 1e6);

        fitGrouped(data, options);
        fitPedigree(data, options);
        double[] groupedTimes = new double[repetitions];
        double[] pedigreeTimes = new double[repetitions];
        BetaMixedModelResult grouped = null;
        BetaMixedModelResult pedigree = null;
        for (int repetition = 0; repetition < repetitions; repetition++) {
            long start = System.nanoTime();
            grouped = fitGrouped(data, options);
            groupedTimes[repetition] = elapsedMillis(start);
            start = System.nanoTime();
            pedigree = fitPedigree(data, options);
            pedigreeTimes[repetition] = elapsedMillis(start);
        }
        System.out.printf("rows=%d groups=%d repetitions=%d%n",
            rows, groups, repetitions);
        report("grouped", grouped, median(groupedTimes));
        report("pedigree", pedigree, median(pedigreeTimes));
        System.out.printf("pedigree_overhead=%.2fx%n",
            median(pedigreeTimes) / median(groupedTimes));
    }

    private static BetaMixedModelResult fitGrouped(
            Data data, BetaMixedModelOptions options) {
        return BetaMixedModel.fit(data.response(), data.fixed(),
            data.response().length, 2,
            List.of(RandomEffectTerm.randomIntercept("group", data.groups())),
            null, null, options, BackendPolicy.CPU);
    }

    private static BetaMixedModelResult fitPedigree(
            Data data, BetaMixedModelOptions options) {
        List<PedigreeIndividual> pedigree = new ArrayList<>(data.groupCount() * 3);
        List<String> observed = new ArrayList<>(data.response().length);
        for (int group = 0; group < data.groupCount(); group++) {
            String sire = "s" + group;
            String dam = "d" + group;
            String child = "c" + group;
            pedigree.add(PedigreeIndividual.founder(sire));
            pedigree.add(PedigreeIndividual.founder(dam));
            pedigree.add(new PedigreeIndividual(child, sire, dam));
        }
        for (String group : data.groups()) observed.add("c" + group);
        PedigreeRandomEffectTerm term = PedigreeRandomEffectTerm.ofUninbred(
            "animal", observed, pedigree);
        return BetaMixedModel.fitPedigree(data.response(), data.fixed(),
            data.response().length, 2, term, null, null, options,
            BackendPolicy.CPU);
    }

    private static void report(String name, BetaMixedModelResult fit,
            double milliseconds) {
        System.out.printf(
            "%s median_ms=%.3f converged=%s beta=%s precision=%.6g variance=%.6g%n",
            name, milliseconds, fit.converged(),
            Arrays.toString(fit.meanCoefficients()), fit.precision(),
            fit.varianceComponents()[0]);
    }

    private static Data data(int rows, int groups) {
        double[] response = new double[rows];
        double[] fixed = new double[rows * 2];
        List<String> labels = new ArrayList<>(rows);
        for (int row = 0; row < rows; row++) {
            int group = row % groups;
            double x = -1.0 + 2.0 * ((row * 37L) % 1009L) / 1008.0;
            double random = 0.55 * Math.sin(1.7 * group)
                + 0.15 * Math.cos(0.31 * group);
            double mu = logistic(-0.3 + 0.8 * x + random);
            double standardDeviation = Math.sqrt(mu * (1.0 - mu) / 19.0);
            double residual = (Math.sin(1.73 * row + 0.2)
                + 0.55 * Math.cos(0.47 * row)) / 1.14;
            response[row] = Math.max(1e-5,
                Math.min(1.0 - 1e-5, mu + standardDeviation * residual));
            fixed[row * 2] = 1.0;
            fixed[row * 2 + 1] = x;
            labels.add(Integer.toString(group));
        }
        return new Data(response, fixed, labels, groups);
    }

    private static double logistic(double value) {
        return value >= 0.0 ? 1.0 / (1.0 + Math.exp(-value))
            : Math.exp(value) / (1.0 + Math.exp(value));
    }
    private static double elapsedMillis(long start) {
        return (System.nanoTime() - start) / 1_000_000.0;
    }
    private static double median(double[] values) {
        double[] sorted = values.clone();
        Arrays.sort(sorted);
        return sorted[sorted.length / 2];
    }
    private record Data(double[] response, double[] fixed,
            List<String> groups, int groupCount) { }
}
