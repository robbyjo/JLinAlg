/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicLong;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.distributional.CorrelatedZeroInflatedRandomEffect;
import org.jlinalg.distributional.SparseZeroInflatedMixedModel;
import org.jlinalg.distributional.ZeroInflatedMixedOptions;
import org.jlinalg.distributional.ZeroInflatedMixedResult;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.pedigree.PedigreeIndividual;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;

/** Reproducible sparse correlated-pedigree ZIP performance and heap gate. */
public final class ZeroInflatedMixedBenchmark {
    private ZeroInflatedMixedBenchmark() { }

    public static void main(String[] arguments) throws InterruptedException {
        int members = Integer.parseInt(System.getProperty(
            "jlinalg.benchmark.members", "500"));
        int repeats = Integer.parseInt(System.getProperty(
            "jlinalg.benchmark.repeats", "3"));
        int rows = members * repeats;
        List<PedigreeIndividual> entries = pedigree(members);
        List<String> ids = new ArrayList<>(rows);
        double[] response = new double[rows];
        double[] count = new double[2 * rows];
        double[] zero = new double[rows];
        Random random = new Random(20260905L);
        for (int row = 0; row < rows; row++) {
            int member = row / repeats;
            double x = repeats == 1 ? 0.0
                : -1.0 + 2.0 * (row % repeats) / (repeats - 1.0);
            ids.add("member-" + member);
            count[2 * row] = 1.0;
            count[2 * row + 1] = x;
            zero[row] = 1.0;
            double mean = Math.exp(0.35 + 0.3 * x
                + 0.35 * Math.sin(member * 0.13));
            double pi = 1.0 / (1.0 + Math.exp(-(-1.2
                + 0.5 * Math.sin(member * 0.13))));
            response[row] = random.nextDouble() < pi
                ? 0.0 : poisson(random, mean);
        }
        PedigreeRandomEffectTerm pedigree = PedigreeRandomEffectTerm.of(
            "pedigree", ids, Pedigree.of(entries));
        CorrelatedZeroInflatedRandomEffect correlated =
            CorrelatedZeroInflatedRandomEffect.pedigree("additive", pedigree);
        long baseline = usedHeap();
        AtomicLong peak = new AtomicLong(baseline);
        AtomicBoolean running = new AtomicBoolean(true);
        Thread sampler = new Thread(() -> {
            while (running.get()) {
                peak.accumulateAndGet(usedHeap(), Math::max);
                try {
                    Thread.sleep(5L);
                } catch (InterruptedException exception) {
                    Thread.currentThread().interrupt();
                    return;
                }
            }
        }, "zero-inflated-heap-sampler");
        sampler.setDaemon(true);
        sampler.start();
        long started = System.nanoTime();
        ZeroInflatedMixedResult fit = SparseZeroInflatedMixedModel.fitPoisson(
            response, count, 2, zero, 1,
            List.of(), List.of(), List.of(), List.of(), List.of(correlated),
            null, new ZeroInflatedMixedOptions(600, 100, 1e-6, 0.35,
                1e-7, 1e3, 1e-5, 1e5, 25.0,
                new double[] {0.3, 0.5}), BackendPolicy.CPU);
        double seconds = (System.nanoTime() - started) / 1e9;
        running.set(false);
        sampler.join();
        peak.accumulateAndGet(usedHeap(), Math::max);
        System.out.println("members,observations,random_coefficients,equation_nonzeros,"
            + "factor_nonzeros,seconds,baseline_heap_bytes,peak_heap_bytes,"
            + "peak_heap_delta_bytes,converged");
        System.out.printf(java.util.Locale.ROOT,
            "%d,%d,%d,%d,%d,%.6f,%d,%d,%d,%s%n", members, rows,
            fit.randomCoefficientCount(), fit.sparseEquationNonzeroCount(),
            fit.factorNonzeroCount(), seconds, baseline, peak.get(),
            peak.get() - baseline, fit.converged());
    }

    private static List<PedigreeIndividual> pedigree(int members) {
        List<PedigreeIndividual> result = new ArrayList<>(members);
        int founders = Math.min(members, Math.max(2, members / 5));
        for (int member = 0; member < founders; member++)
            result.add(PedigreeIndividual.founder("member-" + member));
        for (int member = founders; member < members; member++) {
            int sire = (member - founders) % founders;
            int dam = (sire + Math.max(1, founders / 2)) % founders;
            result.add(new PedigreeIndividual("member-" + member,
                "member-" + sire, "member-" + dam));
        }
        return result;
    }

    private static int poisson(Random random, double mean) {
        double threshold = Math.exp(-mean);
        double product = 1.0;
        int count = -1;
        do {
            count++;
            product *= random.nextDouble();
        } while (product > threshold);
        return count;
    }

    private static long usedHeap() {
        Runtime runtime = Runtime.getRuntime();
        return runtime.totalMemory() - runtime.freeMemory();
    }
}
