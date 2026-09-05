/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.pedigree.PedigreeIndividual;
import org.jlinalg.pedigree.SparsePedigreeReml;
import org.jlinalg.pedigree.SparsePedigreeRemlResult;
import org.jlinalg.reml.RemlOptions;

/** Deterministic sparse diagonal-PEV scale exercise. */
public final class SparsePedigreePevBenchmark {
    private SparsePedigreePevBenchmark() { }

    public static void main(String[] arguments) {
        int members = Integer.parseInt(System.getProperty(
            "jlinalg.benchmark.members", "500"));
        double[] response = new double[members * 2];
        double[][] fixed = new double[response.length][1];
        List<PedigreeIndividual> pedigreeRows = new ArrayList<>(members);
        List<String> observed = new ArrayList<>(response.length);
        for (int member = 0; member < members; member++) {
            String id = "member-" + member;
            pedigreeRows.add(PedigreeIndividual.founder(id));
            for (int repeat = 0; repeat < 2; repeat++) {
                int row = 2 * member + repeat;
                response[row] = member % 17 + (repeat == 0 ? -0.25 : 0.25);
                fixed[row][0] = 1.0;
                observed.add(id);
            }
        }
        long started = System.nanoTime();
        SparsePedigreeRemlResult result = SparsePedigreeReml.fit(response,
            fixed, observed, Pedigree.of(pedigreeRows), RemlOptions.builder()
                .initialVariances(10.0, 1.0).maximumIterations(60).build(),
            BackendPolicy.CPU);
        double seconds = (System.nanoTime() - started) / 1e9;
        long finite = java.util.Arrays.stream(result.predictionErrorVariances())
            .filter(Double::isFinite).count();
        System.out.println("members,observations,pev_batch_size,finite_pev,seconds");
        System.out.printf(java.util.Locale.ROOT, "%d,%d,32,%d,%.6f%n",
            members, response.length, finite, seconds);
    }
}
