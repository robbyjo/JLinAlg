/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.benchmark;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.association.AssociationBatchResult;
import org.jlinalg.association.AssociationEngineOptions;
import org.jlinalg.association.PreparedGamAssociation;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.distributional.DistributionalFamilies;
import org.jlinalg.distributional.DistributionalOptions;
import org.jlinalg.distributional.DistributionalResult;
import org.jlinalg.distributional.PreparedDistributionalModel;
import org.jlinalg.gam.DiscretePSplineBasis;
import org.jlinalg.gam.GaussianSmoothSelectionResult;
import org.jlinalg.gam.GaussianSmoothSelector;
import org.jlinalg.gam.PSplineTerm;
import org.jlinalg.gam.PenalizedPredictor;
import org.jlinalg.gam.QuadraticSmoothTerm;

/** Deterministic cold/prepared/discrete/association additive-model macrobenchmark. */
public final class GamBenchmark {
    private GamBenchmark() { }

    public static void main(String[] arguments) {
        int observations = Integer.getInteger("jlinalg.benchmark.observations", 1200);
        int markers = Integer.getInteger("jlinalg.benchmark.markers", 256);
        double[] x = new double[observations];
        double[] response = new double[observations];
        double[][] intercept = new double[observations][1];
        double[][] genotype = new double[observations][markers];
        List<String> names = new ArrayList<>(markers);
        for (int marker = 0; marker < markers; marker++) names.add("m" + marker);
        for (int row = 0; row < observations; row++) {
            x[row] = (row % 200) / 199.0;
            response[row] = 1.0 + Math.sin(2.0 * Math.PI * x[row])
                + 0.1 * Math.sin(17.0 * row);
            intercept[row][0] = 1.0;
            for (int marker = 0; marker < markers; marker++) {
                genotype[row][marker] =
                    (row * (2 * marker + 1) + row / 7 + marker) % 3;
            }
        }
        long basisStarted = System.nanoTime();
        PSplineTerm term = PSplineTerm.of("s(x)", x, 14);
        long basisNanos = System.nanoTime() - basisStarted;
        DiscretePSplineBasis discrete = DiscretePSplineBasis.compile(term);
        long discreteStarted = System.nanoTime();
        discrete.crossProduct(null);
        discrete.transposeMultiply(response, null);
        long discreteNanos = System.nanoTime() - discreteStarted;

        long fitStarted = System.nanoTime();
        GaussianSmoothSelectionResult nullFit = GaussianSmoothSelector.fitFixed(
            response, intercept, List.of(QuadraticSmoothTerm.from(term)),
            List.of(new double[] {1.0}), BackendPolicy.PREFERRED);
        long fitNanos = System.nanoTime() - fitStarted;

        AssociationBatchResult scan;
        long scanStarted = System.nanoTime();
        try (PreparedGamAssociation prepared = new PreparedGamAssociation(nullFit,
                AssociationEngineOptions.acceleratedSerial())) {
            scan = prepared.scan(genotype, names);
        }
        long scanNanos = System.nanoTime() - scanStarted;

        PenalizedPredictor mean = nullFit.predictor();
        PenalizedPredictor scale = PenalizedPredictor.linear(intercept);
        PreparedDistributionalModel distributional = new PreparedDistributionalModel(
            List.of(mean, scale), DistributionalFamilies.gaussianLocationScale(),
            DistributionalOptions.defaults(), BackendPolicy.PREFERRED);
        DistributionalResult first = distributional.fit(response);
        long warmStarted = System.nanoTime();
        DistributionalResult warm = distributional.refit(first, response);
        long warmNanos = System.nanoTime() - warmStarted;

        System.out.printf(java.util.Locale.ROOT,
            "observations=%d markers=%d basis_ms=%.3f discrete_ms=%.3f "
                + "fixed_fit_ms=%.3f scan_ms=%.3f warm_distributional_ms=%.3f "
                + "unique_ratio=%.5f scan_ok=%s warm_iterations=%d%n",
            observations, markers, basisNanos / 1e6, discreteNanos / 1e6,
            fitNanos / 1e6, scanNanos / 1e6, warmNanos / 1e6,
            discrete.compressionRatio(), scan.successful(), warm.iterations());
    }
}
