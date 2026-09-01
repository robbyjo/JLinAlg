/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import org.jlinalg.gwas.RemlAssociationScanner;

/**
 * GENESIS-style related-sample Gaussian score null backed by one retained
 * REML/P3D projection.
 */
public final class RemlSetTestNullModel implements GaussianSetTestNullModel {
    private final RemlAssociationScanner scanner;
    private final double[] projectedResponse;

    public RemlSetTestNullModel(RemlAssociationScanner scanner) {
        if (scanner == null)
            throw new IllegalArgumentException("REML association scanner is required");
        this.scanner = scanner;
        projectedResponse = scanner.projectedResponse();
    }

    @Override public int observations() { return scanner.observations(); }
    @Override public double degreesOfFreedom() {
        return scanner.associationDegreesOfFreedom();
    }

    @Override
    public SetTestScoreState score(double[][] variantRows) {
        if (variantRows == null || variantRows.length == 0)
            throw new IllegalArgumentException("variant rows are required");
        int variants = variantRows.length;
        double[] sampleByVariant = new double[observations() * variants];
        for (int variant = 0; variant < variants; variant++) {
            if (variantRows[variant] == null
                    || variantRows[variant].length != observations())
                throw new IllegalArgumentException(
                    "variant rows must match null-model observations");
            for (int sample = 0; sample < observations(); sample++) {
                double value = variantRows[variant][sample];
                if (!Double.isFinite(value))
                    throw new IllegalArgumentException(
                        "variant rows must be finite");
                sampleByVariant[sample * variants + variant] = value;
            }
        }
        double[] projected = scanner.project(sampleByVariant, variants);
        double[] scores = new double[variants];
        double[] information = new double[variants * variants];
        for (int left = 0; left < variants; left++) {
            for (int sample = 0; sample < observations(); sample++)
                scores[left] += variantRows[left][sample]
                    * projectedResponse[sample];
            for (int right = 0; right <= left; right++) {
                double value = 0;
                for (int sample = 0; sample < observations(); sample++)
                    value += variantRows[left][sample]
                        * projected[sample * variants + right];
                information[left * variants + right] = value;
                information[right * variants + left] = value;
            }
        }
        return new SetTestScoreState(scores, information, variants);
    }
}
