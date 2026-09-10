/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mixed;

import java.util.List;
import java.util.Objects;

/** Conditional modes and prediction-error variances for one random term. */
public final class RandomEffectEstimates {
    private final String termName;
    private final List<String> coefficientNames;
    private final double variance;
    private final double[] estimates;
    private final double[] predictionErrorVariances;
    private final int predictionErrorBlockSize;
    private final double[] predictionErrorCovariances;

    RandomEffectEstimates(
            String termName,
            List<String> coefficientNames,
            double variance,
            double[] estimates,
            double[] predictionErrorVariances) {
        this(termName,coefficientNames,variance,estimates,
            predictionErrorVariances,1,predictionErrorVariances);
    }

    RandomEffectEstimates(String termName,List<String> coefficientNames,
            double variance,double[] estimates,double[] predictionErrorVariances,
            int predictionErrorBlockSize,double[] predictionErrorCovariances) {
        this.termName = Objects.requireNonNull(termName, "termName");
        this.coefficientNames = List.copyOf(coefficientNames);
        this.variance = variance;
        this.estimates = estimates.clone();
        this.predictionErrorVariances = predictionErrorVariances.clone();
        this.predictionErrorBlockSize=predictionErrorBlockSize;
        this.predictionErrorCovariances=predictionErrorCovariances.clone();
    }

    public String termName() { return termName; }
    public List<String> coefficientNames() { return coefficientNames; }
    public double variance() { return variance; }
    public double[] estimates() { return estimates.clone(); }
    public double[] predictionErrorVariances() {
        return predictionErrorVariances.clone();
    }
    /** Group/block-major PEV matrices, or diagonal 1-by-1 blocks. */
    public double[] predictionErrorCovariances(){return predictionErrorCovariances.clone();}
    public int predictionErrorBlockSize(){return predictionErrorBlockSize;}

    public double estimate(String coefficientName) {
        int index = coefficientNames.indexOf(coefficientName);
        if (index < 0) {
            throw new IllegalArgumentException(
                "unknown random-effect coefficient: " + coefficientName);
        }
        return estimates[index];
    }
}
