/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

/** Score vector and its row-major covariance/information matrix. */
public record SetTestScoreState(
        double[] scores, double[] information, int variants) {
    public SetTestScoreState {
        if (variants < 1 || scores == null || scores.length != variants
                || information == null
                || information.length != variants * variants)
            throw new IllegalArgumentException("set-test score dimensions are invalid");
        scores = scores.clone();
        information = information.clone();
    }
    @Override public double[] scores() { return scores.clone(); }
    @Override public double[] information() { return information.clone(); }
    double[] scoresView() { return scores; }
    double[] informationView() { return information; }
}
