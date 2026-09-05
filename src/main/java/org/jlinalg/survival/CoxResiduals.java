/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

/** Observation- and event-level residuals from a fitted Cox model. */
public record CoxResiduals(
        double[] martingale,
        double[] deviance,
        double[] score,
        double[] dfbeta,
        int observations,
        int coefficients,
        int[] schoenfeldRows,
        double[] schoenfeldTimes,
        double[] schoenfeld) {
    public CoxResiduals {
        martingale = martingale.clone();
        deviance = deviance.clone();
        score = score.clone();
        dfbeta = dfbeta.clone();
        schoenfeldRows = schoenfeldRows.clone();
        schoenfeldTimes = schoenfeldTimes.clone();
        schoenfeld = schoenfeld.clone();
    }

    @Override public double[] martingale() { return martingale.clone(); }
    @Override public double[] deviance() { return deviance.clone(); }
    @Override public double[] score() { return score.clone(); }
    @Override public double[] dfbeta() { return dfbeta.clone(); }
    @Override public int[] schoenfeldRows() { return schoenfeldRows.clone(); }
    @Override public double[] schoenfeldTimes() {
        return schoenfeldTimes.clone();
    }
    @Override public double[] schoenfeld() { return schoenfeld.clone(); }
}
