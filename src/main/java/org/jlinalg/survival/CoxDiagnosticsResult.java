/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

/** Robust covariance, influence residuals, and PH diagnostics. */
public record CoxDiagnosticsResult(
        double[] robustCovariance,
        double[] robustStandardErrors,
        int clusters,
        CoxResiduals residuals,
        CoxProportionalHazardsTest proportionalHazards) {
    public CoxDiagnosticsResult {
        robustCovariance = robustCovariance.clone();
        robustStandardErrors = robustStandardErrors.clone();
    }

    @Override public double[] robustCovariance() {
        return robustCovariance.clone();
    }
    @Override public double[] robustStandardErrors() {
        return robustStandardErrors.clone();
    }
}
