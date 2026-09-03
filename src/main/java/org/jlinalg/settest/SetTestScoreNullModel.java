/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import jdistlib.Normal;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendPolicy;

/** Reusable efficient-score projection for a variant-set null model. */
public interface SetTestScoreNullModel {
    int observations();
    SetTestScoreState score(double[][] variantRows);

    /** Denominator degrees of freedom, or NaN for an asymptotic normal score. */
    default double degreesOfFreedom() { return Double.NaN; }

    /** Two-sided burden-score tail probability for this null model. */
    default double burdenPValue(double statistic) {
        return Math.min(1.0, 2.0 * Normal.cumulative(
            Math.abs(statistic), 0.0, 1.0, false, false));
    }

    default String burdenPValueMethod() { return "normal-score"; }

    /**
     * Returns a retained backend owned by this null model, or {@code null}
     * when set-test operations should acquire one from {@link #backendPolicy()}.
     * Callers must not close the returned backend.
     */
    default ComputeBackend computeBackend() { return null; }

    /** Backend policy used when this null model does not retain a backend. */
    default BackendPolicy backendPolicy() { return BackendPolicy.PREFERRED; }
}
