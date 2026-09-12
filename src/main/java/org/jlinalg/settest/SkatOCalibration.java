/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

/** Calibration algorithms for the minimum-p SKAT-O omnibus statistic. */
public enum SkatOCalibration {
    /** GMMAT-compatible moment matching followed by one-dimensional quadrature. */
    ANALYTIC,
    /** Conditional noncentral Gaussian quadrature with explicit convergence checks. */
    DETERMINISTIC,
    /** Reproducible Monte Carlo calibration from the exact Gaussian score null. */
    PARAMETRIC_SIMULATION
}
