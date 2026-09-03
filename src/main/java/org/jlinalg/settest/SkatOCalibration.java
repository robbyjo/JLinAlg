/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

/** Calibration algorithms for the minimum-p SKAT-O omnibus statistic. */
public enum SkatOCalibration {
    /** GMMAT-compatible moment matching followed by one-dimensional quadrature. */
    ANALYTIC,
    /** Reproducible correlated parametric-null simulation. */
    PARAMETRIC_SIMULATION
}
