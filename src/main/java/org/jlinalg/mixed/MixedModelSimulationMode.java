/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mixed;

/** Whether simulated responses draw new random effects or retain fitted modes. */
public enum MixedModelSimulationMode {
    /** Draw new random effects and residual errors from fitted variances. */
    MARGINAL,
    /** Retain fitted conditional modes and draw only residual errors. */
    CONDITIONAL
}
