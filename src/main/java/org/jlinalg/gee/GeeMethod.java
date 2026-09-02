/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

/** Regression estimating-equation adjustment. */
public enum GeeMethod {
    ORDINARY,
    BIAS_REDUCED,
    BIAS_CORRECTED,
    JEFFREYS,
    ONE_STEP_JEFFREYS,
    HYBRID_JEFFREYS
}
