/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

/** Residual representation returned by a GEE fit. */
public enum GeeResidualType {
    RESPONSE,
    PEARSON,
    DEVIANCE,
    WORKING,
    STANDARDIZED
}
