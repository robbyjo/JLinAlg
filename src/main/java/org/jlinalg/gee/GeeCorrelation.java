/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

/** Working within-cluster association structure. */
public enum GeeCorrelation {
    INDEPENDENCE,
    EXCHANGEABLE,
    AR1,
    M_DEPENDENT,
    TOEPLITZ,
    UNSTRUCTURED,
    FIXED,
    USER_DEFINED
}
