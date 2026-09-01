/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.formula;

/** Coding used when expanding categorical fixed effects. */
public enum ContrastCoding {
    /** R-style treatment coding with one reference level. */
    TREATMENT,
    /** Sum-to-zero coding with the last observed level omitted. */
    SUM
}
