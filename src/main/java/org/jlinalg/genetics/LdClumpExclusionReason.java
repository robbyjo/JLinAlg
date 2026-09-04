/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.genetics;

/** Reason that a candidate was not retained as a clump index variant. */
public enum LdClumpExclusionReason {
    ABSENT_FROM_REFERENCE,
    ABOVE_INDEX_P_VALUE_THRESHOLD,
    IN_LINKAGE_DISEQUILIBRIUM
}
