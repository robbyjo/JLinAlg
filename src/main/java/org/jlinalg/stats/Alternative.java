/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.stats;

/** Alternative hypothesis used by scalar statistical tests. */
public enum Alternative {
    /** The parameter differs from its null value. */
    TWO_SIDED,
    /** The parameter is smaller than its null value. */
    LESS,
    /** The parameter is greater than its null value. */
    GREATER
}
