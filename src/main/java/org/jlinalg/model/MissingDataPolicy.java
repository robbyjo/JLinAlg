/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.model;

/** Policy applied when model inputs contain NaN or infinite values. */
public enum MissingDataPolicy {
    /** Reject the fit immediately. This is the fastest policy and the default. */
    ERROR,

    /** Fit using rows for which response, design, weights, and offset are finite. */
    OMIT
}
