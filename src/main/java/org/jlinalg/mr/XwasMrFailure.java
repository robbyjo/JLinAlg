/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

/** A pair that was analyzable at screening but failed full diagnostics. */
public record XwasMrFailure(
        int exposureIndex,
        int outcomeIndex,
        String exposureId,
        String outcomeId,
        String exceptionType,
        String message) {
}
