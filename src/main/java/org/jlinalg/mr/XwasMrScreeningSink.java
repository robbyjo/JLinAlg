/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

/** Receives successful xWAS screens in exposure-major, outcome-major order. */
@FunctionalInterface
public interface XwasMrScreeningSink {
    /** Accepts one pair after IVW screening and before optional diagnostics. */
    void accept(XwasMrScreeningResult result);
}
