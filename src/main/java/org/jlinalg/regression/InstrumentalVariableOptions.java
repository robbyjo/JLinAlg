/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.regression;

import java.util.Objects;

/** Immutable options for individual-level instrumental-variable regression. */
public record InstrumentalVariableOptions(
        InstrumentalVariableCovariance covariance,
        double confidenceLevel) {

    /** Creates validated options. */
    public InstrumentalVariableOptions {
        Objects.requireNonNull(covariance, "covariance");
        if (!(confidenceLevel > 0.0 && confidenceLevel < 1.0)) {
            throw new IllegalArgumentException(
                "confidenceLevel must be strictly between zero and one");
        }
    }

    /** Conventional homoskedastic 2SLS with 95 percent intervals. */
    public static InstrumentalVariableOptions defaults() {
        return new InstrumentalVariableOptions(
            InstrumentalVariableCovariance.HOMOSKEDASTIC, 0.95);
    }
}
