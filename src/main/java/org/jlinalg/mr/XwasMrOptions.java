/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

import java.util.Objects;

/** Parallelism, blocking, screening, and diagnostic controls for xWAS MR. */
public record XwasMrOptions(
        int parallelism,
        int pairBlockSize,
        XwasMrScreeningMethod screeningMethod,
        XwasMrSignificanceFilter significanceFilter,
        MrOptions diagnosticOptions) {
    /** Validates bounded execution controls. */
    public XwasMrOptions {
        if (parallelism < 1)
            throw new IllegalArgumentException("parallelism must be positive");
        if (pairBlockSize < 1)
            throw new IllegalArgumentException("pairBlockSize must be positive");
        Objects.requireNonNull(screeningMethod, "screeningMethod");
        Objects.requireNonNull(significanceFilter, "significanceFilter");
        Objects.requireNonNull(diagnosticOptions, "diagnosticOptions");
    }

    /** Uses available processors and bounded blocks for the requested filter. */
    public static XwasMrOptions defaults(
            XwasMrSignificanceFilter significanceFilter) {
        int workers = Math.max(1, Runtime.getRuntime().availableProcessors());
        return new XwasMrOptions(workers, Math.max(32, workers * 8),
            XwasMrScreeningMethod.IVW_MULTIPLICATIVE_RANDOM,
            significanceFilter, MrOptions.defaults());
    }
}
