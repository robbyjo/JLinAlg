/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.inference;

/** Percentile interval and empirical moments for one fitted parameter. */
public record BootstrapParameterSummary(
        String name,
        double estimate,
        double bootstrapMean,
        double bias,
        double standardError,
        double lower,
        double upper,
        int successfulReplicates) { }
