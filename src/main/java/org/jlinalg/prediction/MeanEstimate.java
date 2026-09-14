/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.prediction;

/**
 * Expected-response estimate and uncertainty in that fitted mean.
 * This is not an interval for a future outcome.
 */
public record MeanEstimate(
        double estimate,
        double standardError,
        double confidenceLower,
        double confidenceUpper) { }
