/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.prediction;

/** Positive response-mean ratio with delta-method, log-scale inference. */
public record RiskRatio(
        Averaging averaging,
        MeanEstimate firstScenario,
        MeanEstimate secondScenario,
        double covarianceBetweenScenarios,
        double ratio,
        double standardError,
        double logStandardError,
        double confidenceLower,
        double confidenceUpper) { }
