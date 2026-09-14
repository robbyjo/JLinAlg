/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.prediction;

/** Response-scale difference between two covariate scenarios. */
public record ScenarioDifference(
        Averaging averaging,
        MeanEstimate firstScenario,
        MeanEstimate secondScenario,
        double covarianceBetweenScenarios,
        double difference,
        double standardError,
        double confidenceLower,
        double confidenceUpper) { }
