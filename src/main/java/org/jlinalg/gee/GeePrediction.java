/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

/** One new-data marginal prediction with link and response-scale uncertainty. */
public record GeePrediction(
        double linearPredictor,
        double linkStandardError,
        double mean,
        double responseStandardError,
        double confidenceLower,
        double confidenceUpper) { }
