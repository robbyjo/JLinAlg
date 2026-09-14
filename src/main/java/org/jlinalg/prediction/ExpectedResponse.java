/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.prediction;

/** One design-row expected response with link- and response-scale uncertainty. */
public record ExpectedResponse(
        double linearPredictor,
        double linkStandardError,
        MeanEstimate mean) { }
