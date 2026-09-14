/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.prediction;

/** Response-scale derivative per unit of one linear design column. */
public record AverageMarginalEffect(
        Averaging averaging,
        int predictorColumn,
        double estimate,
        double standardError,
        double confidenceLower,
        double confidenceUpper) { }
