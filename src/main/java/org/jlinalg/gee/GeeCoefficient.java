/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

/** One tidy coefficient-inference row. */
public record GeeCoefficient(
        int index,
        double estimate,
        double standardError,
        double statistic,
        double pValue,
        double confidenceLower,
        double confidenceUpper,
        double degreesOfFreedom) { }
