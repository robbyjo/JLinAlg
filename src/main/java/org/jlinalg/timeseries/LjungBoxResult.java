/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.timeseries;

/** Ljung-Box residual-portmanteau statistic and chi-square approximation. */
public record LjungBoxResult(
        double statistic,
        int lags,
        int degreesOfFreedom,
        double pValue) { }
