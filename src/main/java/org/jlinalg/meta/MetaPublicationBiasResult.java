/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

/** Small-study and rank-correlation publication-bias diagnostics. */
public record MetaPublicationBiasResult(
        double eggerIntercept, double eggerStandardError,
        double eggerStatistic, double eggerPValue,
        double rankCorrelation, double rankStatistic, double rankPValue) { }
