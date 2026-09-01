/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

/** One stratum-specific baseline cumulative-hazard step. */
public record BaselineHazardPoint(
        int stratum,
        double time,
        int events,
        double hazardIncrement,
        double cumulativeHazard,
        double baselineSurvival) { }
