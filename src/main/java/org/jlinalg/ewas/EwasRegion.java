/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.ewas;

import java.util.List;

/** One nonoverlapping coordinate-defined EWAS region and its calibrated test. */
public record EwasRegion(
        String genomeBuild,
        String chromosome,
        long start,
        long end,
        List<String> probes,
        double meanEffect,
        double statistic,
        double pValue,
        double adjustedPValue,
        double probesPerKilobase,
        long maximumObservedGap,
        double correlationLength) {
    public EwasRegion { probes = List.copyOf(probes); }
}
