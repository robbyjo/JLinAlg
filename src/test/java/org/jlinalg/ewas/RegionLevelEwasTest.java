/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.ewas;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import java.util.List;
import jdistlib.Normal;
import org.junit.jupiter.api.Test;

final class RegionLevelEwasTest {
    @Test void signedStoufferUsesDeclaredSpatialCovarianceAndCoverage() {
        List<EwasProbe> probes = List.of(
            new EwasProbe("a", "1", 100, 0.2, 0.01),
            new EwasProbe("b", "chr1", 150, 0.3, 0.02),
            new EwasProbe("c", "1", 220, 0.1, 0.05),
            new EwasProbe("d", "2", 1000, -0.2, 0.10),
            new EwasProbe("e", "2", 1060, -0.3, 0.03),
            new EwasProbe("f", "2", 1120, -0.1, 0.04));
        List<EwasRegion> result = RegionLevelEwas.scan(probes, "GRCh38",
            new RegionLevelEwas.Options(100, 3, 200));
        assertEquals(2, result.size());
        EwasRegion first = result.get(0);
        double sumZ = Normal.quantile(0.995, 0, 1, true, false)
            + Normal.quantile(0.99, 0, 1, true, false)
            + Normal.quantile(0.975, 0, 1, true, false);
        long[] positions = {100, 150, 220};
        double variance = 0.0;
        for (long left : positions) for (long right : positions)
            variance += Math.exp(-Math.abs(left - right) / 200.0);
        assertEquals(sumZ / Math.sqrt(variance), first.statistic(), 1e-12);
        assertEquals(3 * 1000.0 / 121.0, first.probesPerKilobase(), 1e-12);
        assertEquals(70, first.maximumObservedGap());
        assertTrue(first.adjustedPValue() >= first.pValue());
    }

    @Test void duplicateCoordinatesAreRejected() {
        List<EwasProbe> probes = List.of(
            new EwasProbe("a", "1", 100, 1, 0.1),
            new EwasProbe("b", "chr1", 100, 1, 0.1));
        assertThrows(IllegalArgumentException.class, () ->
            RegionLevelEwas.scan(probes, "GRCh38",
                RegionLevelEwas.Options.defaults()));
    }
}
