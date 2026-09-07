/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.jlinalg.genetics.LdClumpCandidate;
import org.jlinalg.genetics.LdClumpResult;
import org.jlinalg.genetics.SecondarySignalClumper;
import org.junit.jupiter.api.Test;

class MrWorkflowExtensionsTest {
    private static List<HarmonizedInstrument> instruments() {
        return List.of(
            new HarmonizedInstrument("rs1", "A", "C", .2, .05, .1, .05, .3, .3, false, false),
            new HarmonizedInstrument("rs2", "A", "C", .3, .05, .2, .05, .3, .3, false, false),
            new HarmonizedInstrument("rs3", "A", "C", .4, .05, .1, .05, .3, .3, false, false));
    }

    @Test
    void nativeSvgPlotIsPortable() throws Exception {
        Path file = Files.createTempFile("jlinalg-mr-", ".svg");
        MrPlot.writeSvg(file, instruments(), MendelianRandomization.ivw(instruments(), false, .95));
        assertTrue(Files.readString(file).contains("<svg"));
        Files.deleteIfExists(file);
    }

    @Test
    void secondarySignalsRespectWithinGroupLd() {
        LdClumpResult result = SecondarySignalClumper.clump(
            List.of(new LdClumpCandidate("rs1", .001, "trait"),
                new LdClumpCandidate("rs2", .01, "trait"),
                new LdClumpCandidate("rs3", .02, "trait")),
            new double[][] {{1, .9, .1}, {.9, 1, .1}, {.1, .1, 1}}, .5, 3);
        assertTrue(result.retained().stream().map(LdClumpCandidate::variantId).toList()
            .equals(List.of("rs1", "rs3")));
    }
}
