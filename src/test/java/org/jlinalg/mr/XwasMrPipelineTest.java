/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.inference.PValueScale;
import org.junit.jupiter.api.Test;

class XwasMrPipelineTest {
    @Test
    void scansNestedGridInStableOrderAndReusesClumpedExposures() {
        XwasMrExposure expression = new XwasMrExposure("GENE1", "Gene 1",
            exposure("a", new double[] {0.10, 0.15, 0.20}));
        XwasMrExposure protein = new XwasMrExposure("PROT1", "Protein 1",
            exposure("b", new double[] {0.08, 0.12, 0.17}));
        List<XwasMrOutcome> outcomes = List.of(
            new XwasMrOutcome("CAD", "Coronary artery disease",
                "cardiovascular", combinedOutcome(0.5, 0.0)),
            new XwasMrOutcome("CKD", "Chronic kidney disease",
                "kidney", combinedOutcome(0.0, -0.4)),
            new XwasMrOutcome("COPD", "Chronic obstructive pulmonary disease",
                "lung", List.of(outcome("a1", 0.01, 0.1),
                    outcome("b1", 0.01, 0.1))));
        XwasMrPipeline pipeline = XwasMrPipeline.prepare(
            List.of(expression, protein), outcomes);
        XwasMrSignificanceFilter filter =
            XwasMrSignificanceFilter.pValueAtMost(1e-4);
        MrOptions diagnostics = new MrOptions(0.95, 100, 42L);

        XwasMrBatchResult serial = pipeline.scan(new XwasMrOptions(1, 2,
            XwasMrScreeningMethod.IVW_MULTIPLICATIVE_RANDOM, filter,
            diagnostics));
        XwasMrBatchResult parallel = pipeline.scan(new XwasMrOptions(4, 3,
            XwasMrScreeningMethod.IVW_MULTIPLICATIVE_RANDOM, filter,
            diagnostics));

        assertEquals(6, serial.totalPairs());
        assertEquals(4, serial.analyzablePairs());
        assertEquals(2, serial.belowThresholdPairs());
        assertEquals(2, serial.insufficientInstrumentPairs());
        assertTrue(serial.failures().isEmpty());
        assertEquals(List.of("GENE1->CAD", "PROT1->CKD"), serial.hits().stream()
            .map(hit -> hit.exposureId() + "->" + hit.outcomeId()).toList());
        assertEquals(List.of("cardiovascular", "kidney"), serial.hits().stream()
            .map(XwasMrHit::outcomeCategory).toList());
        assertEquals(serial.hits(), parallel.hits());
        assertEquals(3, expression.clumpedInstruments().size());
        assertEquals(3, protein.clumpedInstruments().size());
    }

    @Test
    void makesPValueScaleDirectionsExplicit() {
        double p = 1e-5;
        XwasMrSignificanceFilter regular =
            XwasMrSignificanceFilter.pValueAtMost(p);
        XwasMrSignificanceFilter log =
            XwasMrSignificanceFilter.log10PAtMost(-5.0);
        XwasMrSignificanceFilter negativeLog =
            XwasMrSignificanceFilter.negativeLog10PAtLeast(5.0);

        assertTrue(regular.includes(p));
        assertTrue(log.includes(p));
        assertTrue(negativeLog.includes(p));
        assertFalse(regular.includes(1e-4));
        assertFalse(log.includes(1e-4));
        assertFalse(negativeLog.includes(1e-4));
        assertEquals(PValueScale.LOG10, log.scale());
        assertThrows(IllegalArgumentException.class,
            () -> XwasMrSignificanceFilter.log10PAtMost(5.0));
        assertThrows(IllegalArgumentException.class,
            () -> XwasMrSignificanceFilter.negativeLog10PAtLeast(-5.0));
    }

    private static List<SummaryAssociation> exposure(String prefix,
            double[] effects) {
        List<SummaryAssociation> result = new ArrayList<>();
        for (int index = 0; index < effects.length; index++)
            result.add(new SummaryAssociation(prefix + (index + 1), "A", "C",
                effects[index], 0.01, 0.2 + 0.05 * index));
        return result;
    }

    private static List<SummaryAssociation> combinedOutcome(
            double firstCausalEffect, double secondCausalEffect) {
        List<SummaryAssociation> result = new ArrayList<>();
        double[] first = {0.10, 0.15, 0.20};
        double[] second = {0.08, 0.12, 0.17};
        for (int index = 0; index < first.length; index++)
            result.add(outcome("a" + (index + 1),
                firstCausalEffect * first[index],
                firstCausalEffect == 0.0 ? 0.1 : 0.01));
        for (int index = 0; index < second.length; index++)
            result.add(outcome("b" + (index + 1),
                secondCausalEffect * second[index],
                secondCausalEffect == 0.0 ? 0.1 : 0.01));
        return result;
    }

    private static SummaryAssociation outcome(String variant, double effect,
            double standardError) {
        return new SummaryAssociation(variant, "A", "C", effect,
            standardError, 0.2);
    }
}
