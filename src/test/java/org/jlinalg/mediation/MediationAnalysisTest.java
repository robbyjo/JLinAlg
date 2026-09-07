/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mediation;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.model.MissingDataPolicy;
import org.jlinalg.ols.Ols;
import org.jlinalg.ols.OlsOptions;
import org.jlinalg.ols.RankDeficiencyStrategy;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.pedigree.PedigreeIndividual;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

class MediationAnalysisTest {
    private static final double[] TREATMENT = {
        -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5
    };
    private static final double[] MEDIATOR = {
        -2.1, -1.1, -0.2, 0.7, 1.1, 2.0,
        2.9, 3.4, 4.7, 5.3, 6.1, 7.2
    };
    private static final double[] OUTCOME = {
        -2.4, -0.4, 0.7, 2.0, 2.8, 4.3,
        5.0, 6.8, 8.1, 9.0, 10.7, 12.0
    };

    @Test
    void estimatesPathsAndIndirectEffectWithoutSampling() {
        MediationResult result = MediationAnalysis.fit(
            OUTCOME, TREATMENT, MEDIATOR);

        OlsResultReference expected = independentFits();
        assertEquals(expected.mediator().coefficients()[1],
            result.aPath().estimate(), 1e-12);
        assertEquals(expected.outcome().coefficients()[2],
            result.bPath().estimate(), 1e-12);
        assertEquals(expected.outcome().coefficients()[1],
            result.directEffect().estimate(), 1e-12);
        assertEquals(expected.total().coefficients()[1],
            result.totalEffect().estimate(), 1e-12);

        double a = result.aPath().estimate();
        double b = result.bPath().estimate();
        double expectedIndirect = a * b;
        double expectedSe = Math.sqrt(
            b * b * expected.mediator().covariance()[1 * 2 + 1]
            + a * a * expected.outcome().covariance()[2 * 3 + 2]);
        assertEquals(expectedIndirect, result.indirectEffect().estimate(), 1e-12);
        assertEquals(expectedSe, result.indirectEffect().standardError(), 1e-12);
        assertEquals(Double.POSITIVE_INFINITY,
            result.indirectEffect().degreesOfFreedom());
        assertEquals(TREATMENT.length, result.observations());
        assertArrayEquals(new int[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
            result.retainedRows());
    }

    @Test
    void adjustsBothPathsForCovariatesAndUsesCommonOmitRows() {
        double[] outcome = OUTCOME.clone();
        double[] treatment = TREATMENT.clone();
        double[] mediator = MEDIATOR.clone();
        double[][] covariates = new double[TREATMENT.length][1];
        for (int row = 0; row < covariates.length; row++) {
            covariates[row][0] = row % 2;
        }
        outcome[2] = Double.NaN;
        covariates[8][0] = Double.NaN;
        OlsOptions options = new OlsOptions(
            RankDeficiencyStrategy.ERROR, 0.90, MissingDataPolicy.OMIT);

        MediationResult result = MediationAnalysis.fit(
            outcome, treatment, mediator, covariates, options,
            BackendPolicy.CPU);

        assertEquals(10, result.observations());
        assertEquals(12, result.originalObservations());
        assertArrayEquals(new int[] {0, 1, 3, 4, 5, 6, 7, 9, 10, 11},
            result.retainedRows());
        assertEquals(3, result.mediatorModel().parameters());
        assertEquals(4, result.outcomeModel().parameters());
        assertEquals(3, result.totalModel().parameters());
    }

    @Test
    void rejectsRankDeficientMediationPaths() {
        double[] constantTreatment = new double[TREATMENT.length];
        OlsOptions options = new OlsOptions(
            RankDeficiencyStrategy.MINIMUM_NORM, 0.95,
            MissingDataPolicy.ERROR);

        assertThrows(IllegalArgumentException.class, () ->
            MediationAnalysis.fit(
                OUTCOME, constantTreatment, MEDIATOR, null, options,
                BackendPolicy.CPU));
    }

    @Test
    void reusesSparseRemlStructureForMixedMediation() {
        List<String> groups = List.of(
            "a", "a", "a", "b", "b", "b", "c", "c", "c",
            "d", "d", "d");
        RandomEffectTerm randomIntercept = RandomEffectTerm.randomIntercept(
            "subject", groups);
        RemlOptions options = RemlOptions.builder()
            .initialVariances(2.0, 1.0)
            .maximumIterations(160)
            .build();

        MediationMixedResult result = MediationAnalysis.fitMixed(
            OUTCOME, TREATMENT, MEDIATOR, null,
            List.of(randomIntercept), options, BackendPolicy.CPU);

        assertEquals(12, result.observations());
        assertTrue(result.converged());
        assertEquals(result.aPath().estimate(),
            result.mediatorModel().beta()[1], 0.0);
        assertEquals(result.bPath().estimate(),
            result.outcomeModel().beta()[2], 0.0);
        assertEquals("subject",
            result.outcomeModel().randomEffects().get(0).termName());
        assertEquals(2, result.outcomeModel().varianceComponents().length);
    }

    @Test
    void includesPedigreePrecisionInAllMediationComponentFits() {
        Pedigree pedigree = Pedigree.of(List.of(
            PedigreeIndividual.founder("A"),
            PedigreeIndividual.founder("B"),
            new PedigreeIndividual("C", "A", "B"),
            new PedigreeIndividual("D", "A", "C")));
        List<String> observed = List.of(
            "A", "A", "A", "B", "B", "B",
            "C", "C", "C", "D", "D", "D");
        PedigreeRandomEffectTerm additive = PedigreeRandomEffectTerm.of(
            "animal", observed, pedigree);
        RemlOptions options = RemlOptions.builder()
            .initialVariances(2.0, 1.0)
            .maximumIterations(160)
            .build();

        MediationMixedResult result = MediationAnalysis.fitPedigree(
            OUTCOME, TREATMENT, MEDIATOR, null,
            List.of(additive), List.of(), options, BackendPolicy.CPU);

        assertTrue(result.converged());
        assertEquals(pedigree.size(), result.outcomeModel()
            .randomEffects("animal").estimates().length);
        assertEquals("animal", result.totalModel().componentNames().get(0));
        assertTrue(result.indirectEffect().standardError() >= 0.0);
    }

    private static OlsResultReference independentFits() {
        double[][] mediatorDesign = new double[TREATMENT.length][2];
        double[][] outcomeDesign = new double[TREATMENT.length][3];
        double[][] totalDesign = new double[TREATMENT.length][2];
        for (int row = 0; row < TREATMENT.length; row++) {
            mediatorDesign[row] = new double[] {1, TREATMENT[row]};
            outcomeDesign[row] = new double[] {
                1, TREATMENT[row], MEDIATOR[row]
            };
            totalDesign[row] = new double[] {1, TREATMENT[row]};
        }
        return new OlsResultReference(
            Ols.fit(MEDIATOR, mediatorDesign),
            Ols.fit(OUTCOME, outcomeDesign),
            Ols.fit(OUTCOME, totalDesign));
    }

    private record OlsResultReference(
            org.jlinalg.ols.OlsResult mediator,
            org.jlinalg.ols.OlsResult outcome,
            org.jlinalg.ols.OlsResult total) {
    }
}
