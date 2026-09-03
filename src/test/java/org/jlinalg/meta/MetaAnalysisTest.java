/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class MetaAnalysisTest {
    private static final List<MetaStudy> STUDIES = List.of(
        new MetaStudy("A", 0.2, 0.1),
        new MetaStudy("B", 0.5, 0.2),
        new MetaStudy("C", 0.1, 0.15),
        new MetaStudy("D", 0.7, 0.25));

    @Test
    void fixedEffectMatchesBaseRInverseVarianceReference() {
        MetaAnalysisResult result = MetaAnalysis.fit(STUDIES,
            MetaAnalysisOptions.fixedEffect(), BackendPolicy.CPU);

        assertEquals(0.25961653684841224, result.pooledEffectSize(), 1e-12);
        assertEquals(0.0734333070, result.standardError(), 1e-9);
        assertEquals(6.035350509286997, result.cochranQ(), 1e-11);
        assertEquals(0.10990295531017164, result.cochranQPValue(), 1e-11);
        assertEquals(50.29286210662165, result.iSquared(), 1e-10);
        assertEquals(2.0117835030956654, result.hSquared(), 1e-11);
        assertEquals(0.0, result.tauSquared());
        assertTrue(Double.isNaN(result.predictionLower()));
        assertEquals(1.0, java.util.Arrays.stream(
            result.normalizedWeights()).sum(), 1e-14);
    }

    @Test
    void randomEffectsDlMatchesBaseRReferenceAndExposesLogP() {
        MetaAnalysisOptions options = MetaAnalysisOptions.builder()
            .tauSquaredEstimator(TauSquaredEstimator.DERSIMONIAN_LAIRD)
            .build();
        MetaAnalysisResult result = MetaAnalysis.fit(
            STUDIES, options, BackendPolicy.CPU);

        assertEquals(0.02614035088, result.tauSquared(), 1e-9);
        assertEquals(0.3088126569134951, result.pooledEffectSize(), 1e-12);
        assertEquals(0.11575603720170137, result.standardError(), 1e-12);
        assertTrue(Double.isFinite(result.predictionLower()));
        assertEquals(Math.log10(result.pValue()), result.log10PValue(), 1e-12);
        assertEquals(-result.log10PValue(), result.negativeLog10PValue(), 1e-12);
    }

    @Test
    void metaRegressionMatchesBaseRGeneralizedDlReference() {
        double[][] moderators = {{-1}, {0}, {1}, {2}};
        MetaAnalysisOptions options = MetaAnalysisOptions.builder()
            .tauSquaredEstimator(TauSquaredEstimator.DERSIMONIAN_LAIRD)
            .build();
        MetaRegressionResult result = MetaRegression.fit(STUDIES, moderators,
            List.of("dose"), true, options, BackendPolicy.CPU);

        assertEquals(0.05264313, result.tauSquared(), 1e-7);
        assertArrayEquals(new double[] {0.3004978, 0.09178625},
            result.beta(), 1e-7);
        assertArrayEquals(new double[] {0.1484654, 0.1295650},
            result.standardErrors(), 1e-7);
        assertEquals(5.21784776902887, result.residualQ(), 1e-11);
        assertEquals(1.0, result.moderatorQDegreesOfFreedom());
        assertTrue(result.moderatorQPValue() >= 0.0
            && result.moderatorQPValue() <= 1.0);
        assertEquals(List.of("(Intercept)", "dose"),
            result.coefficientNames());
    }

    @Test
    void remlAndModifiedKnappHartungRemainFinite() {
        MetaAnalysisOptions options = MetaAnalysisOptions.builder()
            .inferenceMethod(MetaInferenceMethod.MODIFIED_HARTUNG_KNAPP)
            .build();
        MetaAnalysisResult result = MetaAnalysis.fit(
            STUDIES, options, BackendPolicy.CPU);
        assertTrue(result.tauSquared() >= 0.0);
        assertTrue(Double.isFinite(result.standardError()));
        assertTrue(result.pValue() >= 0.0 && result.pValue() <= 1.0);
    }

    @Test
    void preparedBatchMatchesScalarFixedDlAndRemlFits() {
        double[] effects = STUDIES.stream()
            .mapToDouble(MetaStudy::effectSize).toArray();
        double[] errors = STUDIES.stream()
            .mapToDouble(MetaStudy::standardError).toArray();
        PreparedMetaAnalysisBatch prepared = MetaAnalysis.prepareBatch(
            effects, errors, 1, STUDIES.size());

        for (MetaAnalysisOptions options : List.of(
                MetaAnalysisOptions.fixedEffect(),
                MetaAnalysisOptions.builder()
                    .tauSquaredEstimator(TauSquaredEstimator.DERSIMONIAN_LAIRD)
                    .build(),
                MetaAnalysisOptions.builder()
                    .tauSquaredEstimator(TauSquaredEstimator.PAULE_MANDEL)
                    .build(),
                MetaAnalysisOptions.randomEffects())) {
            MetaAnalysisResult scalar = MetaAnalysis.fit(
                STUDIES, options, BackendPolicy.CPU);
            MetaAnalysisBatchResult batch = prepared.fit(options);
            assertEquals(scalar.pooledEffectSize(),
                batch.pooledEffectSizes()[0], 2e-8);
            assertEquals(scalar.standardError(),
                batch.standardErrors()[0], 2e-8);
            assertEquals(scalar.cochranQ(), batch.cochranQ()[0], 1e-12);
            assertEquals(scalar.tauSquared(), batch.tauSquared()[0], 2e-8);
            assertEquals(scalar.pValue(), batch.pValues()[0], 2e-8);
        }
    }

    @Test
    void preparedBatchParallelChunksAreDeterministic() {
        int analyses = 4097;
        double[] effects = new double[analyses * STUDIES.size()];
        double[] errors = new double[effects.length];
        for (int analysis = 0; analysis < analyses; analysis++) {
            for (int study = 0; study < STUDIES.size(); study++) {
                int index = analysis * STUDIES.size() + study;
                effects[index] = STUDIES.get(study).effectSize()
                    + analysis * 1e-7;
                errors[index] = STUDIES.get(study).standardError();
            }
        }
        PreparedMetaAnalysisBatch prepared = MetaAnalysis.prepareBatch(
            effects, errors, analyses, STUDIES.size());
        MetaAnalysisOptions options = MetaAnalysisOptions.builder()
            .tauSquaredEstimator(TauSquaredEstimator.DERSIMONIAN_LAIRD)
            .build();
        MetaAnalysisBatchResult sequential = prepared.fit(options, 1);
        MetaAnalysisBatchResult parallel = prepared.fit(options, 4);
        assertArrayEquals(sequential.pooledEffectSizes(),
            parallel.pooledEffectSizes());
        assertArrayEquals(sequential.pValues(), parallel.pValues());
        assertArrayEquals(sequential.tauSquared(), parallel.tauSquared());
    }

    @Test
    void preparedBatchRejectsInvalidShapesAndValues() {
        assertThrows(IllegalArgumentException.class,
            () -> MetaAnalysis.prepareBatch(
                new double[3], new double[3], 1, 2));
        assertThrows(IllegalArgumentException.class,
            () -> MetaAnalysis.prepareBatch(
                new double[] {0.1, 0.2}, new double[] {0.1, 0.0}, 1, 2));
    }
}
