/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import jdistlib.ChiSquare;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.CpuComputeBackend;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gwas.RemlAssociationScanner;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glmm.GlmmPqlOptions;
import org.jlinalg.ols.OlsOptions;
import org.jlinalg.pipeline.VariantFilterOptions;
import org.jlinalg.pipeline.VariantRecord;
import org.jlinalg.reml.RemlOptions;
import org.jlinalg.reml.VarianceComponent;
import org.junit.jupiter.api.Test;

class SetTestsTest {
    private static final double[][] INTERCEPT = {
        {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}
    };
    private static final double[] FIRST =
        {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2};
    private static final double[] SECOND =
        {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
    private static final double[] RESPONSE =
        {1.0, 1.3, 1.1, 1.4, 2.5, 2.9, 2.6, 3.0, 4.1, 4.6, 4.2, 4.7};

    @Test
    void burdenReturnsEffectUncertaintyAndAlleleOrientation() {
        LinearSetTestNullModel nullModel = nullModel();
        VariantSet alternate = new VariantSet("gene",
            List.of(member("v1", FIRST, EffectAllele.ALTERNATE)));
        VariantSet reference = new VariantSet("gene-ref",
            List.of(member("v1", FIRST, EffectAllele.REFERENCE)));

        SetTestResult alt = SetTests.burden(
            alternate, nullModel, options(100));
        SetTestResult ref = SetTests.burden(
            reference, nullModel, options(100));

        assertTrue(alt.beta() > 1);
        assertTrue(alt.standardError() > 0);
        assertTrue(alt.pValue() < 1e-5);
        assertEquals(-alt.beta(), ref.beta(), 1e-12);
        assertEquals(alt.pValue(), ref.pValue(), 1e-12);
    }

    @Test
    void skatAndSkatOAreFiniteAndSeedReproducible() {
        VariantSet set = new VariantSet("gene", List.of(
            member("v1", FIRST, EffectAllele.ALTERNATE),
            member("v2", SECOND, EffectAllele.ALTERNATE)));
        SetTestOptions options = options(250);

        SetTestResult skat = SetTests.skat(set, nullModel(), options);
        SkatOResult first = SetTests.skatO(set, nullModel(), options);
        SkatOResult second = SetTests.skatO(set, nullModel(), options);

        assertTrue(skat.statistic() > 0);
        assertTrue(skat.pValue() > 0 && skat.pValue() <= 1);
        assertTrue(skat.eigenvalues().length > 0);
        assertEquals(3, first.components().size());
        assertEquals(first.adjustedPValue(), second.adjustedPValue(), 0);
        assertEquals(250, first.simulations());
    }

    @Test
    void quadraticFormUsesExactScaledChiSquareForOneEigenvalue() {
        QuadraticFormDistribution.Tail result =
            QuadraticFormDistribution.survival(3.5, new double[] {2});
        assertEquals(ChiSquare.cumulative(3.5 / 2, 1, false, false),
            result.pValue(), 1e-15);
        assertEquals("exact-scaled-chi-square", result.method());
    }

    @Test
    void analyticSkatOMatchesGmmatReferenceAndUsesRetainedBackend() {
        double[] third =
            {0, 0, 1, 0, 1, 0, 2, 1, 0, 1, 2, 0};
        VariantSet set = new VariantSet("analytic", List.of(
            member("v1", FIRST, EffectAllele.ALTERNATE),
            member("v2", SECOND, EffectAllele.ALTERNATE),
            member("v3", third, EffectAllele.ALTERNATE)));
        SetTestOptions options = new SetTestOptions(
            VariantFilterOptions.defaults(), SetTestMissingPolicy.MEAN_IMPUTE,
            new double[] {0, 0.25, 0.5, 0.75, 1}, 0, 0,
            SkatOCalibration.ANALYTIC);
        FixedScoreNullModel nullModel = new FixedScoreNullModel(
            new double[] {0.4, -0.3, 0.8},
            new double[] {
                1, 0.2, 0.1,
                0.2, 0.8, 0.05,
                0.1, 0.05, 1.2
            });

        SkatOResult result = SetTests.skatO(set, nullModel, options);

        assertEquals(0.7856467311977459,
            result.adjustedPValue(), 1e-4);
        assertEquals(0, result.simulations());
    }

    @Test
    void setMembershipFilteringIsAudited() {
        VariantSet set = new VariantSet("filtered", List.of(
            member("v1", FIRST, EffectAllele.ALTERNATE),
            member("singleton",
                new double[] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                EffectAllele.ALTERNATE)));
        SetTestOptions options = new SetTestOptions(
            VariantFilterOptions.builder().minimumMac(2).build(),
            SetTestMissingPolicy.MEAN_IMPUTE,
            new double[] {0, 1}, 20, 7);

        SetTestResult result = SetTests.skat(set, nullModel(), options);

        assertEquals(1, result.includedVariants());
        assertEquals(1, result.excludedVariants().size());
        assertEquals("singleton",
            result.excludedVariants().get(0).variant().id());
    }

    @Test
    void relatedSampleTestsReuseOneRemlProjection() {
        double[] kinship = new double[RESPONSE.length * RESPONSE.length];
        for (int row = 0; row < RESPONSE.length; row++)
            for (int column = 0; column < RESPONSE.length; column++)
                kinship[row * RESPONSE.length + column] =
                    Math.pow(0.2, Math.abs(row - column));
        RemlAssociationScanner scanner = RemlAssociationScanner.prepare(
            RESPONSE, INTERCEPT,
            List.of(
                new VarianceComponent(
                    "cryptic-relatedness", RESPONSE.length, kinship),
                VarianceComponent.identity("residual", RESPONSE.length)),
            RemlOptions.builder().initialVariances(0.2, 1).build(),
            BackendPolicy.CPU);
        RemlSetTestNullModel nullModel = new RemlSetTestNullModel(scanner);
        VariantSet set = new VariantSet("related-gene", List.of(
            member("v1", FIRST, EffectAllele.ALTERNATE),
            member("v2", SECOND, EffectAllele.ALTERNATE)));

        SetTestResult burden = SetTests.burden(set, nullModel, options(100));
        SetTestResult skat = SetTests.skat(set, nullModel, options(100));
        SkatOResult skatO = SetTests.skatO(set, nullModel, options(100));

        assertTrue(Double.isFinite(burden.beta()));
        assertTrue(burden.standardError() > 0);
        assertTrue(burden.pValue() > 0 && burden.pValue() <= 1);
        assertTrue(skat.statistic() > 0);
        assertTrue(skat.pValue() > 0 && skat.pValue() <= 1);
        assertTrue(skatO.adjustedPValue() > 0
            && skatO.adjustedPValue() <= 1);
    }

    @Test
    void preparedSuiteProjectsOnceAndMatchesStandaloneKernelTests() {
        VariantSet set = new VariantSet("suite", List.of(
            member("v1", FIRST, EffectAllele.ALTERNATE),
            member("v2", SECOND, EffectAllele.ALTERNATE)));
        SetTestOptions options = options(250);
        PreparedVariantSet prepared = SetTests.prepare(
            set, RESPONSE.length, options);
        CountingNullModel counted = new CountingNullModel(nullModel());

        SetTestSuiteResult suite = SetTests.scoreSuite(
            prepared, counted, options);
        assertEquals(1, counted.calls);

        SetTestResult skat = SetTests.skat(prepared, nullModel());
        SkatOResult skatO = SetTests.skatO(prepared, nullModel(), options);
        assertEquals(skat.statistic(), suite.skat().statistic(), 1e-12);
        assertEquals(skat.pValue(), suite.skat().pValue(), 1e-12);
        assertEquals(skatO.adjustedPValue(),
            suite.skatO().adjustedPValue(), 0);
        for (int index = 0; index < skatO.components().size(); index++)
            assertEquals(skatO.components().get(index).result().pValue(),
                suite.skatO().components().get(index).result().pValue(), 1e-10);
    }

    @Test
    void binaryPqlNullSupportsPreparedSetTests() {
        double[] binary = {0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1};
        double[] kinship = new double[binary.length * binary.length];
        for (int row = 0; row < binary.length; row++)
            for (int column = 0; column < binary.length; column++)
                kinship[row * binary.length + column] =
                    Math.pow(0.1, Math.abs(row - column));
        PqlSetTestNullModel nullModel = PqlSetTestNullModel.prepare(
            binary, INTERCEPT, GlmFamilies.binomial(),
            List.of(new VarianceComponent("kinship", binary.length, kinship)),
            GlmmPqlOptions.defaults(), BackendPolicy.CPU);
        VariantSet set = new VariantSet("binary", List.of(
            member("v1", FIRST, EffectAllele.ALTERNATE),
            member("v2", SECOND, EffectAllele.ALTERNATE)));
        PreparedVariantSet prepared = SetTests.prepare(
            set, binary.length, options(50));

        SetTestSuiteResult result = SetTests.scoreSuite(
            prepared, nullModel, options(50));

        assertTrue(result.burden().pValue() > 0
            && result.burden().pValue() <= 1);
        assertTrue(result.skat().pValue() > 0
            && result.skat().pValue() <= 1);
        assertTrue(result.skatO().adjustedPValue() > 0
            && result.skatO().adjustedPValue() <= 1);
    }

    private static LinearSetTestNullModel nullModel() {
        return LinearSetTestNullModel.prepare(
            RESPONSE, INTERCEPT, OlsOptions.defaults(), BackendPolicy.CPU);
    }

    private static SetTestOptions options(int simulations) {
        return new SetTestOptions(VariantFilterOptions.defaults(),
            SetTestMissingPolicy.MEAN_IMPUTE,
            new double[] {0, 0.5, 1}, simulations, 42);
    }

    private static WeightedVariant member(
            String id, double[] dosage, EffectAllele effectAllele) {
        return new WeightedVariant(new VariantRecord(
            id, "1", 1, "A", "G", dosage, Double.NaN),
            effectAllele, 1);
    }

    private static final class CountingNullModel
            implements GaussianSetTestNullModel {
        private final GaussianSetTestNullModel delegate;
        private int calls;

        private CountingNullModel(GaussianSetTestNullModel delegate) {
            this.delegate = delegate;
        }

        @Override public int observations() { return delegate.observations(); }
        @Override public double degreesOfFreedom() {
            return delegate.degreesOfFreedom();
        }
        @Override public SetTestScoreState score(double[][] variantRows) {
            calls++;
            return delegate.score(variantRows);
        }
        @Override public ComputeBackend computeBackend() {
            return delegate.computeBackend();
        }
        @Override public BackendPolicy backendPolicy() {
            return delegate.backendPolicy();
        }
    }

    private static final class FixedScoreNullModel
            implements SetTestScoreNullModel {
        private final SetTestScoreState state;
        private final ComputeBackend backend = new CpuComputeBackend();

        private FixedScoreNullModel(double[] scores, double[] information) {
            state = new SetTestScoreState(
                scores, information, scores.length);
        }

        @Override public int observations() { return RESPONSE.length; }
        @Override public SetTestScoreState score(double[][] variantRows) {
            assertEquals(state.variants(), variantRows.length);
            return state;
        }
        @Override public ComputeBackend computeBackend() { return backend; }
        @Override public BackendPolicy backendPolicy() {
            throw new AssertionError(
                "retained backend should avoid policy selection");
        }
    }
}
