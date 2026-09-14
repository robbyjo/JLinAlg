/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.confounding;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Random;
import org.junit.jupiter.api.Test;

class ConfounderMethodsTest {
    @Test
    void pcaRecoversDominantSampleFactorAndDeterministicOrientation() {
        double[] latent = {-2, -1.5, -1, -.5, 0, .5, 1, 1.5};
        double[][] data = new double[20][latent.length];
        for (int feature = 0; feature < data.length; feature++)
            for (int sample = 0; sample < latent.length; sample++)
                data[feature][sample] = (feature + 1.0) * latent[sample]
                    + 0.001 * Math.sin(feature + sample);
        LatentFactorResult fit = PrincipalComponentConfounders.fit(data,
            PrincipalComponentConfounders.Options.defaults(1));
        assertEquals(latent.length, fit.factors().length);
        assertTrue(Math.abs(correlation(latent, column(fit.factors(), 0))) > 0.999999);
        assertTrue(fit.varianceExplained()[0] > 0.999999);
    }

    @Test
    void standardAndAutoSvaProduceFiniteProtectedFactors() {
        Fixture fixture = fixture(40, 16);
        SvaResult standard = SurrogateVariableAnalysis.fit(fixture.data,
            fixture.full, fixture.nullDesign,
            new SurrogateVariableAnalysis.Options(1, 3));
        assertEquals(16, standard.surrogateVariables().length);
        assertEquals(40, standard.featureWeights().length);
        assertFinite(standard.surrogateVariables());
        SvaResult automatic = AutoSva.fit(fixture.data, fixture.full,
            fixture.nullDesign, AutoSva.Options.fixed(1));
        assertEquals(1, automatic.factorCount());
        assertTrue(automatic.iterations() <= 20);
        assertFinite(automatic.adjusted());
    }

    @Test
    void combatRemovesKnownBatchLocationWhilePreservingCovariate() {
        int n = 12, features = 30;
        String[] batch = new String[n];
        double[][] covariates = new double[n][2];
        double[][] data = new double[features][n];
        for (int sample = 0; sample < n; sample++) {
            batch[sample] = sample < 6 ? "A" : "B";
            covariates[sample][0] = 1.0;
            covariates[sample][1] = sample % 2;
            for (int feature = 0; feature < features; feature++)
                data[feature][sample] = feature * 0.03 + 2.0 * covariates[sample][1]
                    + (sample < 6 ? -3.0 : 4.0) + 0.1 * Math.sin(feature * 2.0 + sample);
        }
        ComBat.Result result = ComBat.adjust(data, batch, covariates, ComBat.Options.defaults());
        double before = Math.abs(groupMean(data, 0, 6) - groupMean(data, 6, 12));
        double after = Math.abs(groupMean(result.adjusted(), 0, 6)
            - groupMean(result.adjusted(), 6, 12));
        assertTrue(after < before * 0.05, "known batch location should be strongly reduced");
        assertEquals(0, result.unadjustedFeatures().length);
    }

    @Test
    void combatReferenceBatchIsBitwiseUnchanged() {
        Fixture fixture = fixture(20, 10);
        String[] batch = new String[10];
        for (int i = 0; i < batch.length; i++) batch[i] = i < 5 ? "one" : "two";
        ComBat.Options options = new ComBat.Options(true, false, "one", 1e-4, 10000);
        ComBat.Result result = ComBat.adjust(fixture.data, batch, null, options);
        for (int feature = 0; feature < fixture.data.length; feature++)
            assertArrayEquals(java.util.Arrays.copyOfRange(fixture.data[feature], 0, 5),
                java.util.Arrays.copyOfRange(result.adjusted()[feature], 0, 5));
    }

    @Test
    void peerReducesResidualVarianceAndReturnsRequestedFactors() {
        Fixture fixture = fixture(24, 14);
        Peer.Options options = new Peer.Options(2, 100, 1e-3, 1e-9,
            true, 42L, 0.001, 0.1, 0.1, 10.0, 100.0);
        PeerResult result = Peer.fit(fixture.data, null, options);
        assertEquals(14, result.factors().length);
        assertEquals(2, result.factors()[0].length);
        assertEquals(24, result.loadings().length);
        assertFinite(result.factors());
        assertTrue(result.residualVariances().get(result.residualVariances().size() - 1)
            < result.residualVariances().get(0));
    }

    private static Fixture fixture(int features, int samples) {
        Random random = new Random(7);
        double[][] data = new double[features][samples];
        double[][] full = new double[samples][2], reduced = new double[samples][1];
        for (int sample = 0; sample < samples; sample++) {
            full[sample][0] = reduced[sample][0] = 1.0;
            full[sample][1] = sample % 2;
            double hidden = Math.sin(sample * 0.7);
            for (int feature = 0; feature < features; feature++)
                data[feature][sample] = (feature % 3 - 1) * hidden
                    + (feature % 5 == 0 ? 1.5 * full[sample][1] : 0.0)
                    + random.nextGaussian() * 0.15;
        }
        return new Fixture(data, full, reduced);
    }

    private static double[] column(double[][] matrix, int column) { double[] result = new double[matrix.length]; for (int i = 0; i < result.length; i++) result[i] = matrix[i][column]; return result; }
    private static double correlation(double[] x, double[] y) { double xm = java.util.Arrays.stream(x).average().orElseThrow(), ym = java.util.Arrays.stream(y).average().orElseThrow(), xy = 0, xx = 0, yy = 0; for (int i = 0; i < x.length; i++) { xy += (x[i]-xm)*(y[i]-ym); xx += (x[i]-xm)*(x[i]-xm); yy += (y[i]-ym)*(y[i]-ym); } return xy / Math.sqrt(xx*yy); }
    private static double groupMean(double[][] data, int from, int to) { double sum = 0; int count = 0; for (double[] row : data) for (int j = from; j < to; j++) { sum += row[j]; count++; } return sum / count; }
    private static void assertFinite(double[][] matrix) { for (double[] row : matrix) for (double value : row) assertTrue(Double.isFinite(value)); }
    private record Fixture(double[][] data, double[][] full, double[][] nullDesign) { }
}
