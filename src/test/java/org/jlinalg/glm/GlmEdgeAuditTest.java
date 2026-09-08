package org.jlinalg.glm;

import static org.junit.jupiter.api.Assertions.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.ols.*;
import org.junit.jupiter.api.Test;

class GlmEdgeAuditTest {
    @Test void gaussianInferenceAndLikelihoodMatchWeightedOffsetOls() {
        double[] y = {1,2,5,7,9}, w = {1,2,3,2,1}, offset = {.1,-.2,.3,.1,-.2};
        double[][] x = {{1,0},{1,1},{1,2},{1,3},{1,4}};
        OlsResult ols = Ols.fit(y, x, w, offset, OlsOptions.defaults(), BackendPolicy.CPU);
        GlmResult fit = Glm.fit(y, x, GlmFamilies.gaussian(), w, offset,
            GlmOptions.defaults(), BackendPolicy.CPU);
        assertTrue(fit.converged(), fit.convergenceMessage());
        assertArrayEquals(ols.coefficients(), fit.coefficients(), 1e-12);
        assertArrayEquals(ols.standardErrors(), fit.standardErrors(), 1e-12);
        assertArrayEquals(ols.pValues(), fit.pValues(), 1e-12);
        assertArrayEquals(ols.pValues(), fit.associationStatistics().pValues(), 1e-12);
        assertArrayEquals(ols.confidenceLower(), fit.confidenceLower(), 1e-12);
        assertEquals(ols.logLikelihood(), fit.logLikelihood(), 1e-12);
        assertEquals(6 - 2*ols.logLikelihood(), fit.aic(), 1e-12);
    }

    @Test void gammaDispersionAndSlopeAreResponseScaleInvariant() {
        GlmResult reference = gamma(1);
        GlmResult tiny = gamma(1e-8);
        assertTrue(tiny.converged(), tiny.convergenceMessage());
        assertEquals(reference.dispersion(), tiny.dispersion(), 1e-12);
        assertArrayEquals(reference.standardErrors(), tiny.standardErrors(), 1e-12);
        assertEquals(reference.coefficients()[0] + Math.log(1e-8), tiny.coefficients()[0], 1e-12);
        assertEquals(reference.logLikelihood() - 4*Math.log(1e-8), tiny.logLikelihood(), 1e-10);
    }

    private static GlmResult gamma(double scale) {
        return Glm.fit(new double[] {scale,2*scale,3*scale,6*scale},
            new double[][] {{1},{1},{1},{1}}, GlmFamilies.gamma(), null, null,
            GlmOptions.defaults(), BackendPolicy.CPU);
    }

    @Test void tinyPoissonMeansAndRankDeficientGaussianInference() {
        assertEquals(Math.exp(-40), GlmFamilies.poisson().inverseLink(-40), 1e-30);
        double[][] x = {{1,0,0},{1,1,1},{1,2,2},{1,3,3},{1,4,4}};
        GlmResult fit = Glm.fit(new double[] {1,2,5,7,9}, x, GlmFamilies.gaussian(), null, null,
            GlmOptions.builder().rankDeficiencyStrategy(RankDeficiencyStrategy.MINIMUM_NORM).build(), BackendPolicy.CPU);
        assertEquals(2, fit.rank());
        assertTrue(Double.isNaN(fit.pValues()[1]));
        assertEquals(6-2*fit.logLikelihood(), fit.aic(), 1e-12);
        assertThrows(IllegalArgumentException.class, () -> fit.testContrast(new double[][] {{0,1,0}}));
    }

    @Test void iterationBudgetDoesNotMasqueradeAsConvergence() {
        GlmResult fit = Glm.fit(new double[] {1,5,3,8}, new double[][] {{1},{1},{1},{1}},
            GlmFamilies.poisson(), null, null,
            GlmOptions.builder().maximumIterations(1).initialCoefficients(-2).build(), BackendPolicy.CPU);
        assertFalse(fit.converged());
    }

    @Test void rareGroupedBinomialRecoversAnalyticProbabilityAndInformation() {
        GlmResult fit=Glm.fit(new double[]{1e-14,2e-14,3e-14},new double[][]{{1},{1},{1}},
            GlmFamilies.binomial(),new double[]{1e14,1e14,1e14},null,
            GlmOptions.defaults(),BackendPolicy.CPU);
        assertTrue(fit.converged(),fit.convergenceMessage());
        assertEquals(Math.log(2e-14/(1-2e-14)),fit.coefficients()[0],1e-10);
        assertEquals(2e-14,fit.fittedMeans()[0],1e-24);
        assertEquals(1/Math.sqrt(6*(1-2e-14)),fit.standardErrors()[0],1e-10);
    }

    @Test void gaussianFastPathRetainsInitialShapeValidation() {
        assertThrows(IllegalArgumentException.class,()->Glm.fit(new double[]{1,2,4},
            new double[][]{{1},{1},{1}},GlmFamilies.gaussian(),null,null,
            GlmOptions.builder().initialCoefficients(0,0).build(),BackendPolicy.CPU));
    }

    @Test void covarianceUsesActualTinyInformationAndExactFitLikelihoodIsUnbounded() {
        var p=Glm.fit(new double[]{1,2,3,4},new double[][]{{1},{1},{1},{1}},GlmFamilies.poisson(),
            new double[]{1e-20,1e-20,1e-20,1e-20},null,GlmOptions.defaults(),BackendPolicy.CPU);
        assertTrue(p.converged());assertEquals(Math.log(2.5),p.coefficients()[0],1e-12);
        assertEquals(1/Math.sqrt(1e-19),p.standardErrors()[0],1e-5);
        var g=Glm.fit(new double[4],new double[][]{{1},{1},{1},{1}},GlmFamilies.gaussian(),
            null,null,GlmOptions.defaults(),BackendPolicy.CPU);
        assertTrue(g.converged());assertEquals(Double.POSITIVE_INFINITY,g.logLikelihood());
        assertEquals(Double.NEGATIVE_INFINITY,g.aic());
    }

    @Test void positiveLogitTailRetainsComplementScoreInformationAndLikelihood() {
        assertEquals(Math.exp(-50)/Math.pow(1+Math.exp(-50),2),
            GlmFamilies.binomial().meanDerivative(50),1e-35);
        var fit=Glm.fit(new double[]{1,0,0,0},new double[][]{{1},{1},{1},{1}},GlmFamilies.binomial(),
            new double[]{1e22,1,1,1},new double[]{50,0,0,0},GlmOptions.defaults(),BackendPolicy.CPU);
        assertTrue(fit.converged(),fit.convergenceMessage());
        double b=fit.coefficients()[0],p=1/(1+Math.exp(-b)),q=1/(1+Math.exp(50+b));
        assertEquals(0,1e22*q-3*p,1e-9);
        assertEquals(.170010356794, b,1e-10);
        double ll=-1e22*Math.log1p(Math.exp(-50-b))-3*Math.log1p(Math.exp(b));
        assertEquals(ll,fit.logLikelihood(),1e-10);
        assertEquals(1/Math.sqrt(1e22*q*(1-q)+3*p*(1-p)),fit.standardErrors()[0],1e-10);
    }
}
