/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

import static org.junit.jupiter.api.Assertions.*;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.formula.*;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

/** Frozen lme4/lmerTest/pbkrtest comparisons, including a nontrivial KR adjustment. */
public final class SparseMixedRReferenceTest {
    record Fixture(CompiledMixedFormula formula, Properties reference) { }
    static Fixture fixture(String name) throws Exception {
        String prefix = "/r-reference/mixed-" + name;
        List<String> lines;
        try (InputStream in = SparseMixedRReferenceTest.class.getResourceAsStream(prefix + ".tsv")) {
            lines = new String(in.readAllBytes(), StandardCharsets.UTF_8).lines().toList();
        }
        int n = lines.size() - 1;
        double[] y = new double[n], x = new double[n], w = new double[n], o = new double[n];
        String[] g = new String[n], h = new String[n];
        for (int i = 0; i < n; i++) {
            String[] f = lines.get(i + 1).split("\t");
            y[i] = Double.parseDouble(f[0]); x[i] = Double.parseDouble(f[1]);
            g[i] = f[2]; h[i] = f[3]; w[i] = Double.parseDouble(f[4]); o[i] = Double.parseDouble(f[5]);
        }
        ModelTable table = ModelTable.builder(n).numeric("y",y).numeric("x",x).numeric("w",w)
            .numeric("o",o).categorical("g",g).categorical("h",h).build();
        CompiledMixedFormula formula = MixedFormula.compile("y~x+offset(o)+(1+x|g)"
            + (name.equals("weighted-crossed") ? "+(1|h)" : name.equals("two-correlated") ? "+(1+x|h)" : ""), table,
            new FormulaOptions(ContrastCoding.TREATMENT,"w"));
        Properties reference = new Properties();
        try (InputStream in = SparseMixedRReferenceTest.class.getResourceAsStream(prefix + ".properties")) {
            reference.load(in);
        }
        return new Fixture(formula, reference);
    }
    static double ref(Properties p, String key) { return Double.parseDouble(p.getProperty(key)); }
    static void near(Properties p, String key, double actual, double relative) {
        assertEquals(ref(p,key),actual,relative * Math.max(1,Math.abs(ref(p,key))),key);
    }
    static void compare(String name, DegreesOfFreedomMethod method) throws Exception {
        Fixture data = fixture(name);
        long started = System.nanoTime();
        var result = data.formula().fitSparseUnstructured(RemlOptions.builder()
            .maximumIterations(300).degreesOfFreedomMethod(method).build(),BackendPolicy.CPU);
        double seconds = (System.nanoTime() - started) / 1e9;
        var fit = result.correlatedFit();
        System.out.printf("%s %s %.4fs ll=%.12f beta=%s covariance=%s SE=%s DF=%s converged=%s eval=%d%n",
            name,method,seconds,fit.logLikelihood(),Arrays.toString(fit.beta()),
            Arrays.toString(fit.randomEffects().get(0).covariance()),Arrays.toString(fit.standardErrors()),
            Arrays.toString(fit.associationStatistics().degreesOfFreedom()),fit.converged(),result.evaluations());
        assertTrue(fit.converged(),"optimizer must be stationary");
        Properties p = data.reference();
        near(p,"logLik",fit.logLikelihood(),1e-8);
        for (int i = 0; i < 2; i++) {
            near(p,"beta"+i,fit.beta()[i],1e-5);
            near(p,(method == DegreesOfFreedomMethod.KENWARD_ROGER ? "krse" : "se")+i,fit.standardErrors()[i],1e-4);
            if (method != DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION)
                near(p,(method == DegreesOfFreedomMethod.KENWARD_ROGER ? "krdf" : "df")+i,
                    fit.associationStatistics().degreesOfFreedom()[i],5e-4);
        }
        double[] covariance = fit.randomEffects().get(0).covariance();
        for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++)
            near(p,"cov"+i+j,covariance[2*i+j],2e-4);
        near(p,"residual",fit.residualVariance(),2e-5);
        if (name.equals("weighted-crossed") || name.equals("two-correlated"))
            near(p,"secondVariance",fit.randomEffects().get(1).covariance()[0],2e-4);
        if (name.equals("two-correlated")) {
            near(p,"secondCovariance",fit.randomEffects().get(1).covariance()[1],2e-4);
            near(p,"secondSlopeVariance",fit.randomEffects().get(1).covariance()[3],2e-4);
        }
        // Sparse equation dimension is q, not n; coefficient and factor nnz are reported.
        assertTrue(result.fit().equationNonzeroCount() < data.formula().fixed().rows() * 12);
    }
    @Test void balancedMatchesR() throws Exception {
        compare("sleepstudy",DegreesOfFreedomMethod.SATTERTHWAITE);
        compare("sleepstudy",DegreesOfFreedomMethod.KENWARD_ROGER);
    }
    @Test void unbalancedMatchesRIncludingKrInflation() throws Exception {
        compare("unbalanced",DegreesOfFreedomMethod.SATTERTHWAITE);
        compare("unbalanced",DegreesOfFreedomMethod.KENWARD_ROGER);
    }
    @Test void weightedCrossedMatchesR() throws Exception {
        compare("weighted-crossed",DegreesOfFreedomMethod.SATTERTHWAITE);
        compare("weighted-crossed",DegreesOfFreedomMethod.KENWARD_ROGER);
    }
    @Test void twoUnstructuredBlocksMatchR() throws Exception {
        compare("two-correlated",DegreesOfFreedomMethod.SATTERTHWAITE);
        compare("two-correlated",DegreesOfFreedomMethod.KENWARD_ROGER);
    }
    @Test void formulaProfileEndpointsMatchR() throws Exception {
        Fixture data = fixture("sleepstudy");
        var ci = data.formula().profileFixedEffect(1,.95,3,18,
            RemlOptions.builder().maximumIterations(300).build(),BackendPolicy.CPU);
        System.out.println("profile="+ci);
        assertTrue(ci.lowerFound() && ci.upperFound());
        near(data.reference(),"profileLower",ci.lower(),1e-4);
        near(data.reference(),"profileUpper",ci.upper(),1e-4);
        var options = RemlOptions.builder().maximumIterations(300).build();
        var residual = data.formula().profileResidualSd(.95,15,40,options,BackendPolicy.CPU);
        var random = data.formula().profileRandomSd(0,0,.95,0,60,options,BackendPolicy.CPU);
        var slope = data.formula().profileRandomSd(0,1,.95,0,15,options,BackendPolicy.CPU);
        System.out.println("residual="+residual+" interceptSD="+random+" slopeSD="+slope);
        near(data.reference(),"residualSdLower",residual.lower(),2e-4);
        near(data.reference(),"residualSdUpper",residual.upper(),2e-4);
        near(data.reference(),"randomSdLower",random.lower(),2e-4);
        near(data.reference(),"randomSdUpper",random.upper(),2e-4);
        near(data.reference(),"slopeSdLower",slope.lower(),2e-4);
        near(data.reference(),"slopeSdUpper",slope.upper(),2e-4);
        var correlation = data.formula().profileCorrelation(0,.95,options,BackendPolicy.CPU);
        System.out.println("correlation="+correlation);
        near(data.reference(),"correlationLower",correlation.lower(),2e-4);
        near(data.reference(),"correlationUpper",correlation.upper(),2e-4);
    }

    @Test void denseAndSparseLikelihoodsAgree() throws Exception {
        Fixture data = fixture("unbalanced");
        var formula = data.formula();
        long started = System.nanoTime();
        var dense = CorrelatedLinearMixedModel.fit(formula.fixed().response(), formula.fixed().design(),
            formula.fixed().rows(), formula.fixed().columns(), formula.correlatedRandomEffects(),
            RemlOptions.builder().maximumIterations(300).build(),BackendPolicy.CPU);
        System.out.printf("dense %.4fs ll=%.12f converged=%s%n",(System.nanoTime()-started)/1e9,dense.logLikelihood(),dense.converged());
        near(data.reference(),"logLik",dense.logLikelihood(),1e-8);
        assertTrue(dense.converged());
        var sparse = formula.fitSparseUnstructured(RemlOptions.defaults(),BackendPolicy.CPU).correlatedFit();
        assertArrayEquals(dense.randomEffects().get(0).covariance(),sparse.randomEffects().get(0).covariance(),.005);
        assertArrayEquals(dense.randomEffects().get(0).modes(),sparse.randomEffects().get(0).modes(),.002);
        assertArrayEquals(dense.fittedValues(),sparse.fittedValues(),.002);
    }

    @Test void independentFiniteDfMatchesBalancedAnova() {
        int groups = 30, size = 8, n = groups * size;
        double[] y = new double[n], x = new double[2*n]; String[] labels = new String[n];
        java.util.Random random = new java.util.Random(817);
        for (int g = 0; g < groups; g++) {
            double b = random.nextGaussian() * 2;
            for (int j = 0; j < size; j++) {
                int r = g * size + j; labels[r] = "g"+g;
                x[2*r] = 1; x[2*r+1] = j - 3.5;
                y[r] = 2 + .8*x[2*r+1] + b + random.nextGaussian();
            }
        }
        for (var method : List.of(DegreesOfFreedomMethod.SATTERTHWAITE,DegreesOfFreedomMethod.KENWARD_ROGER)) {
            var fit = SparseLinearMixedModel.fit(y,x,n,2,List.of(RandomEffectTerm.randomIntercept("g",Arrays.asList(labels))),
                RemlOptions.builder().degreesOfFreedomMethod(method).build(),BackendPolicy.CPU);
            assertEquals(groups-1,fit.associationStatistics().degreesOfFreedom()[0],.005);
            assertEquals(n-groups-1,fit.associationStatistics().degreesOfFreedom()[1],.005);
        }
    }

    @Test void precisionBasesRemainInTheLikelihood() {
        int groups = 9, size = 5, n = groups * size;
        double rho = .4, denominator = 1-rho*rho;
        int[] starts = new int[groups+1], indices = new int[3*groups-2];
        double[] values = new double[indices.length];
        for (int g = 0, k = 0; g < groups; g++) {
            starts[g] = k;
            if (g > 0) { indices[k] = g-1; values[k++] = -rho/denominator; }
            indices[k] = g; values[k++] = (g==0 || g==groups-1 ? 1 : 1+rho*rho)/denominator;
            if (g+1 < groups) { indices[k] = g+1; values[k++] = -rho/denominator; }
        }
        starts[groups] = values.length;
        double[] y = new double[n], fixed = new double[n], design = new double[n*groups];
        String[] labels = new String[n], names = new String[groups];
        for (int g = 0; g < groups; g++) names[g] = "g"+g;
        for (int r = 0; r < n; r++) {
            int g = r/size; labels[r] = names[g]; fixed[r] = 1;
            y[r] = 2 + Math.sin(g)*2 + .4*Math.sin(r*2.7);
            for (int c = 0; c <= g; c++) design[r*groups+c] = Math.pow(rho,g-c)*(c==0 ? 1 : Math.sqrt(denominator));
        }
        var options = RemlOptions.defaults();
        var sparse = SparseLinearMixedModel.fitWithPrecision(y,fixed,n,1,
            List.of(RandomEffectTerm.randomIntercept("g",Arrays.asList(labels))),
            List.of(new SparsePrecisionMatrix(groups,starts,indices,values)),options,BackendPolicy.CPU);
        var dense = LinearMixedModel.fit(y,fixed,n,1,
            List.of(RandomEffectTerm.of("g",design,n,groups,Arrays.asList(names))),options,BackendPolicy.CPU);
        assertEquals(dense.reml().restrictedLogLikelihood(),sparse.logLikelihood(),1e-7);
        assertArrayEquals(dense.reml().varianceComponents(),sparse.varianceComponents(),1e-5);
    }

    @Test void minimumVarianceBoundaryIsAllowedButInteriorFiniteDfIsNotInvented() {
        int n = 40; double[] y = new double[n], x = new double[n];
        String[] group = new String[n]; double[][] random = new double[n][1];
        for (int r = 0; r < n; r++) { y[r] = 2+(r%2==0?-1:1); x[r]=1; group[r]="g"+(r/4); random[r][0]=1; }
        var block = CorrelatedRandomEffectBlock.of("g",Arrays.asList(group),List.of("intercept"),random);
        var fit = SparseUnstructuredCorrelatedModel.fit(y,x,n,1,List.of(block),RemlOptions.defaults(),BackendPolicy.CPU);
        assertTrue(fit.converged());
        double minimum = RemlOptions.defaults().minimumVariance();
        assertTrue(fit.randomEffects().get(0).covariance()[0] >= minimum);
        assertEquals(minimum,fit.randomEffects().get(0).covariance()[0],minimum * .002);
        assertThrows(IllegalArgumentException.class, () -> SparseUnstructuredCorrelatedModel.fit(y,x,n,1,
            List.of(block),RemlOptions.builder().degreesOfFreedomMethod(DegreesOfFreedomMethod.KENWARD_ROGER).build(),BackendPolicy.CPU));
    }

    @Test void invalidProfileRefitCannotProduceAnApparentlyValidCrossing() {
        assertThrows(IllegalStateException.class, () -> MixedModelProfile.interval(
            value -> Math.abs(value)>.5 ? Double.NaN : -.5*value*value,0,0,.95,-4,4,16));
    }
    public static void main(String[] args) throws Exception {
        SparseMixedRReferenceTest test = new SparseMixedRReferenceTest();
        test.balancedMatchesR(); test.unbalancedMatchesRIncludingKrInflation();
        test.weightedCrossedMatchesR(); test.formulaProfileEndpointsMatchR();
        test.denseAndSparseLikelihoodsAgree(); test.independentFiniteDfMatchesBalancedAnova();
        test.precisionBasesRemainInTheLikelihood();
        test.twoUnstructuredBlocksMatchR();
    }
}
