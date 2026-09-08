/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glm.GlmFamily;

/** Also runnable without JUnit for isolated package validation during parallel repairs. */
public final class GlmmQuadratureReferenceCheck {
    private GlmmQuadratureReferenceCheck() { }
    public static void main(String[] args) throws IOException {
        verify(Path.of("src/test/resources/r-reference"));
    }
    static void verify(Path root) throws IOException {
        for (String line : Files.readAllLines(root.resolve("quadrature-integrals.tsv")).subList(1,9)) {
            String[] a = line.split("\t"); int n = Integer.parseInt(a[2]);
            double[] y = new double[n]; double[][] x = new double[n][1];
            double sum = Double.parseDouble(a[3]);
            GlmFamily family = a[1].equals("binomial") ? GlmFamilies.binomial() : GlmFamilies.poisson();
            if (family == GlmFamilies.binomial()) for (int i = 0; i < (int)sum; i++) y[i] = 1;
            else if (a[0].equals("countlarge")) for (int i = 0; i < n; i++) y[i] = i % 2 == 0 ? 90 : 110;
            else if (a[0].equals("countmixed")) { y[3] = 1; y[4] = 2; }
            for (double[] row : x) row[0] = 1;
            GlmmQuadratureEvaluation result = GlmmQuadrature.evaluate(y,x,
                java.util.Collections.nCopies(n,"a"),family,new double[]{Double.parseDouble(a[4])},
                Double.parseDouble(a[5]),GlmmQuadratureOptions.defaults());
            near(result.logLikelihood(),Double.parseDouble(a[6]),2e-9,a[0]);
            require(result.converged(),a[0]+" integration convergence: "+result);
            System.out.println(a[0]+" "+result);
        }
        Properties references = new Properties();
        try (Reader reader = Files.newBufferedReader(root.resolve("quadrature-fits.properties"))) { references.load(reader); }
        for (String kind : List.of("binary","rare","binomial","poisson")) {
            List<String> lines = Files.readAllLines(root.resolve("quadrature-"+kind+".tsv"));
            int n = lines.size()-1; double[] y = new double[n], offset = new double[n], trials = new double[n];
            double[][] x = new double[n][2]; List<String> labels = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                String[] a = lines.get(i+1).split("\t"); y[i]=Double.parseDouble(a[0]);
                x[i][0]=1; x[i][1]=Double.parseDouble(a[1]); labels.add(a[2]);
                offset[i]=Double.parseDouble(a[3]); trials[i]=Double.parseDouble(a[4]);
            }
            long start = System.nanoTime();
            GlmmQuadratureResult result = GlmmQuadrature.fit(y,x,labels,
                kind.equals("poisson")?GlmFamilies.poisson():GlmFamilies.binomial(),trials,offset,GlmmQuadratureOptions.defaults());
            double seconds = (System.nanoTime()-start)*1e-9;
            System.out.println(kind+" beta="+Arrays.toString(result.beta())+" sd="+result.randomStandardDeviation()
                +" ll="+result.logLikelihood()+" se="+Arrays.toString(result.standardErrors())+" score="+result.gradientNorm()
                +" status="+result.status()+" seconds="+seconds);
            require(result.converged(),kind+" optimizer convergence: "+result.status());
            require(result.jointInferenceAvailable(),kind+" joint inference");
            near(result.beta()[0],number(references,kind,"intercept"),3e-4,kind+" intercept");
            near(result.beta()[1],number(references,kind,"slope"),3e-4,kind+" slope");
            near(result.randomStandardDeviation(),number(references,kind,"sd"),3e-4,kind+" sd");
            near(result.logLikelihood(),number(references,kind,"logLikelihood"),2e-6,kind+" LL");
            near(result.standardErrors()[0],number(references,kind,"se0"),3e-4,kind+" se0");
            near(result.standardErrors()[1],number(references,kind,"se1"),3e-4,kind+" se1");
            near(result.fixedEffectCovariance()[0][1],number(references,kind,"cov01"),2e-4,kind+" cov01");
            near(result.parameterCovariance()[0][2],number(references,kind,"cov02"),2e-4,kind+" cov02");
            near(result.parameterCovariance()[1][2],number(references,kind,"cov12"),2e-4,kind+" cov12");
            near(result.parameterCovariance()[2][2],number(references,kind,"varVariance"),2e-4,kind+" variance covariance");
            for(int index:new int[]{0,192,383}) {
                double expected = number(references,kind,"mean"+index);
                near(result.fittedMeans()[index],expected,1e-6*Math.max(1,Math.abs(expected)),kind+" marginal mean "+index);
            }
        }
    }
    private static double number(Properties p,String kind,String key) { return Double.parseDouble(p.getProperty(kind+"."+key)); }
    private static void near(double actual,double expected,double tolerance,String label) {
        require(Double.isFinite(actual) && Math.abs(actual-expected)<=tolerance,label+": "+actual+" expected "+expected);
    }
    private static void require(boolean value,String label) { if (!value) throw new AssertionError(label); }
}
