/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glmm.GlmmQuadrature;
import org.jlinalg.glmm.GlmmQuadratureOptions;
import org.jlinalg.glmm.GlmmQuadratureResult;

/** Full-fit benchmark with an accuracy/convergence gate on every timed fit.
 * Run from repository root after compiling benchmark sources. */
public final class GlmmQuadratureBenchmark {
    private GlmmQuadratureBenchmark() { }
    public static void main(String[] args) throws Exception {
        Path root=Path.of("src/test/resources/r-reference"); Properties reference=new Properties();
        try(Reader reader=Files.newBufferedReader(root.resolve("quadrature-fits.properties"))) { reference.load(reader); }
        System.out.println("family\tmedianSeconds\tchecksum\tfullLLerror\tmaxBetaError\tmaxSEerror");
        for(String kind:List.of("binary","rare","binomial","poisson")) {
            List<String> lines=Files.readAllLines(root.resolve("quadrature-"+kind+".tsv"));
            int n=lines.size()-1; double[] y=new double[n], trials=new double[n], offset=new double[n];
            double[][] x=new double[n][2]; List<String> groups=new ArrayList<>();
            for(int i=0;i<n;i++) {
                String[] a=lines.get(i+1).split("\t"); y[i]=Double.parseDouble(a[0]); x[i][0]=1;
                x[i][1]=Double.parseDouble(a[1]); groups.add(a[2]); offset[i]=Double.parseDouble(a[3]);trials[i]=Double.parseDouble(a[4]);
            }
            double[] times=new double[7]; double checksum=0, llError=0,betaError=0,seError=0;
            for(int iteration=-3;iteration<times.length;iteration++) {
                long start=System.nanoTime();
                GlmmQuadratureResult result=GlmmQuadrature.fit(y,x,groups,
                    kind.equals("poisson")?GlmFamilies.poisson():GlmFamilies.binomial(),trials,offset,GlmmQuadratureOptions.defaults());
                double seconds=(System.nanoTime()-start)*1e-9;
                llError=Math.abs(result.logLikelihood()-number(reference,kind,"logLikelihood"));
                betaError=Math.max(Math.abs(result.beta()[0]-number(reference,kind,"intercept")),Math.abs(result.beta()[1]-number(reference,kind,"slope")));
                seError=Math.max(Math.abs(result.standardErrors()[0]-number(reference,kind,"se0")),Math.abs(result.standardErrors()[1]-number(reference,kind,"se1")));
                if(!result.converged() || !result.jointInferenceAvailable() || llError>2e-6 || betaError>1e-5 || seError>1e-5
                        || Math.abs(result.randomStandardDeviation()-number(reference,kind,"sd"))>1e-5)
                    throw new AssertionError(kind+" accuracy gate failed: "+result.status());
                if(iteration>=0) { times[iteration]=seconds; checksum+=result.logLikelihood()+result.beta()[0]+result.standardErrors()[0]; }
            }
            Arrays.sort(times);
            System.out.printf(java.util.Locale.ROOT,"%s\t%.9f\t%.12f\t%.3g\t%.3g\t%.3g%n",kind,times[3],checksum,llError,betaError,seError);
        }
    }
    private static double number(Properties p,String kind,String key) { return Double.parseDouble(p.getProperty(kind+"."+key)); }
}
