/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;
import java.nio.file.*;
import java.util.*;
import org.jlinalg.distributional.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.glmm.GlmmLaplaceOptions;

/** Warm timings are emitted only after likelihood, coefficient and convergence gates. */
public final class DistributionalAuditBenchmark {
    public static void main(String[] args) throws Exception {
        Path root=Path.of("src/test/resources/r-reference");
        List<String> data=Files.readAllLines(root.resolve("distributional-audit-beta.tsv"));int n=data.size()-1;
        double[] y=new double[n],x=new double[2*n],z=new double[2*n];
        for(int r=0;r<n;r++){String[] f=data.get(r+1).split("\t");y[r]=Double.parseDouble(f[0]);x[2*r]=z[2*r]=1;x[2*r+1]=Double.parseDouble(f[1]);z[2*r+1]=Double.parseDouble(f[2]);}
        Properties p=new Properties();try(var in=Files.newInputStream(root.resolve("distributional-audit-beta.properties"))){p.load(in);}
        double[] elapsed=new double[5];double checksum=0,maxError=0;
        for(int run=-3;run<5;run++) {
            long start=System.nanoTime();
            var fit=BetaRegression.fit(y,x,2,z,2,BetaRegressionOptions.variablePrecisionDefaults(),BackendPolicy.CPU);
            double seconds=(System.nanoTime()-start)/1e9;
            if(!fit.converged())throw new AssertionError(fit.convergenceMessage());
            if(Math.abs(fit.logLikelihood()-Double.parseDouble(p.getProperty("logLik")))>1e-7)throw new AssertionError("likelihood gate");
            double[] coefficients={fit.meanCoefficients()[0],fit.meanCoefficients()[1],fit.precisionCoefficients()[0],fit.precisionCoefficients()[1]};
            for(int c=0;c<4;c++){double error=Math.abs(coefficients[c]-Double.parseDouble(p.getProperty("coefficient"+c)));maxError=Math.max(maxError,error);if(error>1e-5)throw new AssertionError("coefficient gate");}
            if(run>=0){elapsed[run]=seconds;checksum+=fit.logLikelihood()+Arrays.stream(coefficients).sum();}
        }
        Arrays.sort(elapsed);
        System.out.printf("beta n=%d JavaMedianSeconds=%.6f RMedianSeconds=%s maxCoefficientError=%.3g checksum=%.12f RChecksumPerFit=%s%n",n,elapsed[2],p.getProperty("seconds"),maxError,checksum,p.getProperty("checksum"));
        mixed(root);
    }
    private static void mixed(Path root) throws Exception {
        Properties reference=new Properties();
        try(var in=Files.newInputStream(Path.of("src/benchmark/resources/distributional-audit/r-timings.properties"))){reference.load(in);}
        for(String name:List.of("betaMixed","zip","zinb")) {
            boolean beta=name.equals("betaMixed"),nb=name.equals("zinb");
            List<String> rows=Files.readAllLines(root.resolve(beta?"beta-mixed-glmmtmb.tsv":"zero-inflated-mixed-data.tsv"));
            int n=rows.size()-1;double[] y=new double[n],x=new double[n*2],one=new double[n];List<String> groups=new ArrayList<>();
            for(int r=0;r<n;r++){String[] f=rows.get(r+1).split("\t");y[r]=Double.parseDouble(f[nb?1:0]);x[r*2]=one[r]=1;x[r*2+1]=Double.parseDouble(f[beta?1:2]);groups.add(f[beta?2:3]);}
            var term=RandomEffectTerm.randomIntercept("group",groups);
            var opts=new ZeroInflatedMixedOptions(900,120,1e-7,.35,1e-8,1e4,1e-5,1e5,25,null);
            double[] times=new double[5];double checksum=0,error=0;
            for(int run=-3;run<5;run++) {
                long start=System.nanoTime();double ll,sum;boolean converged;
                if(beta){var fit=BetaMixedModel.fit(y,x,n,2,List.of(term),null,null,
                    new BetaMixedModelOptions(new GlmmLaplaceOptions(35,100,1e-7,.75,1e-8,100,null),10,1e-6,1e8),BackendPolicy.CPU);
                    ll=fit.marginalLogLikelihood();sum=Arrays.stream(fit.meanCoefficients()).sum();converged=fit.converged();
                } else {
                    ZeroInflatedMixedResult fit;
                    if(nb)fit=SparseZeroInflatedMixedModel.fitNegativeBinomial(y,x,2,x,2,one,1,List.of(term),null,List.of(),null,List.of(),null,opts,BackendPolicy.CPU);
                    else try(var prepared=SparseZeroInflatedMixedModel.preparePoisson(n,List.of(term),null,List.of(term),null,List.of(),opts,BackendPolicy.CPU)){
                        fit=prepared.fitWithInference(y,x,2,x,2,null,0,null);
                    }
                    ll=fit.marginalLogLikelihood();sum=Arrays.stream(fit.countCoefficients()).sum()+Arrays.stream(fit.zeroCoefficients()).sum();converged=fit.converged();
                }
                double seconds=(System.nanoTime()-start)/1e9;
                error=Math.max(error,Math.abs(ll-Double.parseDouble(reference.getProperty(name+".logLik"))));
                if(!converged||error>1e-6)throw new AssertionError(name+" convergence/likelihood gate "+error);
                if(run>=0){times[run]=seconds;checksum+=ll+sum;}
            }
            System.out.println(name+" n="+n+" warmSeconds="+Arrays.toString(times));Arrays.sort(times);
            System.out.printf("%s JavaMedianSeconds=%.6f RMedianSeconds=%s maxLikelihoodError=%.3g checksum=%.12f RChecksum=%s%n",name,times[2],reference.getProperty(name+".seconds"),error,checksum,reference.getProperty(name+".checksum"));
        }
    }
}
