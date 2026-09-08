/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;
import java.nio.file.*;
import java.util.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glmm.*;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.reml.*;

/** Full-fit timings with frozen same-likelihood accuracy gates on every run. */
public final class RemainingMixedBenchmark {
    private RemainingMixedBenchmark() { }
    private record Data(double[] y,double[] x,double[] offsets,List<String> groups) { }
    private static final Path ROOT=Path.of("src/test/resources/r-reference/remaining-mixed");
    private static double[] expected(Properties p,String key){return Arrays.stream(p.getProperty(key).split(",")).mapToDouble(Double::parseDouble).toArray();}
    private static Data data(String name) throws Exception {
        var lines=Files.readAllLines(ROOT.resolve(name+".tsv"));int n=lines.size()-1;
        double[] y=new double[n],x=new double[2*n],offsets=new double[n];List<String> groups=new ArrayList<>();
        for(int i=0;i<n;i++){var a=lines.get(i+1).split("\t");y[i]=Double.parseDouble(a[0]);x[2*i]=1;x[2*i+1]=Double.parseDouble(a[1]);groups.add(a[2]);if(a.length>3)offsets[i]=Double.parseDouble(a[3]);}
        return new Data(y,x,offsets,groups);
    }
    private static void check(double[] got,double[] want,double tolerance) {
        for(int i=0;i<got.length;i++)if(!Double.isFinite(got[i])||Math.abs(got[i]-want[i])>tolerance)
            throw new AssertionError("accuracy gate: got="+Arrays.toString(got)+" expected="+Arrays.toString(want));
    }
    public static void main(String[] args) throws Exception {
        var reference=new Properties();try(var reader=Files.newBufferedReader(ROOT.resolve("reference.properties"))){reference.load(reader);}
        System.out.println("engine\tcase\trun\tseconds\tchecksum\tLLerror");
        for(String kind:List.of("rare","poisson"))for(boolean dense:new boolean[]{false,true}) {
            Data d=data(kind);var family=kind.equals("rare")?GlmFamilies.binomial():GlmFamilies.poisson();
            var terms=List.of(RandomEffectTerm.randomIntercept("g",d.groups));var bases=List.of(VarianceComponent.randomIntercept("g",d.groups));
            var options=new GlmmLaplaceOptions(100,100,1e-8,1,1e-8,100,null);double[] times=new double[7];double checksum=0;
            for(int run=-3;run<7;run++) {
                long start=System.nanoTime();
                var fit=dense?GlmmLaplace.fit(d.y,d.x,d.y.length,2,family,bases,null,d.offsets,options,BackendPolicy.CPU)
                    :SparseGlmmLaplace.fit(d.y,d.x,d.y.length,2,family,terms,null,d.offsets,options,BackendPolicy.CPU);
                double seconds=(System.nanoTime()-start)*1e-9;
                if(!fit.converged())throw new AssertionError("nonconverged fit");
                check(fit.beta(),expected(reference,kind+".beta"),2e-6);check(fit.varianceComponents(),expected(reference,kind+".variance"),3e-6);
                check(fit.fixedEffectCovariance(),expected(reference,kind+".covariance"),3e-6);
                double error=Math.abs(fit.marginalLogLikelihood()-expected(reference,kind+".LL")[0]);if(error>2e-8)throw new AssertionError("LL error "+error);
                checksum+=fit.marginalLogLikelihood()+Arrays.stream(fit.beta()).sum()+Arrays.stream(fit.standardErrors()).sum();
                if(run>=0){times[run]=seconds;System.out.printf(Locale.ROOT,"%s\t%s\t%d\t%.9f\t%.14g\t%.3g%n",dense?"java-dense":"java-sparse",kind,run,seconds,checksum,error);}
            }
            Arrays.sort(times);System.out.printf(Locale.ROOT,"%s\t%s\tmedian\t%.9f\t%.14g\t0%n",dense?"java-dense":"java-sparse",kind,times[3],checksum);
        }
        Data d=data("reml");var bases=List.of(VarianceComponent.randomIntercept("g",d.groups),VarianceComponent.identity("residual",d.y.length));
        for(var method:VarianceEstimation.values()) {
            var options=RemlOptions.builder().varianceEstimation(method).build();double[] times=new double[7];double checksum=0;String prefix="reml."+method;
            for(int run=-3;run<7;run++) {
                long start=System.nanoTime();var fit=Reml.fit(d.y,d.x,d.y.length,2,bases,options,BackendPolicy.CPU);double seconds=(System.nanoTime()-start)*1e-9;
                if(!fit.converged())throw new AssertionError(fit.convergenceMessage());
                check(fit.beta(),expected(reference,prefix+".beta"),2e-6);check(fit.varianceComponents(),expected(reference,prefix+".variance"),2e-6);
                check(fit.fixedEffectCovariance(),expected(reference,prefix+".covariance"),2e-6);
                double error=Math.abs(fit.restrictedLogLikelihood()-expected(reference,prefix+".LL")[0]);if(error>2e-8)throw new AssertionError("LL error "+error);
                checksum+=fit.restrictedLogLikelihood()+Arrays.stream(fit.beta()).sum()+Arrays.stream(fit.standardErrors()).sum();
                if(run>=0){times[run]=seconds;System.out.printf(Locale.ROOT,"java-dense\t%s\t%d\t%.9f\t%.14g\t%.3g%n",method,run,seconds,checksum,error);}
            }
            Arrays.sort(times);System.out.printf(Locale.ROOT,"java-dense\t%s\tmedian\t%.9f\t%.14g\t0%n",method,times[3],checksum);
        }
    }
}
