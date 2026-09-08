/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;
import java.util.function.Supplier;
import org.jlinalg.sem.*;

/** Accuracy-gated warm timings for the checked-in joint SEM / lavaan fixtures.
 * Run from the repository root; no Gradle task or data download is required. */
public final class SemJointBenchmark {
    private static volatile double checksum;
    private record Fit(List<SemParameterEstimate> parameters,double[] covariance,double ll,boolean converged) { }
    private static Fit plain(SemFitResult f){return new Fit(f.parameters(),f.parameterCovariance(),f.logLikelihood(),f.converged());}
    public static void main(String[] args) throws Exception {
        Path root=Path.of("src/test/resources/r-reference/sem-joint");
        double[][] continuous=read(root.resolve("continuous.csv")),missing=read(root.resolve("missing.csv"));
        int[][] ordinal=Arrays.stream(read(root.resolve("ordinal.csv"))).map(r->Arrays.stream(r).mapToInt(v->(int)v).toArray()).toArray(int[][]::new);
        SemModel.Builder builder=SemModel.builder("x1","x2","x3","y1","y2","y3").latent("f","g")
            .fixedLoading("x1","f",1).loading("l2","x2","f",.8).loading("l3","x3","f",1.1)
            .fixedLoading("y1","g",1).loading("l5","y2","g",.9).loading("l6","y3","g",.7)
            .regression("b","g","f",.5).variance("vf","f",1).variance("vg","g",.6);
        SemModel latent=builder.build(),means=builder.meanStructure().build();
        SemModel ordered=SemModel.builder("z1","z2","z3","z4").latent("f").fixedVariance("f",1)
            .loading("a","z1","f",.8).loading("b","z2","f",.8).loading("c","z3","f",.8).loading("d","z4","f",.8)
            .fixedVariance("z1",1).fixedVariance("z2",1).fixedVariance("z3",1).fixedVariance("z4",1).build();
        System.out.println("case\tmaxEstimateError\tmaxCovarianceError\tlikelihoodError\tjavaMedianSeconds\tRMedianSeconds\tRoverJava");
        run(root,"latent",()->plain(Sem.fit(continuous,latent)));
        run(root,"means",()->plain(Sem.fit(continuous,means)));
        run(root,"fiml",()->plain(SemFiml.fit(missing,means).fit()));
        run(root,"robust",()->{SemFitResult f=Sem.fit(continuous,means);return new Fit(f.parameters(),SemInference.robust(f).parameterCovariance(),f.logLikelihood(),f.converged());});
        run(root,"ordinal",()->{SemOrdinal.Result f=SemOrdinal.fit(ordinal,new int[]{3,3,3,3},ordered);return new Fit(f.parameters(),f.parameterCovariance(),f.pairwiseLogLikelihood(),f.converged());});
        System.out.println("consumedChecksum="+checksum);
    }
    private static void run(Path root,String name,Supplier<Fit> supplier) throws Exception {
        Properties p=new Properties();try(var in=Files.newInputStream(root.resolve(name+".properties"))){p.load(in);}
        List<String> labels=List.of(p.getProperty("labels").split(","));int k=labels.size();
        double[] estimates=numbers(p,"estimates"),cov=numbers(p,"covariance");Fit fit=supplier.get();
        double maxEstimate=0,maxCovariance=0;
        for(int i=0;i<k;i++) {
            int a=labels.indexOf(fit.parameters.get(i).label());if(a<0)throw new AssertionError("unknown parameter");
            maxEstimate=Math.max(maxEstimate,Math.abs(estimates[a]-fit.parameters.get(i).estimate()));
            for(int j=0;j<k;j++){int b=labels.indexOf(fit.parameters.get(j).label());maxCovariance=Math.max(maxCovariance,Math.abs(cov[a*k+b]-fit.covariance[i*k+j]));}
        }
        double likelihoodError=Math.abs(fit.ll-numbers(p,name.equals("ordinal")?"pairwiseLogLikelihood":"logLikelihood")[0]);
        if(!fit.converged || !(maxEstimate<3e-5 && maxCovariance<3e-5 && likelihoodError<2e-5))throw new AssertionError(name+" accuracy gate: "+maxEstimate+", "+maxCovariance+", "+likelihoodError);
        for(int i=0;i<10;i++)consume(supplier.get());
        double[] seconds=new double[5];
        for(int i=0;i<seconds.length;i++){long start=System.nanoTime();for(int j=0;j<5;j++)consume(supplier.get());seconds[i]=(System.nanoTime()-start)*1e-9/5;}
        Arrays.sort(seconds);double r=numbers(p,"R.medianSeconds")[0];
        System.out.printf(java.util.Locale.ROOT,"%s\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.3f%n",name,maxEstimate,maxCovariance,likelihoodError,seconds[2],r,r/seconds[2]);
    }
    private static void consume(Fit fit){if(!fit.converged)throw new AssertionError("timed fit did not converge");checksum+=fit.ll+fit.parameters.get(0).estimate()+fit.covariance[0];}
    private static double[][] read(Path file)throws Exception{return Files.readAllLines(file).stream().skip(1).map(l->Arrays.stream(l.split(",")).mapToDouble(Double::parseDouble).toArray()).toArray(double[][]::new);}
    private static double[] numbers(Properties p,String key){return Arrays.stream(p.getProperty(key).split(",")).mapToDouble(Double::parseDouble).toArray();}
}
