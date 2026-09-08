/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import org.jlinalg.regression.QuantileLinearProgram;

/** Accuracy-gated, consumed exact quantile LP timings on the same matrices as rq.fit.br/fnb. */
public final class ExactQuantileBenchmark {
    private ExactQuantileBenchmark() { }
    public static void main(String[] args) throws Exception {
        int repeats=args.length==0?100:Integer.parseInt(args[0]);
        if(repeats<1)throw new IllegalArgumentException("positive repeat count required");
        Path root=Path.of("src/benchmark/resources/exact-quantile-benchmark");
        Map<String,List<String[]>> references=new LinkedHashMap<>();
        List<String> refLines=Files.readAllLines(root.resolve("reference.csv"));
        for(String line:refLines.subList(1,refLines.size())){String[] a=line.split(",");references.computeIfAbsent(a[0],ignored->new ArrayList<>()).add(a);}
        System.out.println("case,method,repeats,milliseconds_per_fit,objective_difference,max_beta_difference,iterations,relative_duality_gap,checksum");
        for(var entry:references.entrySet()){
            String name=entry.getKey();List<String> lines=Files.readAllLines(root.resolve(name+".csv"));
            int n=lines.size(),p=lines.get(0).split(",").length-1;double[] y=new double[n];double[][] x=new double[n][p];
            for(int i=0;i<n;i++){String[] a=lines.get(i).split(",");y[i]=Double.parseDouble(a[0]);for(int j=0;j<p;j++)x[i][j]=Double.parseDouble(a[j+1]);}
            double tau=Double.parseDouble(entry.getValue().get(0)[1]),target=Double.parseDouble(entry.getValue().get(0)[2]);
            var options=QuantileLinearProgram.Options.defaults();
            var result=QuantileLinearProgram.solve(y,x,tau,options);
            double error=Math.abs(result.fit().objective()-target),betaError=0;
            for(String[] a:entry.getValue())betaError=Math.max(betaError,Math.abs(result.fit().beta()[Integer.parseInt(a[3])]-Double.parseDouble(a[4])));
            if(!result.fit().converged()||error>1e-6||betaError>1e-5)throw new AssertionError(name+" accuracy gate: objective="+error+" beta="+betaError+" "+result.certificate());
            for(int k=0;k<10;k++)QuantileLinearProgram.solve(y,x,tau,options);
            double checksum=0;long start=System.nanoTime();
            for(int k=0;k<repeats;k++){
                var fit=QuantileLinearProgram.solve(y,x,tau,options);
                if(!fit.fit().converged())throw new AssertionError("nonconverged timed LP");
                for(double b:fit.fit().beta())checksum+=b;
            }
            double millis=(System.nanoTime()-start)/1e6/repeats;
            System.out.printf(Locale.ROOT,"%s,primal_dual_lp,%d,%.9f,%.12g,%.12g,%d,%.12g,%.12g%n",name,repeats,millis,error,betaError,
                result.fit().iterations(),result.certificate().relativeGap(),checksum);
        }
    }
}
