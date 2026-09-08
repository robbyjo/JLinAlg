/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Properties;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.formula.*;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.reml.RemlOptions;

/** Accuracy-gated ML/REML-equivalent R fixture timings; run scalability with -Xmx96m. */
public final class SparseMixedParityBenchmark {
    public static void main(String[] args) throws Exception {
        String[] cases = args.length == 0 ? new String[] {"sleepstudy","unbalanced","weighted-crossed"} : args;
        double checksum = 0;
        for (String name : cases) {
            Path root = Path.of(name.equals("scalability") ? "src/benchmark/resources/r-reference" : "src/test/resources/r-reference");
            var lines = Files.readAllLines(root.resolve("mixed-"+name+".tsv"));
            int n = lines.size()-1;
            double[] y=new double[n],x=new double[n],w=new double[n],o=new double[n];
            String[] g=new String[n],h=new String[n];
            for (int i=0;i<n;i++) {
                String[] f=lines.get(i+1).split("\t");
                y[i]=Double.parseDouble(f[0]); x[i]=Double.parseDouble(f[1]); g[i]=f[2]; h[i]=f[3];
                w[i]=Double.parseDouble(f[4]); o[i]=Double.parseDouble(f[5]);
            }
            var table=ModelTable.builder(n).numeric("y",y).numeric("x",x).numeric("w",w).numeric("o",o)
                .categorical("g",g).categorical("h",h).build();
            var formula=MixedFormula.compile("y~x+offset(o)+(1+x|g)"+(name.equals("weighted-crossed")?"+(1|h)":""),
                table,new FormulaOptions(ContrastCoding.TREATMENT,"w"));
            Properties reference=new Properties();
            try(var in=Files.newInputStream(root.resolve("mixed-"+name+".properties"))) { reference.load(in); }
            for (var method : new DegreesOfFreedomMethod[] {DegreesOfFreedomMethod.SATTERTHWAITE,DegreesOfFreedomMethod.KENWARD_ROGER}) {
                double[] times=new double[3];
                for(int run=-1;run<3;run++) {
                    long started=System.nanoTime();
                    var result=formula.fitSparseUnstructured(RemlOptions.builder().maximumIterations(300)
                        .degreesOfFreedomMethod(method).build(),BackendPolicy.CPU);
                    double elapsed=(System.nanoTime()-started)/1e9;
                    var fit=result.correlatedFit();
                    if(!fit.converged() || Math.abs(fit.logLikelihood()-Double.parseDouble(reference.getProperty("logLik")))>1e-5)
                        throw new AssertionError("likelihood/convergence accuracy gate failed");
                    for(int c=0;c<2;c++) {
                        String key=(method==DegreesOfFreedomMethod.KENWARD_ROGER?"krdf":"df")+c;
                        if(reference.containsKey(key)) {
                            double expected=Double.parseDouble(reference.getProperty(key));
                            if(Math.abs(fit.associationStatistics().degreesOfFreedom()[c]-expected)>1e-3*Math.max(1,expected))
                                throw new AssertionError("DF accuracy gate failed: "+key);
                        }
                        if (name.equals("scalability") && method == DegreesOfFreedomMethod.KENWARD_ROGER
                                && Math.abs(fit.associationStatistics().degreesOfFreedom()[c] - 1499) > .2)
                            throw new AssertionError("balanced-model KR DF must equal groups-1");
                    }
                    checksum+=fit.logLikelihood()+Arrays.stream(fit.beta()).sum()+Arrays.stream(fit.standardErrors()).sum()
                        +Arrays.stream(fit.randomEffects().get(0).covariance()).sum();
                    if(run>=0)times[run]=elapsed;
                    if(run==2)System.out.printf("%s %s n=%d q=%d equationNnz=%d factorNnz=%d backend=%s DF=%s ",name,method,n,
                        result.fit().randomCoefficientCount(),result.fit().equationNonzeroCount(),result.fit().factorNonzeroCount(),fit.backend().selectedBackend(),
                        Arrays.toString(fit.associationStatistics().degreesOfFreedom()));
                }
                Arrays.sort(times);
                System.out.printf("medianSeconds=%.6f RfitSeconds=%s RkrSeconds=%s checksum=%.12f%n",
                    times[1],reference.getProperty("fitSeconds"),reference.getProperty("krSeconds","not-run"),checksum);
            }
        }
    }
}
