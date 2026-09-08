/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.penalized.PolyhedralSelectiveInference;

/** Full fixed-lambda selection plus conditional inference, accuracy-gated on
 * the independent n160 R fixture. Uses tight inversion, not R's native CI grid. */
public final class SelectiveInferenceBenchmark {
    private SelectiveInferenceBenchmark() { }
    private static volatile double consumed;
    public static void main(String[] args) throws Exception {
        int n=160;double[] y=new double[n];double[][] x=new double[n][4];
        for(int i=0;i<n;i++) {
            double z=-3+6.0*i/(n-1);
            x[i]=new double[]{Math.cos(i*.7)+.2*z,Math.sin(i*.23)-.1*z,Math.cos(i*.33),Math.sin(i*.41)};
            y[i]=.5+.9*x[i][0]-.6*x[i][1]+.5*Math.sin(i*1.13);
        }
        var lines=Files.readAllLines(Path.of("src/benchmark/resources/r-reference/selective-n160.csv"));
        var validated=fit(y,x);
        if(validated.effects().size()!=lines.size()-1)throw new IllegalStateException("active-set count mismatch");
        double maxError=0;
        for(int i=0;i<validated.effects().size();i++) {
            var effect=validated.effects().get(i);String[] fields=lines.get(i+1).split(",");
            if(effect.predictorIndex()!=Integer.parseInt(fields[0]))throw new IllegalStateException("active-set identity mismatch");
            double[] actual=values(effect);
            for(int j=0;j<actual.length;j++) {
                double expected=Double.parseDouble(fields[j+1]),error=Math.abs(actual[j]-expected);
                maxError=Math.max(maxError,error);
                double tolerance=j==4?1e-8*Math.abs(expected):1e-8;
                if(!Double.isFinite(actual[j])||error>tolerance)throw new IllegalStateException("R accuracy gate failed: row="+i+" field="+j+" expected="+expected+" actual="+actual[j]);
            }
        }
        int repetitions=Integer.getInteger("jlinalg.benchmark.selective.repetitions",200);
        if(repetitions<1)throw new IllegalArgumentException("repetitions must be positive");
        for(int j=0;j<100;j++)consumed=checksum(fit(y,x));
        double[] elapsed=new double[3],checks=new double[3];
        for(int batch=0;batch<3;batch++) {
            long start=System.nanoTime();
            for(int j=0;j<repetitions;j++)checks[batch]+=checksum(fit(y,x));
            elapsed[batch]=(System.nanoTime()-start)/1e9;consumed=checks[batch];
        }
        Arrays.sort(elapsed);
        System.out.printf(java.util.Locale.ROOT,"Java_tight_inversion repetitions=%d median_seconds=%.9f checksum=%.17g max_reference_absolute_error=%.9g%n",repetitions,elapsed[1],consumed,maxError);
    }
    private static PolyhedralSelectiveInference.Result fit(double[] y,double[][] x) {
        var result=PolyhedralSelectiveInference.fit(y,x,.12,1,.5,true,.95,BackendPolicy.CPU);
        if(!result.selectionFit().converged())throw new IllegalStateException("selection did not converge");
        return result;
    }
    private static double[] values(PolyhedralSelectiveInference.Effect e) {
        return new double[]{e.estimate(),e.untruncatedStandardError(),e.truncationLower(),e.truncationUpper(),e.pValue(),e.confidenceLower(),e.confidenceUpper()};
    }
    private static double checksum(PolyhedralSelectiveInference.Result fit) {
        return fit.effects().stream().mapToDouble(e->Arrays.stream(values(e)).sum()).sum();
    }
}
