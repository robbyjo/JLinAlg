/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.util.*;
import java.util.function.DoubleSupplier;
import org.jlinalg.regression.ProductKernelRegression;
import org.jlinalg.xwas.PredictedOmics;

/** Deterministic CPU macrobenchmarks; one warm-up, median of three measured runs.
 * Run with gradlew benchmarkSummaryKernelAudit. Timings are diagnostic, not JMH. */
public final class SummaryKernelAuditBenchmark {
    private SummaryKernelAuditBenchmark() { }
    private static volatile double sink;
    public static void main(String[] args) {
        for(int n:new int[]{100,500,1000}) {
            int k=30;double[] ld=new double[n*n],z=new double[n],sd=new double[n];
            double[][] weights=new double[k][n];Random random=new Random(42);Arrays.fill(sd,1);
            for(int i=0;i<n;i++) {
                z[i]=random.nextGaussian();
                for(int j=0;j<n;j++)ld[i*n+j]=Math.pow(.5,Math.abs(i-j));
                for(int j=0;j<k;j++)weights[j][i]=random.nextGaussian();
            }
            measure("joint-omics n="+n+" models="+k,()->PredictedOmics.joint(z,weights,sd,ld).chiSquare());
        }
        for(int n:new int[]{1000,5000,10000,20000}) {
            double[][] x=new double[n][1];double[] y=new double[n];Random random=new Random(39);
            for(int i=0;i<n;i++){x[i][0]=i%2;y[i]=random.nextGaussian();}
            double mean=0,ss=0;int m=n/2;
            for(int i=0;i<n;i+=2)mean+=y[i]/m;
            for(int i=0;i<n;i+=2)ss+=Math.pow(y[i]-mean,2);
            double expectedMean=mean,expectedSe=Math.sqrt(ss)/(m-1);
            measure("categorical-kernel n="+n,()->{
                var fit=ProductKernelRegression.infer(x,y,new double[]{0},
                    new ProductKernelRegression.Type[]{ProductKernelRegression.Type.UNORDERED},.95);
                if(Math.abs(fit.estimate()-expectedMean)>1e-12||Math.abs(fit.standardError()-expectedSe)>1e-12)
                    throw new AssertionError("categorical HC3 reference mismatch");
                return fit.standardError();
            });
        }
    }
    private static void measure(String name,DoubleSupplier action) {
        sink=action.getAsDouble();double[] seconds=new double[3];
        for(int i=0;i<3;i++){long start=System.nanoTime();sink=action.getAsDouble();seconds[i]=(System.nanoTime()-start)*1e-9;}
        Arrays.sort(seconds);System.out.printf(Locale.ROOT,"%s median_seconds=%.6f result=%.12g%n",name,seconds[1],sink);
    }
}
