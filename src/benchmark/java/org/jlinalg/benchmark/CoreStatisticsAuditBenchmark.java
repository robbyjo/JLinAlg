/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;
import java.util.Arrays;
import org.jlinalg.stats.Alternative;
import org.jlinalg.stats.CorrelationMethod;
import org.jlinalg.stats.StatisticalTests;

/** Paired deterministic workload from src/test/R/core-statistics-audit.R. */
public final class CoreStatisticsAuditBenchmark {
    private CoreStatisticsAuditBenchmark() { }
    public static void main(String[] args) {
        System.out.println("n,median_seconds,checksum");
        for(int n:new int[]{5000,20000}){
            double[]x=new double[n],y=new double[n];
            for(int i=0;i<n;i++){x[i]=i%101;y[i]=(i*37)%97;}
            for(int i=0;i<10;i++) fit(x,y);
            double[] times=new double[7];double checksum=0;
            for(int r=0;r<times.length;r++){
                long start=System.nanoTime();
                for(int k=0;k<10;k++) checksum=fit(x,y);
                times[r]=(System.nanoTime()-start)/1e10;
            }
            Arrays.sort(times);
            System.out.println(n+","+times[times.length/2]+","+checksum);
        }
    }
    private static double fit(double[]x,double[]y){
        double tau=StatisticalTests.correlation(x,y,CorrelationMethod.KENDALL,
            Alternative.TWO_SIDED,.95,false,false).estimates().get("tau");
        double reference = x.length == 5000
            ? 0.0009037736641225347 : 0.0001043595718463651;
        if(!Double.isFinite(tau) || Math.abs(tau-reference) > 1e-16)
            throw new AssertionError("Kendall estimate differs from independent R fixture");
        return tau;
    }
}
