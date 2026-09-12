/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.lang.management.ManagementFactory;
import java.util.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.meta.*;

/** Synthetic missing-cohort macrobenchmark; includes preparation and consumes results. */
public final class MetaArrayBenchmark {
    private static volatile double consumed;
    private MetaArrayBenchmark() { }
    public static void main(String[] args) {
        int rows=args.length==0?5000:Integer.parseInt(args[0]),cohorts=8;
        double[] y=new double[rows*cohorts],se=new double[y.length];
        for(int i=0;i<rows;i++)for(int j=0;j<cohorts;j++) {
            y[i*cohorts+j]=.2*Math.sin(i*.13+j)+.15*Math.cos(j*2.1);
            se[i*cohorts+j]=.05+.01*j;
            if(j==5&&i%2==0) { y[i*cohorts+j]=Double.NaN;se[i*cohorts+j]=Double.NaN; }
        }
        var bean=(com.sun.management.ThreadMXBean)ManagementFactory.getThreadMXBean();
        if(bean.isThreadAllocatedMemorySupported()) bean.setThreadAllocatedMemoryEnabled(true);
        System.out.println("model,engine,rows,median_ms,median_allocated_bytes,max_beta_error,checksum");
        for(var options:List.of(MetaAnalysisOptions.fixedEffect(),MetaAnalysisOptions.randomEffects())) {
            double[] reference=scalar(y,se,rows,cohorts,options),actual=batch(y,se,rows,cohorts,options);
            double error=0;for(int i=0;i<rows;i++)error=Math.max(error,Math.abs(reference[i]-actual[i]));
            if(error>3e-7)throw new AssertionError("scalar/batch discrepancy: "+error);
            for(boolean arrays:new boolean[]{false,true}) {
                long[] nanos=new long[3],bytes=new long[3];
                for(int run=0;run<3;run++) {
                    long before=allocated(bean),start=System.nanoTime();
                    double[] result=arrays?batch(y,se,rows,cohorts,options):scalar(y,se,rows,cohorts,options);
                    nanos[run]=System.nanoTime()-start;bytes[run]=allocated(bean)-before;
                    consumed=Arrays.stream(result).sum();
                }
                Arrays.sort(nanos);Arrays.sort(bytes);
                System.out.printf(Locale.ROOT,"%s,%s,%d,%.3f,%d,%.6g,%.12g%n",options.method(),
                    arrays?"array_batch":"study_list",rows,nanos[1]/1e6,bytes[1],error,consumed);
            }
        }
    }
    private static long allocated(com.sun.management.ThreadMXBean bean) {
        return bean.isThreadAllocatedMemorySupported()?bean.getThreadAllocatedBytes(Thread.currentThread().getId()):0;
    }
    private static double[] batch(double[] y,double[] se,int rows,int k,MetaAnalysisOptions options) {
        return MetaAnalysis.prepareBatch(y,se,rows,k,1).fit(options).pooledEffectSizes();
    }
    private static double[] scalar(double[] y,double[] se,int rows,int k,MetaAnalysisOptions options) {
        double[] result=new double[rows];
        for(int i=0;i<rows;i++) {
            List<MetaStudy> studies=new ArrayList<>(k);
            for(int j=0;j<k;j++)if(!Double.isNaN(y[i*k+j]))studies.add(new MetaStudy("c"+j,y[i*k+j],se[i*k+j]));
            result[i]=MetaAnalysis.fit(studies,options,BackendPolicy.CPU).pooledEffectSize();
        }
        return result;
    }
}
