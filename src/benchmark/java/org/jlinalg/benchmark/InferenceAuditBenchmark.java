/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.util.Arrays;
import java.util.Locale;
import java.util.Random;
import org.jlinalg.confounding.AutoSva;
import org.jlinalg.confounding.SurrogateVariableAnalysis;
import org.jlinalg.differential.EmpiricalBayesDifferential;

/** Deterministic CPU scaling probes, median of three after one warm-up.
 * Run with gradlew benchmarkInferenceAudit. Timings are diagnostic, not JMH. */
public final class InferenceAuditBenchmark {
    private InferenceAuditBenchmark() { }
    public static void main(String[] args) {
        double[][] design = new double[12][2];
        for(int i=0;i<12;i++){design[i][0]=1;design[i][1]=i<6?0:1;}
        for(int features:new int[]{1000,4000,20000}) {
            Random random=new Random(42);double[][] counts=new double[features][12];
            for(var row:counts)for(int j=0;j<12;j++)row[j]=random.nextInt(500)+1;
            measure("voom",features,12,()->EmpiricalBayesDifferential.fitVoom(counts,design,new double[]{0,1}));
        }
        double[][] full=new double[32][2],reduced=new double[32][1];
        for(int i=0;i<32;i++){full[i][0]=reduced[i][0]=1;full[i][1]=i<16?0:1;}
        for(int features:new int[]{2000,16000}) {
            Random random=new Random(42);double[] latent=new double[32];
            for(int j=0;j<32;j++)latent[j]=random.nextGaussian();
            double[][] data=new double[features][32];
            for(int i=0;i<features;i++)for(int j=0;j<32;j++)
                data[i][j]=random.nextGaussian()+((i%3)-1)*latent[j]+(i%5==0?full[j][1]:0);
            measure("SVA-5-iterations",features,32,()->SurrogateVariableAnalysis.fit(data,full,reduced,new SurrogateVariableAnalysis.Options(2,5)));
            measure("AutoSVA-20-iterations",features,32,()->AutoSva.fit(data,full,reduced,AutoSva.Options.fixed(2)));
        }
    }
    private static void measure(String method,int features,int samples,Runnable action) {
        action.run();double[] seconds=new double[3];
        for(int i=0;i<3;i++){long start=System.nanoTime();action.run();seconds[i]=(System.nanoTime()-start)/1e9;}
        Arrays.sort(seconds);
        System.out.printf(Locale.ROOT,"%s features=%d samples=%d median_seconds=%.6f%n",method,features,samples,seconds[1]);
    }
}
