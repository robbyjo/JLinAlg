/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.util.*;
import org.jlinalg.xwas.*;

/** Reproducible synthetic workload timings; excludes input parsing and JVM startup. */
public final class XwasFollowupBenchmark {
    private static volatile double checksum;
    private XwasFollowupBenchmark() { }
    public static void main(String[] args) {
        Random random=new Random(421);
        int n=20_000,t=4;
        double[] ld=new double[n],wld=new double[n];double[][] z=new double[n][t],ns=new double[n][t];int[] blocks=new int[n];
        for(int i=0;i<n;i++) {
            ld[i]=1+(i*17%97)/5.0;wld[i]=1+ld[i]/2;blocks[i]=i/100;
            double shared=random.nextGaussian();
            for(int j=0;j<t;j++){ns[i][j]=10000;z[i][j]=Math.sqrt(.1*ld[i])*(.4+.1*j)*shared+Math.sqrt(1+.02*ld[i])*random.nextGaussian();}
        }
        time("ldsc_20000_variants_4_traits_200_blocks",()->checksum+=LdScoreRegression.fit(ld,wld,z,ns,blocks,100000).geneticCovariance()[0]);
        int p=100;double[] r=new double[p*p],sd=new double[p],gz=new double[p],weights=new double[p];
        for(int i=0;i<p;i++){sd[i]=.6;gz[i]=random.nextGaussian();weights[i]=random.nextGaussian()*.1;for(int j=0;j<p;j++)r[i*p+j]=Math.pow(.5,Math.abs(i-j));}
        time("twas_100_variants",()->checksum+=PredictedOmics.test(gz,weights,sd,r).z());
        double[] load={.5,.6,.7,.8},s=new double[16],v=new double[100];
        for(int i=0;i<4;i++)for(int j=0;j<4;j++)s[i*4+j]=load[i]*load[j]+(i==j?.2:0);
        for(int i=0;i<10;i++)v[i*10+i]=.0001;
        time("factor_4_traits",()->checksum+=GenomicFactor.fit(s,4,v).loadings()[0]);
        double[][] x=new double[2000][50];double[] y=new double[2000];
        for(int i=0;i<y.length;i++){for(int j=0;j<50;j++)x[i][j]=random.nextGaussian();y[i]=x[i][0]-.5*x[i][1]+random.nextGaussian();}
        time("ridge_score_2000_samples_50_features",()->checksum+=PredictionScores.train(y,x,new double[]{.1},0,5,42).weights()[0]);
        System.out.println("checksum="+checksum);
    }
    private static void time(String name,Runnable action) {
        for(int i=0;i<2;i++)action.run();double[] ms=new double[5];
        for(int i=0;i<ms.length;i++){long start=System.nanoTime();action.run();ms[i]=(System.nanoTime()-start)/1e6;}
        Arrays.sort(ms);System.out.printf(Locale.ROOT,"%s\tmedian_ms=%.3f\tmin_ms=%.3f\tmax_ms=%.3f%n",name,ms[2],ms[0],ms[4]);
    }
}
