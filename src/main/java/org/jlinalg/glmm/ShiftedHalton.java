/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

import java.util.SplittableRandom;
import java.util.function.ToDoubleFunction;
import jdistlib.Normal;

/** Mode-centred Gaussian importance integration with eight randomized QMC replicates.
 * The fixed seed gives common random numbers across likelihood evaluations. */
final class ShiftedHalton {
    private ShiftedHalton() { }
    record Estimate(double logIntegral,double error) { }
    static Estimate integrate(int dimension,int points,double[] mode,double[] lower,
            ToDoubleFunction<double[]> logKernel,double logNormalizer) {
        return integrate(dimension,points,mode,lower,logKernel,logNormalizer,null,null,0);
    }
    /** Defensive mixture with the prior prevents an overly narrow mode proposal
     * from having infinite importance-weight variance. The supported likelihoods
     * are probabilities bounded by one, so target/prior is bounded. */
    static Estimate integrate(int dimension,int points,double[] mode,double[] lower,
            ToDoubleFunction<double[]> logKernel,double logNormalizer,double[] priorLower,
            double[] precision,double logPriorDeterminant) {
        int[] primes=new int[dimension];int candidate=2;
        for(int j=0;j<dimension;j++) {
            while(!prime(candidate))candidate++;
            primes[j]=candidate++;
        }
        SplittableRandom random=new SplittableRandom(0x4a4c696e416c67L);
        double[] estimates=new double[8];
        for(int r=0;r<8;r++) {
            double[] shifts=new double[dimension];
            for(int j=0;j<dimension;j++)shifts[j]=random.nextDouble();
            double sum=Double.NEGATIVE_INFINITY;
            for(int i=1;i<=points;i++) {
                boolean usePrior=priorLower!=null&&(i%2==1);
                double[] z=new double[dimension],b=usePrior?new double[dimension]:mode.clone();double squares=0;
                double[] transform=usePrior?priorLower:lower;
                for(int j=0;j<dimension;j++) {
                    double u=radical(priorLower==null?i:(i+1)/2,primes[j])+shifts[j];u-=Math.floor(u);
                    z[j]=Normal.quantile(Math.max(1e-15,Math.min(1-1e-15,u)),0,1,true,false);
                    squares+=z[j]*z[j];
                    for(int k=0;k<=j;k++)b[j]+=transform[j*dimension+k]*z[k];
                }
                double value=logKernel.applyAsDouble(b)+.5*squares;
                if(priorLower!=null) {
                    double[] standardized=new double[dimension];double modeLog=0,priorQuadratic=0;
                    for(int j=0;j<dimension;j++) {
                        double residual=b[j]-mode[j];for(int k=0;k<j;k++)residual-=lower[j*dimension+k]*standardized[k];
                        standardized[j]=residual/lower[j*dimension+j];
                        modeLog-=Math.log(lower[j*dimension+j])+.5*standardized[j]*standardized[j];
                        for(int k=0;k<dimension;k++)priorQuadratic+=b[j]*precision[j*dimension+k]*b[k];
                    }
                    double priorFraction=(points-points/2)/(double)points;
                    double mixture=logAdd(Math.log1p(-priorFraction)+modeLog,
                        Math.log(priorFraction)+.5*(logPriorDeterminant-priorQuadratic));
                    value=logKernel.applyAsDouble(b)+.5*logPriorDeterminant-mixture;
                }
                if(Double.isNaN(value)||value==Double.POSITIVE_INFINITY)
                    return new Estimate(Double.NEGATIVE_INFINITY,Double.POSITIVE_INFINITY);
                sum=logAdd(sum,value);
            }
            estimates[r]=sum-Math.log(points)+(priorLower==null?logNormalizer:0);
        }
        double max=java.util.Arrays.stream(estimates).max().orElseThrow();
        if(!Double.isFinite(max))return new Estimate(max,Double.POSITIVE_INFINITY);
        double mean=0,ss=0;
        for(double x:estimates)mean+=Math.exp(x-max)/8;
        for(double x:estimates)ss+=Math.pow(Math.exp(x-max)-mean,2);
        return new Estimate(max+Math.log(mean),3*Math.sqrt(ss/(8*7))/mean);
    }
    private static double logAdd(double a,double b) {
        if(a==Double.NEGATIVE_INFINITY)return b;
        if(b==Double.NEGATIVE_INFINITY)return a;
        double m=Math.max(a,b);return m+Math.log1p(Math.exp(Math.min(a,b)-m));
    }
    private static double radical(int value,int base) {
        double result=0,factor=1.0/base;
        while(value>0){result+=(value%base)*factor;value/=base;factor/=base;}
        return result;
    }
    private static boolean prime(int n) {
        for(int i=2;i<=n/i;i++)if(n%i==0)return false;
        return true;
    }
}
