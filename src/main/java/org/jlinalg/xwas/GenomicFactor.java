/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.xwas;

import java.util.Arrays;
import jdistlib.ChiSquare;
import org.jlinalg.genetics.GeneticCovarianceValidation;

/** Full-WLS single genetic factor with factor variance fixed to one and positive residual variances. */
public final class GenomicFactor {
    private GenomicFactor() { }
    public record Fit(double[] loadings,double[] residualVariances,double[] parameterCovariance,
            double chiSquare,int degreesOfFreedom,double pValue,int iterations) {
        public Fit {loadings=loadings.clone();residualVariances=residualVariances.clone();parameterCovariance=parameterCovariance.clone();}
        @Override public double[] loadings(){return loadings.clone();}
        @Override public double[] residualVariances(){return residualVariances.clone();}
        @Override public double[] parameterCovariance(){return parameterCovariance.clone();}
    }
    public record Association(double beta,double standardError,double z,double pValue,
            double heterogeneity,int heterogeneityDf,double heterogeneityP) { }
    /** V indexes the upper triangle of S by rows, including the diagonal. */
    public static Fit fit(double[] s,int traits,double[] v) {
        int t=traits,m=t*(t+1)/2,p=2*t;
        if(t<3||t>16)throw new IllegalArgumentException("a genetic factor requires 3 to 16 traits");
        GeneticCovarianceValidation.requirePositiveSemidefinite(s,t);
        double[] precision=SummaryMath.inverse(v,m),observed=new double[m];
        int k=0;for(int i=0;i<t;i++)for(int j=i;j<t;j++)observed[k++]=s[i*t+j];
        double[] theta=new double[p];
        for(int i=0;i<t;i++) {
            if(!(s[i*t+i]>0))throw new IllegalArgumentException("positive genetic variances are required");
            theta[i]=Math.sqrt(s[i*t+i]*.5)*(i==0||s[i]>=0?1:-1);
            theta[t+i]=Math.log(s[i*t+i]*.5);
        }
        int iteration;double[] covariance=null;double objective=Double.POSITIVE_INFINITY;
        for(iteration=0;iteration<500;iteration++) {
            State state=state(theta,observed,precision,t);
            objective=state.q;
            covariance=SummaryMath.inverse(state.information,p);
            double[] step=SummaryMath.multiply(covariance,state.score);
            double decrement=SummaryMath.dot(step,state.score);
            if(decrement<1e-12*(1+objective))break;
            boolean accepted=false;
            for(double scale=1;scale>=Math.scalb(1.0,-30);scale*=.5) {
                double[] next=theta.clone();for(int i=0;i<p;i++)next[i]+=scale*step[i];
                if(Arrays.stream(next).anyMatch(x->!Double.isFinite(x))
                        ||Arrays.stream(next,t,p).anyMatch(x->x < -50 || x > 50))continue;
                if(state(next,observed,precision,t).q<objective){theta=next;accepted=true;break;}
            }
            if(!accepted)throw new IllegalArgumentException("genetic factor line search failed; no inference returned");
        }
        if(iteration==500)throw new IllegalArgumentException("genetic factor did not converge");
        double[] load=Arrays.copyOf(theta,t),res=new double[t];
        double sign=load[0]<0?-1:1;
        for(int i=0;i<t;i++) {
            load[i]*=sign;res[i]=Math.exp(theta[t+i]);
            if(res[i]<1e-8*s[i*t+i])throw new IllegalArgumentException("boundary residual variance; ordinary factor inference is unavailable");
        }
        // Delta method from (loadings, log residual variances) to reported parameters.
        for(int i=0;i<p;i++)for(int j=0;j<p;j++)
            covariance[i*p+j]*=(i<t?sign:res[i-t])*(j<t?sign:res[j-t]);
        return new Fit(load,res,covariance,objective,m-p,m==p?Double.NaN:ChiSquare.cumulative(objective,m-p,false,false),iteration);
    }
    private record State(double q,double[] score,double[] information) { }
    private static State state(double[] theta,double[] observed,double[] w,int t) {
        int m=observed.length,p=theta.length,k=0;
        double[] residual=new double[m];double[][] jac=new double[m][p];
        for(int i=0;i<t;i++)for(int j=i;j<t;j++,k++) {
            residual[k]=observed[k]-theta[i]*theta[j]-(i==j?Math.exp(theta[t+i]):0);
            jac[k][i]+=theta[j];jac[k][j]+=theta[i];
            if(i==j)jac[k][t+i]=Math.exp(theta[t+i]);
        }
        double[] wr=SummaryMath.multiply(w,residual),score=new double[p],info=new double[p*p];
        double[][] wj=new double[m][p];
        for(int i=0;i<m;i++)for(int j=0;j<m;j++)for(int c=0;c<p;c++)wj[i][c]+=w[i*m+j]*jac[j][c];
        for(int i=0;i<m;i++)for(int a=0;a<p;a++) {
            score[a]+=jac[i][a]*wr[i];
            for(int b=0;b<p;b++)info[a*p+b]+=jac[i][a]*wj[i][b];
        }
        return new State(SummaryMath.dot(residual,wr),score,info);
    }
    /** GLS common-factor SNP test, conditional on fixed measurement loadings. */
    public static Association associate(double[] effects,double[] samplingCovariance,double[] loadings) {
        SummaryMath.finite(effects);SummaryMath.finite(loadings);
        int n=effects.length;
        if(n<2||loadings.length!=n)throw new IllegalArgumentException("invalid factor association dimensions");
        double[] w=SummaryMath.inverse(samplingCovariance,n),wl=SummaryMath.multiply(w,loadings);
        double information=SummaryMath.dot(loadings,wl);
        if(!(information>0))throw new IllegalArgumentException("factor loadings have no information");
        double beta=SummaryMath.dot(effects,wl)/information,se=1/Math.sqrt(information),z=beta/se;
        double[] residual=effects.clone();for(int i=0;i<n;i++)residual[i]-=beta*loadings[i];
        double q=Math.max(0,SummaryMath.dot(residual,SummaryMath.multiply(w,residual)));
        return new Association(beta,se,z,SummaryMath.p(z),q,n-1,ChiSquare.cumulative(q,n-1,false,false));
    }
}
