/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

import java.util.Arrays;
import jdistlib.Normal;
import org.jlinalg.compute.*;
import org.jlinalg.internal.LeastSquaresSolver;

/** Location-scale ML for censored Gaussian and log-time AFT responses.
 * Censor codes: -1 = left censored, 0 = observed, 1 = right censored.
 * The supplied response is the observed value or known censoring limit.
 * Censoring is assumed independent conditional on the supplied design. */
public final class CensoredRegression {
    private CensoredRegression() { }
    public enum Distribution { GAUSSIAN, LOGNORMAL, WEIBULL, EXPONENTIAL }

    public static Result fit(double[] response,double[][] design,int[] censor,Distribution distribution) {
        if(response==null||response.length<3||design==null||design.length!=response.length||design[0]==null
                ||censor==null||censor.length!=response.length||distribution==null)
            throw new IllegalArgumentException("matching nonempty response, design and censor codes required");
        int n=response.length,p=design[0].length,k=p+(distribution==Distribution.EXPONENTIAL?0:1);
        if(p<1||n<=k)throw new IllegalArgumentException("more observations than parameters required");
        double[] x=new double[n*p],y=response.clone();int exact=0;
        for(int i=0;i<n;i++) {
            if(design[i]==null||design[i].length!=p||!Double.isFinite(y[i])||censor[i]<-1||censor[i]>1)
                throw new IllegalArgumentException("invalid design, response, or censor code");
            if(censor[i]==0)exact++;
            if(distribution!=Distribution.GAUSSIAN){if(!(y[i]>0))throw new IllegalArgumentException("AFT times and limits must be positive");y[i]=Math.log(y[i]);}
            for(int j=0;j<p;j++){x[i*p+j]=design[i][j];if(!Double.isFinite(x[i*p+j]))throw new IllegalArgumentException("finite design required");}
        }
        if(exact<2)throw new IllegalArgumentException("at least two uncensored observations required");
        double[] scales=RegressionOptimizer.scaleColumns(x,n,p),base=new double[k];
        try(var context=BackendContext.select(BackendPolicy.CPU)) {
            var initial=LeastSquaresSolver.solve(x,y,n,p,false,context.backend());
            System.arraycopy(initial.coefficients(),0,base,0,p);
            if(k>p){double sse=0;for(int i=0;i<n;i++){double r=y[i];for(int j=0;j<p;j++)r-=x[i*p+j]*base[j];sse+=r*r;}
                if(!(sse>0))throw new IllegalArgumentException("zero initial residual scale");base[p]=Math.log(Math.sqrt(sse/n));}
        }
        RegressionOptimizer.Objective objective=delta->{double[] theta=add(base,delta);return evaluate(y,x,censor,n,p,distribution,theta);};
        var optimized=RegressionOptimizer.minimize(objective,k,2000,2e-8,1);
        if(!optimized.converged())throw new IllegalArgumentException("censored likelihood failed score convergence");
        double[] point=optimized.point(),theta=add(base,point);
        double[] covariance=LikelihoodInference.covariance(objective,point,n);
        for(int j=0;j<p;j++)theta[j]/=scales[j];
        for(int j=0;j<k;j++)for(int l=0;l<k;l++)covariance[j*k+l]/=(j<p?scales[j]:1)*(l<p?scales[l]:1);
        double likelihood=-objective.evaluate(point).value()*n;
        // Densities for uncensored AFT observations include the log-time Jacobian.
        if(distribution!=Distribution.GAUSSIAN)for(int i=0;i<n;i++)if(censor[i]==0)likelihood-=y[i];
        return new Result(Arrays.copyOf(theta,p),k>p?Math.exp(theta[p]):1,covariance,likelihood,distribution,optimized.iterations());
    }

    private static RegressionOptimizer.Evaluation evaluate(double[] y,double[] x,int[] censor,int n,int p,Distribution distribution,double[] theta) {
        int k=theta.length;double scale=k>p?Math.exp(theta[p]):1,value=0;double[] gradient=new double[k];
        boolean normal=distribution==Distribution.GAUSSIAN||distribution==Distribution.LOGNORMAL;
        for(int i=0;i<n;i++) {
            double eta=0;for(int j=0;j<p;j++)eta+=x[i*p+j]*theta[j];
            double z=(y[i]-eta)/scale,loss,gEta,gScale;
            if(normal) {
                if(censor[i]==0){loss=.5*z*z+Math.log(scale)+.5*Math.log(2*Math.PI);gEta=-z/scale;gScale=1-z*z;}
                else {
                    boolean left=censor[i]<0;double logTail=Normal.cumulative(z,0,1,left,true);
                    double hazard=Math.exp(Normal.density(z,0,1,true)-logTail);
                    loss=-logTail;gEta=(left?hazard:-hazard)/scale;gScale=(left?hazard:-hazard)*z;
                }
            } else {
                double ez=Math.exp(z);
                if(censor[i]==0){loss=ez-z+Math.log(scale);gEta=(1-ez)/scale;gScale=1+z*(1-ez);}
                else if(censor[i]>0){loss=ez;gEta=-ez/scale;gScale=-z*ez;}
                else {double logCdf=z < -35 ? z : Math.log(-Math.expm1(-ez));
                    double hazard=Math.exp(z-ez-logCdf);loss=-logCdf;gEta=hazard/scale;gScale=hazard*z;}
            }
            value+=loss/n;for(int j=0;j<p;j++)gradient[j]+=gEta*x[i*p+j]/n;if(k>p)gradient[p]+=gScale/n;
        }
        return new RegressionOptimizer.Evaluation(value,gradient);
    }
    private static double[] add(double[] a,double[] b){double[] out=a.clone();for(int i=0;i<out.length;i++)out[i]+=b[i];return out;}

    /** Covariance order is beta, then log(scale), except exponential has only beta. */
    public record Result(double[] coefficients,double scale,double[] covariance,double logLikelihood,Distribution distribution,int iterations) {
        public Result {coefficients=coefficients.clone();covariance=covariance.clone();}
        @Override public double[] coefficients(){return coefficients.clone();}
        @Override public double[] covariance(){return covariance.clone();}
        private double eta(double[] row){if(row==null||row.length!=coefficients.length)throw new IllegalArgumentException("prediction row width differs");double v=0;
            for(int j=0;j<row.length;j++){if(!Double.isFinite(row[j]))throw new IllegalArgumentException("finite prediction row required");v+=row[j]*coefficients[j];}
            if(!Double.isFinite(v))throw new IllegalArgumentException("prediction exceeds numerical range");return v;}
        /** Gaussian latent mean or AFT log-time location. */
        public double location(double[] row){return eta(row);}
        /** Probability that a Gaussian latent response is at/below a limit. */
        public double leftCensoringProbability(double[] row,double limit){requireGaussian();if(!Double.isFinite(limit))throw new IllegalArgumentException("finite limit required");return Normal.cumulative(limit,eta(row),scale,true,false);}
        /** Mean of a Gaussian response clipped at known lower/upper limits. */
        public double observedMean(double[] row,double lower,double upper){
            requireGaussian();if(!(lower<upper))throw new IllegalArgumentException("ordered censoring limits required");double mean=eta(row);
            double a=(lower-mean)/scale,b=(upper-mean)/scale;
            double pa=Normal.cumulative(a,0,1,true,false),qb=Normal.cumulative(b,0,1,false,false);
            double middle=a>0?Normal.cumulative(a,0,1,false,false)-qb:Normal.cumulative(b,0,1,true,false)-pa;
            return (Double.isFinite(lower)?lower*pa:0)+(Double.isFinite(upper)?upper*qb:0)+mean*middle
                +scale*(Normal.density(a,0,1,false)-Normal.density(b,0,1,false));
        }
        /** AFT survival probability at positive time. */
        public double survival(double[] row,double time){requireAft();if(!(time>0)||!Double.isFinite(time))throw new IllegalArgumentException("positive finite time required");
            double z=(Math.log(time)-eta(row))/scale;return distribution==Distribution.LOGNORMAL?Normal.cumulative(z,0,1,false,false):Math.exp(-Math.exp(z));}
        /** AFT time quantile; this is not an observation prediction interval. */
        public double timeQuantile(double[] row,double probability){requireAft();if(!(probability>0&&probability<1))throw new IllegalArgumentException("probability in (0,1) required");
            double z=distribution==Distribution.LOGNORMAL?Normal.quantile(probability,0,1,true,false):Math.log(-Math.log1p(-probability));
            double value=Math.exp(eta(row)+scale*z);if(!Double.isFinite(value)||value==0)throw new IllegalArgumentException("time quantile exceeds numerical range");return value;}
        /** Confidence interval for an estimated AFT quantile, including scale uncertainty. */
        public double[] timeQuantileInterval(double[] row,double probability,double level){
            double value=timeQuantile(row,probability);if(!(level>0&&level<1))throw new IllegalArgumentException("confidence level in (0,1) required");
            int p=coefficients.length,k=distribution==Distribution.EXPONENTIAL?p:p+1;double[] gradient=Arrays.copyOf(row,k);
            if(k>p)gradient[p]=Math.log(value)-eta(row);double variance=0;
            for(int j=0;j<k;j++)for(int l=0;l<k;l++)variance+=gradient[j]*covariance[j*k+l]*gradient[l];
            if(variance<0||!Double.isFinite(variance))throw new IllegalArgumentException("invalid quantile variance");
            double width=Normal.quantile((1+level)/2,0,1,true,false)*Math.sqrt(variance);
            double lower=value*Math.exp(-width),upper=value*Math.exp(width);
            if(!Double.isFinite(upper)||lower==0)throw new IllegalArgumentException("quantile interval exceeds numerical range");return new double[]{value,lower,upper};
        }
        /** Time ratio exp((B-A) beta), with asymptotic log-scale Wald limits. */
        public double[] timeRatio(double[] a,double[] b,double level){requireAft();if(!(level>0&&level<1))throw new IllegalArgumentException("confidence level in (0,1) required");
            double difference=eta(b)-eta(a),variance=0;int p=coefficients.length,k=distribution==Distribution.EXPONENTIAL?p:p+1;
            for(int j=0;j<p;j++)for(int l=0;l<p;l++)variance+=(b[j]-a[j])*covariance[j*k+l]*(b[l]-a[l]);
            if(variance<0||!Double.isFinite(variance))throw new IllegalArgumentException("invalid time-ratio variance");double width=Normal.quantile((1+level)/2,0,1,true,false)*Math.sqrt(variance);
            double value=Math.exp(difference),lower=Math.exp(difference-width),upper=Math.exp(difference+width);
            if(lower==0||!Double.isFinite(upper))throw new IllegalArgumentException("time ratio exceeds numerical range");
            return new double[]{value,lower,upper};}
        private void requireGaussian(){if(distribution!=Distribution.GAUSSIAN)throw new IllegalArgumentException("Gaussian model required");}
        private void requireAft(){if(distribution==Distribution.GAUSSIAN)throw new IllegalArgumentException("AFT model required");}
    }
}
