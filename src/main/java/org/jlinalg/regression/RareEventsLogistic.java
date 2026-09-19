/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

import org.jlinalg.compute.*;
import org.jlinalg.glm.*;
import org.jlinalg.internal.LeastSquaresSolver;

/** King-Zeng first-order coefficient correction and separate prior-prevalence correction.
 * Requires independent binary observations and an explicit first-column intercept.
 * Covariance is the underlying ML Fisher covariance, a first-order approximation;
 * externally supplied population prevalence is treated as known. */
public final class RareEventsLogistic {
    private RareEventsLogistic() { }
    public static Result fit(double[] y,double[][] x,boolean biasCorrection,Double populationPrevalence) {
        if(y==null||x==null||y.length==0||x.length!=y.length||x[0]==null||x[0].length==0)
            throw new IllegalArgumentException("matching binary response and intercept design required");
        int n=y.length,p=x[0].length;double sample=0;
        for(int i=0;i<n;i++) {
            if(y[i]!=0&&y[i]!=1||x[i]==null||x[i].length!=p||x[i][0]!=1)
                throw new IllegalArgumentException("binary responses and first-column unit intercept required");sample+=y[i]/n;
        }
        if(!(sample>0&&sample<1)||populationPrevalence!=null&&!(populationPrevalence>0&&populationPrevalence<1))
            throw new IllegalArgumentException("both outcome classes and prevalence strictly inside (0,1) required");
        var fit=Glm.fit(y,x,GlmFamilies.binomial(),null,null,GlmOptions.defaults(),BackendPolicy.CPU);
        if(!fit.converged()||fit.rank()!=p)throw new IllegalArgumentException("underlying logit fit failed or separated: "+fit.convergenceMessage());
        double[] beta=fit.coefficients(),bias=new double[p];
        if(biasCorrection) {
            double w1=populationPrevalence==null?1:populationPrevalence/sample,
                w0=populationPrevalence==null?1:(1-populationPrevalence)/(1-sample);
            double[] mu=fit.fittedMeans(),eta=fit.linearPredictor(),wx=new double[n*p],w=new double[n];
            for(int i=0;i<n;i++) {
                w[i]=GlmFamilies.binomial().meanDerivative(eta[i])*(y[i]==1?w1:w0);
                for(int j=0;j<p;j++)wx[i*p+j]=Math.sqrt(w[i])*x[i][j];
            }
            try(var context=BackendContext.select(BackendPolicy.CPU)) {
                var solution=LeastSquaresSolver.solve(wx,new double[n],n,p,false,context.backend());
                double[] covariance=solution.unscaledCovariance(),target=new double[n];
                for(int i=0;i<n;i++) {
                    double q=0;for(int j=0;j<p;j++)for(int k=0;k<p;k++)q+=x[i][j]*covariance[j*p+k]*x[i][k];
                    target[i]=Math.sqrt(w[i])*.5*q*((1+w1)*mu[i]-w1);
                }
                bias=LeastSquaresSolver.solve(wx,target,n,p,false,context.backend()).coefficients();
                for(int j=0;j<p;j++)beta[j]-=bias[j];
            }
        }
        double prior=populationPrevalence==null?0:Math.log(populationPrevalence)-Math.log1p(-populationPrevalence)-Math.log(sample)+Math.log1p(-sample);
        beta[0]+=prior;
        for(double b:beta)if(!Double.isFinite(b))throw new IllegalArgumentException("nonfinite corrected coefficient");
        return new Result(beta,fit.covariance(),bias,prior,sample,biasCorrection,populationPrevalence);
    }
    public record Result(double[] coefficients,double[] covariance,double[] coefficientBias,double interceptPriorCorrection,
            double samplePrevalence,boolean biasCorrected,Double populationPrevalence) {
        public Result{coefficients=coefficients.clone();covariance=covariance.clone();coefficientBias=coefficientBias.clone();}
        @Override public double[] coefficients(){return coefficients.clone();}
        @Override public double[] covariance(){return covariance.clone();}
        @Override public double[] coefficientBias(){return coefficientBias.clone();}
    }
}
