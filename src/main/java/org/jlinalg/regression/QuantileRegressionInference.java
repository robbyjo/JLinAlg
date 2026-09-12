/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

import java.util.Arrays;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.LeastSquaresSolver;

/** Asymptotic covariance for exact, nonsmoothed linear quantile regression.
 * Assumes independent observations, a correctly specified conditional quantile,
 * and positive continuous response densities at that quantile. The supplied
 * density path permits heteroskedasticity; the pooled kernel path assumes iid
 * errors. Neither is valid for response mass points or dependent observations.
 * Predictors explicitly include an intercept when wanted. */
public final class QuantileRegressionInference {
    private QuantileRegressionInference() { }

    /** Automatic normal-reference residual-density bandwidth for iid continuous
     * errors. h = 1.06 s n^(-1/5) satisfies h -> 0, n h -> infinity for a stable
     * positive residual scale. Requires a fixed-dimensional regular design. */
    public static Result fitExactIidKernel(double[] response,double[][] predictors,double quantile) {
        QuantileRegressionResult fit=QuantileRegression.fitExact(response,predictors,quantile);
        double[] residuals=fit.residuals();double mean=Arrays.stream(residuals).average().orElseThrow(),ss=0;
        for(double r:residuals)ss+=Math.pow(r-mean,2);
        double bandwidth=1.06*Math.sqrt(ss/(residuals.length-1))*Math.pow(residuals.length,-.2);
        if(!(bandwidth>0)||!Double.isFinite(bandwidth))throw new IllegalArgumentException("positive continuous residual scale required");
        return fitExactIidKernel(response,predictors,quantile,bandwidth);
    }

    /** CR0 cluster-robust covariance with supplied consistent conditional
     * densities. Clusters must be independent, with many non-dominant clusters;
     * continuous conditional responses and a regular conditional quantile are
     * still required. Within-cluster dependence is unrestricted. */
    public static Result fitExactCluster(double[] response,double[][] predictors,double quantile,
            double[] conditionalDensities,int[] clusters) {
        validateDensities(conditionalDensities,response==null?0:response.length);
        if(clusters==null||clusters.length!=response.length||Arrays.stream(clusters).distinct().count()<2)
            throw new IllegalArgumentException("at least two independent clusters required");
        var fit=QuantileRegression.fitExact(response,predictors,quantile);
        if(!fit.converged())throw new IllegalStateException("cluster inference requires a converged exact fit");
        int n=response.length,p=fit.coefficients().length;
        double[] weighted=new double[n*p];
        for(int i=0;i<n;i++)for(int j=0;j<p;j++)weighted[i*p+j]=Math.sqrt(conditionalDensities[i])*predictors[i][j];
        double[] covariance=new double[p*p];
        try(var context=BackendContext.select(BackendPolicy.CPU)) {
            double[] bread=LeastSquaresSolver.solve(weighted,new double[n],n,p,false,context.backend()).unscaledCovariance();
            java.util.Map<Integer,double[]> scores=new java.util.LinkedHashMap<>();double[] residuals=fit.residuals();
            for(int i=0;i<n;i++) {
                double[] score=scores.computeIfAbsent(clusters[i],unused->new double[p]);
                double psi=quantile-(residuals[i]<0?1:0);
                for(int j=0;j<p;j++)score[j]+=predictors[i][j]*psi;
            }
            for(double[] score:scores.values()) {
                double[] influence=new double[p];for(int j=0;j<p;j++)for(int k=0;k<p;k++)influence[j]+=bread[j*p+k]*score[k];
                for(int j=0;j<p;j++)for(int k=0;k<p;k++)covariance[j*p+k]+=influence[j]*influence[k];
            }
        }
        double[] se=new double[p];for(int j=0;j<p;j++)se[j]=Math.sqrt(covariance[j*p+j]);
        return new Result(fit,covariance,se,conditionalDensities,"cluster-CR0-supplied-density",Double.NaN);
    }

    /** Distribution-free iid marginal quantile interval obtained by binomial
     * order-statistic inversion. Valid conservatively at response mass points;
     * this is a marginal quantile method, not a regression covariance. Infinite
     * endpoints explicitly represent insufficient tail information. */
    public static QuantileInterval marginalInterval(double[] response,double quantile,double level) {
        if(response==null||response.length<1||!(quantile>0&&quantile<1)||!(level>0&&level<1))
            throw new IllegalArgumentException("responses, quantile and level required");
        double[] sorted=response.clone();for(double v:sorted)if(!Double.isFinite(v))throw new IllegalArgumentException("finite responses required");
        Arrays.sort(sorted);int n=sorted.length,lo=0,hi=n+1;double tail=(1-level)/2;
        for(int k=1;k<=n;k++)if(jdistlib.Binomial.cumulative(k-1,n,quantile,true,false)<=tail)lo=k;
        for(int k=n;k>=1;k--)if(jdistlib.Binomial.cumulative(k-1,n,quantile,false,false)<=tail)hi=k;
        return new QuantileInterval(sorted[Math.max(0,(int)Math.ceil(n*quantile)-1)],
            lo==0?Double.NEGATIVE_INFINITY:sorted[lo-1],hi==n+1?Double.POSITIVE_INFINITY:sorted[hi-1],level);
    }
    public record QuantileInterval(double estimate,double lower,double upper,double confidenceLevel) { }

    /** Bartlett HAC covariance for ordered, weakly dependent continuous-response
     * quantile scores, using supplied consistent marginal conditional densities.
     * Row order defines time. The lag is caller-selected: for general dependence
     * it must increase slowly with n, with adequate mixing/moment conditions.
     * This is not valid for unit-root designs or response mass points. */
    public static Result fitExactHac(double[] response,double[][] predictors,double quantile,
            double[] conditionalDensities,int maximumLag) {
        validateDensities(conditionalDensities,response==null?0:response.length);
        if(maximumLag<0||maximumLag>=response.length)throw new IllegalArgumentException("HAC lag must be in 0..n-1");
        var fit=QuantileRegression.fitExact(response,predictors,quantile);
        if(!fit.converged())throw new IllegalStateException("HAC inference requires a converged exact fit");
        int n=response.length,p=fit.coefficients().length;double[] weighted=new double[n*p],residuals=fit.residuals();
        for(int i=0;i<n;i++)for(int j=0;j<p;j++)weighted[i*p+j]=Math.sqrt(conditionalDensities[i])*predictors[i][j];
        double[] covariance=new double[p*p];double[][] influences=new double[n][p];
        try(var context=BackendContext.select(BackendPolicy.CPU)) {
            double[] bread=LeastSquaresSolver.solve(weighted,new double[n],n,p,false,context.backend()).unscaledCovariance();
            for(int i=0;i<n;i++)for(int j=0;j<p;j++)for(int k=0;k<p;k++)influences[i][j]+=bread[j*p+k]*predictors[i][k]*(quantile-(residuals[i]<0?1:0));
        }
        double[] mean=new double[p];for(double[] row:influences)for(int j=0;j<p;j++)mean[j]+=row[j]/n;
        for(double[] row:influences)for(int j=0;j<p;j++)row[j]-=mean[j];
        for(int i=0;i<n;i++)for(int j=0;j<p;j++)for(int k=0;k<p;k++)covariance[j*p+k]+=influences[i][j]*influences[i][k];
        for(int lag=1;lag<=maximumLag;lag++) {
            double weight=1-(double)lag/(maximumLag+1);
            for(int i=lag;i<n;i++)for(int j=0;j<p;j++)for(int k=0;k<p;k++)
                covariance[j*p+k]+=weight*(influences[i][j]*influences[i-lag][k]+influences[i-lag][j]*influences[i][k]);
        }
        double[] se=new double[p];for(int j=0;j<p;j++) {
            if(!(covariance[j*p+j]>0)||!Double.isFinite(covariance[j*p+j]))throw new IllegalArgumentException("unresolved HAC covariance");
            se[j]=Math.sqrt(covariance[j*p+j]);
        }
        return new Result(fit,covariance,se,conditionalDensities,"Bartlett-HAC-lag-"+maximumLag,Double.NaN);
    }

    /** Fit and sandwich covariance with externally supplied conditional densities
     * f_i(x_i' beta). These must be known or consistently estimated; they are
     * densities in inverse response units, not weights or sparsities (1/f_i).
     * Covariance is tau(1-tau) A^-1 X'X A^-1, A = X' diag(f_i) X. */
    public static Result fitExact(double[] response,double[][] predictors,double quantile,
            double[] conditionalDensities) {
        return fitExact(response,predictors,quantile,conditionalDensities,
            QuantileLinearProgram.Options.defaults());
    }

    public static Result fitExact(double[] response,double[][] predictors,double quantile,
            double[] conditionalDensities,QuantileLinearProgram.Options options) {
        validateDensities(conditionalDensities,response==null?0:response.length);
        QuantileRegressionResult fit=QuantileRegression.fitExact(response,predictors,quantile,options);
        return inference(fit,predictors,conditionalDensities,"supplied-conditional-density",Double.NaN);
    }

    /** Exact fit with a pooled Gaussian residual density at zero and a caller-set
     * bandwidth in response units. Requires iid errors. Bandwidth is not chosen
     * automatically; asymptotic validity requires h -> 0 and n*h -> infinity,
     * plus the usual density smoothness and quantile regularity conditions.
     * A finite positive bandwidth alone does not establish valid inference. */
    public static Result fitExactIidKernel(double[] response,double[][] predictors,
            double quantile,double bandwidth) {
        if(!(bandwidth>0) || !Double.isFinite(bandwidth))
            throw new IllegalArgumentException("finite positive density bandwidth required");
        QuantileRegressionResult fit=QuantileRegression.fitExact(response,predictors,quantile);
        double[] residuals=fit.residuals();double density=0;
        for(double residual:residuals) {
            double z=residual/bandwidth;
            density+=Math.exp(-.5*z*z)/residuals.length;
        }
        density=density/Math.sqrt(2*Math.PI)/bandwidth;
        double[] densities=new double[residuals.length];Arrays.fill(densities,density);
        return inference(fit,predictors,densities,"iid-gaussian-residual-kernel",bandwidth);
    }

    private static Result inference(QuantileRegressionResult fit,double[][] x,
            double[] densities,String method,double bandwidth) {
        if(!fit.converged())throw new IllegalStateException("quantile inference requires a converged exact fit");
        int n=x.length,p=fit.coefficients().length;
        if(n<=p)throw new IllegalArgumentException("quantile inference requires more rows than coefficients");
        validateDensities(densities,n);
        double densityScale=Arrays.stream(densities).max().orElseThrow();
        double[] scales=new double[p],weighted=new double[n*p];
        for(double[] row:x)for(int j=0;j<p;j++)scales[j]=Math.max(scales[j],Math.abs(row[j]));
        for(int j=0;j<p;j++)if(!(scales[j]>0))throw new IllegalArgumentException("quantile design has a zero column");
        for(int i=0;i<n;i++)for(int j=0;j<p;j++)
            weighted[i*p+j]=Math.sqrt(densities[i]/densityScale)*(x[i][j]/scales[j]);
        double[] covariance=new double[p*p];
        try(var context=BackendContext.select(BackendPolicy.CPU)) {
            double[] bread=LeastSquaresSolver.solve(weighted,new double[n],n,p,false,
                context.backend()).unscaledCovariance();
            double multiplier=Math.sqrt(fit.quantile()*(1-fit.quantile()));
            // Outer products of observation influences preserve positive variances.
            for(double[] row:x) {
                double[] influence=new double[p];
                for(int j=0;j<p;j++) {
                    for(int k=0;k<p;k++)influence[j]+=bread[j*p+k]*(row[k]/scales[k]);
                    influence[j]=influence[j]/densityScale/scales[j]*multiplier;
                }
                for(int j=0;j<p;j++)for(int k=0;k<=j;k++)
                    covariance[j*p+k]+=influence[j]*influence[k];
            }
        }
        double[] se=new double[p];
        for(int j=0;j<p;j++) {
            if(!(covariance[j*p+j]>0) || !Double.isFinite(covariance[j*p+j]))
                throw new IllegalArgumentException("quantile covariance is numerically unresolved");
            se[j]=Math.sqrt(covariance[j*p+j]);
            for(int k=0;k<j;k++)covariance[k*p+j]=covariance[j*p+k];
        }
        return new Result(fit,covariance,se,densities,method,bandwidth);
    }

    private static void validateDensities(double[] densities,int n) {
        if(densities==null || n==0 || densities.length!=n)
            throw new IllegalArgumentException("one conditional density per response required");
        for(double density:densities)if(!(density>0) || !Double.isFinite(density))
            throw new IllegalArgumentException("conditional densities must be finite and positive");
    }

    /** Row-major covariance and normal-asymptotic standard errors in coefficient
     * order. No finite-sample Student-t or automatic bandwidth guarantee is made. */
    public record Result(QuantileRegressionResult fit,double[] coefficientCovariance,
            double[] standardErrors,double[] conditionalDensities,String method,double bandwidth) {
        public Result {
            coefficientCovariance=coefficientCovariance.clone();standardErrors=standardErrors.clone();
            conditionalDensities=conditionalDensities.clone();
        }
        public double[] coefficientCovariance(){return coefficientCovariance.clone();}
        public double[] standardErrors(){return standardErrors.clone();}
        public double[] conditionalDensities(){return conditionalDensities.clone();}
    }
}
