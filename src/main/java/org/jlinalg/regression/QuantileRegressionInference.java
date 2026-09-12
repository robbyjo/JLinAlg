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
