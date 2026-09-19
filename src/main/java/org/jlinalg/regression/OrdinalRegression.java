/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

import java.util.Arrays;
import jdistlib.Normal;
import org.jlinalg.compute.*;
import org.jlinalg.internal.LeastSquaresSolver;

/** Cumulative-link ordinal ML: {@code P(Y <= k | x) = F(threshold[k] - x beta)}.
 * Predictors exclude an intercept; thresholds supply the location parameter. */
public final class OrdinalRegression {
    private OrdinalRegression() { }
    public enum Link { LOGIT, PROBIT }
    public static Result fit(int[] response,double[][] predictors,int categories,Link link) {
        if(response==null||predictors==null||response.length!=predictors.length||response.length==0||predictors[0]==null||categories<2||link==null)
            throw new IllegalArgumentException("matching ordinal data and link required");
        int n=response.length,p=predictors[0].length,m=categories-1,k=p+m;
        if(n<=k)throw new IllegalArgumentException("more observations than parameters required");
        int[] counts=new int[categories];double[] x=new double[n*p],augmented=new double[n*(p+1)];
        for(int i=0;i<n;i++) {
            if(response[i]<0||response[i]>=categories||predictors[i]==null||predictors[i].length!=p)
                throw new IllegalArgumentException("categories must be integers from zero to K-1 and designs rectangular");
            counts[response[i]]++;augmented[i*(p+1)]=1;
            for(int j=0;j<p;j++){double value=predictors[i][j];if(!Double.isFinite(value))throw new IllegalArgumentException("finite predictors required");x[i*p+j]=value;augmented[i*(p+1)+j+1]=value;}
        }
        for(int count:counts)if(count==0)throw new IllegalArgumentException("all prespecified response categories need observations");
        try(var context=BackendContext.select(BackendPolicy.CPU)){
            LeastSquaresSolver.solve(augmented,new double[n],n,p+1,false,context.backend());
        }
        double[] scales=RegressionOptimizer.scaleColumns(x,n,p),base=new double[k];int cumulative=0;double previous=0;
        for(int j=0;j<m;j++) {
            cumulative+=counts[j];double prob=(double)cumulative/n;
            double threshold=link==Link.PROBIT?Normal.quantile(prob,0,1,true,false):Math.log(prob)-Math.log1p(-prob);
            base[p+j]=j==0?threshold:Math.log(threshold-previous);previous=threshold;
        }
        RegressionOptimizer.Objective objective=delta->{double[] theta=base.clone();for(int j=0;j<k;j++)theta[j]+=delta[j];return evaluate(response,x,n,p,m,link,theta);};
        var optimized=RegressionOptimizer.minimize(objective,k,2000,2e-8,1);
        if(!optimized.converged())throw new IllegalArgumentException("ordinal likelihood failed score convergence (possible separation)");
        double[] theta=base.clone();for(int j=0;j<k;j++)theta[j]+=optimized.point()[j];
        double[] threshold=thresholds(theta,p,m),cov=LikelihoodInference.covariance(objective,optimized.point(),n),jacobian=new double[k*k];
        for(int j=0;j<p;j++)jacobian[j*k+j]=1/scales[j];
        for(int j=0;j<m;j++){jacobian[(p+j)*k+p]=1;for(int l=1;l<=j;l++)jacobian[(p+j)*k+p+l]=Math.exp(theta[p+l]);}
        double[] transformed=new double[k*k],temporary=new double[k*k];
        for(int i=0;i<k;i++)for(int j=0;j<k;j++)for(int l=0;l<k;l++)temporary[i*k+j]+=jacobian[i*k+l]*cov[l*k+j];
        for(int i=0;i<k;i++)for(int j=0;j<k;j++)for(int l=0;l<k;l++)transformed[i*k+j]+=temporary[i*k+l]*jacobian[j*k+l];
        double[] beta=Arrays.copyOf(theta,p);for(int j=0;j<p;j++)beta[j]/=scales[j];
        return new Result(beta,threshold,transformed,-objective.evaluate(optimized.point()).value()*n,link,optimized.iterations());
    }
    private static RegressionOptimizer.Evaluation evaluate(int[] y,double[] x,int n,int p,int m,Link link,double[] theta) {
        double[] cuts=thresholds(theta,p,m),gradient=new double[p+m],cutGradient=new double[m];double value=0;
        for(int i=0;i<n;i++) {
            double eta=0;for(int j=0;j<p;j++)eta+=x[i*p+j]*theta[j];int category=y[i];
            double lo=category==0?Double.NEGATIVE_INFINITY:cuts[category-1]-eta,
                hi=category==m?Double.POSITIVE_INFINITY:cuts[category]-eta,logP=logProbability(lo,hi,link);
            double lower=category==0?0:Math.exp(logDensity(lo,link)-logP),upper=category==m?0:Math.exp(logDensity(hi,link)-logP);
            value-=logP/n;for(int j=0;j<p;j++)gradient[j]+=(upper-lower)*x[i*p+j]/n;
            if(category>0)cutGradient[category-1]+=lower/n;if(category<m)cutGradient[category]-=upper/n;
        }
        double sum=0;for(int j=m-1;j>=0;j--){sum+=cutGradient[j];gradient[p+j]=sum*(j==0?1:Math.exp(theta[p+j]));}
        return new RegressionOptimizer.Evaluation(value,gradient);
    }
    private static double[] thresholds(double[] theta,int p,int m){double[] cuts=new double[m];cuts[0]=theta[p];for(int j=1;j<m;j++)cuts[j]=cuts[j-1]+Math.exp(theta[p+j]);return cuts;}
    private static double logDensity(double x,Link link){return link==Link.PROBIT?Normal.density(x,0,1,true):-Math.abs(x)-2*Math.log1p(Math.exp(-Math.abs(x)));}
    private static double logCdf(double x,Link link){return link==Link.PROBIT?Normal.cumulative(x,0,1,true,true):x<0?x-Math.log1p(Math.exp(x)):-Math.log1p(Math.exp(-x));}
    private static double logProbability(double lo,double hi,Link link){
        if(lo==Double.NEGATIVE_INFINITY)return logCdf(hi,link);
        if(hi==Double.POSITIVE_INFINITY)return logCdf(-lo,link);
        double a=lo>0?logCdf(-lo,link):logCdf(hi,link),b=lo>0?logCdf(-hi,link):logCdf(lo,link);
        return a+Math.log(-Math.expm1(b-a));
    }
    /** Covariance order is slopes followed by thresholds in their original coordinates. */
    public record Result(double[] coefficients,double[] thresholds,double[] covariance,double logLikelihood,Link link,int iterations) {
        public Result {coefficients=coefficients.clone();thresholds=thresholds.clone();covariance=covariance.clone();}
        @Override public double[] coefficients(){return coefficients.clone();}
        @Override public double[] thresholds(){return thresholds.clone();}
        @Override public double[] covariance(){return covariance.clone();}
        public double[] probabilities(double[] row){
            if(row==null||row.length!=coefficients.length)throw new IllegalArgumentException("prediction row width differs");
            double eta=0;for(int j=0;j<row.length;j++){if(!Double.isFinite(row[j]))throw new IllegalArgumentException("finite prediction rows required");eta+=row[j]*coefficients[j];}
            if(!Double.isFinite(eta))throw new IllegalArgumentException("prediction exceeds numerical range");
            double[] result=new double[thresholds.length+1];
            for(int j=0;j<result.length;j++)result[j]=Math.exp(logProbability(j==0?Double.NEGATIVE_INFINITY:thresholds[j-1]-eta,
                j==thresholds.length?Double.POSITIVE_INFINITY:thresholds[j]-eta,link));return result;
        }
    }
}
