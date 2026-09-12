/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;

import java.util.Arrays;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.LeastSquaresSolver;
import org.jlinalg.internal.MatrixOps;

/** Gaussian regression with exact, possibly seasonal and integrated ARIMA errors.
 * Design columns explicitly include any intercept or trend. Missing responses
 * are marginalized; every design row, including future rows, must be finite.
 * ArimaOptions includeMean/includeDrift are not used: the design defines these. */
public final class ArimaRegression {
    private ArimaRegression() { }
    public static Result fit(double[] response,double[][] design,ArimaOrder order,ArimaOptions options) {
        if(response==null||order==null)throw new IllegalArgumentException("response and order required");
        MatrixOps.rowMajor(design,response.length);
        double[] y=response.clone();double[][] x=Arrays.stream(design).map(double[]::clone).toArray(double[][]::new);
        for(double v:y)if(Double.isInfinite(v))throw new IllegalArgumentException("use NaN for missing responses");
        ArimaOptions controls=options==null?ArimaOptions.defaults():options;
        var seasonal=controls.seasonalOrder();
        int dynamic=order.autoregressive()+order.movingAverage()+seasonal.autoregressive()+seasonal.movingAverage();
        double[] difference=ArimaMath.differencingPolynomial(order,seasonal);
        double[] initial=new double[dynamic],lower=new double[dynamic],upper=new double[dynamic];
        Arrays.fill(lower,-3.8);Arrays.fill(upper,3.8);
        try(var context=BackendContext.select(BackendPolicy.CPU)) {
            jdistlib.math.MultivariableFunction objective=point->{
                try{return profile(y,x,state(point,order,seasonal,difference),context.backend()).likelihood.negativeLogLikelihood();}
                catch(IllegalArgumentException e){return Double.MAX_VALUE/8;}
            };
            var optimum=Arima.multiStart(initial,lower,upper,objective,dynamic,
                controls.maximumFunctionEvaluations(),controls.optimizationTolerance(),controls.optimizationStarts());
            var coefficients=ArimaMath.decode(optimum.parameters(),order,seasonal,false,false);
            var model=state(optimum.parameters(),order,seasonal,difference);
            Profile fitted=profile(y,x,model,context.backend());
            if(!Double.isFinite(fitted.likelihood.negativeLogLikelihood()))throw new IllegalArgumentException("unidentified or degenerate regression ARIMA");
            double[] covariance=fitted.covariance.clone();
            for(int i=0;i<covariance.length;i++)covariance[i]=optimum.converged()?covariance[i]*fitted.likelihood.variance():Double.NaN;
            return new Result(fitted.beta,covariance,coefficients.effectiveAr(),coefficients.effectiveMa(),
                difference,fitted.residuals,fitted.likelihood.variance(),-fitted.likelihood.negativeLogLikelihood(),
                optimum.converged(),fitted.likelihood.state());
        }
    }
    private static ArimaStateSpace state(double[] point,ArimaOrder order,SeasonalArimaOrder seasonal,double[] difference) {
        var c=ArimaMath.decode(point,order,seasonal,false,false);return new ArimaStateSpace(c.effectiveAr(),c.effectiveMa(),difference);
    }
    private static Profile profile(double[] y,double[][] x,ArimaStateSpace model,jdistlib.accelerator.ComputeBackend backend) {
        var base=model.filter(y,0,0,true);int count=base.observations(),p=x[0].length;
        if(count<=p)throw new IllegalArgumentException("too few finite innovations for design");
        double[] transformed=new double[count*p],target=new double[count];
        for(int t=0,i=0;t<y.length;t++)if(Double.isFinite(base.innovations()[t]))target[i++]=base.innovations()[t];
        for(int j=0;j<p;j++) {
            double[] column=new double[y.length];for(int t=0;t<y.length;t++)column[t]=Double.isNaN(y[t])?Double.NaN:x[t][j];
            double[] innovation=model.filter(column,0,0,true).innovations();
            for(int t=0,i=0;t<y.length;t++)if(Double.isFinite(base.innovations()[t]))transformed[i++*p+j]=innovation[t];
        }
        var solution=LeastSquaresSolver.solve(transformed,target,count,p,false,backend);
        double[] beta=solution.coefficients(),residuals=y.clone();
        for(int t=0;t<y.length;t++)for(int j=0;j<p;j++)residuals[t]-=x[t][j]*beta[j];
        return new Profile(beta,solution.unscaledCovariance(),residuals,model.filter(residuals,0,0,true));
    }
    private record Profile(double[] beta,double[] covariance,double[] residuals,ArimaStateSpace.Evaluation likelihood) { }
    public static final class Result {
        private final double[] beta,covariance,ar,ma,difference,residuals;
        private final double variance,logLikelihood;private final boolean converged;
        private final ArimaStateSpace.ForecastState state;
        private Result(double[] beta,double[] covariance,double[] ar,double[] ma,double[] difference,
                double[] residuals,double variance,double logLikelihood,boolean converged,ArimaStateSpace.ForecastState state) {
            this.beta=beta.clone();this.covariance=covariance.clone();this.ar=ar.clone();this.ma=ma.clone();
            this.difference=difference.clone();this.residuals=residuals.clone();this.variance=variance;
            this.logLikelihood=logLikelihood;this.converged=converged;this.state=state;
        }
        public double[] coefficients(){return beta.clone();}
        public double[] effectiveAutoregressive(){return ar.clone();}
        public double[] effectiveMovingAverage(){return ma.clone();}
        /** GLS covariance conditional on fitted dynamics; excludes dynamic estimation uncertainty. */
        public double[] conditionalCoefficientCovariance(){return covariance.clone();}
        public double innovationVariance(){return variance;}
        public double logLikelihood(){return logLikelihood;}
        public boolean converged(){return converged;}
        /** Historical ARIMA error states, conditional on estimated regression coefficients. */
        public ArimaSmoothing.Result smoothErrors(){return ArimaSmoothing.smooth(residuals,ar,ma,difference,variance);}
        public ArimaForecast forecast(double[][] futureDesign,double level) {
            MatrixOps.rowMajor(futureDesign,futureDesign==null?0:futureDesign.length);
            if(futureDesign[0].length!=beta.length)throw new IllegalArgumentException("future design columns differ");
            ArimaForecast f=state.forecast(futureDesign.length,level,variance);
            double[] mean=f.means(),lo=f.lowerBounds(),hi=f.upperBounds();
            for(int t=0;t<mean.length;t++)for(int j=0;j<beta.length;j++) {
                double value=futureDesign[t][j]*beta[j];mean[t]+=value;lo[t]+=value;hi[t]+=value;
            }
            return new ArimaForecast(mean,f.standardErrors(),lo,hi,level);
        }
    }
}
