/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;

import java.util.Arrays;
import jdistlib.Normal;
import jdistlib.math.MultivariableFunction;
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
                optimum.converged(),fitted.likelihood.state(),y,x,order,seasonal,
                optimum.parameters());
        }
    }

    private static JointInference jointInference(double[] y,double[][] x,ArimaOrder order,
            SeasonalArimaOrder seasonal,double[] difference,double[] dynamicPoint,
            boolean converged,double[] fittedBeta,double[] conditionalCovariance) {
        int dynamic=dynamicPoint.length,p=fittedBeta.length,k=dynamic+p;
        double[] point=new double[k],scales=new double[p];System.arraycopy(dynamicPoint,0,point,0,dynamic);
        if(!converged)return JointInference.unavailable(point,scales,k);
        for(double value:dynamicPoint)if(Math.abs(value)>=3.8-1e-4)return JointInference.unavailable(point,scales,k);
        for(int j=0;j<p;j++) {
            double variance=conditionalCovariance[j*p+j];
            if(!(variance>0)||!Double.isFinite(variance))return JointInference.unavailable(point,scales,k);
            scales[j]=Math.sqrt(variance);
        }
        MultivariableFunction objective=raw->{
            try {
                double[] beta=beta(raw,dynamic,fittedBeta,scales);
                double[] residual=residual(y,x,beta);
                var coefficients=ArimaMath.decode(Arrays.copyOf(raw,dynamic),order,seasonal,false,false);
                double value=new ArimaStateSpace(coefficients.effectiveAr(),coefficients.effectiveMa(),difference)
                    .filter(residual,0,0,false).negativeLogLikelihood();
                return Double.isFinite(value)?value:Double.NaN;
            } catch(IllegalArgumentException exception){return Double.NaN;}
        };
        try(BackendContext context=BackendContext.select(BackendPolicy.CPU)) {
            double[] rawCovariance=ExactArma.inverseHessian(point,objective,context.backend());
            if(Arrays.stream(rawCovariance).anyMatch(v->!Double.isFinite(v)))
                return JointInference.unavailable(point,scales,k);
            double[] jacobian=new double[k*k];
            for(int j=0;j<dynamic;j++) {
                double step=1e-5*(1+Math.abs(point[j]));double[] plus=point.clone(),minus=point.clone();
                plus[j]+=step;minus[j]-=step;double[] a=reported(plus,dynamic,fittedBeta,scales,order,seasonal);
                double[] b=reported(minus,dynamic,fittedBeta,scales,order,seasonal);
                for(int i=0;i<k;i++)jacobian[i*k+j]=(a[i]-b[i])/(2*step);
            }
            for(int j=0;j<p;j++)jacobian[j*k+dynamic+j]=scales[j];
            double[] reportedCovariance=covarianceTransform(jacobian,rawCovariance,k);
            for(int i=0;i<k;i++)if(!(reportedCovariance[i*k+i]>0)||!Double.isFinite(reportedCovariance[i*k+i]))
                return JointInference.unavailable(point,scales,k);
            return new JointInference(point,rawCovariance,reportedCovariance,scales,true);
        }
    }

    private static double[] reported(double[] raw,int dynamic,double[] fittedBeta,double[] scales,
            ArimaOrder order,SeasonalArimaOrder seasonal) {
        double[] result=new double[raw.length];double[] beta=beta(raw,dynamic,fittedBeta,scales);
        System.arraycopy(beta,0,result,0,beta.length);
        var c=ArimaMath.decode(Arrays.copyOf(raw,dynamic),order,seasonal,false,false);int at=beta.length;
        for(double[] values:new double[][]{c.ar(),c.ma(),c.seasonalAr(),c.seasonalMa()}) {
            System.arraycopy(values,0,result,at,values.length);at+=values.length;
        }
        return result;
    }

    private static double[] beta(double[] raw,int dynamic,double[] fitted,double[] scales) {
        double[] result=fitted.clone();for(int j=0;j<result.length;j++)result[j]+=raw[dynamic+j]*scales[j];return result;
    }
    private static double[] residual(double[] y,double[][] x,double[] beta) {
        double[] result=y.clone();for(int i=0;i<result.length;i++)for(int j=0;j<beta.length;j++)result[i]-=x[i][j]*beta[j];return result;
    }
    private static double[] covarianceTransform(double[] jacobian,double[] covariance,int n) {
        double[] temporary=new double[n*n],result=new double[n*n];
        for(int i=0;i<n;i++)for(int k=0;k<n;k++)for(int j=0;j<n;j++)temporary[i*n+j]+=jacobian[i*n+k]*covariance[k*n+j];
        for(int i=0;i<n;i++)for(int k=0;k<n;k++)for(int j=0;j<n;j++)result[i*n+j]+=temporary[i*n+k]*jacobian[j*n+k];
        return result;
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
    private record JointInference(double[] point,double[] rawCovariance,double[] reportedCovariance,
            double[] betaScales,boolean available) {
        private static JointInference unavailable(double[] point,double[] scales,int size) {
            double[] raw=new double[size*size],reported=new double[size*size];
            Arrays.fill(raw,Double.NaN);Arrays.fill(reported,Double.NaN);
            return new JointInference(point.clone(),raw,reported,scales.clone(),false);
        }
    }
    public static final class Result {
        private final double[] beta,covariance,ar,ma,difference,residuals;
        private final double variance,logLikelihood;private final boolean converged;
        private final ArimaStateSpace.ForecastState state;
        private final double[] response;private final double[][] design;
        private final ArimaOrder order;private final SeasonalArimaOrder seasonal;
        private final double[] dynamicPoint;
        private volatile JointInference inference;
        private Result(double[] beta,double[] covariance,double[] ar,double[] ma,double[] difference,
                double[] residuals,double variance,double logLikelihood,boolean converged,ArimaStateSpace.ForecastState state,
                double[] response,double[][] design,ArimaOrder order,SeasonalArimaOrder seasonal,
                double[] dynamicPoint) {
            this.beta=beta.clone();this.covariance=covariance.clone();this.ar=ar.clone();this.ma=ma.clone();
            this.difference=difference.clone();this.residuals=residuals.clone();this.variance=variance;
            this.logLikelihood=logLikelihood;this.converged=converged;this.state=state;
            this.response=response.clone();this.design=Arrays.stream(design).map(double[]::clone).toArray(double[][]::new);
            this.order=order;this.seasonal=seasonal;this.dynamicPoint=dynamicPoint.clone();
        }
        public double[] coefficients(){return beta.clone();}
        public double[] effectiveAutoregressive(){return ar.clone();}
        public double[] effectiveMovingAverage(){return ma.clone();}
        /** GLS covariance conditional on fitted dynamics; excludes dynamic estimation uncertainty. */
        public double[] conditionalCoefficientCovariance(){return covariance.clone();}
        /** Joint covariance ordered as regression coefficients, AR, MA, seasonal AR, seasonal MA.
         * Computed and cached on first request. Innovation-variance uncertainty
         * is not included. */
        public double[] jointParameterCovariance(){return jointInference().reportedCovariance.clone();}
        /** Estimates matching {@link #jointParameterCovariance()}: regression
         * coefficients, AR, MA, seasonal AR, then seasonal MA parameters. */
        public double[] jointParameterEstimates() {
            int dynamic=dynamicPoint.length;
            double[] point=new double[dynamic+beta.length];
            System.arraycopy(dynamicPoint,0,point,0,dynamic);
            return reported(point,dynamic,beta,new double[beta.length],order,seasonal);
        }
        public boolean jointParameterInferenceAvailable(){return jointInference().available;}
        public double innovationVariance(){return variance;}
        public double logLikelihood(){return logLikelihood;}
        public boolean converged(){return converged;}
        /** Historical ARIMA error states conditional on fitted regression,
         * dynamics, and innovation variance; excludes parameter uncertainty. */
        public ArimaSmoothing.Result smoothErrors(){return ArimaSmoothing.smooth(residuals,ar,ma,difference,variance);}
        /** Process-noise forecast conditional on fitted regression and dynamics. */
        public ArimaForecast forecast(double[][] futureDesign,double level) {
            return forecastConditional(futureDesign,level);
        }
        public ArimaForecast forecastConditional(double[][] futureDesign,double level) {
            MatrixOps.rowMajor(futureDesign,futureDesign==null?0:futureDesign.length);
            if(futureDesign[0].length!=beta.length)throw new IllegalArgumentException("future design columns differ");
            ArimaForecast f=state.forecast(futureDesign.length,level,variance);
            double[] mean=f.means(),lo=f.lowerBounds(),hi=f.upperBounds();
            for(int t=0;t<mean.length;t++)for(int j=0;j<beta.length;j++) {
                double value=futureDesign[t][j]*beta[j];mean[t]+=value;lo[t]+=value;hi[t]+=value;
            }
            return new ArimaForecast(mean,f.standardErrors(),lo,hi,level);
        }

        /** Forecast including delta-method uncertainty from jointly estimated regression and dynamics. */
        public ArimaForecast forecastWithParameterUncertainty(double[][] futureDesign,double level) {
            ArimaForecast conditional=forecastConditional(futureDesign,level);
            JointInference joint=jointInference();
            if(!joint.available)throw new IllegalStateException("joint ARIMA regression parameter inference is unavailable");
            int horizon=futureDesign.length,k=joint.point.length;double[][] gradient=new double[horizon][k];
            for(int j=0;j<k;j++) {
                double step=1e-5*(1+Math.abs(joint.point[j]));
                double[] plus=joint.point.clone(),minus=joint.point.clone();plus[j]+=step;minus[j]-=step;
                double[] a=forecastMeans(plus,futureDesign,joint),b=forecastMeans(minus,futureDesign,joint);
                for(int h=0;h<horizon;h++)gradient[h][j]=(a[h]-b[h])/(2*step);
            }
            double[] means=conditional.means(),standardErrors=conditional.standardErrors();
            double[] lower=new double[horizon],upper=new double[horizon];
            double critical=Normal.quantile(.5+level/2,0,1,true,false);
            for(int h=0;h<horizon;h++) {
                double parameterVariance=0;
                for(int i=0;i<k;i++)for(int j=0;j<k;j++)parameterVariance+=gradient[h][i]*joint.rawCovariance[i*k+j]*gradient[h][j];
                if(parameterVariance < -1e-10*(1+standardErrors[h]*standardErrors[h]))
                    throw new IllegalStateException("forecast parameter variance lost precision");
                standardErrors[h]=Math.sqrt(standardErrors[h]*standardErrors[h]+Math.max(0,parameterVariance));
                lower[h]=means[h]-critical*standardErrors[h];upper[h]=means[h]+critical*standardErrors[h];
            }
            return new ArimaForecast(means,standardErrors,lower,upper,level);
        }

        private double[] forecastMeans(double[] point,double[][] futureDesign,
                JointInference joint) {
            int dynamic=point.length-beta.length;double[] candidateBeta=beta(point,dynamic,beta,joint.betaScales);
            var c=ArimaMath.decode(Arrays.copyOf(point,dynamic),order,seasonal,false,false);
            ArimaStateSpace model=new ArimaStateSpace(c.effectiveAr(),c.effectiveMa(),difference);
            var evaluation=model.filter(residual(response,design,candidateBeta),0,0,true);
            if(!Double.isFinite(evaluation.negativeLogLikelihood()))throw new IllegalStateException("forecast perturbation is unidentified");
            double[] means=evaluation.state().forecast(futureDesign.length,.95,variance).means();
            for(int h=0;h<means.length;h++)for(int j=0;j<candidateBeta.length;j++)means[h]+=futureDesign[h][j]*candidateBeta[j];
            return means;
        }

        private JointInference jointInference() {
            JointInference current=inference;
            if(current!=null)return current;
            synchronized(this) {
                current=inference;
                if(current==null) {
                    current=ArimaRegression.jointInference(response,design,order,
                        seasonal,difference,dynamicPoint,converged,beta,covariance);
                    inference=current;
                }
            }
            return current;
        }

        boolean jointParameterInferenceComputed(){return inference!=null;}
    }
}
