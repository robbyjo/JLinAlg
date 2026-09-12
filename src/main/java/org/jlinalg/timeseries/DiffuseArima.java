/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;

import java.util.Arrays;
import jdistlib.math.MultivariableFunction;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;

/** Exact diffuse Gaussian likelihood for integrated (including seasonal) ARIMA. */
public final class DiffuseArima {
    private DiffuseArima() { }

    /**
     * Fits the original observation sequence, marginalizing NaN observations.
     * The stationary ARMA block has its unconditional covariance; integrated
     * states have symbolic infinite variance. As in R's arima ML convention,
     * diffuse observations contribute neither a quadratic nor a determinant.
     */
    public static Result fit(double[] series, ArimaOrder order, ArimaOptions options) {
        ArimaOptions controls = options == null ? ArimaOptions.defaults() : options;
        SeasonalArimaOrder seasonal = controls.seasonalOrder();
        if (series == null || order == null
                || order.differences() + seasonal.differences() < 1)
            throw new IllegalArgumentException("an integrated ARIMA order is required");
        double[] values = series.clone();
        int observed = 0;
        for (double value : values) {
            if (Double.isInfinite(value)) throw new IllegalArgumentException("infinite observations are invalid; use NaN for missing");
            if (!Double.isNaN(value)) observed++;
        }
        if (controls.includeDrift() && order.differences() + seasonal.differences() != 1)
            throw new IllegalArgumentException("drift requires exactly one differencing operator");
        double[] difference = ArimaMath.differencingPolynomial(order, seasonal);
        int diffuse = difference.length - 1;
        int dynamic = order.autoregressive() + order.movingAverage()
            + seasonal.autoregressive() + seasonal.movingAverage();
        int parameters = dynamic + (controls.includeDrift() ? 1 : 0);
        if (values.length <= diffuse || observed - diffuse <= parameters + 1)
            throw new IllegalArgumentException("too few observed values for the requested diffuse ARIMA");
        double[] initial = new double[dynamic], lower = new double[dynamic], upper = new double[dynamic];
        Arrays.fill(lower, -3.8); Arrays.fill(upper, 3.8);
        // location is the mean of the differenced process. A seasonal difference
        // of period s maps a linear trend with slope location/s to that mean.
        double driftDivisor = seasonal.differences() > 0 ? seasonal.period() : 1;
        MultivariableFunction objective = point -> {
            try {
                ArimaMath.Coefficients c = ArimaMath.decode(point, order, seasonal, false, controls.includeDrift());
                ArimaStateSpace model = new ArimaStateSpace(c.effectiveAr(), c.effectiveMa(), difference);
                double nll = (controls.includeDrift() ? model.profileSlope(values, false)
                    : model.filter(values, 0, 0, false)).negativeLogLikelihood();
                return Double.isFinite(nll) ? nll : Double.MAX_VALUE / 8;
            } catch (IllegalArgumentException exception) { return Double.MAX_VALUE / 8; }
        };
        BoundedOptimizer.Result optimized = Arima.multiStart(initial, lower, upper, objective, dynamic,
            controls.maximumFunctionEvaluations(), controls.optimizationTolerance(), controls.optimizationStarts());
        ArimaMath.Coefficients c = ArimaMath.decode(optimized.parameters(), order, seasonal,
            false, controls.includeDrift());
        ArimaStateSpace model = new ArimaStateSpace(c.effectiveAr(), c.effectiveMa(), difference);
        ArimaStateSpace.Evaluation likelihood = controls.includeDrift() ? model.profileSlope(values, true)
            : model.filter(values, 0, 0, true);
        if (!Double.isFinite(likelihood.negativeLogLikelihood()))
            throw new IllegalArgumentException("diffuse states are unidentified or the innovation variance is degenerate");
        int count = likelihood.observations(), k = parameters + 1;
        double ll = -likelihood.negativeLogLikelihood(), aic = -2 * ll + 2 * k;
        double aicc = count > k + 1 ? aic + 2.0 * k * (k + 1) / (count - k - 1) : Double.POSITIVE_INFINITY;
        ArimaResult result = new ArimaResult(order, seasonal, c.ar(), c.ma(), c.seasonalAr(), c.seasonalMa(),
            likelihood.fittedSlope() * driftDivisor, c.drift(), likelihood.variance(), likelihood.innovations(),
            ArimaMath.difference(values, order, seasonal), values, ll, aic, aicc, -2 * ll + Math.log(count) * k,
            count, optimized.evaluations(), optimized.converged(), optimized.converged()
                ? "exact diffuse likelihood optimized" : "optimizer stopped before convergence",
            likelihood.state());
        double[] covariance=coefficientCovariance(values,order,seasonal,difference,
            controls.includeDrift(),driftDivisor,optimized,likelihood);
        return new Result(result, diffuse, count, true, covariance);
    }

    private static double[] coefficientCovariance(double[] values,ArimaOrder order,
            SeasonalArimaOrder seasonal,double[] difference,boolean drift,double divisor,
            BoundedOptimizer.Result optimized,ArimaStateSpace.Evaluation likelihood) {
        int dynamic=optimized.parameters().length,k=dynamic+(drift?1:0);
        double[] unavailable=new double[k*k];Arrays.fill(unavailable,Double.NaN);
        if(!optimized.converged())return unavailable;
        for(double parameter:optimized.parameters())if(Math.abs(parameter)>=3.8-1e-4)return unavailable;
        double[] point=Arrays.copyOf(optimized.parameters(),k),centered=values.clone();
        double scale=Math.sqrt(likelihood.variance());
        // Fit profiles drift for speed. Inference uses the joint drift/dynamic
        // information, profiling only innovation variance, retaining cross terms.
        // A centered, noise-scaled drift coordinate makes finite differences
        // insensitive to the units and magnitude of the fitted drift.
        if(drift)for(int t=0;t<centered.length;t++)centered[t]-=likelihood.fittedSlope()*(t+1.0);
        MultivariableFunction joint=x->{
            try {
                ArimaMath.Coefficients c=ArimaMath.decode(x,order,seasonal,false,false);
                double nll=new ArimaStateSpace(c.effectiveAr(),c.effectiveMa(),difference)
                    .filter(centered,0,drift?x[dynamic]*scale/divisor:0,false).negativeLogLikelihood();
                return Double.isFinite(nll)?nll:Double.NaN;
            } catch(IllegalArgumentException exception){return Double.NaN;}
        };
        try(var context=BackendContext.select(BackendPolicy.CPU)) {
            double[] raw=ExactArma.inverseHessian(point,joint,context.backend());
            double[] covariance=ExactArma.deltaCovariance(point,raw,order,seasonal,drift,context.backend());
            if(drift)for(int i=0;i<k;i++) {
                covariance[i*k+dynamic]*=scale;covariance[dynamic*k+i]*=scale;
            }
            for(int i=0;i<k;i++)if(!(covariance[i*k+i]>0))return unavailable;
            for(double value:covariance)if(!Double.isFinite(value))return unavailable;
            return covariance;
        }
    }

    /** Covariance order: AR, MA, seasonal AR, seasonal MA, then drift location
     * when fitted. Drift location is the mean of the differenced process;
     * for seasonal differencing it is period times the slope per time step.
     * Innovation variance is profiled out. These are observed-information,
     * asymptotic Gaussian-likelihood errors, not robust errors or forecast bands. */
    public record Result(ArimaResult fit, int diffuseStateCount, int likelihoodObservations,
                         boolean diffuseLikelihood,double[] coefficientCovariance) {
        public Result { coefficientCovariance=coefficientCovariance.clone(); }
        /** Compatibility constructor for results without computed information. */
        public Result(ArimaResult fit,int diffuseStateCount,int likelihoodObservations,boolean diffuseLikelihood) {
            this(fit,diffuseStateCount,likelihoodObservations,diffuseLikelihood,unavailable(fit));
        }
        private static double[] unavailable(ArimaResult fit) {
            int k=fit.autoregressive().length+fit.movingAverage().length
                +fit.seasonalAutoregressive().length+fit.seasonalMovingAverage().length+(fit.drift()?1:0);
            double[] result=new double[k*k];Arrays.fill(result,Double.NaN);return result;
        }
        public double[] coefficientCovariance(){return coefficientCovariance.clone();}
        public double[] standardErrors() {
            int k=(int)Math.sqrt(coefficientCovariance.length);double[] result=new double[k];
            for(int i=0;i<k;i++)result[i]=Math.sqrt(coefficientCovariance[i*k+i]);return result;
        }
        /** False for nonconvergence, a transform-bound solution, or information
         * that is singular or unresolved across two differentiation step sizes. */
        public boolean coefficientInferenceAvailable() {
            return fit.converged() && Arrays.stream(coefficientCovariance).allMatch(Double::isFinite);
        }
    }
}
