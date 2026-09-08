/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.penalized;

import java.util.ArrayList;
import java.util.List;
import jdistlib.Normal;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;

/** Exact Gaussian inference conditional on the active set and signs at a fixed
 * LASSO/elastic-net penalty. Noise SD must be known independently; selecting
 * lambda using these responses invalidates this conditioning contract.
 * Coefficient targets are selected-model least-squares projections, not the
 * shrunken penalized estimates. Predictors are not standardized.
 * Based on Lee, Sun, Sun and Taylor (2016), doi:10.1214/15-AOS1371;
 * elastic-net constraints follow by retaining the ridge term in the KKT system. */
public final class PolyhedralSelectiveInference {
    private PolyhedralSelectiveInference() { }

    public static Result fit(double[] response,double[][] predictors,double lambda,double alpha,
                              double noiseStandardDeviation,boolean includeIntercept,
                              double confidenceLevel,BackendPolicy policy) {
        if(response==null||predictors==null||response.length<3||predictors.length!=response.length
                ||predictors[0]==null||!(lambda>0)||!Double.isFinite(lambda)||!(alpha>0&&alpha<=1)
                ||!(noiseStandardDeviation>0)||!Double.isFinite(noiseStandardDeviation)
                ||!(confidenceLevel>0&&confidenceLevel<1)||policy==null)throw new IllegalArgumentException("invalid selective inference inputs");
        int n=response.length,p=predictors[0].length;
        double[][] x=new double[n][p];double[] y=response.clone(),means=new double[p];double meanY=0;
        for(int i=0;i<n;i++) {
            if(predictors[i]==null||predictors[i].length!=p||!Double.isFinite(y[i]))throw new IllegalArgumentException("finite rectangular data required");
            meanY+=y[i]/n;
            for(int j=0;j<p;j++){double v=predictors[i][j];if(!Double.isFinite(v))throw new IllegalArgumentException("finite predictors required");means[j]+=v/n;x[i][j]=v;}
        }
        if(includeIntercept)for(int i=0;i<n;i++){y[i]-=meanY;for(int j=0;j<p;j++)x[i][j]-=means[j];}
        var selection=PenalizedRegression.fit(response,predictors,lambda,ElasticNetOptions.builder().alpha(alpha)
            .fitIntercept(includeIntercept).standardize(false).relativeTolerance(1e-12).maximumIterations(200000).build());
        if(!selection.converged())throw new IllegalStateException("selection optimization did not converge");
        double[] selected=selection.coefficients();List<Integer> active=new ArrayList<>();
        for(int j=0;j<p;j++)if(selected[j]!=0)active.add(j);
        int k=active.size();if(k==0)return new Result(selection,List.of(),confidenceLevel);
        if(k>=n-(includeIntercept?1:0))throw new IllegalArgumentException("selected model needs full column rank and residual observations");
        double[] gram=new double[k*k],cross=new double[k*n],signs=new double[k];
        for(int a=0;a<k;a++) {
            signs[a]=Math.signum(selected[active.get(a)]);
            for(int i=0;i<n;i++)cross[a*n+i]=x[i][active.get(a)];
            for(int b=0;b<k;b++)for(int i=0;i<n;i++)gram[a*k+b]+=x[i][active.get(a)]*x[i][active.get(b)];
        }
        double[] ridgeGram=gram.clone();for(int a=0;a<k;a++)ridgeGram[a*k+a]+=n*lambda*(1-alpha);
        double[] w,c,etaRows;
        try(var context=BackendContext.select(policy)) {
            var ridge=context.backend().dpotrf(ridgeGram,k);
            w=ridge.solve(cross,n);c=ridge.solve(signs);
            etaRows=context.backend().dpotrf(gram,k).solve(cross,n);
        }
        for(int a=0;a<k;a++)c[a]*=n*lambda*alpha;
        List<double[]> constraints=new ArrayList<>();List<Double> bounds=new ArrayList<>();
        for(int a=0;a<k;a++) {
            double[] row=new double[n];for(int i=0;i<n;i++)row[i]=-signs[a]*w[a*n+i];
            constraints.add(row);bounds.add(-signs[a]*c[a]);
        }
        for(int j=0;j<p;j++)if(selected[j]==0) {
            double[] product=new double[k];for(int a=0;a<k;a++)for(int i=0;i<n;i++)product[a]+=x[i][j]*x[i][active.get(a)];
            double shift=0;for(int a=0;a<k;a++)shift+=product[a]*c[a];
            double[] plus=new double[n],minus=new double[n];
            for(int i=0;i<n;i++){plus[i]=x[i][j];for(int a=0;a<k;a++)plus[i]-=product[a]*w[a*n+i];minus[i]=-plus[i];}
            constraints.add(plus);bounds.add(n*lambda*alpha-shift);
            constraints.add(minus);bounds.add(n*lambda*alpha+shift);
        }
        // Inequalities may be multiplied by any positive constant. Normalize
        // them before checking feasibility; never use an absolute unit floor.
        for(int r=0;r<constraints.size();r++) {
            double[] row=constraints.get(r);double bound=bounds.get(r),scale=Math.abs(bound);
            for(double value:row)scale=Math.max(scale,Math.abs(value));
            if(!Double.isFinite(scale))throw new IllegalStateException("nonfinite KKT constraint");
            if(scale>0) {for(int i=0;i<n;i++)row[i]/=scale;bound/=scale;bounds.set(r,bound);}
            double magnitude=Math.abs(bound)+absoluteProducts(row,y);
            if(dot(row,y)-bound>1e-9*magnitude+32*Math.ulp(magnitude))
                throw new IllegalStateException("selection KKT constraints not satisfied to relative tolerance");
        }
        List<Effect> effects=new ArrayList<>();
        for(int a=0;a<k;a++) {
            double[] eta=new double[n];System.arraycopy(etaRows,a*n,eta,0,n);
            double variance=dot(eta,eta),estimate=dot(eta,y),sd=noiseStandardDeviation*Math.sqrt(variance);
            if(!(variance>0)||!Double.isFinite(variance)||!Double.isFinite(estimate)
                    ||!(sd>0)||!Double.isFinite(sd))throw new IllegalStateException("nonfinite selective target or variance");
            double lower=Double.NEGATIVE_INFINITY,upper=Double.POSITIVE_INFINITY;
            for(int r=0;r<constraints.size();r++) {
                double[] row=constraints.get(r);double projection=dot(row,eta);
                // Detect cancellation at dot-product precision, not by the
                // units of the coefficient or the size of the ridge penalty.
                if(Math.abs(projection)<=64*Math.ulp(1.0)*absoluteProducts(row,eta))continue;
                double slope=projection/variance;
                double remainder=bounds.get(r)-dot(row,y)+slope*estimate;
                if(slope==0||!Double.isFinite(slope)||!Double.isFinite(remainder))
                    throw new IllegalStateException("selection boundary exceeds numerical range");
                if(slope>0)upper=Math.min(upper,remainder/slope);
                else lower=Math.max(lower,remainder/slope);
            }
            if(!(lower<estimate&&estimate<upper))throw new IllegalStateException("estimate is on a selection boundary");
            double[] tails=logTruncatedTails(estimate,0,sd,lower,upper);
            double pValue=Math.min(1,2*Math.exp(Math.min(tails[0],tails[1])));
            double tail=(1-confidenceLevel)/2;
            double logOdds=Math.log1p(-tail)-Math.log(tail);
            double ciLower=invert(estimate,sd,lower,upper,logOdds),ciUpper=invert(estimate,sd,lower,upper,-logOdds);
            if(!(ciLower<ciUpper))throw new IllegalStateException("selective confidence interval is not numerically resolved");
            effects.add(new Effect(active.get(a),estimate,sd,lower,upper,pValue,ciLower,ciUpper));
        }
        return new Result(selection,effects,confidenceLevel);
    }
    private static double invert(double t,double sd,double lower,double upper,double goal) {
        double lo=-1,hi=1;
        int attempts=0;
        while(pivot(t,Math.fma(lo,sd,t),sd,lower,upper)<goal&&attempts++<1024)lo*=2;
        attempts=0;
        while(pivot(t,Math.fma(hi,sd,t),sd,lower,upper)>goal&&attempts++<1024)hi*=2;
        if(!(pivot(t,Math.fma(lo,sd,t),sd,lower,upper)>=goal
                &&pivot(t,Math.fma(hi,sd,t),sd,lower,upper)<=goal))
            throw new IllegalStateException("could not bracket selective confidence endpoint");
        double result=Double.NaN;
        for(int i=0;i<200;i++) {
            double mid=.5*lo+.5*hi;result=Math.fma(mid,sd,t);
            double discrepancy=pivot(t,result,sd,lower,upper)-goal;
            if(Math.abs(discrepancy)<1e-11)return result;
            if(mid==lo||mid==hi)break;
            if(discrepancy>0)lo=mid;else hi=mid;
        }
        if(!Double.isFinite(result)||Math.abs(pivot(t,result,sd,lower,upper)-goal)>1e-8)
            throw new IllegalStateException("selective confidence endpoint failed attained-tail check");
        return result;
    }

    private static double pivot(double t,double mean,double sd,double lower,double upper) {
        double[] tails=logTruncatedTails(t,mean,sd,lower,upper);
        double value=tails[0]-tails[1];
        if(Double.isNaN(value))throw new IllegalStateException("selective pivot is numerically unresolved");
        return value;
    }

    static double truncatedCdf(double value,double mean,double sd,double lower,double upper) {
        return Math.exp(logTruncatedTails(value,mean,sd,lower,upper)[0]);
    }

    /** Log lower and upper probabilities, preserving small gaps before shifting
     * by a possibly enormous mean. Never subtract two enormous log CDFs. */
    private static double[] logTruncatedTails(double value,double mean,double sd,double lower,double upper) {
        if(!Double.isFinite(mean)||!Double.isFinite(value)||!(sd>0)||!Double.isFinite(sd)||!(lower<upper))
            throw new IllegalStateException("invalid truncated-normal numerical range");
        if(value<=lower)return new double[]{Double.NEGATIVE_INFINITY,0};
        if(value>=upper)return new double[]{0,Double.NEGATIVE_INFINITY};
        if(mean>=upper) {
            double[] reflected=logTruncatedTails(-value,-mean,sd,-upper,-lower);
            return new double[]{reflected[1],reflected[0]};
        }
        double left,right;
        if(mean<=lower) {
            double a=(lower-mean)/sd,d=(value-lower)/sd,width=(upper-lower)/sd;
            if(!Double.isFinite(a)||!(d>0))throw new IllegalStateException("truncation distances exceed numerical range");
            double ratio=logTailRatio(a,d);
            left=logOneMinusExp(ratio);
            right=ratio+logOneMinusExp(logTailRatio(a+d,(upper-value)/sd));
            // Both masses have the same Q(a) divisor, so normalization cancels
            // it exactly even far into the tail or in a very narrow interval.
            if(!(width>0))throw new IllegalStateException("unresolved truncation width");
        } else {
            left=logInterval((lower-mean)/sd,(value-mean)/sd);
            right=logInterval((value-mean)/sd,(upper-mean)/sd);
        }
        double normalization=logAdd(left,right);
        if(!Double.isFinite(normalization))throw new IllegalStateException("truncated probability mass is numerically unresolved");
        return new double[]{left-normalization,right-normalization};
    }

    /** log(Q(a+d)/Q(a)), a>=0. Factored exponent and Mills ratio avoid
     * cancellation when a is huge but d is a representably tiny boundary gap. */
    private static double logTailRatio(double a,double d) {
        if(d==Double.POSITIVE_INFINITY)return Double.NEGATIVE_INFINITY;
        if(d==0)return 0;
        if(!(a>=0)||!(d>0)||!Double.isFinite(a))throw new IllegalStateException("invalid tail-ratio arguments");
        if(d<1e-4) {
            // Three-point Gauss-Legendre integration of the normal hazard;
            // this avoids subtracting nearly identical moderate log tails.
            double mid=a+.5*d,offset=.5*d*Math.sqrt(3.0/5.0);
            return -d*(5.0/18.0*hazard(mid-offset)+4.0/9.0*hazard(mid)+5.0/18.0*hazard(mid+offset));
        }
        if(a>=10) return -d*(a+.5*d)-Math.log1p(d/a)
            +logMillsCorrection(a+d)-logMillsCorrection(a);
        return Normal.cumulative(-(a+d),0,1,true,true)-Normal.cumulative(-a,0,1,true,true);
    }

    private static double hazard(double a) {
        if(a>=10)return a/Math.exp(logMillsCorrection(a));
        return Math.exp(-.5*a*a-.5*Math.log(2*Math.PI)-Normal.cumulative(-a,0,1,true,true));
    }

    private static double logMillsCorrection(double a) {
        double inverse=1/a,term=1,sum=1;
        for(int k=1;k<100;k++) {
            double next=-term*(2*k-1)*inverse*inverse;
            if(Math.abs(next)>=Math.abs(term))break;
            sum+=next;term=next;
            if(Math.abs(term)<Math.ulp(sum)*.25)break;
        }
        return Math.log(sum);
    }

    private static double logInterval(double a,double b) {
        if(a>=b)return Double.NEGATIVE_INFINITY;
        if(a>=0)return Normal.cumulative(-a,0,1,true,true)+logOneMinusExp(logTailRatio(a,b-a));
        if(b<=0)return logInterval(-b,-a);
        return Math.log(.5)+logAdd(logOneMinusExp(logTailRatio(0,-a)),logOneMinusExp(logTailRatio(0,b)));
    }
    private static double logOneMinusExp(double logValue) {
        return logValue<-.6931471805599453?Math.log1p(-Math.exp(logValue)):Math.log(-Math.expm1(logValue));
    }
    private static double logAdd(double a,double b) {
        double high=Math.max(a,b),low=Math.min(a,b);
        return high==Double.NEGATIVE_INFINITY?high:high+Math.log1p(Math.exp(low-high));
    }
    private static double absoluteProducts(double[] a,double[] b) {double v=0;for(int i=0;i<a.length;i++)v+=Math.abs(a[i]*b[i]);return v;}
    private static double dot(double[] a,double[] b){double v=0,error=0;for(int i=0;i<a.length;i++){double term=a[i]*b[i]-error,next=v+term;error=(next-v)-term;v=next;}return v;}
    public record Effect(int predictorIndex,double estimate,double untruncatedStandardError,
                         double truncationLower,double truncationUpper,double pValue,
                         double confidenceLower,double confidenceUpper) { }
    public record Result(PenalizedRegressionResult selectionFit,List<Effect> effects,double confidenceLevel) {
        public Result{effects=List.copyOf(effects);}
    }
}
