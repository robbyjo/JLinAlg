/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.math.MathFunctions;
import jdistlib.math.opt.Bobyqa;
import jdistlib.math.opt.OptimizationResult;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparsePrecisionMatrix;

/**
 * Low-dimensional adaptive Gauss-Hermite GLMM likelihood for arbitrary crossed
 * and sparse-precision Gaussian random effects. Tensor growth is explicit and
 * bounded; use sparse Laplace/PQL for structures exceeding the node budget.
 */
public final class MultidimensionalGlmmQuadrature {
    private MultidimensionalGlmmQuadrature(){ }

    public static MultidimensionalQuadratureResult fit(double[] response,
            double[][] fixed,List<RandomEffectTerm> randomEffects,
            List<SparsePrecisionMatrix> precisionBases,GlmFamily family,
            MultidimensionalQuadratureOptions options,BackendPolicy policy){
        try(BackendContext context=BackendContext.select(policy)){
            Data data=new Data(response,fixed,randomEffects,precisionBases,family,
                options,context.backend());
            int dimensions=data.p+data.terms;
            double[] initial=new double[dimensions],lower=new double[dimensions],upper=new double[dimensions];
            for(int i=0;i<data.p;i++){lower[i]=-options.maximumAbsoluteCoefficient();upper[i]=options.maximumAbsoluteCoefficient();}
            for(int i=data.p;i<dimensions;i++){
                initial[i]=Math.log(.5);lower[i]=Math.log(options.minimumVariance());upper[i]=Math.log(options.maximumVariance());
            }
            int[] calls={0};
            java.util.function.ToDoubleFunction<double[]> objective=point->{
                calls[0]++;double[] beta=Arrays.copyOf(point,data.p),variance=new double[data.terms];
                for(int i=0;i<variance.length;i++)variance[i]=Math.exp(point[data.p+i]);
                MultidimensionalQuadratureEvaluation value=data.evaluate(beta,variance);
                return value.converged()?-value.logLikelihood():1e100;
            };
            OptimizationResult optimum;
            try{optimum=Bobyqa.bobyqa(initial,lower,upper,objective::applyAsDouble,
                Math.min(2*dimensions+1,(dimensions+1)*(dimensions+2)/2),.5,
                options.relativeTolerance(),options.maximumOuterEvaluations(),true);}
            catch(ArithmeticException|ArrayIndexOutOfBoundsException failure){
                optimum=coordinate(initial,lower,upper,objective,options);
            }
            double[] point=optimum.mX==null?initial:optimum.mX;
            double[] beta=Arrays.copyOf(point,data.p),variance=new double[data.terms];
            for(int i=0;i<variance.length;i++)variance[i]=Math.exp(point[data.p+i]);
            MultidimensionalQuadratureEvaluation value=data.evaluate(beta,variance);
            boolean stationary=stationary(point,lower,upper,objective,
                Math.max(1e-5,Math.sqrt(options.relativeTolerance())));
            boolean converged=stationary&&value.converged();
            return new MultidimensionalQuadratureResult(family.name(),beta,variance,
                value.logLikelihood(),value.order(),value.nodes(),value.estimatedError(),
                value.converged(),calls[0],converged,!value.converged()
                    ?"quadrature tolerance or node budget not met"
                    :!stationary?"marginal likelihood score tolerance not met":"converged");
        }
    }

    public static MultidimensionalQuadratureEvaluation evaluate(double[] response,
            double[][] fixed,List<RandomEffectTerm> randomEffects,
            List<SparsePrecisionMatrix> precisionBases,GlmFamily family,
            double[] beta,double[] variances,
            MultidimensionalQuadratureOptions options,BackendPolicy policy){
        try(BackendContext context=BackendContext.select(policy)){
            return new Data(response,fixed,randomEffects,precisionBases,family,
                options,context.backend()).evaluate(beta,variances);
        }
    }

    private static OptimizationResult coordinate(double[] initial,double[] lower,
            double[] upper,java.util.function.ToDoubleFunction<double[]> objective,
            MultidimensionalQuadratureOptions options){
        double[] point=initial.clone();double best=objective.applyAsDouble(point);int calls=1;double step=.5;
        while(calls<options.maximumOuterEvaluations()&&step>options.relativeTolerance()){
            boolean improved=false;
            for(int i=0;i<point.length&&calls<options.maximumOuterEvaluations();i++){
                double original=point[i];
                for(int sign:new int[]{-1,1}){
                    point[i]=Math.max(lower[i],Math.min(upper[i],original+sign*step));
                    double candidate=objective.applyAsDouble(point);calls++;
                    if(candidate<best){best=candidate;original=point[i];improved=true;}else point[i]=original;
                }
            }
            if(!improved)step*=.5;
        }
        return new OptimizationResult(point,best,calls,calls<options.maximumOuterEvaluations());
    }
    private static boolean stationary(double[] point,double[] lower,double[] upper,
            java.util.function.ToDoubleFunction<double[]> objective,double tolerance){
        double center=objective.applyAsDouble(point);
        for(int i=0;i<point.length;i++){
            double h=1e-5*(1+Math.abs(point[i]));double lo=Math.max(lower[i],point[i]-h),hi=Math.min(upper[i],point[i]+h);
            double[] trial=point.clone();trial[i]=lo;double below=objective.applyAsDouble(trial);trial[i]=hi;double above=objective.applyAsDouble(trial);
            double gradient=(above-below)/(hi-lo);
            if(point[i]<=lower[i]&&gradient>0||point[i]>=upper[i]&&gradient<0)gradient=0;
            if(!Double.isFinite(gradient)||Math.abs(gradient)>tolerance)return false;
        }
        return center<1e100;
    }

    private static final class Data{
        final double[] y,constant;final double[][] x,z;final double[][] bases;
        final int n,p,q,terms;final boolean binomial;final int[] starts;
        final MultidimensionalQuadratureOptions options;final ComputeBackend backend;
        Data(double[] response,double[][] fixed,List<RandomEffectTerm> effects,
                List<SparsePrecisionMatrix> precision,GlmFamily family,
                MultidimensionalQuadratureOptions options,ComputeBackend backend){
            if(family!=GlmFamilies.binomial()&&family!=GlmFamilies.poisson())
                throw new IllegalArgumentException("multidimensional quadrature supports binomial(logit) and Poisson(log)");
            if(response==null||fixed==null||effects==null||effects.isEmpty()
                    ||options==null||backend==null||fixed.length!=response.length||response.length==0)
                throw new IllegalArgumentException("invalid multidimensional quadrature inputs");
            this.options=options;this.backend=backend;this.binomial=family==GlmFamilies.binomial();
            n=response.length;p=fixed[0].length;y=response.clone();x=new double[n][p];constant=new double[n];
            for(int i=0;i<n;i++){
                if(fixed[i]==null||fixed[i].length!=p||!Double.isFinite(y[i])||y[i]<0||y[i]!=Math.rint(y[i])||(binomial&&y[i]>1))
                    throw new IllegalArgumentException("invalid response or fixed design row");
                constant[i]=-MathFunctions.lgammafn(y[i]+1);
                for(int j=0;j<p;j++){if(!Double.isFinite(fixed[i][j]))throw new IllegalArgumentException("nonfinite fixed design");x[i][j]=fixed[i][j];}
            }
            terms=effects.size();starts=new int[terms+1];int columns=0;HashSet<String> names=new HashSet<>();
            for(int term=0;term<terms;term++){
                RandomEffectTerm value=effects.get(term);
                if(value==null||value.observations()!=n||!names.add(value.name()))throw new IllegalArgumentException("random terms must match rows and have unique names");
                if(value.coefficients()>Integer.MAX_VALUE-columns)
                    throw new IllegalArgumentException("random-effect dimension exceeds numerical range");
                starts[term]=columns;columns+=value.coefficients();
            }
            starts[terms]=columns;q=columns;
            if(q<2)throw new IllegalArgumentException("use GlmmQuadrature for a scalar random intercept");
            // Reject before allocating dense random-effect designs/precisions
            // or entering the optimizer, not only when refining the grid.
            long initialNodes=power(options.initialOrder(),q);
            if(initialNodes<0||initialNodes>options.maximumTotalNodes())
                throw new IllegalArgumentException("initial quadrature grid exceeds maximumTotalNodes");
            z=new double[n][q];
            for(int term=0;term<terms;term++){
                double[] design=effects.get(term).design();int width=effects.get(term).coefficients();
                for(int i=0;i<n;i++)System.arraycopy(design,i*width,z[i],starts[term],width);
            }
            if(precision!=null&&precision.size()!=terms)throw new IllegalArgumentException("precision bases must match random terms");
            bases=new double[terms][];
            for(int term=0;term<terms;term++){
                SparsePrecisionMatrix value=precision==null?SparsePrecisionMatrix.identity(starts[term+1]-starts[term]):precision.get(term);
                if(value.dimension()!=starts[term+1]-starts[term])throw new IllegalArgumentException("precision dimensions must match random terms");
                bases[term]=dense(value);
                backend.dpotrf(bases[term].clone(),value.dimension());
            }
            double[] gram=new double[p*p];for(double[] row:x)for(int a=0;a<p;a++)for(int b=0;b<p;b++)gram[a*p+b]+=row[a]*row[b];
            backend.dpotrf(gram,p);
        }
        MultidimensionalQuadratureEvaluation evaluate(double[] beta,double[] variance){
            if(beta==null||beta.length!=p||variance==null||variance.length!=terms)
                throw new IllegalArgumentException("coefficient and variance dimensions are invalid");
            double[] precision=new double[q*q];double logDet=0;
            for(int term=0;term<terms;term++){
                if(!(variance[term]>0)||!Double.isFinite(variance[term]))throw new IllegalArgumentException("variances must be finite and positive");
                int width=starts[term+1]-starts[term];
                CholeskyFactor base=backend.dpotrf(bases[term].clone(),width);
                logDet+=base.logDeterminant()-width*Math.log(variance[term]);
                for(int i=0;i<width;i++)for(int j=0;j<width;j++)precision[(starts[term]+i)*q+starts[term]+j]=bases[term][i*width+j]/variance[term];
            }
            double[] eta=new double[n];
            for(int i=0;i<n;i++)for(int j=0;j<p;j++){if(!Double.isFinite(beta[j]))throw new IllegalArgumentException("nonfinite coefficient");eta[i]+=x[i][j]*beta[j];}
            Mode mode=mode(eta,precision);if(mode==null)return new MultidimensionalQuadratureEvaluation(Double.NEGATIVE_INFINITY,0,0,Double.POSITIVE_INFINITY,false,new double[q]);
            int order=options.initialOrder();long nodes=power(order,q);double previous=integrate(eta,precision,logDet,mode,order),error=Double.POSITIVE_INFINITY;int stable=0;
            while(order<options.maximumOrder()){
                int next=Math.min(options.maximumOrder(),order+Math.max(2,order/2));long nextNodes=power(next,q);
                if(nextNodes<0||nextNodes>options.maximumTotalNodes())break;
                double value=integrate(eta,precision,logDet,mode,next);error=Math.abs(value-previous);stable=error<=options.quadratureTolerance()?stable+1:0;
                previous=value;order=next;nodes=nextNodes;if(stable>=2)break;
            }
            return new MultidimensionalQuadratureEvaluation(previous,order,nodes,error,stable>=2&&Double.isFinite(previous),mode.value());
        }
        private Mode mode(double[] eta,double[] precision){
            double[] b=new double[q];
            for(int iteration=0;iteration<100;iteration++){
                Derivatives d=derivatives(eta,b,precision);double[] step;
                try{step=backend.dpotrf(d.hessian(),q).solve(d.score());}catch(RuntimeException failure){return null;}
                double norm=Arrays.stream(step).map(Math::abs).max().orElse(0);if(norm<=1e-10*(1+Arrays.stream(b).map(Math::abs).max().orElse(0)))
                    return new Mode(b,d.hessian());
                double before=kernel(eta,b,precision),scale=1;boolean accepted=false;
                while(scale>1e-8){double[] candidate=b.clone();for(int i=0;i<q;i++)candidate[i]+=scale*step[i];double value=kernel(eta,candidate,precision);
                    if(value>=before){b=candidate;accepted=true;break;}scale*=.5;}
                if(!accepted)return null;
            }
            return null;
        }
        private Derivatives derivatives(double[] eta,double[] b,double[] precision){
            double[] score=new double[q],hessian=precision.clone();
            for(int i=0;i<n;i++){
                double predictor=eta[i];for(int j=0;j<q;j++)predictor+=z[i][j]*b[j];
                double mean=binomial?logistic(predictor):Math.exp(predictor);
                double residual=y[i]-mean,weight=binomial?mean*(1-mean):mean;
                for(int a=0;a<q;a++){score[a]+=z[i][a]*residual;for(int c=0;c<q;c++)hessian[a*q+c]+=weight*z[i][a]*z[i][c];}
            }
            for(int a=0;a<q;a++)for(int c=0;c<q;c++)score[a]-=precision[a*q+c]*b[c];
            return new Derivatives(score,hessian);
        }
        private double kernel(double[] eta,double[] b,double[] precision){
            double value=0;
            for(int i=0;i<n;i++){
                double predictor=eta[i];for(int j=0;j<q;j++)predictor+=z[i][j]*b[j];
                value+=binomial?(y[i]*predictor-softplus(predictor))
                    :constant[i]+(y[i]==0?0:y[i]*predictor)-Math.exp(predictor);
            }
            for(int i=0;i<q;i++)for(int j=0;j<q;j++)value-=.5*b[i]*precision[i*q+j]*b[j];
            return value;
        }
        private double integrate(double[] eta,double[] precision,double logDet,Mode mode,int order){
            double[] inverse=backend.dpotrf(mode.hessian(),q).solve(MatrixOps.identity(q),q);
            CholeskyFactor covariance=backend.dpotrf(inverse.clone(),q);
            double[] lower=lower(inverse,q);
            AdaptiveHermiteRule rule=AdaptiveHermiteRule.of(order);Accumulator sum=new Accumulator();
            enumerate(0,new double[q],0,rule,lower,eta,precision,mode.value(),sum);
            return sum.logSum+covariance.logDeterminant()/2-.5*q*Math.log(Math.PI)+.5*logDet;
        }
        private void enumerate(int dimension,double[] node,double logWeight,AdaptiveHermiteRule rule,
                double[] lower,double[] eta,double[] precision,double[] mode,Accumulator sum){
            if(dimension==q){double[] b=mode.clone();double squares=0;for(int i=0;i<q;i++){squares+=node[i]*node[i];for(int j=0;j<=i;j++)b[i]+=Math.sqrt(2)*lower[i*q+j]*node[j];}
                sum.add(logWeight+kernel(eta,b,precision)+squares);return;}
            for(int i=0;i<rule.nodes.length;i++){node[dimension]=rule.nodes[i];enumerate(dimension+1,node,logWeight+rule.logWeights[i],rule,lower,eta,precision,mode,sum);}
        }
        private double[] lower(double[] inverse,int dimension){
            double[] result=new double[dimension*dimension];
            for(int i=0;i<dimension;i++)for(int j=0;j<=i;j++){
                double value=inverse[i*dimension+j];for(int k=0;k<j;k++)value-=result[i*dimension+k]*result[j*dimension+k];
                result[i*dimension+j]=i==j?Math.sqrt(value):value/result[j*dimension+j];
            }
            return result;
        }
        private static double[] dense(SparsePrecisionMatrix matrix){
            int n=matrix.dimension();double[] result=new double[n*n],values=matrix.values();int[] starts=matrix.rowStarts(),columns=matrix.columnIndices();
            for(int i=0;i<n;i++)for(int k=starts[i];k<starts[i+1];k++)result[i*n+columns[k]]=values[k];return result;
        }
        // A negative sentinel distinguishes overflow from an exact user budget.
        private static long power(int base,int exponent){long result=1;for(int i=0;i<exponent;i++){if(result>Long.MAX_VALUE/base)return -1;result*=base;}return result;}
    }
    private record Mode(double[] value,double[] hessian){ }
    private record Derivatives(double[] score,double[] hessian){ }
    private static final class Accumulator{double logSum=Double.NEGATIVE_INFINITY;void add(double x){if(logSum==Double.NEGATIVE_INFINITY)logSum=x;else{double m=Math.max(logSum,x);logSum=m+Math.log(Math.exp(logSum-m)+Math.exp(x-m));}}}
    private static double logistic(double x){return x>=0?1/(1+Math.exp(-x)):Math.exp(x)/(1+Math.exp(x));}
    private static double softplus(double x){return Math.max(x,0)+Math.log1p(Math.exp(-Math.abs(x)));}
}
