/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.Arrays;
import java.util.List;
import java.util.function.DoubleUnaryOperator;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTriangle;
import jdistlib.accelerator.SparseCholeskyFactor;
import jdistlib.accelerator.SparseOrdering;
import jdistlib.matrix.CsrMatrix;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.mixed.SparsePrecisionMatrix;

/**
 * Large-study GLS meta-analysis from a known sparse sampling precision.
 * Independent heterogeneity is estimated without materializing the sampling
 * covariance: (Q^-1 + tau I)^-1 = (I + tau Q)^-1 Q.
 */
public final class SparseMetaAnalysis {
    private SparseMetaAnalysis() { }

    public static MetaCorrelatedResult fit(List<MetaStudy> studies,
            SparsePrecisionMatrix samplingPrecision,
            MetaAnalysisOptions options, BackendPolicy policy) {
        MetaMath.Data data=MetaMath.data(studies); int n=studies.size();
        if(samplingPrecision==null||samplingPrecision.dimension()!=n
                ||options==null||policy==null)
            throw new IllegalArgumentException("sampling precision, options and backend must match studies");
        if(options.tauSquaredEstimator()==TauSquaredEstimator.DERSIMONIAN_LAIRD)
            throw new IllegalArgumentException("sparse precision meta-analysis supports REML or Paule-Mandel heterogeneity; DL requires an observation-scale trace");
        try(BackendContext context=BackendContext.select(policy)){
            Engine engine=new Engine(samplingPrecision,context.backend());
            double[] x=new double[n];Arrays.fill(x,1);
            Fit zero=engine.fit(data.effects(),x,1,0);
            double tau=options.method()==MetaAnalysisMethod.FIXED_EFFECT?0
                : estimateTau(engine,data.effects(),x,options,zero);
            Fit fitted=engine.fit(data.effects(),x,1,tau);
            double variance=fitted.bread()[0]*MetaAnalysis.inferenceScale(options,fitted.q(),n-1);
            double margin=MetaAnalysis.critical(options,n-1)*Math.sqrt(variance);
            return new MetaCorrelatedResult(studies.stream().map(MetaStudy::name).toList(),
                MetaAnalysis.statistics(options,fitted.beta(),new double[]{Math.sqrt(variance)},n-1),
                new double[]{variance},tau,zero.q(),n-1,MetaAnalysis.normalize(fitted.wx()),
                context.provenance(),fitted.beta()[0]-margin,fitted.beta()[0]+margin);
        }
    }

    private static double estimateTau(Engine engine,double[] y,double[] x,
            MetaAnalysisOptions options,Fit zero){
        int n=y.length,p=1;double mean=Arrays.stream(y).average().orElseThrow(),sumSquares=0;
        for(double value:y)sumSquares+=(value-mean)*(value-mean);
        // Optimize a dimensionless variance multiplier. An absolute variance
        // floor or stopping tolerance changes the fit when effect units change.
        double scale=sumSquares/(n-1);
        if(scale==0)return 0;
        if(!Double.isFinite(scale))throw new ArithmeticException("effect variance exceeds numerical range");
        double upper=1;
        if(options.tauSquaredEstimator()==TauSquaredEstimator.PAULE_MANDEL){
            if(zero.q()<=n-p)return 0;
            while(engine.fit(y,x,p,scale*upper).q()>n-p){upper*=4;if(!Double.isFinite(upper))throw new ArithmeticException("cannot bracket sparse heterogeneity");}
            double lower=0;
            for(int i=0;i<options.maximumIterations();i++){
                double middle=(lower+upper)/2;
                if(engine.fit(y,x,p,scale*middle).q()>n-p)lower=middle;else upper=middle;
                if(upper-lower<=options.tolerance()*Math.max(1,upper))return scale*((lower+upper)/2);
            }
            throw new ArithmeticException("sparse Paule-Mandel did not converge");
        }
        DoubleUnaryOperator objective=t->-engine.fit(y,x,p,scale*t).reml();
        double before=objective.applyAsDouble(0),current=objective.applyAsDouble(upper);
        while(current<before){before=current;upper*=4;if(!Double.isFinite(upper))throw new ArithmeticException("cannot bracket sparse REML");current=objective.applyAsDouble(upper);}
        double lower=0,ratio=(Math.sqrt(5)-1)/2;
        double a=upper-ratio*upper,b=ratio*upper,fa=objective.applyAsDouble(a),fb=objective.applyAsDouble(b);
        for(int i=0;i<options.maximumIterations();i++){
            if(fa<fb){upper=b;b=a;fb=fa;a=upper-ratio*(upper-lower);fa=objective.applyAsDouble(a);}
            else{lower=a;a=b;fa=fb;b=lower+ratio*(upper-lower);fb=objective.applyAsDouble(b);}
            if(upper-lower<=options.tolerance()*Math.max(1,upper)){
                double candidate=(lower+upper)/2;
                candidate=refineReml(objective,candidate);
                return objective.applyAsDouble(0)<=objective.applyAsDouble(candidate)?0:scale*candidate;
            }
        }
        throw new ArithmeticException("sparse REML did not converge");
    }

    private static double refineReml(DoubleUnaryOperator objective,double candidate){
        // Objective comparisons flatten at roughly sqrt(machine epsilon) near
        // the optimum. A local fourth-order score estimate resolves that plateau
        // without dense inverse/trace calculations in this sparse fitter.
        double h=Math.min(candidate/2,Math.pow(Math.ulp(1.0),.2)*(1+candidate));
        if(!(h>0))return candidate;
        double center=objective.applyAsDouble(candidate);
        double below=objective.applyAsDouble(candidate-h),above=objective.applyAsDouble(candidate+h);
        double halfBelow=objective.applyAsDouble(candidate-h/2),halfAbove=objective.applyAsDouble(candidate+h/2);
        double score=(8*(halfAbove-halfBelow)-(above-below))/(6*h);
        double curvature=((above-center)+(below-center))/(h*h);
        double step=score/curvature;
        if(!(curvature>0)||!Double.isFinite(step)||Math.abs(step)>h)return candidate;
        double refined=candidate-step;
        return objective.applyAsDouble(refined)<=center+16*Math.ulp(Math.max(1,Math.abs(center)))
            ?refined:candidate;
    }

    private record Fit(double[] beta,double[] bread,double[] wx,double q,
            double ml,double reml){ }

    private static final class Engine{
        private final int n;private final int[] starts,columns;private final double[] values;
        private final ComputeBackend backend;private final double logDetQ;
        Engine(SparsePrecisionMatrix q,ComputeBackend backend){
            n=q.dimension();starts=q.rowStarts();columns=q.columnIndices();values=q.values();this.backend=backend;
            validateSymmetric();logDetQ=factor(0,true).logDeterminant();
        }
        Fit fit(double[] y,double[] x,int p,double tau){
            if(!(tau>=0)||!Double.isFinite(tau))throw new IllegalArgumentException("heterogeneity must be finite and nonnegative");
            SparseCholeskyFactor a=factor(tau,false);
            double[] qx=multiply(x,p),qy=multiply(y,1);
            double[] wx=a.solve(qx,p),wy=a.solve(qy);
            double[] info=new double[p*p],rhs=new double[p],xtx=new double[p*p];
            for(int i=0;i<n;i++)for(int left=0;left<p;left++){
                rhs[left]+=x[i*p+left]*wy[i];
                for(int right=0;right<p;right++){
                    info[left*p+right]+=x[i*p+left]*wx[i*p+right];
                    xtx[left*p+right]+=x[i*p+left]*x[i*p+right];
                }
            }
            CholeskyFactor information=backend.dpotrf(info,p);
            double[] beta=information.solve(rhs),bread=information.solve(MatrixOps.identity(p),p),residual=y.clone();
            for(int i=0;i<n;i++)for(int j=0;j<p;j++)residual[i]-=x[i*p+j]*beta[j];
            double[] wr=a.solve(multiply(residual,1));double quadratic=0;
            for(int i=0;i<n;i++)quadratic+=residual[i]*wr[i];
            double logDet=a.logDeterminant()-logDetQ;
            double ml=-.5*(n*Math.log(2*Math.PI)+logDet+quadratic);
            double reml=-.5*((n-p)*Math.log(2*Math.PI)+logDet+information.logDeterminant()
                -backend.dpotrf(xtx,p).logDeterminant()+quadratic);
            return new Fit(beta,bread,wx,quadratic,ml,reml);
        }
        private double[] multiply(double[] matrix,int p){
            double[] result=new double[n*p];
            for(int row=0;row<n;row++)for(int index=starts[row];index<starts[row+1];index++)
                for(int j=0;j<p;j++)result[row*p+j]+=values[index]*matrix[columns[index]*p+j];
            return result;
        }
        private SparseCholeskyFactor factor(double tau,boolean precision){
            int count=0;for(int row=0;row<n;row++)for(int k=starts[row];k<starts[row+1];k++)if(columns[k]<=row)count++;
            double[] numeric=new double[count];int[] col=new int[count],rowStart=new int[n+1];int position=0;
            for(int row=0;row<n;row++){
                rowStart[row]=position+1;boolean diagonal=false;
                for(int k=starts[row];k<starts[row+1];k++)if(columns[k]<=row){
                    col[position]=columns[k]+1;
                    numeric[position]=precision?values[k]:(columns[k]==row?1+tau*values[k]:tau*values[k]);
                    diagonal|=columns[k]==row;position++;
                }
                if(!diagonal)throw new IllegalArgumentException("sampling precision requires an explicit positive diagonal");
            }
            rowStart[n]=position+1;
            return backend.dcsrpotrf(new CsrMatrix(n,n,numeric,col,rowStart),
                MatrixTriangle.LOWER,SparseOrdering.MINIMUM_DEGREE);
        }
        private void validateSymmetric(){
            for(int row=0;row<n;row++)for(int k=starts[row];k<starts[row+1];k++){
                int column=columns[k];int reverse=Arrays.binarySearch(columns,starts[column],starts[column+1],row);
                if(reverse<0||Math.abs(values[k]-values[reverse])>1e-12*Math.max(1,Math.abs(values[k])))
                    throw new IllegalArgumentException("sampling precision must be full-CSR symmetric");
            }
        }
    }
}
