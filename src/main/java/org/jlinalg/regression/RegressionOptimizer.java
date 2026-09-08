/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

/** Small dense BFGS solver for scaled regression coordinates. */
final class RegressionOptimizer {
    private RegressionOptimizer() { }
    record Evaluation(double value, double[] gradient) { }
    record Result(double[] point, int iterations, boolean converged) { }
    interface Objective { Evaluation evaluate(double[] point); }

    static Result minimize(Objective objective, int dimension, int maximumIterations,
                           double tolerance, double initialStep) {
        double[] point = new double[dimension], inverse = identity(dimension);
        Evaluation current = objective.evaluate(point);
        for (int iteration = 0; iteration < maximumIterations; iteration++) {
            if (maxAbs(current.gradient) <= tolerance)
                return new Result(point, iteration, true);
            double[] direction = new double[dimension];
            for (int i=0;i<dimension;i++) for (int j=0;j<dimension;j++)
                direction[i] -= inverse[i*dimension+j]*current.gradient[j];
            double slope = dot(direction,current.gradient);
            if (!(slope < 0) || !Double.isFinite(slope)) {
                inverse=identity(dimension);
                for (int i=0;i<dimension;i++) direction[i]=-current.gradient[i];
                slope=-dot(current.gradient,current.gradient);
            }
            double step=iteration==0?initialStep:1;
            Evaluation trial=null; double[] next=new double[dimension]; boolean accepted=false;
            for(int attempt=0;attempt<64;attempt++) {
                for(int i=0;i<dimension;i++) next[i]=point[i]+step*direction[i];
                trial=objective.evaluate(next);
                if(Double.isFinite(trial.value) && trial.value<=current.value+1e-4*step*slope
                        && Double.isFinite(maxAbs(trial.gradient))) { accepted=true;break; }
                step*=.5;
            }
            // Stalling is not convergence: only the score can certify convergence.
            if(!accepted) return new Result(point,iteration+1,false);
            double[] s=new double[dimension], y=new double[dimension], hy=new double[dimension];
            for(int i=0;i<dimension;i++){s[i]=next[i]-point[i];y[i]=trial.gradient[i]-current.gradient[i];}
            double ys=dot(y,s);
            if(ys>1e-12*Math.sqrt(dot(y,y)*dot(s,s))) {
                for(int i=0;i<dimension;i++)for(int j=0;j<dimension;j++)hy[i]+=inverse[i*dimension+j]*y[j];
                double coefficient=(ys+dot(y,hy))/(ys*ys);
                for(int i=0;i<dimension;i++)for(int j=0;j<dimension;j++)
                    inverse[i*dimension+j]+=coefficient*s[i]*s[j]-(hy[i]*s[j]+s[i]*hy[j])/ys;
            }
            point=next;current=trial;
        }
        return new Result(point,maximumIterations,maxAbs(current.gradient)<=tolerance);
    }
    static double[] scaleColumns(double[] x,int rows,int columns) {
        double[] scales=new double[columns];
        for(int j=0;j<columns;j++) {
            for(int i=0;i<rows;i++)scales[j]=Math.max(scales[j],Math.abs(x[i*columns+j]));
            if(scales[j]==0)scales[j]=1;
            for(int i=0;i<rows;i++)x[i*columns+j]/=scales[j];
        }
        return scales;
    }
    private static double[] identity(int n){double[] a=new double[n*n];for(int i=0;i<n;i++)a[i*n+i]=1;return a;}
    private static double dot(double[] a,double[] b){double sum=0;for(int i=0;i<a.length;i++)sum+=a[i]*b[i];return sum;}
    private static double maxAbs(double[] x){double max=0;for(double v:x)max=Math.max(max,Math.abs(v));return max;}
}
