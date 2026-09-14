/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;

import java.util.Arrays;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;

/** Exact Gaussian historical state smoothing with symbolic diffuse levels.
 * Diffuse states are unrestricted regression effects estimated through the
 * finite-state innovations. A backward information recursion supplies state
 * moments without an observation-by-observation covariance matrix.
 * Returned moments condition on the supplied dynamics and innovation variance;
 * parameter-estimation uncertainty is not integrated. */
public final class ArimaSmoothing {
    private static final int MAXIMUM_STATES = 128;
    private ArimaSmoothing() { }

    public static Result smooth(double[] series,double[] ar,double[] ma,
            double[] difference,double innovationVariance) {
        if(series==null||series.length==0||ar==null||ma==null||difference==null
                ||difference.length==0||Math.max(ar.length,ma.length+1)+difference.length-1>MAXIMUM_STATES
                ||difference[0]!=1||!(innovationVariance>0)||!Double.isFinite(innovationVariance))
            throw new IllegalArgumentException("finite variance, positive dates, and at most 128 states required");
        for(double[] values:new double[][]{ar,ma,difference})for(double value:values)
            if(!Double.isFinite(value))throw new IllegalArgumentException("nonfinite dynamics");
        boolean[] observed=new boolean[series.length];int observationCount=0;
        for(int t=0;t<series.length;t++) {
            if(Double.isInfinite(series[t]))throw new IllegalArgumentException("use NaN for missing dates");
            observed[t]=Double.isFinite(series[t]);if(observed[t])observationCount++;
        }
        ArimaStateSpace model=new ArimaStateSpace(ar,ma,difference);
        int states=model.stateSize(),diffuse=states-model.stationarySize();
        if(observationCount<=diffuse)throw new IllegalArgumentException("insufficient observations to identify diffuse states");
        FilterPattern pattern=new FilterPattern(model,observed);
        double[][] diffuseStates=diffuseStates(model,series.length,diffuse);
        double[][] diffuseObservations=new double[diffuse][series.length];double[] z=model.observationVector();
        for(int t=0;t<series.length;t++)for(int j=0;j<diffuse;j++)
            diffuseObservations[j][t]=dotColumn(z,diffuseStates[t],states,diffuse,j);

        double[] yInnovation=pattern.innovations(series);double[][] xInnovation=new double[diffuse][];
        for(int j=0;j<diffuse;j++)xInnovation[j]=pattern.innovations(diffuseObservations[j]);
        double[] levels=new double[diffuse],informationInverse=new double[diffuse*diffuse];
        if(diffuse>0) {
            double[] information=new double[diffuse*diffuse],right=new double[diffuse];
            for(int t=0;t<series.length;t++)if(observed[t]) {
                double inverse=1/pattern.innovationVariances[t];
                for(int i=0;i<diffuse;i++) {
                    double x=xInnovation[i][t];right[i]+=x*yInnovation[t]*inverse;
                    for(int j=0;j<=i;j++)information[i*diffuse+j]+=x*xInnovation[j][t]*inverse;
                }
            }
            for(int i=0;i<diffuse;i++)for(int j=0;j<i;j++)information[j*diffuse+i]=information[i*diffuse+j];
            try(BackendContext context=BackendContext.select(BackendPolicy.CPU)) {
                var factor=context.backend().dpotrf(information,diffuse);levels=factor.solve(right);
                informationInverse=factor.solve(MatrixOps.identity(diffuse),diffuse);
            } catch(IllegalArgumentException|IllegalStateException exception) {
                throw new IllegalArgumentException("diffuse initial states are not identified by the observed dates",exception);
            }
        }
        double[] residual=series.clone();
        for(int t=0;t<series.length;t++)if(observed[t])for(int j=0;j<diffuse;j++)
            residual[t]-=diffuseObservations[j][t]*levels[j];
        Smoothed finite=pattern.smooth(residual,true);
        double[][] jacobians=new double[series.length][states*diffuse];
        for(int j=0;j<diffuse;j++) {
            Smoothed response=pattern.smooth(diffuseObservations[j],false);
            for(int t=0;t<series.length;t++)for(int i=0;i<states;i++)
                jacobians[t][i*diffuse+j]=diffuseStates[t][i*diffuse+j]-response.means[t][i];
        }
        double[][] means=finite.means,covariances=finite.covariances;
        for(int t=0;t<series.length;t++) {
            for(int i=0;i<states;i++)for(int j=0;j<diffuse;j++)means[t][i]+=diffuseStates[t][i*diffuse+j]*levels[j];
            double[] extra=diffuse==0?new double[states*states]:covarianceProduct(jacobians[t],states,diffuse,informationInverse);
            for(int i=0;i<covariances[t].length;i++)covariances[t][i]=innovationVariance*(covariances[t][i]+extra[i]);
            stabilizeCovariance(covariances[t],states,innovationVariance);
        }
        double[] signal=new double[series.length],variance=new double[series.length];
        for(int t=0;t<series.length;t++){signal[t]=dot(z,means[t]);variance[t]=Math.max(0,quadratic(z,covariances[t],states));}
        return new Result(means,covariances,signal,variance);
    }

    private static double[][] diffuseStates(ArimaStateSpace model,int dates,int diffuse) {
        int states=model.stateSize();double[][] result=new double[dates][states*diffuse];
        for(int j=0;j<diffuse;j++)result[0][(model.stationarySize()+j)*diffuse+j]=1;
        double[] transition=model.transitionMatrix();
        for(int t=1;t<dates;t++)result[t]=multiply(transition,states,states,result[t-1],diffuse);
        return result;
    }

    private static double[] covarianceProduct(double[] jacobian,int rows,int columns,double[] covariance) {
        double[] temporary=multiply(jacobian,rows,columns,covariance,columns);
        return multiply(temporary,rows,columns,transpose(jacobian,rows,columns),rows);
    }

    static void stabilizeCovariance(double[] covariance,int states,double scale) {
        for(int i=0;i<states;i++) {
            for(int j=0;j<i;j++) {
                double value=.5*(covariance[i*states+j]+covariance[j*states+i]);
                if(!Double.isFinite(value))
                    throw new IllegalArgumentException("smoothing covariance is nonfinite");
                covariance[i*states+j]=covariance[j*states+i]=value;
            }
            double variance=covariance[i*states+i];
            if(!Double.isFinite(variance)
                    || variance < -1e-8*(scale+Math.abs(variance)))
                throw new IllegalArgumentException("smoothing covariance lost precision");
            covariance[i*states+i]=Math.max(0,variance);
        }
    }

    /** Finite-state filter pattern and backward information smoother. */
    private static final class FilterPattern {
        private final int dates,states;
        private final boolean[] observed;
        private final double[] transition,noise,observation;
        private final double[][] predictedCovariances,updateGains;
        private final double[] innovationVariances;

        private FilterPattern(ArimaStateSpace model,boolean[] observed) {
            dates=observed.length;states=model.stateSize();this.observed=observed.clone();
            transition=model.transitionMatrix();noise=model.noiseMatrix();observation=model.observationVector();
            predictedCovariances=new double[dates][states*states];updateGains=new double[dates][states];
            innovationVariances=new double[dates];double[] covariance=model.initialCovariance();
            for(int t=0;t<dates;t++) {
                predictedCovariances[t]=covariance.clone();double[] filtered=covariance;
                if(observed[t]) {
                    double[] projected=matrixVector(covariance,observation,states);double variance=dot(observation,projected);
                    if(!(variance>0)||!Double.isFinite(variance))
                        throw new IllegalArgumentException("nonpositive finite-state innovation variance");
                    innovationVariances[t]=variance;
                    for(int i=0;i<states;i++)updateGains[t][i]=projected[i]/variance;
                    filtered=covariance.clone();
                    for(int i=0;i<states;i++)for(int j=0;j<=i;j++) {
                        double value=covariance[i*states+j]-projected[i]*projected[j]/variance;
                        filtered[i*states+j]=filtered[j*states+i]=value;
                    }
                }
                if(t+1<dates)covariance=predictCovariance(filtered);
            }
        }

        private double[] innovations(double[] values){return forward(values,false).innovations;}

        private Smoothed smooth(double[] values,boolean covariance) {
            Forward pass=forward(values,true);double[][] means=new double[dates][states];
            double[][] covariances=covariance?new double[dates][states*states]:new double[0][];
            double[] nextR=new double[states],nextN=covariance?new double[states*states]:null;
            for(int t=dates-1;t>=0;t--) {
                double[] l=transition.clone();
                if(observed[t]) {
                    double[] predictedGain=matrixVector(transition,updateGains[t],states);
                    for(int i=0;i<states;i++)for(int j=0;j<states;j++)
                        l[i*states+j]-=predictedGain[i]*observation[j];
                }
                double[] previousR=transposeVector(l,nextR,states);
                if(observed[t]) {
                    double standardized=pass.innovations[t]/innovationVariances[t];
                    for(int i=0;i<states;i++)previousR[i]+=observation[i]*standardized;
                }
                double[] correction=matrixVector(predictedCovariances[t],previousR,states);
                for(int i=0;i<states;i++)means[t][i]=pass.predictedMeans[t][i]+correction[i];
                if(covariance) {
                    double[] previousN=sandwichTranspose(l,nextN,states);
                    if(observed[t]) {
                        double inverse=1/innovationVariances[t];
                        for(int i=0;i<states;i++)for(int j=0;j<states;j++)
                            previousN[i*states+j]+=observation[i]*observation[j]*inverse;
                    }
                    double[] reduction=sandwich(predictedCovariances[t],previousN,states);
                    for(int i=0;i<reduction.length;i++)covariances[t][i]=predictedCovariances[t][i]-reduction[i];
                    nextN=previousN;
                }
                nextR=previousR;
            }
            return new Smoothed(means,covariances);
        }

        private Forward forward(double[] values,boolean retainMeans) {
            if(values==null||values.length!=dates)throw new IllegalArgumentException("smoothing values must match the date count");
            double[] state=new double[states];double[][] predictedMeans=retainMeans?new double[dates][states]:new double[0][];
            double[] innovations=new double[dates];Arrays.fill(innovations,Double.NaN);
            for(int t=0;t<dates;t++) {
                if(retainMeans)predictedMeans[t]=state.clone();
                if(observed[t]) {
                    if(!Double.isFinite(values[t]))throw new IllegalArgumentException("every observed-pattern date needs a finite value");
                    double innovation=values[t]-dot(observation,state);innovations[t]=innovation;
                    for(int i=0;i<states;i++)state[i]+=updateGains[t][i]*innovation;
                }
                if(t+1<dates)state=matrixVector(transition,state,states);
            }
            return new Forward(predictedMeans,innovations);
        }

        private double[] predictCovariance(double[] covariance) {
            double[] temporary=multiply(transition,states,states,covariance,states);
            double[] result=multiply(temporary,states,states,transpose(transition,states,states),states);
            for(int i=0;i<result.length;i++)result[i]+=noise[i];return result;
        }
    }

    private static double[] sandwichTranspose(double[] left,double[] center,int dimension) {
        double[] temporary=multiply(center,dimension,dimension,left,dimension);
        return multiply(transpose(left,dimension,dimension),dimension,dimension,temporary,dimension);
    }

    private static double[] sandwich(double[] covariance,double[] information,int dimension) {
        double[] temporary=multiply(covariance,dimension,dimension,information,dimension);
        return multiply(temporary,dimension,dimension,covariance,dimension);
    }

    private static double[] matrixVector(double[] matrix,double[] vector,int dimension) {
        double[] result=new double[dimension];
        for(int i=0;i<dimension;i++)for(int j=0;j<dimension;j++)result[i]+=matrix[i*dimension+j]*vector[j];
        return result;
    }

    private static double[] transposeVector(double[] matrix,double[] vector,int dimension) {
        double[] result=new double[dimension];
        for(int i=0;i<dimension;i++)for(int j=0;j<dimension;j++)result[j]+=matrix[i*dimension+j]*vector[i];
        return result;
    }

    private static double dotColumn(double[] vector,double[] matrix,int rows,int columns,int column) {
        double result=0;for(int row=0;row<rows;row++)result+=vector[row]*matrix[row*columns+column];return result;
    }
    private static double dot(double[] first,double[] second) {
        double result=0;for(int i=0;i<first.length;i++)result+=first[i]*second[i];return result;
    }
    private static double quadratic(double[] vector,double[] matrix,int dimension) {
        return dot(vector,matrixVector(matrix,vector,dimension));
    }
    private static double[] multiply(double[] first,int rows,int shared,double[] second,int columns) {
        double[] result=new double[rows*columns];
        for(int i=0;i<rows;i++)for(int k=0;k<shared;k++) {
            double value=first[i*shared+k];if(value==0)continue;
            for(int j=0;j<columns;j++)result[i*columns+j]+=value*second[k*columns+j];
        }
        return result;
    }
    private static double[] transpose(double[] matrix,int rows,int columns) {
        double[] result=new double[matrix.length];
        for(int i=0;i<rows;i++)for(int j=0;j<columns;j++)result[j*rows+i]=matrix[i*columns+j];
        return result;
    }

    private record Forward(double[][] predictedMeans,double[] innovations) { }
    private record Smoothed(double[][] means,double[][] covariances) { }

    public record Result(double[][] states,double[][] stateCovariances,double[] signal,double[] signalVariance) {
        public Result {states=copy(states);stateCovariances=copy(stateCovariances);signal=signal.clone();signalVariance=signalVariance.clone();}
        public double[][] states(){return copy(states);}
        public double[][] stateCovariances(){return copy(stateCovariances);}
        public double[] signal(){return signal.clone();}
        public double[] signalVariance(){return signalVariance.clone();}
        private static double[][] copy(double[][] values){return Arrays.stream(values).map(double[]::clone).toArray(double[][]::new);}
    }
}
