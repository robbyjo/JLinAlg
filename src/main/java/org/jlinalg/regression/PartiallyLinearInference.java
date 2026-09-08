/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.LeastSquaresSolver;

/** Robinson partialling-out estimator with smoother-aware slope inference.
 * Identification requires variation in the linear predictors conditional on z.
 * The bandwidth is fixed, not selected using the reported inference.
 * HC3-style inference uses original-observation influence and diagonal residual corrections;
 * it is asymptotic, not finite-sample exact, and can be conservative at high
 * smoothing leverage. Homoskedastic covariance additionally assumes the mean
 * is reproduced by the fitted smoother (otherwise smoothing bias remains). */
public final class PartiallyLinearInference {
    private PartiallyLinearInference() { }
    public static Result fit(double[] y,double[][] x,double[] z,double bandwidth,BackendPolicy policy) {
        if(y==null||x==null||z==null||x.length!=y.length||z.length!=y.length||y.length<4||x[0]==null
                ||!(bandwidth>0)||!Double.isFinite(bandwidth)||policy==null)throw new IllegalArgumentException("invalid partially linear data");
        int n=y.length,p=x[0].length,intercept=-1;
        if(p==0)throw new IllegalArgumentException("linear predictors are required");
        for(int i=0;i<n;i++) {
            if(x[i]==null||x[i].length!=p||!Double.isFinite(y[i])||!Double.isFinite(z[i]))throw new IllegalArgumentException("finite rectangular data required");
            for(double value:x[i])if(!Double.isFinite(value))throw new IllegalArgumentException("finite design required");
        }
        List<Integer> slopes=new ArrayList<>();
        for(int j=0;j<p;j++) {
            boolean constant=true;for(int i=1;i<n;i++)constant&=x[i][j]==x[0][j];
            if(constant) {
                if(x[0][j]==0||intercept>=0)throw new IllegalArgumentException("linear design has redundant constants");
                intercept=j;
            } else slopes.add(j);
        }
        int k=slopes.size();if(k==0||n<=k+1)throw new IllegalArgumentException("identifiable slopes and residual degrees of freedom required");
        double[] meanY=KernelRegression.predict(z,y,z,bandwidth),partialY=new double[n],partialX=new double[n*k];
        for(int i=0;i<n;i++)partialY[i]=y[i]-meanY[i];
        for(int j=0;j<k;j++) {
            double[] column=new double[n];for(int i=0;i<n;i++)column[i]=x[i][slopes.get(j)];
            double[] mean=KernelRegression.predict(z,column,z,bandwidth);
            for(int i=0;i<n;i++)partialX[i*k+j]=column[i]-mean[i];
        }
        double[] beta=new double[p],covariance,homoskedasticCovariance;
        double residualNoiseDegreesOfFreedom=0;
        try(var context=BackendContext.select(policy)) {
            var solved=LeastSquaresSolver.solve(partialX,partialY,n,k,false,context.backend());
            double[] coefficients=solved.coefficients(),bread=solved.unscaledCovariance();
            for(int j=0;j<k;j++)beta[slopes.get(j)]=coefficients[j];
            // beta = L y, L = (Xt'Xt)^-1 Xt' (I-S).  Smoothing induces
            // correlated transformed errors, so Xt rows are not independent scores.
            double[] influence=new double[k*n],partialInfluence=new double[k*n];
            for(int j=0;j<k;j++)for(int i=0;i<n;i++)for(int l=0;l<k;l++)
                partialInfluence[j*n+i]+=bread[j*k+l]*partialX[i*k+l];
            System.arraycopy(partialInfluence,0,influence,0,influence.length);
            for(int i=0;i<n;i++) {
                double[] weights=smootherRow(z,i,bandwidth);
                for(int j=0;j<k;j++)for(int r=0;r<n;r++)
                    influence[j*n+r]-=partialInfluence[j*n+i]*weights[r];
            }
            covariance=new double[k*k];homoskedasticCovariance=new double[k*k];
            double sumResidualSquares=0;
            for(int i=0;i<n;i++) {
                double residual=partialY[i];
                for(int j=0;j<k;j++)residual-=partialX[i*k+j]*coefficients[j];
                double[] weights=smootherRow(z,i,bandwidth);double residualDiagonal=0;
                // B = (I-Ht)(I-S) = I-S-Xt L maps original noise to residuals.
                for(int r=0;r<n;r++) {
                    double transfer=(i==r?1:0)-weights[r];
                    for(int j=0;j<k;j++)transfer-=partialX[i*k+j]*influence[j*n+r];
                    residualNoiseDegreesOfFreedom+=transfer*transfer;
                    if(i==r)residualDiagonal=transfer;
                }
                if(!(residualDiagonal>1e-12))throw new IllegalArgumentException("smoother has nonpositive residual leverage");
                double adjusted=residual/residualDiagonal;sumResidualSquares+=residual*residual;
                for(int j=0;j<k;j++)for(int l=0;l<k;l++) {
                    double crossInfluence=influence[j*n+i]*influence[l*n+i];
                    covariance[j*k+l]+=crossInfluence*adjusted*adjusted;
                    homoskedasticCovariance[j*k+l]+=crossInfluence;
                }
            }
            if(!(residualNoiseDegreesOfFreedom>1e-12))throw new IllegalArgumentException("smoother has no residual noise degrees of freedom");
            double variance=sumResidualSquares/residualNoiseDegreesOfFreedom;
            for(int j=0;j<k*k;j++)homoskedasticCovariance[j]*=variance;
        }
        double[] residualLinear=new double[n];
        for(int i=0;i<n;i++){residualLinear[i]=y[i];for(int j:slopes)residualLinear[i]-=x[i][j]*beta[j];}
        double[] smooth=KernelRegression.predict(z,residualLinear,z,bandwidth);
        if(intercept>=0) {
            double mean=0;for(double value:smooth)mean+=value/n;
            beta[intercept]=mean/x[0][intercept];for(int i=0;i<n;i++)smooth[i]-=mean;
        }
        double[] fitted=new double[n],residuals=new double[n],se=new double[k];
        for(int i=0;i<n;i++){fitted[i]=smooth[i];for(int j=0;j<p;j++)fitted[i]+=x[i][j]*beta[j];residuals[i]=y[i]-fitted[i];}
        for(int j=0;j<k;j++)se[j]=Math.sqrt(covariance[j*k+j]);
        return new Result(new PartiallyLinearRegression.Result(beta,smooth,fitted,residuals,bandwidth,1,true),
            slopes.stream().mapToInt(Integer::intValue).toArray(),covariance,se,
            homoskedasticCovariance,residualNoiseDegreesOfFreedom);
    }
    private static double[] smootherRow(double[] z,int row,double bandwidth) {
        double[] weights=new double[z.length];double total=0;
        for(int i=0;i<z.length;i++){double d=(z[row]-z[i])/bandwidth;weights[i]=Math.exp(-.5*d*d);total+=weights[i];}
        for(int i=0;i<z.length;i++)weights[i]/=total;
        return weights;
    }
    /** Covariance and generalized HC3 standard errors concern slopes only, not
     * the normalization intercept. residualNoiseDegreesOfFreedom is tr(BB'),
     * not an exact Student-t denominator. */
    public record Result(PartiallyLinearRegression.Result fit,int[] slopeIndices,double[] slopeCovariance,double[] slopeStandardErrors,
                         double[] homoskedasticSlopeCovariance,double residualNoiseDegreesOfFreedom) {
        public Result{slopeIndices=slopeIndices.clone();slopeCovariance=slopeCovariance.clone();slopeStandardErrors=slopeStandardErrors.clone();homoskedasticSlopeCovariance=homoskedasticSlopeCovariance.clone();}
        public int[] slopeIndices(){return slopeIndices.clone();}
        public double[] slopeCovariance(){return slopeCovariance.clone();}
        public double[] slopeStandardErrors(){return slopeStandardErrors.clone();}
        public double[] homoskedasticSlopeCovariance(){return homoskedasticSlopeCovariance.clone();}
    }
}
