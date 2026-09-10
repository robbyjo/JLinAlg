/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

import jdistlib.F;
import jdistlib.accelerator.CholeskyFactor;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;

/** Joint multi-degree-of-freedom Kenward-Roger fixed-effect tests. */
public final class KenwardRogerTest {
    private KenwardRogerTest() { }

    public record Result(double statistic, int numeratorDegreesOfFreedom,
            double denominatorDegreesOfFreedom, double scaling,
            double pValue, double unscaledStatistic) { }

    /** Tests L beta = rhs; contrast is row-major with {@code restrictions} rows. */
    public static Result test(SparseLinearMixedModelResult fit,
            double[] contrast, int restrictions, double[] rhs,
            BackendPolicy policy) {
        if(fit==null||contrast==null||rhs==null||policy==null
                ||restrictions<1||rhs.length!=restrictions
                ||contrast.length!=restrictions*fit.beta().length)
            throw new IllegalArgumentException("fit, contrast, null values and backend are required");
        SparseFiniteDf.JointState state=fit.jointState();
        if(state==null)throw new IllegalArgumentException("joint KR testing requires a REML fit requested with KENWARD_ROGER inference");
        int p=fit.beta().length,q=restrictions,m=state.covarianceGradient().length;
        try(BackendContext context=BackendContext.select(policy)){
            double[] unadjusted=contrastCovariance(contrast,q,p,state.unadjustedCovariance());
            CholeskyFactor restrictionFactor=context.backend().dpotrf(unadjusted,q);
            double[] inverse=restrictionFactor.solve(identity(q),q);
            double[] theta=new double[p*p];
            for(int a=0;a<p;a++)for(int b=0;b<p;b++)
                for(int i=0;i<q;i++)for(int j=0;j<q;j++)
                    theta[a*p+b]+=contrast[i*p+a]*inverse[i*q+j]*contrast[j*p+b];
            double[] traces=new double[m];double[][] products=new double[m][];
            for(int i=0;i<m;i++){
                products[i]=multiply(theta,state.covarianceGradient()[i],p);
                for(int a=0;a<p;a++)traces[i]+=products[i][a*p+a];
            }
            double a1=0,a2=0;double[] parameterCovariance=state.varianceParameterCovariance();
            for(int i=0;i<m;i++)for(int j=0;j<m;j++){
                a1+=parameterCovariance[i*m+j]*traces[i]*traces[j];
                double trace=0;double[] product=multiply(products[i],products[j],p);
                for(int a=0;a<p;a++)trace+=product[a*p+a];
                a2+=parameterCovariance[i*m+j]*trace;
            }
            double b=(a1+6*a2)/(2*q);
            double g=((q+1)*a1-(q+4)*a2)/((q+2)*a2);
            double denominator=3*q+2*(1-g);
            double c1=g/denominator,c2=(q-g)/denominator,c3=(q+2-g)/denominator;
            double v0=1+c1*b,v1=1-c2*b,v2=1-c3*b;
            if(Math.abs(v0)<1e-10)v0=0;
            double rho=(1.0/q)*Math.pow((1-a2/q)/v1,2)*v0/v2;
            double ddf=4+(q+2)/(q*rho-1);
            double scaling=Math.abs(ddf-2)<1e-2?1:ddf*(1-a2/q)/(ddf-2);
            if(!(ddf>0)||!Double.isFinite(ddf)||!(scaling>0)||!Double.isFinite(scaling))
                throw new IllegalArgumentException("joint KR moments do not define an F approximation");
            double[] difference=new double[q],beta=fit.beta();
            for(int i=0;i<q;i++){
                difference[i]=-rhs[i];
                for(int j=0;j<p;j++)difference[i]+=contrast[i*p+j]*beta[j];
            }
            double[] adjusted=contrastCovariance(contrast,q,p,fit.fixedEffectCovariance());
            double[] solved=context.backend().dpotrf(adjusted,q).solve(difference);
            double wald=0;for(int i=0;i<q;i++)wald+=difference[i]*solved[i];
            double unscaled=wald/q,statistic=scaling*unscaled;
            return new Result(statistic,q,ddf,scaling,
                F.cumulative(statistic,q,ddf,false,false),unscaled);
        }
    }
    private static double[] contrastCovariance(double[] l,int q,int p,double[] covariance){
        double[] result=new double[q*q];
        for(int i=0;i<q;i++)for(int j=0;j<q;j++)for(int a=0;a<p;a++)for(int b=0;b<p;b++)
            result[i*q+j]+=l[i*p+a]*covariance[a*p+b]*l[j*p+b];
        return result;
    }
    private static double[] multiply(double[] a,double[] b,int n){
        double[] result=new double[n*n];
        for(int i=0;i<n;i++)for(int k=0;k<n;k++)for(int j=0;j<n;j++)result[i*n+j]+=a[i*n+k]*b[k*n+j];
        return result;
    }
    private static double[] identity(int n){double[] x=new double[n*n];for(int i=0;i<n;i++)x[i*n+i]=1;return x;}
}
