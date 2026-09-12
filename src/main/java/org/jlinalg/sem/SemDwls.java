/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;

/** Diagonally weighted least squares with full asymptotic moment covariance
 * for sandwich errors and a mean/variance adjusted, scaled-and-shifted WLSMV
 * test. Supports continuous/ordinal mixed moments. Input Gamma is the covariance
 * of asymptotically scaled sample statistics, NOT their finite-sample covariance.
 * The conventional n-1 finite-sample multiplier is used for tests and errors. */
public final class SemDwls {
    private SemDwls() { }
    /** Statistic order: each variable's mean (continuous) or standardized
     * thresholds (ordinal), followed by row-wise lower covariance triangle,
     * omitting ordinal diagonals. Ordinal latent responses have unit variance;
     * continuous-ordinal covariances are in continuous response units.
     * Gamma must include every threshold/covariance cross-covariance. */
    public record Moments(int observations,int[] categoryCounts,double[] statistics,double[] gamma) {
        public Moments {
            if(observations<3||categoryCounts==null||statistics==null||gamma==null)
                throw new IllegalArgumentException("moments and asymptotic covariance required");
            categoryCounts=categoryCounts.clone();statistics=statistics.clone();gamma=gamma.clone();
            int m=dimension(categoryCounts);
            if(statistics.length!=m||gamma.length!=m*m)throw new IllegalArgumentException("moment dimensions differ");
            for(double v:statistics)if(!Double.isFinite(v))throw new IllegalArgumentException("nonfinite statistic");
            for(int i=0;i<m;i++)for(int j=0;j<m;j++)if(!Double.isFinite(gamma[i*m+j])||Math.abs(gamma[i*m+j]-gamma[j*m+i])>1e-10*(1+Math.abs(gamma[i*m+j])))
                throw new IllegalArgumentException("Gamma must be finite symmetric");
            RamFit.chol(gamma,m);
        }
        public int[] categoryCounts(){return categoryCounts.clone();}
        public double[] statistics(){return statistics.clone();}
        public double[] gamma(){return gamma.clone();}
    }
    /** MAR-capable moment preparation by saturated joint mixed likelihood.
     * This ML moment estimator differs from lavaan's default two-stage ordinal
     * sample statistics. The resulting DWLS/WLSMV algebra is the same. */
    public static Moments prepare(double[][] data,int[] categories,String... variables) {
        SemModel.Builder b=SemModel.builder(variables).meanStructure();
        if(categories==null||categories.length!=variables.length)throw new IllegalArgumentException("category dimensions differ");
        for(int j=0;j<variables.length;j++) {
            if(categories[j]>0)b.fixedVariance(variables[j],1).fixedIntercept(variables[j],0);
            for(int k=0;k<j;k++)b.covariance(variables[j],variables[k],0);
        }
        var fit=SemMixed.fit(data,categories,b.build());
        if(!fit.informationAvailable())throw new IllegalArgumentException("saturated mixed moments did not converge with identified information");
        Mapping mapping=new Mapping(categories,fit.engine.model);
        double[] statistics=mapping.values(fit.point);
        // Empirical case-score sandwich also propagates uncertainty in thresholds.
        int n=fit.engine.data.length,k=fit.point.length;double[][] scores=new double[n][k];
        for(int j=0;j<k;j++) {
            double h=2e-5*(1+Math.abs(fit.point[j]));double[] x=fit.point.clone();x[j]+=h;
            var above=RamFit.distribution(fit.engine.model,Arrays.copyOf(x,fit.engine.model.freeParameterCount()));double[][] ta=fit.engine.thresholds(x);
            x[j]-=2*h;var below=RamFit.distribution(fit.engine.model,Arrays.copyOf(x,fit.engine.model.freeParameterCount()));double[][] tb=fit.engine.thresholds(x);
            double average=0;
            for(int i=0;i<n;i++){scores[i][j]=(fit.engine.logDensity(fit.engine.data[i],above,ta)-fit.engine.logDensity(fit.engine.data[i],below,tb))/(2*h);average+=scores[i][j]/n;}
            for(int i=0;i<n;i++)scores[i][j]-=average;
        }
        double[] raw=SemInference.sandwich(scores,java.util.stream.IntStream.range(0,n).toArray(),fit.rawCovariance);
        double[] gamma=SemMixed.transformedCovariance(fit.point,raw,mapping::values);
        for(int i=0;i<gamma.length;i++)gamma[i]*=n-1;
        return new Moments(n,categories,statistics,gamma);
    }
    public static Result fit(Moments moments,SemModel model) {return fit(moments,model,SemOptions.defaults());}
    public static Result fit(Moments moments,SemModel model,SemOptions options) {
        Mapping mapping=new Mapping(moments.categoryCounts,model);int m=moments.statistics.length,k=mapping.start.length;
        if(k>m)throw new IllegalArgumentException("more parameters than sample moments");
        double[] weights=new double[m];for(int i=0;i<m;i++)weights[i]=1/moments.gamma[i*m+i];
        var optimum=SemOptimizer.minimize(x->objective(x,mapping,moments.statistics,weights),mapping.start,options.maximumEvaluations(),options.tolerance());
        double[] implied=mapping.values(optimum.point()),residual=new double[m];double chi=0;
        // lavaan/Mplus conventional finite-sample multiplier for WLS tests.
        // The same finite-sample convention is used for the parameter sandwich.
        for(int i=0;i<m;i++){residual[i]=moments.statistics[i]-implied[i];chi+=(moments.observations-1)*weights[i]*residual[i]*residual[i];}
        double[] covariance=SemInformation.unavailable(k);double scale=Double.NaN,shift=Double.NaN,adjusted=Double.NaN;int df=m-k;
        if(optimum.converged())try {
            double[][] jac=jacobian(mapping,optimum.point());double[] information=new double[k*k];
            for(int a=0;a<k;a++)for(int b=0;b<k;b++)for(int i=0;i<m;i++)information[a*k+b]+=jac[i][a]*weights[i]*jac[i][b];
            double[] bread=RamFit.informationInverse(information,k),transfer=new double[k*m];
            for(int a=0;a<k;a++)for(int i=0;i<m;i++)for(int b=0;b<k;b++)transfer[a*m+i]+=bread[a*k+b]*jac[i][b]*weights[i];
            covariance=new double[k*k];
            for(int a=0;a<k;a++)for(int b=0;b<k;b++)for(int i=0;i<m;i++)for(int j=0;j<m;j++)
                covariance[a*k+b]+=transfer[a*m+i]*moments.gamma[i*m+j]*transfer[b*m+j]/(moments.observations-1);
            double[] u=new double[m*m];
            for(int i=0;i<m;i++)for(int j=0;j<m;j++) {
                u[i*m+j]=i==j?weights[i]:0;
                for(int a=0;a<k;a++)u[i*m+j]-=weights[i]*jac[i][a]*transfer[a*m+j];
            }
            double[] ug=RamFit.mm(u,moments.gamma,m);double trace=0,trace2=0;
            for(int i=0;i<m;i++){trace+=ug[i*m+i];for(int j=0;j<m;j++)trace2+=ug[i*m+j]*ug[j*m+i];}
            if(df>0&&trace>0&&trace2>0){scale=Math.sqrt(trace2/df);shift=df-trace/scale;adjusted=chi/scale+shift;}
        }catch(IllegalArgumentException unavailable){df=-1;}
        if(!optimum.converged()||!SemInformation.finiteCovariance(covariance,k)) {
            df=-1;scale=Double.NaN;shift=Double.NaN;adjusted=Double.NaN;
        }
        return new Result(mapping.labels,mapping.natural(optimum.point()),SemMixed.transformedCovariance(optimum.point(),covariance,mapping::natural),
            chi,adjusted,df,scale,shift,df>0?jdistlib.ChiSquare.cumulative(adjusted,df,false,false):Double.NaN,optimum.converged());
    }
    private static SemOptimizer.Value objective(double[] x,Mapping mapping,double[] observed,double[] weights) {
        try {
            double[] fitted=mapping.values(x);double[][] jac=jacobian(mapping,x);double value=0;double[] gradient=new double[x.length];
            for(int i=0;i<fitted.length;i++){double e=fitted[i]-observed[i];value+=.5*weights[i]*e*e;for(int j=0;j<x.length;j++)gradient[j]+=jac[i][j]*weights[i]*e;}
            return new SemOptimizer.Value(value,gradient);
        }catch(IllegalArgumentException failure){return new SemOptimizer.Value(Double.POSITIVE_INFINITY,null);}
    }
    static double[][] jacobian(Mapping mapping,double[] x) {
        double[][] jac=new double[dimension(mapping.categories)][x.length];
        for(int j=0;j<x.length;j++) {
            double h=1e-5*(1+Math.abs(x[j]));double[] point=x.clone();point[j]+=h;double[] above=mapping.values(point);point[j]-=2*h;double[] below=mapping.values(point);
            for(int i=0;i<above.length;i++)jac[i][j]=(above[i]-below[i])/(2*h);
        }return jac;
    }
    static int dimension(int[] categories) {
        int m=0;for(int j=0;j<categories.length;j++) {
            if(categories[j]!=0&&categories[j]<2)throw new IllegalArgumentException("invalid category count");
            m+=categories[j]==0?2:categories[j]-1;m+=j;
        }return m;
    }
    static final class Mapping {
        final int[] categories,offset;final SemModel model;final double[] start;final List<String> labels;
        Mapping(int[] categories,SemModel model) {
            if(model==null||categories.length!=model.variables().size())throw new IllegalArgumentException("model dimension differs");
            this.categories=categories.clone();this.model=model;offset=new int[categories.length];
            int k=model.freeParameterCount();List<String> names=new ArrayList<>(model.freeParameterLabels());
            for(int j=0;j<categories.length;j++) {
                final int index=j;offset[j]=k;
                if(categories[j]>0) {
                    if(model.elements().stream().noneMatch(e->e.kind()==SemModel.Kind.VARIANCE&&e.first()==index&&e.fixed()&&e.start()>0)
                            ||model.elements().stream().anyMatch(e->e.kind()==SemModel.Kind.INTERCEPT&&e.first()==index&&!e.fixed()))
                        throw new IllegalArgumentException("fix ordinal residual scales and intercepts");
                    for(int c=1;c<categories[j];c++){names.add(model.variables().get(j)+"|t"+c);k++;}
                }
            }
            labels=List.copyOf(names);start=Arrays.copyOf(RamFit.initial(model),k);
            for(int j=0;j<categories.length;j++)if(categories[j]>0)start[offset[j]]=-.5*(categories[j]-2);
        }
        double[] natural(double[] point) {
            double[] result=point.clone();
            for(int j=0;j<model.freeParameterCount();j++)if(model.varianceParameter(j))result[j]=Math.exp(point[j]);
            for(int j=0;j<categories.length;j++)for(int c=1;c<categories[j]-1;c++)result[offset[j]+c]=result[offset[j]+c-1]+Math.exp(point[offset[j]+c]);
            return result;
        }
        double[] values(double[] point) {
            var d=RamFit.distribution(model,Arrays.copyOf(point,model.freeParameterCount()));double[] natural=natural(point);
            int p=categories.length,index=0;double[] values=new double[dimension(categories)],scale=new double[p];
            for(int j=0;j<p;j++) {
                scale[j]=categories[j]>0?Math.sqrt(d.covariance()[j*p+j]):1;
                if(categories[j]==0)values[index++]=d.mean()[j];
                else for(int c=0;c<categories[j]-1;c++)values[index++]=(natural[offset[j]+c]-d.mean()[j])/scale[j];
            }
            for(int i=0;i<p;i++)for(int j=0;j<=i;j++)if(i!=j||categories[i]==0)values[index++]=d.covariance()[i*p+j]/scale[i]/scale[j];
            return values;
        }
    }
    public record Result(List<String> parameterLabels,double[] estimates,double[] parameterCovariance,
            double dwlsStatistic,double wlsmvStatistic,int degreesOfFreedom,double scalingFactor,
            double shift,double pValue,boolean converged) {
        public Result{parameterLabels=List.copyOf(parameterLabels);estimates=estimates.clone();parameterCovariance=parameterCovariance.clone();}
        public double[] estimates(){return estimates.clone();}
        public double[] parameterCovariance(){return parameterCovariance.clone();}
    }
}
