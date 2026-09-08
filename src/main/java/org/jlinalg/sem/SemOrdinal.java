/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import jdistlib.Normal;

/** Joint ordinal probit RAM estimation by pairwise maximum likelihood (PML).
 * Thresholds are estimated jointly with the structural parameters. Inference uses
 * the case/cluster composite-score sandwich, not inverse composite information alone.
 * Observed residual variances must be fixed for response-scale identification. */
public final class SemOrdinal {
    private SemOrdinal() { }
    public static Result fit(int[][] data,int[] categoryCounts,SemModel model) {
        return fit(data,categoryCounts,model,SemOptions.defaults(),null);
    }
    /** Complete ordinal rows, coded 0..categoryCounts[j]-1. */
    public static Result fit(int[][] data,int[] categoryCounts,SemModel model,SemOptions options,int[] clusters) {
        Engine engine=new Engine(data,categoryCounts,model);
        if(options==null)throw new IllegalArgumentException("ordinal options required");
        int n=data.length,k=engine.start.length,ram=model.freeParameterCount();
        int[] ids=clusters==null?java.util.stream.IntStream.range(0,n).toArray():clusters.clone();
        if(ids.length!=n || Arrays.stream(ids).distinct().count()<2)throw new IllegalArgumentException("at least two independent ordinal clusters required");
        SemOptimizer.Optimum optimum=SemOptimizer.minimize(engine::at,engine.start,options.maximumEvaluations(),options.tolerance());
        double[] x=optimum.point();Evaluation evaluated=engine.evaluate(x,true);
        double[] bread=SemInformation.unavailable(k);
        if(engine.identified(x))try { bread=RamFit.informationInverse(RamFit.hessian(engine::at,x,n),k); }
        catch(IllegalArgumentException unavailable) { /* Preserve the likelihood; all inference stays unavailable. */ }
        double[] raw=SemInference.sandwich(evaluated.scores,ids,bread),jac=new double[k*k];
        double[] values=x.clone();
        for(int h=0;h<ram;h++){values[h]=model.varianceParameter(h)?Math.exp(x[h]):x[h];jac[h*k+h]=model.varianceParameter(h)?values[h]:1;}
        double[][] thresholds=engine.thresholds(x);
        List<String> labels=new ArrayList<>(model.freeParameterLabels());
        for(int j=0;j<categoryCounts.length;j++)for(int c=0;c<categoryCounts[j]-1;c++) {
            int h=engine.offset[j]+c;values[h]=thresholds[j][c];labels.add(model.variables().get(j)+"|t"+(c+1));
            jac[h*k+engine.offset[j]]=1;
            for(int a=1;a<=c;a++)jac[h*k+engine.offset[j]+a]=Math.exp(x[engine.offset[j]+a]);
        }
        double[] covariance=RamFit.mm(RamFit.mm(jac,raw,k),RamFit.transpose(jac,k),k);
        if(!SemInformation.finiteCovariance(covariance,k))covariance=SemInformation.unavailable(k);
        List<SemParameterEstimate> parameters=new ArrayList<>();
        for(int h=0;h<k;h++){double se=Math.sqrt(covariance[h*k+h]),z=values[h]/se;parameters.add(new SemParameterEstimate(labels.get(h),values[h],se,z,2*Normal.cumulative(-Math.abs(z),0,1,true,false)));}
        double[] implied=RamFit.distribution(model,Arrays.copyOf(x,ram)).covariance();int p=categoryCounts.length;
        double[] correlation=new double[p*p];for(int i=0;i<p;i++)for(int j=0;j<p;j++)correlation[i*p+j]=implied[i*p+j]/Math.sqrt(implied[i*p+i]*implied[j*p+j]);
        return new Result(parameters,covariance,thresholds,correlation,-n*evaluated.value.value(),n,optimum.evaluations(),optimum.converged(),optimum.scoreNorm());
    }

    public record Result(List<SemParameterEstimate> parameters,double[] parameterCovariance,
                         double[][] thresholds,double[] impliedCorrelation,double pairwiseLogLikelihood,
                         int observations,int functionEvaluations,boolean converged,double scoreNorm) {
        public Result { parameters=List.copyOf(parameters);parameterCovariance=parameterCovariance.clone();thresholds=copy(thresholds);impliedCorrelation=impliedCorrelation.clone(); }
        public double[] parameterCovariance(){return parameterCovariance.clone();}
        public double[][] thresholds(){return copy(thresholds);}
        public double[] impliedCorrelation(){return impliedCorrelation.clone();}
        /** False suppresses inference for every parameter, including thresholds. */
        public boolean informationAvailable(){return SemInformation.finiteCovariance(parameterCovariance,parameters.size());}
        public SemParameterEstimate parameter(String label){return parameters.stream().filter(v->v.label().equals(label)).findFirst().orElseThrow(()->new IllegalArgumentException("unknown ordinal parameter: "+label));}
        private static double[][] copy(double[][] a){return Arrays.stream(a).map(double[]::clone).toArray(double[][]::new);}
    }
    private record Pair(int first,int second,int[][] counts) { }
    private record Evaluation(SemOptimizer.Value value,double[][] scores) { }
    private static final class Engine {
        final int[][] data;final int[] categories,offset;final SemModel model;final double[] start;
        final List<Pair> pairs=new ArrayList<>();
        Engine(int[][] rows,int[] counts,SemModel specification) {
            if(rows==null || rows.length<3 || counts==null || specification==null || specification.hasMeanStructure())
                throw new IllegalArgumentException("ordinal SEM needs complete categories and a model without intercepts (thresholds determine location)");
            model=specification;categories=counts.clone();int p=model.variables().size();
            if(counts.length!=p)throw new IllegalArgumentException("ordinal categories must match observed variables");
            for(int j=0;j<p;j++) {
                final int index=j;
                if(counts[j]<2 || model.elements().stream().noneMatch(e->e.kind()==SemModel.Kind.VARIANCE && e.first()==index && e.fixed() && e.start()>0))
                    throw new IllegalArgumentException("fix every observed residual variance to a positive value to identify ordinal response scales");
            }
            data=Arrays.stream(rows).map(r->{if(r==null||r.length!=p)throw new IllegalArgumentException("invalid ordinal row");return r.clone();}).toArray(int[][]::new);
            int[][] marginal=new int[p][];offset=new int[p];int k=model.freeParameterCount();
            for(int j=0;j<p;j++){marginal[j]=new int[counts[j]];offset[j]=k;k+=counts[j]-1;}
            if(model.freeParameterCount()>p*(p-1)/2)
                throw new IllegalArgumentException("ordinal model has more free structural parameters than correlations; it is unidentified");
            for(int[] row:data)for(int j=0;j<p;j++){if(row[j]<0||row[j]>=counts[j])throw new IllegalArgumentException("ordinal category outside range");marginal[j][row[j]]++;}
            start=Arrays.copyOf(RamFit.initial(model),k);double[] covariance=RamFit.distribution(model,RamFit.initial(model)).covariance();
            for(int j=0;j<p;j++) {
                int cumulative=0;double previous=0;
                for(int c=0;c<counts[j];c++) {
                    if(marginal[j][c]==0)throw new IllegalArgumentException("empty marginal ordinal category; collapse or remove it explicitly");
                    cumulative+=marginal[j][c];if(c==counts[j]-1)break;
                    double threshold=Math.sqrt(covariance[j*p+j])*Normal.quantile((double)cumulative/data.length,0,1,true,false);
                    start[offset[j]+c]=c==0?threshold:Math.log(threshold-previous);previous=threshold;
                }
            }
            for(int i=0;i<p;i++)for(int j=0;j<i;j++) {
                int[][] table=new int[counts[i]][counts[j]];for(int[] row:data)table[row[i]][row[j]]++;
                pairs.add(new Pair(i,j,table));
            }
        }
        double[][] thresholds(double[] x) {
            double[][] t=new double[categories.length][];
            for(int j=0;j<t.length;j++){t[j]=new double[categories[j]-1];t[j][0]=x[offset[j]];for(int c=1;c<t[j].length;c++)t[j][c]=t[j][c-1]+Math.exp(x[offset[j]+c]);}
            return t;
        }
        SemOptimizer.Value at(double[] x) {
            try{return evaluate(x,false).value;}catch(IllegalArgumentException e){return new SemOptimizer.Value(Double.POSITIVE_INFINITY,null);}
        }
        boolean identified(double[] x) {
            int p=categories.length,k=x.length,ram=model.freeParameterCount();
            RamFit.Distribution d=RamFit.distribution(model,Arrays.copyOf(x,ram));double[][] t=thresholds(x);
            double[][] columns=new double[k][k-ram+p*(p-1)/2];int row=0;
            for(int j=0;j<p;j++)for(int c=0;c<t[j].length;c++,row++) {
                double variance=d.covariance()[j*p+j],sd=Math.sqrt(variance);
                for(int h=0;h<ram;h++)columns[h][row]=-.5*t[j][c]/sd*d.dCovariance()[h][j*p+j]/variance;
                columns[offset[j]][row]=1/sd;
                for(int h=1;h<=c;h++)columns[offset[j]+h][row]=Math.exp(x[offset[j]+h])/sd;
            }
            for(int i=0;i<p;i++)for(int j=0;j<i;j++,row++) {
                double vi=d.covariance()[i*p+i],vj=d.covariance()[j*p+j],scale=Math.sqrt(vi*vj),rho=d.covariance()[i*p+j]/scale;
                for(int h=0;h<ram;h++)columns[h][row]=d.dCovariance()[h][i*p+j]/scale-.5*rho*(d.dCovariance()[h][i*p+i]/vi+d.dCovariance()[h][j*p+j]/vj);
            }
            return SemInformation.fullColumnRank(columns);
        }
        Evaluation evaluate(double[] x,boolean caseScores) {
            int p=categories.length,k=x.length,ram=model.freeParameterCount(),n=data.length;
            RamFit.Distribution d=RamFit.distribution(model,Arrays.copyOf(x,ram));double[][] t=thresholds(x);double[] sd=new double[p];
            for(int j=0;j<p;j++)sd[j]=Math.sqrt(d.covariance()[j*p+j]);
            double value=0;double[] gradient=new double[k];double[][] scores=caseScores?new double[n][k]:null;
            for(Pair pair:pairs) {
                int i=pair.first,j=pair.second;double rho=d.covariance()[i*p+j]/(sd[i]*sd[j]);
                if(Math.abs(rho)>=.9999)throw new IllegalArgumentException("ordinal correlation too near its boundary for reliable quadrature");
                double[] dr=new double[ram];
                for(int h=0;h<ram;h++)dr[h]=d.dCovariance()[h][i*p+j]/(sd[i]*sd[j])-.5*rho*(d.dCovariance()[h][i*p+i]/(sd[i]*sd[i])+d.dCovariance()[h][j*p+j]/(sd[j]*sd[j]));
                double[][][] cellScores=caseScores?new double[categories[i]][categories[j]][]:null;
                for(int a=0;a<categories[i];a++)for(int b=0;b<categories[j];b++) {
                    if(pair.counts[a][b]==0)continue;
                    double loI=a==0?Double.NEGATIVE_INFINITY:t[i][a-1]/sd[i],hiI=a==t[i].length?Double.POSITIVE_INFINITY:t[i][a]/sd[i];
                    double loJ=b==0?Double.NEGATIVE_INFINITY:t[j][b-1]/sd[j],hiJ=b==t[j].length?Double.POSITIVE_INFINITY:t[j][b]/sd[j];
                    double prob=SemNormalRectangle.rectangle(loI,hiI,loJ,hiJ,rho);
                    if(!(prob>0) || !Double.isFinite(prob))throw new IllegalArgumentException("ordinal cell probability lost numerical precision");
                    double[] dp=new double[k];double rhoDerivative=density(hiI,hiJ,rho)-density(loI,hiJ,rho)-density(hiI,loJ,rho)+density(loI,loJ,rho);
                    for(int h=0;h<ram;h++)dp[h]=rhoDerivative*dr[h];
                    boundary(dp,x,d,i,a,hiI,conditional(hiI,loJ,hiJ,rho),sd[i]);
                    boundary(dp,x,d,i,a-1,loI,-conditional(loI,loJ,hiJ,rho),sd[i]);
                    boundary(dp,x,d,j,b,hiJ,conditional(hiJ,loI,hiI,rho),sd[j]);
                    boundary(dp,x,d,j,b-1,loJ,-conditional(loJ,loI,hiI,rho),sd[j]);
                    double weight=(double)pair.counts[a][b]/n;value-=weight*Math.log(prob);
                    for(int h=0;h<k;h++){dp[h]/=prob;gradient[h]-=weight*dp[h];}
                    if(caseScores)cellScores[a][b]=dp;
                }
                if(caseScores)for(int r=0;r<n;r++)for(int h=0;h<k;h++)scores[r][h]+=cellScores[data[r][i]][data[r][j]][h];
            }
            return new Evaluation(new SemOptimizer.Value(value,gradient),scores);
        }
        void boundary(double[] dp,double[] x,RamFit.Distribution d,int j,int c,double z,double derivative,double sd) {
            if(c<0||c>=categories[j]-1)return;
            for(int h=0;h<model.freeParameterCount();h++)dp[h]-=.5*derivative*z*d.dCovariance()[h][j*categories.length+j]/(sd*sd);
            dp[offset[j]]+=derivative/sd;for(int h=1;h<=c;h++)dp[offset[j]+h]+=derivative*Math.exp(x[offset[j]+h])/sd;
        }
    }
    private static double phi(double x){return Double.isFinite(x)?Math.exp(-.5*x*x)/Math.sqrt(2*Math.PI):0;}
    private static double cdf(double x){return Normal.cumulative(x,0,1,true,false);}
    private static double conditional(double h,double lo,double hi,double r){return Double.isFinite(h)?phi(h)*SemNormalRectangle.interval((lo-r*h)/Math.sqrt(1-r*r),(hi-r*h)/Math.sqrt(1-r*r)):0;}
    private static double density(double h,double k,double r){return Double.isFinite(h)&&Double.isFinite(k)?Math.exp(-.5*((h-r*k)*(h-r*k)/((1-r)*(1+r))+k*k))/(2*Math.PI*Math.sqrt((1-r)*(1+r))):0;}
    /** Error-controlled Plackett integral after a sine change of variable. */
    static double bvn(double h,double k,double r) {
        return SemNormalRectangle.cdf(h,k,r);
    }
}
