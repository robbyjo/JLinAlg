/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import jdistlib.Normal;

/** Joint Gaussian/ordinal-probit observed-data likelihood. Arbitrary ignorable
 * MAR patterns are marginalized, including singleton rows. Category count zero
 * denotes a continuous response; ordinal responses use integer codes 0..K-1.
 * NaN alone denotes missingness. Up to four ordinal responses may be observed
 * in any row. Numerical rectangle and information failures remain explicit. */
public final class SemMixed {
    private SemMixed() { }
    public static Result fit(double[][] data,int[] categoryCounts,SemModel model) {
        return fit(data,categoryCounts,model,SemOptions.defaults());
    }
    public static Result fit(double[][] data,int[] categoryCounts,SemModel model,SemOptions options) {
        Engine engine=new Engine(data,categoryCounts,model);
        var optimum=SemOptimizer.minimize(engine::at,engine.start,options.maximumEvaluations(),options.tolerance());
        double[] covariance=SemInformation.unavailable(engine.start.length);
        if(optimum.converged())try {
            covariance=RamFit.informationInverse(RamFit.hessian(engine::at,optimum.point(),engine.data.length),engine.start.length);
        }catch(IllegalArgumentException unavailable){ /* Preserve converged estimates, suppress information. */ }
        return new Result(engine,optimum.point(),covariance,-engine.data.length*engine.value(optimum.point()),optimum.converged());
    }
    /** Efficient likelihood score tests for omitted/fixed structural parameters,
     * with thresholds and all fitted paths projected out as nuisance parameters. */
    public static List<SemInference.ModificationIndex> modificationIndices(Result fit,SemInference.Modification... candidates) {
        if(!fit.converged||!fit.informationAvailable())throw new IllegalArgumentException("identified converged mixed fit required");
        List<SemInference.ModificationIndex> result=new ArrayList<>();Engine old=fit.engine;
        for(var candidate:candidates) {
            SemModel.Kind kind=switch(candidate.kind()) {
                case "regression"->SemModel.Kind.REGRESSION;case "covariance"->SemModel.Kind.COVARIANCE;
                default->throw new IllegalArgumentException("mixed modification supports paths and covariances");
            };
            String label="__mixed_modification__";while(old.labels.contains(label))label+="_";
            SemModel expanded=old.model.freeElement(kind,candidate.first(),candidate.second(),label);
            Engine engine=new Engine(old.data,old.categories,expanded);double[] point=engine.start.clone();
            for(int j=0;j<old.labels.size();j++)point[engine.labels.indexOf(old.labels.get(j))]=fit.point[j];
            int k=point.length,index=engine.labels.indexOf(label),n=engine.data.length;
            double[] information=RamFit.hessian(engine::at,point,n),score=engine.at(point).gradient();
            int[] nuisance=java.util.stream.IntStream.range(0,k).filter(i->i!=index).toArray();
            double[] inverse=RamFit.informationInverse(RamFit.sub(information,nuisance,k),k-1),cross=new double[k-1],s=new double[k-1];
            for(int j=0;j<k-1;j++){cross[j]=information[index*k+nuisance[j]];s[j]=-n*score[nuisance[j]];}
            double efficient=information[index*k+index]-SemOptimizer.dot(cross,RamFit.mv(inverse,cross));
            if(!(efficient>1e-9*Math.abs(information[index*k+index])))throw new IllegalArgumentException("unidentified modification");
            double scoreValue=-n*score[index]-SemOptimizer.dot(cross,RamFit.mv(inverse,s));
            double statistic=scoreValue*scoreValue/efficient;
            result.add(new SemInference.ModificationIndex(candidate,statistic,scoreValue/efficient,jdistlib.ChiSquare.cumulative(statistic,1,false,false)));
        }
        return List.copyOf(result);
    }
    static final class Engine {
        final double[][] data;final int[] categories,offset;final SemModel model;
        final double[] start;final List<String> labels;
        final double[][] uniqueRows;final int[] frequencies;
        Engine(double[][] rows,int[] counts,SemModel specification) {
            if(rows==null||rows.length<3||counts==null||specification==null||counts.length!=specification.variables().size())
                throw new IllegalArgumentException("mixed SEM requires data, categories and a matching model");
            categories=counts.clone();model=specification;int p=counts.length,k=model.freeParameterCount();
            offset=new int[p];List<String> names=new ArrayList<>(model.freeParameterLabels());
            for(int j=0;j<p;j++) {
                if(counts[j]!=0&&counts[j]<2)throw new IllegalArgumentException("category counts must be zero or at least two");
                offset[j]=k;k+=Math.max(0,counts[j]-1);
                if(counts[j]>0) {
                    final int col=j;
                    if(model.elements().stream().noneMatch(e->e.kind()==SemModel.Kind.VARIANCE&&e.first()==col&&e.fixed()&&e.start()>0))
                        throw new IllegalArgumentException("ordinal residual variances must be fixed positive for scale identification");
                    if(model.elements().stream().anyMatch(e->e.kind()==SemModel.Kind.INTERCEPT&&e.first()==col&&!e.fixed()))
                        throw new IllegalArgumentException("ordinal intercepts must be fixed; thresholds determine location");
                    for(int c=1;c<counts[j];c++)names.add(model.variables().get(j)+"|t"+c);
                }
            }
            if(names.stream().distinct().count()!=names.size())throw new IllegalArgumentException("structural and threshold labels collide");
            labels=List.copyOf(names);List<double[]> retained=new ArrayList<>();
            int[][] marginal=new int[p][];for(int j=0;j<p;j++)marginal[j]=new int[Math.max(1,counts[j])];
            for(double[] row:rows) {
                if(row==null||row.length!=p)throw new IllegalArgumentException("invalid mixed response row");
                int seen=0,ordinal=0;
                for(int j=0;j<p;j++) {
                    double value=row[j];if(Double.isNaN(value))continue;
                    if(!Double.isFinite(value)||(counts[j]>0&&(value!=Math.rint(value)||value<0||value>=counts[j])))
                        throw new IllegalArgumentException("invalid mixed response; use NaN for missing");
                    marginal[j][counts[j]>0?(int)value:0]++;seen++;if(counts[j]>0)ordinal++;
                }
                if(ordinal>4)throw new IllegalArgumentException("at most four observed ordinal responses per row supported");
                if(seen>0)retained.add(row.clone());
            }
            data=retained.toArray(double[][]::new);if(data.length<3)throw new IllegalArgumentException("too few informative rows");
            java.util.Map<String,Integer> ids=new java.util.LinkedHashMap<>();List<double[]> unique=new ArrayList<>();List<Integer> frequency=new ArrayList<>();
            for(double[] row:data) {
                String key=Arrays.toString(row);Integer id=ids.get(key);
                if(id==null){ids.put(key,unique.size());unique.add(row);frequency.add(1);}
                else frequency.set(id,frequency.get(id)+1);
            }
            uniqueRows=unique.toArray(double[][]::new);frequencies=frequency.stream().mapToInt(Integer::intValue).toArray();
            start=Arrays.copyOf(RamFit.initial(model),k);var distribution=RamFit.distribution(model,RamFit.initial(model));
            for(int j=0;j<p;j++) {
                int total=Arrays.stream(marginal[j]).sum(),cumulative=0;double previous=0;
                if(total<2)throw new IllegalArgumentException("each response needs observations");
                for(int c=0;c<counts[j];c++) {
                    if(marginal[j][c]==0)throw new IllegalArgumentException("empty ordinal category; collapse explicitly");
                    cumulative+=marginal[j][c];if(c==counts[j]-1)continue;
                    double t=distribution.mean()[j]+Math.sqrt(distribution.covariance()[j*p+j])*Normal.quantile((double)cumulative/total,0,1,true,false);
                    start[offset[j]+c]=c==0?t:Math.log(t-previous);previous=t;
                }
            }
        }
        double[][] thresholds(double[] point) {
            double[][] t=new double[categories.length][];
            for(int j=0;j<t.length;j++) {
                t[j]=new double[Math.max(0,categories[j]-1)];
                for(int c=0;c<t[j].length;c++)t[j][c]=c==0?point[offset[j]]:t[j][c-1]+Math.exp(point[offset[j]+c]);
            }return t;
        }
        SemOptimizer.Value at(double[] point) {
            try {
                double center=value(point);double[] gradient=new double[point.length];
                for(int j=0;j<point.length;j++) {
                    double h=2e-5*(1+Math.abs(point[j]));double[] x=point.clone();
                    x[j]+=h;double above=value(x);x[j]-=2*h;double below=value(x);gradient[j]=(above-below)/(2*h);
                }
                return new SemOptimizer.Value(center,gradient);
            }catch(IllegalArgumentException failure){return new SemOptimizer.Value(Double.POSITIVE_INFINITY,null);}
        }
        double value(double[] point) {
            var d=RamFit.distribution(model,Arrays.copyOf(point,model.freeParameterCount()));double[][] t=thresholds(point);
            double total=0;for(int i=0;i<uniqueRows.length;i++)total-=frequencies[i]*logDensity(uniqueRows[i],d,t)/data.length;return total;
        }
        double logDensity(double[] row,RamFit.Distribution d,double[][] thresholds) {
            int p=categories.length;
            int[] continuous=java.util.stream.IntStream.range(0,p).filter(j->categories[j]==0&&Double.isFinite(row[j])).toArray();
            int[] ordinal=java.util.stream.IntStream.range(0,p).filter(j->categories[j]>0&&Double.isFinite(row[j])).toArray();
            int c=continuous.length,o=ordinal.length;double log=0;
            double[] delta=new double[c],inverse=new double[c*c];
            if(c>0) {
                double[] l=RamFit.chol(RamFit.sub(d.covariance(),continuous,p),c);inverse=RamFit.cholInverse(l,c);
                for(int j=0;j<c;j++){delta[j]=row[continuous[j]]-d.mean()[continuous[j]];log-=Math.log(l[j*c+j])+.5*Math.log(2*Math.PI);}
                log-=.5*SemOptimizer.dot(delta,RamFit.mv(inverse,delta));
            }
            double[] lo=new double[o],hi=new double[o],cov=RamFit.sub(d.covariance(),ordinal,p),adjusted=RamFit.mv(inverse,delta);
            for(int a=0;a<o;a++) {
                int j=ordinal[a],category=(int)row[j];double mean=d.mean()[j];
                for(int k=0;k<c;k++)mean+=d.covariance()[j*p+continuous[k]]*adjusted[k];
                lo[a]=(category==0?Double.NEGATIVE_INFINITY:thresholds[j][category-1])-mean;
                hi[a]=(category==categories[j]-1?Double.POSITIVE_INFINITY:thresholds[j][category])-mean;
                for(int b=0;b<o;b++)for(int k=0;k<c;k++)for(int l=0;l<c;l++)
                    cov[a*o+b]-=d.covariance()[j*p+continuous[k]]*inverse[k*c+l]*d.covariance()[continuous[l]*p+ordinal[b]];
            }
            double probability=SemRectangle.probability(lo,hi,cov);
            if(!(probability>0)||!Double.isFinite(probability))throw new IllegalArgumentException("ordinal probability unresolved");
            return log+Math.log(probability);
        }
        double[] natural(double[] point) {
            double[] values=point.clone();for(int j=0;j<model.freeParameterCount();j++)if(model.varianceParameter(j))values[j]=Math.exp(point[j]);
            double[][] t=thresholds(point);for(int j=0;j<t.length;j++)System.arraycopy(t[j],0,values,offset[j],t[j].length);return values;
        }
    }
    static double[] transformedCovariance(double[] point,double[] covariance,java.util.function.Function<double[],double[]> transform) {
        double[] values=transform.apply(point);int k=point.length,m=values.length;double[][] jacobian=new double[m][k];
        for(int j=0;j<k;j++) {
            double h=1e-5*(1+Math.abs(point[j]));double[] x=point.clone();x[j]+=h;double[] above=transform.apply(x);x[j]-=2*h;double[] below=transform.apply(x);
            for(int i=0;i<m;i++)jacobian[i][j]=(above[i]-below[i])/(2*h);
        }
        double[] result=new double[m*m];for(int i=0;i<m;i++)for(int j=0;j<m;j++)for(int a=0;a<k;a++)for(int b=0;b<k;b++)result[i*m+j]+=jacobian[i][a]*covariance[a*k+b]*jacobian[j][b];return result;
    }
    public static final class Result {
        final Engine engine;final double[] point,rawCovariance;
        private final double logLikelihood;private final boolean converged;
        Result(Engine engine,double[] point,double[] covariance,double ll,boolean converged) {
            this.engine=engine;this.point=point.clone();rawCovariance=covariance.clone();logLikelihood=ll;this.converged=converged;
        }
        public List<String> parameterLabels(){return engine.labels;}
        public double[] estimates(){return engine.natural(point);}
        public double[] parameterCovariance(){return transformedCovariance(point,rawCovariance,engine::natural);}
        public double[][] thresholds(){return engine.thresholds(point);}
        public double logLikelihood(){return logLikelihood;}
        public boolean converged(){return converged;}
        public boolean informationAvailable(){return converged&&SemInformation.finiteCovariance(rawCovariance,point.length);}
    }
}
