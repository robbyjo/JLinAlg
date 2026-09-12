/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;

/** Independent-group joint likelihood with explicit cross-group equality labels.
 * Empty shared labels fit configural models. Share loading labels for metric
 * invariance, continuous intercept/ordinal threshold labels for scalar
 * invariance, and residual variance labels for strict invariance. Identification
 * and comparable ordinal latent scales must be specified in the group models.
 * Partial invariance is expressed by sharing only selected labels. */
public final class SemMultigroup {
    private SemMultigroup() { }
    public static Result fit(List<double[][]> groups,int[] categoryCounts,List<SemModel> models,
            Set<String> sharedLabels,SemOptions options) {
        if(groups==null||models==null||sharedLabels==null||options==null||groups.size()<2||groups.size()!=models.size())
            throw new IllegalArgumentException("at least two independent groups and matching models required");
        List<SemMixed.Engine> engines=new ArrayList<>();
        for(int g=0;g<groups.size();g++) {
            if(!models.get(g).variables().equals(models.get(0).variables()))throw new IllegalArgumentException("group response names/order differ");
            engines.add(new SemMixed.Engine(groups.get(g),categoryCounts,models.get(g)));
        }
        for(String label:sharedLabels)if(engines.stream().filter(e->e.labels.contains(label)).count()<2)
            throw new IllegalArgumentException("shared label must occur in at least two groups: "+label);
        LinkedHashMap<String,Integer> indices=new LinkedHashMap<>();List<Boolean> variance=new ArrayList<>();
        List<Double> starts=new ArrayList<>();int[][] map=new int[groups.size()][];
        for(int g=0;g<groups.size();g++) {
            var engine=engines.get(g);map[g]=new int[engine.start.length];double[] natural=engine.natural(engine.start);
            for(int j=0;j<engine.start.length;j++) {
                String label=engine.labels.get(j),key=sharedLabels.contains(label)?label:"group"+g+":"+label;
                boolean v=j<engine.model.freeParameterCount()&&engine.model.varianceParameter(j);
                if(!indices.containsKey(key)){indices.put(key,indices.size());variance.add(v);starts.add(v?Math.log(natural[j]):natural[j]);}
                int index=indices.get(key);if(variance.get(index)!=v)throw new IllegalArgumentException("shared label mixes variance and nonvariance parameters");
                map[g][j]=index;
            }
        }
        double[] start=starts.stream().mapToDouble(Double::doubleValue).toArray();int n=engines.stream().mapToInt(e->e.data.length).sum();
        java.util.function.ToDoubleFunction<double[]> value=point->{
            double total=0;for(int g=0;g<engines.size();g++)total+=engines.get(g).data.length*engines.get(g).value(local(point,map[g],engines.get(g)))/n;return total;
        };
        SemOptimizer.Objective objective=point->{
            try {
                double center=value.applyAsDouble(point);double[] gradient=new double[point.length];
                for(int j=0;j<point.length;j++) {
                    double h=2e-5*(1+Math.abs(point[j]));double[] x=point.clone();x[j]+=h;double above=value.applyAsDouble(x);x[j]-=2*h;
                    gradient[j]=(above-value.applyAsDouble(x))/(2*h);
                }return new SemOptimizer.Value(center,gradient);
            }catch(IllegalArgumentException failure){return new SemOptimizer.Value(Double.POSITIVE_INFINITY,null);}
        };
        var optimum=SemOptimizer.minimize(objective,start,options.maximumEvaluations(),options.tolerance());
        double[] covariance=SemInformation.unavailable(start.length);
        if(optimum.converged())try{covariance=RamFit.informationInverse(RamFit.hessian(objective,optimum.point(),n),start.length);}
        catch(IllegalArgumentException unavailable){ /* Suppress unresolved covariance. */ }
        double[] estimates=optimum.point().clone();
        for(int i=0;i<estimates.length;i++)if(variance.get(i))estimates[i]=Math.exp(estimates[i]);
        for(int i=0;i<estimates.length;i++)for(int j=0;j<estimates.length;j++)covariance[i*estimates.length+j]*=(variance.get(i)?estimates[i]:1)*(variance.get(j)?estimates[j]:1);
        return new Result(engines,Set.copyOf(sharedLabels),List.copyOf(indices.keySet()),estimates,covariance,
            -n*value.applyAsDouble(optimum.point()),optimum.converged());
    }
    private static double[] local(double[] point,int[] map,SemMixed.Engine engine) {
        double[] result=new double[map.length];for(int j=0;j<map.length;j++)result[j]=point[map[j]];
        // Global thresholds are actual thresholds, not logarithmic increments:
        // sharing t2 must not accidentally constrain t2-t1 under partial invariance.
        for(int j=0;j<engine.categories.length;j++)for(int c=engine.categories[j]-2;c>=1;c--) {
            int index=engine.offset[j]+c;double gap=result[index]-result[index-1];
            if(!(gap>0))throw new IllegalArgumentException("thresholds must remain strictly ordered");result[index]=Math.log(gap);
        }
        return result;
    }
    /** Conventional likelihood-ratio test for nested equality sets. Requires the
     * same data and group models, regular interior identification and a correct
     * joint distribution; no naive subtraction of WLSMV statistics is used. */
    public static Comparison compare(Result unconstrained,Result constrained) {
        if(!unconstrained.informationAvailable()||!constrained.informationAvailable()
                ||!constrained.shared.containsAll(unconstrained.shared)||unconstrained.engines.size()!=constrained.engines.size())
            throw new IllegalArgumentException("identified converged nested fits required");
        for(int g=0;g<unconstrained.engines.size();g++) {
            var a=unconstrained.engines.get(g);var b=constrained.engines.get(g);
            if(!Arrays.deepEquals(a.data,b.data)||!Arrays.equals(a.categories,b.categories)||!a.model.elements().equals(b.model.elements()))
                throw new IllegalArgumentException("comparison must use identical data and group models");
        }
        int df=unconstrained.estimates.length-constrained.estimates.length;
        double statistic=2*(unconstrained.logLikelihood-constrained.logLikelihood);
        if(df<=0||statistic< -1e-6)throw new IllegalArgumentException("invalid nesting or unresolved optimization");
        statistic=Math.max(0,statistic);return new Comparison(statistic,df,jdistlib.ChiSquare.cumulative(statistic,df,false,false));
    }
    public record Comparison(double chiSquare,int degreesOfFreedom,double pValue) { }
    public static final class Result {
        private final List<SemMixed.Engine> engines;private final Set<String> shared;private final List<String> labels;
        private final double[] estimates,covariance;private final double logLikelihood;private final boolean converged;
        private Result(List<SemMixed.Engine> engines,Set<String> shared,List<String> labels,double[] estimates,double[] covariance,double ll,boolean converged) {
            this.engines=List.copyOf(engines);this.shared=shared;this.labels=labels;this.estimates=estimates.clone();this.covariance=covariance.clone();logLikelihood=ll;this.converged=converged;
        }
        public List<String> parameterLabels(){return labels;}
        public double[] estimates(){return estimates.clone();}
        public double[] parameterCovariance(){return covariance.clone();}
        public double logLikelihood(){return logLikelihood;}
        public boolean converged(){return converged;}
        public boolean informationAvailable(){return converged&&SemInformation.finiteCovariance(covariance,estimates.length);}
    }
}
