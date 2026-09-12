/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import jdistlib.Normal;

/** Deterministic SEM inference utilities that operate on exported scores. */
public final class SemInference {
    private SemInference() { }

    /** Product delta inference using every covariance term from the fitted model. */
    public static IndirectEffect indirect(SemFitResult fit, double confidenceLevel, String... paths) {
        return indirect(fit, fit.parameterCovariance(), confidenceLevel, paths);
    }
    /** Supply a robust covariance here to obtain robust product delta inference. */
    public static IndirectEffect indirect(SemFitResult fit, double[] covariance,
            double confidenceLevel, String... paths) {
        if(paths==null || paths.length<2 || !(confidenceLevel>0 && confidenceLevel<1))
            throw new IllegalArgumentException("at least two paths and a valid confidence level required");
        int k=fit.parameters().size();if(covariance.length!=k*k)throw new IllegalArgumentException("invalid parameter covariance");
        int[] indices=new int[paths.length];double[] values=new double[paths.length];double product=1;
        for(int i=0;i<paths.length;i++) {
            indices[i]=fit.parameters().indexOf(fit.parameter(paths[i]));values[i]=fit.parameter(paths[i]).estimate();product*=values[i];
        }
        double[] gradient=new double[k];
        for(int i=0;i<paths.length;i++){double v=1;for(int j=0;j<paths.length;j++)if(i!=j)v*=values[j];gradient[indices[i]]+=v;}
        double variance=SemOptimizer.dot(gradient,RamFit.mv(covariance,gradient));
        double se=Math.sqrt(variance),z=product/se,critical=Normal.quantile(.5+confidenceLevel/2,0,1,true,false);
        return new IndirectEffect(product,se,z,product-critical*se,product+critical*se,2*Normal.cumulative(-Math.abs(z),0,1,true,false));
    }

    /** Model-derived case-score sandwich; expected bread for complete ML, observed bread for FIML. */
    public static RobustResult robust(SemFitResult fit) { return robust(fit,null); }
    /** CR0 cluster sandwich. Cluster IDs follow the informative rows retained by the fit. */
    public static RobustResult robust(SemFitResult fit,int[] clusters) {
        if(!fit.converged()||!fit.informationAvailable())throw new IllegalArgumentException("robust inference requires a converged identified SEM fit");
        RamFit.State s=fit.state();
        if(s.data()==null)throw new IllegalArgumentException("case data are required for robust inference");
        int n=s.observations(),k=s.point().length;
        if(k==0)throw new IllegalArgumentException("robust inference requires free parameters");
        int[] ids=clusters==null?java.util.stream.IntStream.range(0,n).toArray():clusters.clone();
        if(ids.length!=n)throw new IllegalArgumentException("cluster IDs must match retained rows");
        int groups=(int)java.util.Arrays.stream(ids).distinct().count();
        if(groups<2)throw new IllegalArgumentException("at least two independent clusters required");
        RamFit.Distribution d=RamFit.distribution(s.model(),s.point());
        double[] information=s.missing()?RamFit.hessian(x->RamFit.evaluate(s.model(),x,s.patterns(),n),s.point(),n)
            :RamFit.expectedInformation(d,s.patterns());
        double[] bread=RamFit.informationInverse(information,k);
        double[] raw=sandwich(RamFit.scores(d,s.data()),ids,bread);
        double scaling;
        try{scaling=scaling(s,d,ids,bread,fit.degreesOfFreedom());}
        catch(IllegalArgumentException unavailable){scaling=Double.NaN;}
        if(!(scaling>0)||!Double.isFinite(scaling))scaling=Double.NaN;
        double chi=fit.chiSquare()/scaling;
        return new RobustResult(RamFit.naturalCovariance(raw,s.model(),s.point()),scaling,chi,
            fit.degreesOfFreedom(),fit.degreesOfFreedom()>0?jdistlib.ChiSquare.cumulative(chi,fit.degreesOfFreedom(),false,false):Double.NaN,groups);
    }

    private static double scaling(RamFit.State s,RamFit.Distribution d,int[] ids,double[] bread,int df) {
        if(df<=0)return Double.NaN;
        int p=d.mean().length,k=s.point().length,q=p*(p+1)/2+((s.model().hasMeanStructure()||s.missing())?p:0);
        double[][] dm=new double[q][p],dc=new double[q][p*p];double[][] jacobian=new double[q][k];int h=0;
        if(s.model().hasMeanStructure()||s.missing())for(int i=0;i<p;i++,h++) {
            dm[h][i]=1;for(int a=0;a<k;a++)jacobian[h][a]=d.dMean()[a][i];
        }
        for(int i=0;i<p;i++)for(int j=0;j<=i;j++,h++) {
            dc[h][i*p+j]=dc[h][j*p+i]=1;
            for(int a=0;a<k;a++)jacobian[h][a]=d.dCovariance()[a][i*p+j];
        }
        RamFit.Distribution saturated=new RamFit.Distribution(d.mean(),d.covariance(),dm,dc);
        // Under MAR, pattern membership can depend on observed responses.
        // Use the observed H1 information, including all saturated means, in
        // the LR curvature projection. Expected complete-data information is
        // not an appropriate replacement for this missing-data curvature.
        double[] h1=s.missing()?RamFit.hessian(point->{
            double[] mean=d.mean().clone(),cov=d.covariance().clone();
            for(int a=0;a<point.length;a++) {
                for(int i=0;i<p;i++)mean[i]+=dm[a][i]*point[a];
                for(int i=0;i<p*p;i++)cov[i]+=dc[a][i]*point[a];
            }
            return RamFit.evaluate(new RamFit.Distribution(mean,cov,dm,dc),s.patterns(),s.observations());
        },new double[q],s.observations()):RamFit.expectedInformation(saturated,s.patterns());
        double[] u=RamFit.informationInverse(h1,q);
        for(int i=0;i<q;i++)for(int j=0;j<q;j++)for(int a=0;a<k;a++)for(int b=0;b<k;b++)
            u[i*q+j]-=jacobian[i][a]*bread[a*k+b]*jacobian[j][b];
        double[][] scores=RamFit.scores(saturated,s.data());double[] mean=new double[q];
        for(double[] row:scores)for(int i=0;i<q;i++)mean[i]+=row[i]/scores.length;
        for(double[] row:scores)for(int i=0;i<q;i++)row[i]-=mean[i];
        java.util.Map<Integer,double[]> sums=new java.util.LinkedHashMap<>();
        for(int r=0;r<scores.length;r++){double[] sum=sums.computeIfAbsent(ids[r],unused->new double[q]);for(int i=0;i<q;i++)sum[i]+=scores[r][i];}
        double trace=0;for(double[] sum:sums.values())trace+=SemOptimizer.dot(sum,RamFit.mv(u,sum));
        return trace/df;
    }

    /** Efficient one-df score tests after projecting out all fitted nuisance parameters. */
    public static java.util.List<ModificationIndex> modificationIndices(SemFitResult fit) {
        SemModel model=fit.state().model();java.util.List<Modification> candidates=new java.util.ArrayList<>();
        // The default scan covers omitted observed disturbance covariances. Directed
        // paths, fixed paths and latent covariances can be supplied explicitly.
        for(int i=0;i<model.variables().size();i++)for(int j=0;j<i;j++) {
            final int a=i,b=j;
            if(model.elements().stream().noneMatch(e->e.kind()==SemModel.Kind.COVARIANCE &&
                    ((e.first()==a && e.second()==b)||(e.first()==b && e.second()==a))))
                candidates.add(Modification.covariance(model.variables().get(j),model.variables().get(i)));
        }
        return modificationIndices(fit,candidates.toArray(Modification[]::new));
    }
    /** Candidate-specific score tests; fixed factor-scale releases can be unidentified. */
    public static java.util.List<ModificationIndex> modificationIndices(SemFitResult fit, Modification... candidates) {
        RamFit.State s=fit.state();java.util.List<ModificationIndex> result=new java.util.ArrayList<>();
        for(Modification c:candidates) {
            SemModel.Kind kind=switch(c.kind()) {
                case "regression" -> SemModel.Kind.REGRESSION;
                case "covariance" -> SemModel.Kind.COVARIANCE;
                case "intercept" -> SemModel.Kind.INTERCEPT;
                default -> throw new IllegalArgumentException("unknown modification kind");
            };
            String label="__modification__";while(s.model().freeParameterLabels().contains(label))label+="_";
            SemModel expanded=s.model().freeElement(kind,c.first(),c.second(),label);
            double[] x=RamFit.initial(expanded);int k=x.length,j=expanded.freeIndex(label);
            for(String name:s.model().freeParameterLabels())x[expanded.freeIndex(name)]=s.point()[s.model().freeIndex(name)];
            double[] grad=RamFit.evaluate(expanded,x,s.patterns(),s.observations()).gradient();
            double[] info=RamFit.expectedInformation(RamFit.distribution(expanded,x),s.patterns());
            int[] nuisance=java.util.stream.IntStream.range(0,k).filter(i->i!=j).toArray();
            double[] nuisanceInverse=RamFit.informationInverse(RamFit.sub(info,nuisance,k),k-1),cross=new double[k-1],score=new double[k-1];
            for(int a=0;a<k-1;a++){cross[a]=info[nuisance[a]*k+j];score[a]=-s.observations()*grad[nuisance[a]];}
            double efficientInfo=info[j*k+j]-SemOptimizer.dot(cross,RamFit.mv(nuisanceInverse,cross));
            double efficientScore=-s.observations()*grad[j]-SemOptimizer.dot(cross,RamFit.mv(nuisanceInverse,score));
            boolean identified=efficientInfo>1e-10*Math.abs(info[j*k+j]);
            double mi=identified?efficientScore*efficientScore/efficientInfo:Double.NaN;
            result.add(new ModificationIndex(c,mi,identified?efficientScore/efficientInfo:Double.NaN,jdistlib.ChiSquare.cumulative(mi,1,false,false)));
        }
        return java.util.List.copyOf(result);
    }
    public record Modification(String kind,String first,String second) {
        public static Modification regression(String outcome,String predictor){return new Modification("regression",outcome,predictor);}
        public static Modification covariance(String first,String second){return new Modification("covariance",first,second);}
        public static Modification intercept(String variable){return new Modification("intercept",variable,variable);}
    }
    public record ModificationIndex(Modification parameter,double chiSquare,double expectedChange,double pValue) { }
    public record RobustResult(double[] parameterCovariance,double scalingFactor,double chiSquare,
                               int degreesOfFreedom,double pValue,int clusters) {
        public RobustResult { parameterCovariance=parameterCovariance.clone(); }
        public double[] parameterCovariance(){return parameterCovariance.clone();}
    }

    /** Sobel/product-of-coefficients delta-method inference for an indirect effect. */
    public static IndirectEffect indirect(double a, double b, double varianceA,
                                          double varianceB, double covarianceAB,
                                          double confidenceLevel) {
        validateVariance(varianceA); validateVariance(varianceB);
        if (!Double.isFinite(a) || !Double.isFinite(b) || !Double.isFinite(covarianceAB)
                || !(confidenceLevel > 0.0) || !(confidenceLevel < 1.0))
            throw new IllegalArgumentException("indirect-effect inputs are invalid");
        if(Math.abs(covarianceAB)>Math.sqrt(varianceA)*Math.sqrt(varianceB)*(1+1e-12))
            throw new IllegalArgumentException("indirect-effect covariance is not positive semidefinite");
        double estimate = a * b;
        double variance = b * b * varianceA + a * a * varianceB
            + 2.0 * a * b * covarianceAB;
        double standardError = Math.sqrt(Math.max(0.0, variance));
        double critical = Normal.quantile(0.5 + confidenceLevel / 2.0, 0.0, 1.0, true, false);
        double z = estimate / standardError;
        return new IndirectEffect(estimate, standardError, z,
            estimate - critical * standardError, estimate + critical * standardError,
            2.0 * Normal.cumulative(-Math.abs(z), 0.0, 1.0, true, false));
    }

    /** Sandwich covariance from case or cluster score rows and an inverse bread. */
    public static double[] sandwich(double[][] scores, int[] clusters, double[] bread) {
        if (scores == null || scores.length == 0 || clusters == null
                || clusters.length != scores.length || bread == null
                || bread.length == 0) {
            throw new IllegalArgumentException("robust SEM dimensions are invalid");
        }
        int parameters = bread.length;
        int dimension = (int) Math.sqrt(parameters);
        if (dimension * dimension != parameters) throw new IllegalArgumentException("bread must be square");
        double[] meat = new double[parameters];
        java.util.Map<Integer, double[]> grouped = new java.util.LinkedHashMap<>();
        for (int row = 0; row < scores.length; row++) {
            if (scores[row] == null || scores[row].length != dimension)
                throw new IllegalArgumentException("score rows must match bread");
            double[] sum = grouped.computeIfAbsent(clusters[row], ignored -> new double[dimension]);
            for (int column = 0; column < dimension; column++) {
                if(!Double.isFinite(scores[row][column]))throw new IllegalArgumentException("nonfinite case score");
                sum[column] += scores[row][column];
            }
        }
        if(grouped.size()<2)throw new IllegalArgumentException("at least two independent score groups required");
        // Form outer products of influence vectors directly. B J B' assembled
        // through two dense products can acquire negative variances by cancellation.
        for (double[] sum : grouped.values()) {
            double[] influence=RamFit.mv(bread,sum);
            for(int row=0;row<dimension;row++)for(int column=0;column<=row;column++)
                meat[row*dimension+column]+=influence[row]*influence[column];
        }
        for(int row=0;row<dimension;row++)for(int column=0;column<row;column++)
            meat[column*dimension+row]=meat[row*dimension+column];
        return meat;
    }

    /** One-score-per-parameter modification index approximation. */
    public static double[] modificationIndices(double[] score, double[] information) {
        if (score == null || information == null || score.length != information.length)
            throw new IllegalArgumentException("score and information dimensions must match");
        double[] result = new double[score.length];
        for (int i = 0; i < result.length; i++) {
            if (!(information[i] > 0.0) || !Double.isFinite(information[i])) result[i] = Double.NaN;
            else result[i] = score[i] * score[i] / information[i];
        }
        return result;
    }

    /** Estimates normal latent-response thresholds from ordinal category counts. */
    public static double[] ordinalThresholds(int[] categories, int categoryCount) {
        if (categories == null || categories.length < 2 || categoryCount < 2)
            throw new IllegalArgumentException("ordinal categories are invalid");
        int[] counts = new int[categoryCount];
        for (int category : categories) { if (category < 0 || category >= categoryCount) throw new IllegalArgumentException("ordinal category is outside its range"); counts[category]++; }
        double[] result = new double[categoryCount - 1]; int cumulative = 0;
        for (int index = 0; index < result.length; index++) {
            cumulative += counts[index]; double probability = (double) cumulative / categories.length;
            result[index] = Normal.quantile(Math.max(1e-12, Math.min(1.0 - 1e-12, probability)), 0.0, 1.0, true, false);
        }
        return result;
    }

    private static void validateVariance(double value) { if (!(value >= 0.0) || !Double.isFinite(value)) throw new IllegalArgumentException("variance must be finite and nonnegative"); }
    private static double[] multiply(double[] left, double[] right, int dimension) { double[] result = new double[left.length]; for (int row = 0; row < dimension; row++) for (int column = 0; column < dimension; column++) for (int k = 0; k < dimension; k++) result[row * dimension + column] += left[row * dimension + k] * right[k * dimension + column]; return result; }
    private static void symmetrize(double[] matrix, int dimension) { for (int row = 0; row < dimension; row++) for (int column = row + 1; column < dimension; column++) { double value = 0.5 * (matrix[row * dimension + column] + matrix[column * dimension + row]); matrix[row * dimension + column] = value; matrix[column * dimension + row] = value; } }

    public record IndirectEffect(double estimate, double standardError, double statistic,
                                 double confidenceLower, double confidenceUpper, double pValue) { }
}
