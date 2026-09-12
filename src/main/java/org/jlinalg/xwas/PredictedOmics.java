/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.xwas;

import jdistlib.ChiSquare;
import org.jlinalg.genetics.GeneticCovarianceValidation;

/** Genetically predicted molecular association using raw-dosage weights and reference LD. */
public final class PredictedOmics {
    private PredictedOmics() { }
    public record Association(double z, double pValue, double logPValue, double predictedVariance) { }
    public record JointAssociation(double chiSquare, int degreesOfFreedom, double pValue) { }
    public static Association test(double[] z, double[] dosageWeights, double[] genotypeSd, double[] ld) {
        double[] a=standardized(dosageWeights,genotypeSd);
        SummaryMath.finite(z);
        if(z.length!=a.length)throw new IllegalArgumentException("GWAS and model dimensions differ");
        correlation(ld,a.length);
        double variance=SummaryMath.dot(a,SummaryMath.multiply(ld,a));
        if(!(variance>1e-12*SummaryMath.dot(a,a)))
            throw new IllegalArgumentException("predicted molecular variance is zero or numerically unresolved");
        double score=SummaryMath.dot(a,z)/Math.sqrt(variance);
        if(!Double.isFinite(score))throw new IllegalArgumentException("nonfinite molecular score");
        return new Association(score,SummaryMath.p(score),SummaryMath.logP(score),variance);
    }
    /** Joint tissue/model test on the same complete SNP set; dependent models reject. */
    public static JointAssociation joint(double[] z,double[][] dosageWeights,double[] genotypeSd,double[] ld) {
        int k=dosageWeights.length,n=z.length;
        if(k<1)throw new IllegalArgumentException("at least one model is required");
        double[][] a=new double[k][];double[] zs=new double[k],cov=new double[k*k],var=new double[k];
        for(int i=0;i<k;i++) {
            Association fit=test(z,dosageWeights[i],genotypeSd,ld);
            zs[i]=fit.z();var[i]=fit.predictedVariance();a[i]=standardized(dosageWeights[i],genotypeSd);
        }
        for(int i=0;i<k;i++)for(int j=0;j<k;j++)
            cov[i*k+j]=SummaryMath.dot(a[i],SummaryMath.multiply(ld,a[j]))/Math.sqrt(var[i]*var[j]);
        double q=SummaryMath.dot(zs,SummaryMath.multiply(SummaryMath.inverse(cov,k),zs));
        return new JointAssociation(q,k,ChiSquare.cumulative(q,k,false,false));
    }
    private static double[] standardized(double[] w,double[] sd) {
        SummaryMath.finite(w);SummaryMath.finite(sd);
        if(w.length<1||w.length!=sd.length)throw new IllegalArgumentException("invalid model dimensions");
        double[] a=w.clone();
        for(int i=0;i<a.length;i++) {
            if(!(sd[i]>0))throw new IllegalArgumentException("genotype SD must be positive");
            a[i]*=sd[i];
        }
        SummaryMath.finite(a);return a;
    }
    public static void correlation(double[] matrix,int n) {
        GeneticCovarianceValidation.requirePositiveSemidefinite(matrix,n);
        for(int i=0;i<n;i++)if(Math.abs(matrix[i*n+i]-1)>1e-8)
            throw new IllegalArgumentException("correlation matrix diagonal must equal one");
    }
    /** Exact forward-strand allele match or swap. Strand complements are never guessed. */
    public static int alleleSign(String ea,String oa,String referenceEa,String referenceOa) {
        if(ea==null||oa==null||referenceEa==null||referenceOa==null||ea.isBlank()||oa.isBlank()
                ||ea.equals(oa)||referenceEa.equals(referenceOa))throw new IllegalArgumentException("invalid alleles");
        for(String allele:new String[]{ea,oa,referenceEa,referenceOa})
            if(!allele.matches("[ACGT]+|<[^\\s<>\"]+>"))throw new IllegalArgumentException("alleles must be uppercase forward-strand bases or symbolic alleles");
        if(ea.equals(referenceEa)&&oa.equals(referenceOa))return 1;
        if(ea.equals(referenceOa)&&oa.equals(referenceEa))return -1;
        throw new IllegalArgumentException("allele mismatch: "+ea+"/"+oa+" versus "+referenceEa+"/"+referenceOa);
    }
}
