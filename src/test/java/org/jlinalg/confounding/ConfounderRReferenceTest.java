/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.confounding;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Test;

class ConfounderRReferenceTest {
    private static final String ROOT = "/r-reference/confounding/";

    @Test
    void pcaMatchesRPrcompSampleSubspace() throws Exception {
        double[][] data=table("input-data.tsv",1);
        double[][] expected=table("pca-factors.tsv",1);
        double[][] actual=PrincipalComponentConfounders.fit(data,
            PrincipalComponentConfounders.Options.defaults(2)).factors();
        assertTrue(projectionRmse(actual,expected)<1e-10);
    }

    @Test
    void standardSvaTracksBioconductorReferenceSubspace() throws Exception {
        double[][] data=table("input-data.tsv",1),pheno=caseValues();
        double[][] full=design(pheno),reduced=intercept(pheno.length);
        double[][] expected=table("sva-factors.tsv",1);
        SvaResult actual=SurrogateVariableAnalysis.fit(data,full,reduced,
            new SurrogateVariableAnalysis.Options(2,5));
        assertTrue(projectionRmse(actual.surrogateVariables(),expected)<0.12);
        double[][] weights=table("sva-weights.tsv",1);
        double heterogeneityRmse=rmse(actual.probabilityHeterogeneity(),column(weights,0));
        double biologicalCorrelation=correlation(actual.probabilityBiological(),column(weights,1));
        assertTrue(heterogeneityRmse<0.2,"heterogeneity RMSE="+heterogeneityRmse);
        assertTrue(biologicalCorrelation>0.85,"biological correlation="+biologicalCorrelation);
    }

    @Test
    void autoSvaTracksSourceReferenceSubspace() throws Exception {
        double[][] data=table("input-data.tsv",1),pheno=caseValues();
        double[][] expected=table("autosva-factors.tsv",1);
        SvaResult actual=AutoSva.fit(data,design(pheno),intercept(pheno.length),
            AutoSva.Options.fixed(2));
        assertTrue(projectionRmse(actual.surrogateVariables(),expected)<0.12);
        double[][] weights=table("autosva-weights.tsv",1);
        assertTrue(correlation(actual.featureWeights(),column(weights,0))>0.85);
    }

    @Test
    void combatMatchesBioconductorParametricAdjustment() throws Exception {
        double[][] data=table("input-data.tsv",1),pheno=caseValues();
        String[] batch=strings("input-pheno.tsv",2);
        double[][] expected=table("combat-adjusted.tsv",1);
        ComBat.Result actual=ComBat.adjust(data,batch,design(pheno),ComBat.Options.defaults());
        double maximum=0.0;
        for(int i=0;i<expected.length;i++)for(int j=0;j<expected[i].length;j++)maximum=Math.max(maximum,Math.abs(expected[i][j]-actual.adjusted()[i][j]));
        assertTrue(maximum<1e-8,"maximum absolute error="+maximum);
    }

    private static double[][] design(double[][] pheno){double[][] result=new double[pheno.length][2];for(int i=0;i<result.length;i++){result[i][0]=1;result[i][1]=pheno[i][0];}return result;}
    private static double[][] caseValues() throws Exception {String[] values=strings("input-pheno.tsv",1);double[][] result=new double[values.length][1];for(int i=0;i<values.length;i++)result[i][0]=Double.parseDouble(values[i]);return result;}
    private static double[][] intercept(int n){double[][] result=new double[n][1];for(double[] row:result)row[0]=1;return result;}
    private static double projectionRmse(double[][] first,double[][] second){double[][] a=orthonormal(first),b=orthonormal(second);double sum=0;int n=a.length;for(int i=0;i<n;i++)for(int j=0;j<n;j++){double pa=0,pb=0;for(int k=0;k<a[0].length;k++)pa+=a[i][k]*a[j][k];for(int k=0;k<b[0].length;k++)pb+=b[i][k]*b[j][k];sum+=(pa-pb)*(pa-pb);}return Math.sqrt(sum/(n*n));}
    private static double[] column(double[][] matrix,int column){double[] result=new double[matrix.length];for(int i=0;i<result.length;i++)result[i]=matrix[i][column];return result;}
    private static double correlation(double[] x,double[] y){double xm=java.util.Arrays.stream(x).average().orElseThrow(),ym=java.util.Arrays.stream(y).average().orElseThrow(),xy=0,xx=0,yy=0;for(int i=0;i<x.length;i++){xy+=(x[i]-xm)*(y[i]-ym);xx+=(x[i]-xm)*(x[i]-xm);yy+=(y[i]-ym)*(y[i]-ym);}return xy/Math.sqrt(xx*yy);}
    private static double rmse(double[] x,double[] y){double sum=0;for(int i=0;i<x.length;i++)sum+=(x[i]-y[i])*(x[i]-y[i]);return Math.sqrt(sum/x.length);}
    private static double[][] orthonormal(double[][] x){int n=x.length,p=x[0].length;double[][] q=new double[n][p];for(int c=0;c<p;c++){for(int i=0;i<n;i++)q[i][c]=x[i][c];for(int k=0;k<c;k++){double z=0;for(int i=0;i<n;i++)z+=q[i][k]*q[i][c];for(int i=0;i<n;i++)q[i][c]-=z*q[i][k];}double norm=0;for(int i=0;i<n;i++)norm+=q[i][c]*q[i][c];norm=Math.sqrt(norm);for(int i=0;i<n;i++)q[i][c]/=norm;}return q;}
    private static double[][] table(String name,int skip)throws Exception{List<double[]> rows=new ArrayList<>();try(BufferedReader input=new BufferedReader(new InputStreamReader(ConfounderRReferenceTest.class.getResourceAsStream(ROOT+name),StandardCharsets.UTF_8))){input.readLine();for(String line;(line=input.readLine())!=null;){String[] values=line.split("\t",-1);double[] row=new double[values.length-skip];for(int i=skip;i<values.length;i++)row[i-skip]=Double.parseDouble(values[i]);rows.add(row);}}return rows.toArray(double[][]::new);}
    private static String[] strings(String name,int column)throws Exception{List<String> result=new ArrayList<>();try(BufferedReader input=new BufferedReader(new InputStreamReader(ConfounderRReferenceTest.class.getResourceAsStream(ROOT+name),StandardCharsets.UTF_8))){input.readLine();for(String line;(line=input.readLine())!=null;)result.add(line.split("\t",-1)[column]);}return result.toArray(String[]::new);}
}
