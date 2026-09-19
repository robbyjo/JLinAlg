/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.nio.file.*;
import java.util.*;
import java.util.function.DoubleSupplier;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.*;
import org.jlinalg.regression.*;

/** In-memory CPU timings on the checked-in 240-row independently checked fixtures. */
public final class RegressionInferenceBenchmark {
    private RegressionInferenceBenchmark() { }
    private static volatile double sink;
    public static void main(String[] args) throws Exception {
        double[][] censored=read("censored"),ordinal=read("ordinal"),survey=read("survey");
        System.out.println("method\trows\tmedian_ms\tchecksum");
        for(var dist:CensoredRegression.Distribution.values()) {
            boolean gaussian=dist==CensoredRegression.Distribution.GAUSSIAN;double[] cy=column(censored,gaussian?0:4);
            int[] censor=Arrays.stream(censored).mapToInt(r->(int)r[gaussian?3:5]).toArray();double[][] cx=design(censored,1,2);
            time("censored_"+dist.name().toLowerCase(Locale.ROOT),()->CensoredRegression.fit(cy,cx,censor,dist).logLikelihood());
        }
        int[] oy=Arrays.stream(ordinal).mapToInt(r->(int)r[0]).toArray();double[][] ox=matrix(ordinal,1,2);
        for(var link:OrdinalRegression.Link.values())time("ordinal_"+link.name().toLowerCase(Locale.ROOT),()->OrdinalRegression.fit(oy,ox,4,link).logLikelihood());
        time("rare_events",()->RareEventsLogistic.fit(column(survey,3),design(survey,1,2),true,.04).coefficients()[0]);
        var sampling=new SurveyRegression.Design(column(survey,5),Arrays.stream(survey).map(r->Double.toString(r[6])).toArray(String[]::new),Arrays.stream(survey).map(r->Double.toString(r[7])).toArray(String[]::new));
        time("survey_gaussian",()->SurveyRegression.fit(column(survey,0),design(survey,1,2),GlmFamilies.gaussian(),sampling).coefficients()[0]);
    }
    private static void time(String name,DoubleSupplier action){
        for(int i=0;i<3;i++)sink=action.getAsDouble();double[] elapsed=new double[7];
        for(int i=0;i<elapsed.length;i++){long start=System.nanoTime();sink=action.getAsDouble();elapsed[i]=(System.nanoTime()-start)/1e6;}
        Arrays.sort(elapsed);System.out.printf(Locale.ROOT,"%s\t240\t%.6f\t%.12g%n",name,elapsed[3],sink);
    }
    private static double[][] read(String name)throws Exception{return Files.readAllLines(Path.of("src/test/resources/r-reference/regression-inference/"+name+".tsv")).stream().skip(1).map(r->Arrays.stream(r.split("\t")).mapToDouble(Double::parseDouble).toArray()).toArray(double[][]::new);}
    private static double[] column(double[][] data,int j){return Arrays.stream(data).mapToDouble(r->r[j]).toArray();}
    private static double[][] matrix(double[][] data,int... cols){double[][] x=new double[data.length][cols.length];for(int i=0;i<x.length;i++)for(int j=0;j<cols.length;j++)x[i][j]=data[i][cols[j]];return x;}
    private static double[][] design(double[][] data,int... cols){double[][] x=new double[data.length][cols.length+1];for(int i=0;i<x.length;i++){x[i][0]=1;for(int j=0;j<cols.length;j++)x[i][j+1]=data[i][cols[j]];}return x;}
}
