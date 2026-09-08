/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.util.ArrayList;
import java.util.List;
import java.util.Properties;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.function.DoubleSupplier;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.genetics.ConditionalAssociation;
import org.jlinalg.genetics.ConditionalAssociationModel;
import org.jlinalg.mr.CorrelatedMendelianRandomization;
import org.jlinalg.mr.HarmonizedInstrument;
import org.jlinalg.mr.MrEstimate;
import org.jlinalg.mr.OverlapAwareMendelianRandomization;

/** Same five-SNP estimands as generate-mr-conditional-reference.R; consumes inference. */
public final class MrConditionalBenchmark {
    private MrConditionalBenchmark() { }
    public static void main(String[] args) throws Exception {
        Path root=Path.of("src/test/resources/r-reference");
        List<HarmonizedInstrument> values=new ArrayList<>();
        for(String line:Files.readAllLines(root.resolve("mr-conditional-input.tsv")).subList(1,6)) {
            String[] row=line.split("\t");
            values.add(new HarmonizedInstrument(row[0],"A","C",Double.parseDouble(row[1]),Double.parseDouble(row[2]),
                Double.parseDouble(row[3]),Double.parseDouble(row[4]),.3,.3,false,false));
        }
        double[][] ld=new double[5][5];
        for(int i=0;i<5;i++) for(int j=0;j<5;j++) ld[i][j]=Math.pow(.35,Math.abs(i-j));
        var associations=values.stream().map(v->new ConditionalAssociation(v.variantId(),v.exposureEffect(),v.exposureStandardError(),"locus")).toList();
        double[] covariance=values.stream().mapToDouble(v->.2*v.exposureStandardError()*v.outcomeStandardError()).toArray();
        Properties ref=new Properties();
        try(var input=Files.newInputStream(root.resolve("mr-conditional.properties"))) {ref.load(input);}
        DoubleSupplier conditional=()->new ConditionalAssociationModel(associations,ld).condition(new int[]{0,2},.9).stream()
            .mapToDouble(f->f.effect()+f.standardError()+f.pValue()+f.confidenceLower()+f.confidenceUpper()).sum();
        DoubleSupplier generalized=()->sum(CorrelatedMendelianRandomization.ivw(values,ld,false,.9,BackendPolicy.CPU).estimate());
        DoubleSupplier egger=()->{
            var f=CorrelatedMendelianRandomization.egger(values,ld,.9,BackendPolicy.CPU).estimate();
            return sum(f.slope())+f.intercept()+f.interceptStandardError()+f.interceptPValue()+f.interceptConfidenceLower()+f.interceptConfidenceUpper();
        };
        DoubleSupplier overlap=()->{
            var f=OverlapAwareMendelianRandomization.ivw(values,covariance,.9);
            if(!f.converged()) throw new IllegalStateException("overlap did not converge");
            return sum(f.estimate());
        };
        benchmark("conditional",conditional,referenceSum(ref,"conditional.1","conditional.2","conditional.3","conditional.4","conditional.5"));
        benchmark("generalized",generalized,referenceSum(ref,"ivw.fixed"));
        benchmark("egger",egger,referenceSum(ref,"egger.slope","egger.intercept"));
        benchmark("overlap",overlap,referenceSum(ref,"overlap"));
    }
    private static double sum(MrEstimate f) {return f.estimate()+f.standardError()+f.pValue()+f.confidenceLower()+f.confidenceUpper();}
    private static double referenceSum(Properties ref,String...keys) {
        double total=0;
        for(String key:keys) for(String field:List.of("beta","se","p","lower","upper")) total+=Double.parseDouble(ref.getProperty(key+"."+field));
        return total;
    }
    private static void benchmark(String name,DoubleSupplier f,double expected) {
        double actual=f.getAsDouble();
        if(Math.abs(actual-expected)>1e-10) throw new IllegalStateException(name+" reference mismatch: "+actual+" vs "+expected);
        for(int i=0;i<2000;i++) f.getAsDouble();
        int repetitions=20000; double checksum=0; long start=System.nanoTime();
        for(int i=0;i<repetitions;i++) checksum+=f.getAsDouble();
        System.out.printf(java.util.Locale.ROOT,"%s Java repetitions=%d seconds=%.6f checksum=%.12g%n",name,repetitions,(System.nanoTime()-start)/1e9,checksum);
    }
}
