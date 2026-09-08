/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.function.Supplier;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.meta.*;

/** Small dense meta-analysis macrobenchmark, with metafor accuracy gates and consumed results.
 * Run directly after compiling; does not need a new Gradle task. */
public final class MetaParityBenchmark {
    private MetaParityBenchmark() { }
    public static void main(String[] args) throws Exception {
        int repeats = args.length == 0 ? 100 : Integer.parseInt(args[0]);
        if (repeats < 1) throw new IllegalArgumentException("positive repeat count required");
        Path root = Path.of("src/test/resources/meta");
        List<String> input = Files.readAllLines(root.resolve("hierarchical.csv")), sampling = Files.readAllLines(root.resolve("sampling.csv"));
        int n = input.size() - 1;
        List<MetaStudy> studies = new ArrayList<>(), bias = new ArrayList<>();
        double[][] mods = new double[n][1], random = new double[n][2], v = new double[n][n]; String[] groups = new String[n], ids = new String[n];
        for (int i = 0; i < n; i++) {
            String[] row = input.get(i + 1).split(","), cov = sampling.get(i).split(",");
            studies.add(new MetaStudy("e" + i, Double.parseDouble(row[0]), Math.sqrt(Double.parseDouble(row[1]))));
            mods[i][0] = Double.parseDouble(row[2]); random[i][0] = 1; random[i][1] = mods[i][0]; groups[i] = row[3]; ids[i] = row[4];
            for (int j = 0; j < n; j++) v[i][j] = Double.parseDouble(cov[j]);
        }
        List<String> biasLines = Files.readAllLines(root.resolve("bias.csv"));
        for (int i = 1; i < biasLines.size(); i++) {
            String[] row = biasLines.get(i).split(","); bias.add(new MetaStudy("b" + i, Double.parseDouble(row[0]), Double.parseDouble(row[1])));
        }
        Map<String, Double> ref = new HashMap<>(); List<String> references = Files.readAllLines(root.resolve("metafor-reference.csv"));
        for (String row : references.subList(1,references.size())) { String[] a = row.split(","); ref.put(a[0],Double.parseDouble(a[1])); }
        MetaAnalysisOptions options = MetaAnalysisOptions.builder().maximumIterations(1000).tolerance(1e-11).build();
        MetaAnalysisOptions fixedT = MetaAnalysisOptions.builder().method(MetaAnalysisMethod.FIXED_EFFECT).inferenceMethod(MetaInferenceMethod.STUDENT_T).build();
        List<MetaRandomEffect> unstructured = List.of(new MetaRandomEffect("study",groups,random,MetaRandomEffect.Structure.UNSTRUCTURED));
        List<MetaRandomEffect> nested = List.of(MetaRandomEffect.intercept("study",groups),MetaRandomEffect.intercept("effect",ids));
        List<Scenario> cases = new ArrayList<>();
        cases.add(new Scenario("known_gls", () -> {
            var f = MetaMultilevelRegression.fit(studies,mods,List.of("x"),v,true,BackendPolicy.CPU);
            return join(f.beta(),f.covariance(),new double[] {f.logLikelihood()});
        }, expected(ref,"known.beta",2,"known.cov",4,"known.ll",1),1e-9));
        for (boolean useNested : new boolean[] {false,true}) {
            String name = useNested ? "nested_reml" : "random_moderators_reml", key = useNested ? "nested.REML" : "hier.REML";
            cases.add(new Scenario(name, () -> {
                var f = MetaMultilevelRegression.fit(studies,mods,List.of("x"),v,true,useNested ? nested : unstructured,
                    MetaMultilevelRegression.Estimation.REML,options,BackendPolicy.CPU);
                if (!f.converged()) throw new AssertionError("nonconverged hierarchical benchmark");
                double[] covariance = useNested ? new double[] {f.randomCovariances().get(0)[0],f.randomCovariances().get(1)[0]} : f.randomCovariances().get(0);
                return join(f.beta(),covariance,new double[] {f.logLikelihood()});
            },expected(ref,key+".beta",2,key+(useNested ? ".variance" : ".G"),useNested ? 2 : 4,key+".ll",1),2e-5));
        }
        cases.add(new Scenario("cr2_correlated", () -> {
            var f = MetaClusterRobust.fit(studies,groups,mods,List.of("x"),v,fixedT,MetaRobustCorrection.CR2,BackendPolicy.CPU);
            return join(f.covariance(),f.associationStatistics().degreesOfFreedom(),f.associationStatistics().pValues());
        },expected(ref,"CR2.cov",4,"CR2.df",2,"CR2.p",2),1e-9));
        cases.add(new Scenario("trimfill_fe_l0", () -> {
            var f = MetaBiasCorrections.trimAndFill(bias,MetaAnalysisOptions.fixedEffect(),MetaBiasCorrections.Side.LEFT,
                MetaBiasCorrections.TrimEstimator.L0,BackendPolicy.CPU);
            return new double[] {f.adjustedEffect(),f.adjustedFit().standardError(),f.imputedStudyCount()};
        },expected(ref,"trim.FE.L0.left.beta",1,"trim.FE.L0.left.se",1,"trim.FE.L0.left.k0",1),1e-10));
        cases.add(new Scenario("pet", () -> {
            var f = MetaBiasCorrections.pet(bias); return join(new double[] {f.intercept(),f.slope()},f.covariance(),f.pValues());
        },expected(ref,"pet.beta",2,"pet.cov",4,"pet.p",2),1e-10));
        System.out.println("method,repeats,milliseconds_per_fit,max_absolute_reference_error,checksum");
        for (Scenario scenario : cases) {
            double error = 0;
            double[] actual = scenario.run.get();
            for (int i = 0; i < actual.length; i++) error = Math.max(error,Math.abs(actual[i]-scenario.expected[i]));
            if (!(error <= scenario.tolerance)) throw new AssertionError(scenario.name+" reference error "+error);
            for (int i = 0; i < 10; i++) scenario.run.get();
            double checksum = 0; long start = System.nanoTime();
            for (int i = 0; i < repeats; i++) for (double value : scenario.run.get()) checksum += value;
            double millis = (System.nanoTime()-start)/1e6/repeats;
            System.out.printf(Locale.ROOT,"%s,%d,%.9f,%.12g,%.12g%n",scenario.name,repeats,millis,error,checksum);
        }
    }
    private static double[] expected(Map<String,Double> ref,String a,int na,String b,int nb,String c,int nc) {
        double[] out = new double[na+nb+nc];
        for (int i=0;i<na;i++) out[i]=ref.get(a+"."+i);
        for (int i=0;i<nb;i++) out[na+i]=ref.get(b+"."+i);
        for (int i=0;i<nc;i++) out[na+nb+i]=ref.get(c+"."+i);
        return out;
    }
    private static double[] join(double[] a,double[] b,double[] c) {
        double[] out = new double[a.length+b.length+c.length];
        System.arraycopy(a,0,out,0,a.length); System.arraycopy(b,0,out,a.length,b.length); System.arraycopy(c,0,out,a.length+b.length,c.length); return out;
    }
    private record Scenario(String name,Supplier<double[]> run,double[] expected,double tolerance) { }
}
