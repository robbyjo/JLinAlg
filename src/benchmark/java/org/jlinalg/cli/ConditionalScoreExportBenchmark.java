/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.Properties;

/** End-to-end CLI timings include parsing, null fitting, scans, FDR and disk
 * output. Repeated fixture genotypes measure a full 64-variant covariance block;
 * every score and matrix entry is checked against the corresponding R value. */
public final class ConditionalScoreExportBenchmark {
    private ConditionalScoreExportBenchmark() { }
    public static void main(String[] ignored) throws Exception {
        Path fixture=Path.of("src/test/resources/conditional-score-export");
        Path directory=Path.of("build/conditional-score-benchmark");Files.createDirectories(directory);
        List<String> input=Files.readAllLines(fixture.resolve("variants.vcf"));
        List<String> headers=input.stream().filter(s->s.startsWith("#")).toList();
        List<String> records=input.stream().filter(s->!s.startsWith("#")).toList();
        List<String> variants=new ArrayList<>(headers);
        for(int j=0;j<64;j++) {
            String[] row=records.get(j%4).split("\t");row[1]=Integer.toString(100*(j+1));row[2]="bench"+j;
            variants.add(String.join("\t",row));
        }
        Path genotypes=directory.resolve("variants.vcf");Files.write(genotypes,variants);
        System.out.println("case,standard_median_ms,export_median_ms,max_score_error,max_covariance_error");
        for(String model:List.of("binomial","poisson","gaussian","cox-efron","cox-breslow")) {
            Properties ref=new Properties();try(var reader=Files.newBufferedReader(fixture.resolve(model+".properties"))){ref.load(reader);}
            double[] u=numbers(ref.getProperty("u")),v=numbers(ref.getProperty("v"));
            double[] times=new double[2],errors=new double[2];
            for(int export=0;export<2;export++) {
                boolean cox=model.startsWith("cox");Path out=directory.resolve(model+"-"+export+".tsv");
                List<String> args=new ArrayList<>(List.of("--omics",genotypes.toString(),"--pheno",fixture.resolve("phenotype.tsv").toString(),
                    "--id","IID","--model",cox?"cox":"glm","--formula",
                    (cox?"Surv(start,stop,event)":model.equals("binomial")?"binary":"count")+"~x+offset(o)+<omics>",
                    "--threads","1","--backend","cpu","--block-size","64","--no-log","--overwrite","--out",out.toString()));
                if(cox)args.addAll(List.of("--ties",model.substring(4)));else args.addAll(List.of("--family",model));
                if(export==1)args.addAll(List.of("--conditional-gwas-summary","--score-genome-build","GRCh38","--score-block-size","64"));
                run(args);
                if(export==1) {
                    List<String> summary=Files.readAllLines(out),header=List.of(summary.get(0).split("\t"));
                    if(summary.size()!=65)throw new AssertionError("wrong summary row count");
                    for(int j=0;j<64;j++) {
                        String[] row=summary.get(j+1).split("\t",-1);
                        errors[0]=Math.max(errors[0],Math.abs(Double.parseDouble(row[header.indexOf("score_u")])-u[j%4]));
                    }
                    List<String> covariance=Files.readAllLines(Path.of(out+".score-cov.tsv"));
                    if(covariance.size()!=2081)throw new AssertionError("incomplete covariance block");
                    int line=1;
                    for(int j=0;j<64;j++)for(int k=j;k<64;k++) {
                        double actual=Double.parseDouble(covariance.get(line++).split("\t")[2]);
                        errors[1]=Math.max(errors[1],Math.abs(actual-v[(j%4)*4+k%4]));
                    }
                    if(errors[0]>2e-6 || errors[1]>2e-6)throw new AssertionError("R accuracy gate failed: "+Arrays.toString(errors));
                }
                int warmups=Integer.getInteger("jlinalg.benchmark.warmups",3),measurements=Integer.getInteger("jlinalg.benchmark.measurements",7);
                if(warmups<0 || measurements<1)throw new IllegalArgumentException("invalid benchmark repetitions");
                double[] elapsed=new double[measurements];
                for(int i=-warmups;i<measurements;i++) {
                    long start=System.nanoTime();run(args);if(i>=0)elapsed[i]=(System.nanoTime()-start)/1e6;
                }
                Arrays.sort(elapsed);times[export]=elapsed[measurements/2];
            }
            System.out.printf(Locale.ROOT,"%s-n179-m64,%.6f,%.6f,%.9g,%.9g%n",model,times[0],times[1],errors[0],errors[1]);
        }
    }
    private static void run(List<String> args) {
        ByteArrayOutputStream bytes=new ByteArrayOutputStream();PrintStream stream=new PrintStream(bytes);
        if(JLinAlgCli.run(args.toArray(String[]::new),stream,stream)!=0)throw new AssertionError(bytes.toString());
    }
    private static double[] numbers(String value){return Arrays.stream(value.split(",")).mapToDouble(Double::parseDouble).toArray();}
}
