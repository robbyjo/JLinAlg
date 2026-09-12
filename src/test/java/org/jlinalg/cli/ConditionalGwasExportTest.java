/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class ConditionalGwasExportTest {
    @TempDir Path temporary;
    private static final Path FIXTURE=Path.of("src/test/resources/conditional-score-export");

    @Test void scoresAndAllCovariancesMatchIndependentRNullRefits() throws Exception {
        for(String model:List.of("binomial","poisson","gaussian","cox-efron","cox-breslow","cox-right-efron","cox-right-breslow"))
            for(boolean conditional:List.of(false,true)) {
                String name=model+(conditional?"-conditional":"");
                Path out=temporary.resolve(name+".tsv");
                List<String> args=arguments(model,out,true);
                if(conditional)args.addAll(List.of("--condition-on","g1"));
                Result result=run(args);assertEquals(0,result.code,result.console);
                Properties expected=new Properties();
                try(var reader=Files.newBufferedReader(FIXTURE.resolve(name+".properties"))){expected.load(reader);}
                double[] u=numbers(expected.getProperty("u")),v=numbers(expected.getProperty("v"));
                DelimitedData summary=DelimitedData.read(out);
                assertEquals(u.length,summary.rows().size());
                String modelId=null;
                for(int j=0;j<u.length;j++) {
                    String[] row=summary.rows().get(j);
                    assertEquals("ok",row[summary.column("status")]);
                    assertEquals(u[j],value(summary,row,"score_u"),2e-6,name+" U "+j);
                    assertEquals(v[j*u.length+j],value(summary,row,"score_variance"),2e-6,name+" V "+j);
                    assertEquals(u[j]/v[j*u.length+j],value(summary,row,"beta_score"),2e-7);
                    assertEquals(179,value(summary,row,"n_analyzed"));
                    assertEquals("normal_only",row[summary.column("calibration_status")]);
                    assertEquals("",row[summary.column("p_score_calibrated")]);
                    assertEquals(conditional?"1:"+((j+2)*100)+":A:G":"1:"+((j+1)*100)+":A:G",row[summary.column("score_variant_key")]);
                    if(modelId==null)modelId=row[summary.column("null_model_id")];
                    assertEquals(modelId,row[summary.column("null_model_id")]);
                    assertEquals(!conditional,row[summary.column("conditioning_set_id")].equals("none"));
                    if(model.equals("binomial")) {
                        assertEquals(179,value(summary,row,"n_cases")+value(summary,row,"n_controls"));
                        assertEquals(numbers(expected.getProperty("nCases"))[0],value(summary,row,"n_cases"));
                        assertEquals(numbers(expected.getProperty("nControls"))[0],value(summary,row,"n_controls"));
                        assertEquals(numbers(expected.getProperty("calledCases"))[j],value(summary,row,"n_called_cases"));
                        assertEquals(numbers(expected.getProperty("calledControls"))[j],value(summary,row,"n_called_controls"));
                        assertEquals(numbers(expected.getProperty("acCases"))[j],value(summary,row,"effect_ac_cases"),1e-12);
                        assertEquals(numbers(expected.getProperty("acControls"))[j],value(summary,row,"effect_ac_controls"),1e-12);
                        int called=(int)(value(summary,row,"n_called_cases")+value(summary,row,"n_called_controls"));
                        assertTrue(called<=179);if(row[summary.column("score_variant_key")].contains(":300:"))assertTrue(called<179);
                    } else if(model.startsWith("cox"))assertTrue(value(summary,row,"n_events")>0);
                }
                DelimitedData cov=DelimitedData.read(Path.of(out+".score-cov.tsv"));
                assertEquals(u.length*(u.length+1)/2,cov.rows().size());
                int index=0;
                for(int j=0;j<u.length;j++)for(int k=j;k<u.length;k++) {
                    String[] row=cov.rows().get(index++);
                    assertEquals(v[j*u.length+k],value(cov,row,"score_covariance"),2e-6,name+" covariance "+j+","+k);
                    assertEquals(modelId,row[cov.column("null_model_id")]);
                    assertEquals(summary.rows().get(j)[summary.column("score_variant_key")],row[cov.column("variant_i")]);
                    assertEquals(summary.rows().get(k)[summary.column("score_variant_key")],row[cov.column("variant_j")]);
                }
                String manifest=Files.readString(Path.of(out+".score-manifest.json"));
                assertTrue(manifest.contains("\"status\": \"complete\""));
                assertTrue(manifest.contains("cross-block covariance unavailable, not zero"));
                if(!model.equals("gaussian"))assertTrue(Files.readString(Path.of(out+".log")).contains("Conditional-GWAS export enabled"));
            }
    }

    @Test void optInPreservesOrdinaryColumnsAndScopesConsoleAndLogNotice() throws Exception {
        for(String model:List.of("binomial","poisson","cox-efron","gaussian")) {
            Path ordinary=temporary.resolve(model+"-ordinary.tsv"),extended=temporary.resolve(model+"-extended.tsv");
            Result standard=run(arguments(model,ordinary,false));assertEquals(0,standard.code,standard.console);
            Result extra=run(arguments(model,extended,true));assertEquals(0,extra.code,extra.console);
            boolean warn=!model.equals("gaussian");
            assertEquals(warn,standard.console.contains("current standard output is insufficient"));
            assertEquals(warn,Files.readString(Path.of(ordinary+".log")).contains("current standard output is insufficient"));
            assertFalse(Files.exists(Path.of(ordinary+".score-cov.tsv")));
            DelimitedData a=DelimitedData.read(ordinary),b=DelimitedData.read(extended);
            assertFalse(a.header().contains("score_u"));
            assertEquals(a.rows().size(),b.rows().size());
            for(int i=0;i<a.rows().size();i++)for(String column:a.header())
                assertEquals(a.rows().get(i)[a.column(column)],b.rows().get(i)[b.column(column)],column);
        }
    }

    @Test void blockCoverageIsExplicitAndConditionAliasesAndMissingKeysAreValidated() throws Exception {
        Path out=temporary.resolve("blocks.tsv");List<String> args=arguments("binomial",out,true);
        args.addAll(List.of("--score-block-size","2","--block-size","1","--no-log"));
        Result result=run(args);assertEquals(0,result.code,result.console);
        assertFalse(Files.exists(Path.of(out+".log")));
        assertEquals(6,DelimitedData.read(Path.of(out+".score-cov.tsv")).rows().size());
        assertTrue(Files.readString(Path.of(out+".score-manifest.json")).contains("\"covariance_blocks\": \"2\""));
        assertTrue(Files.readString(Path.of(out+".score-manifest.json")).contains("\"case_value\": \"1\""));
        assertTrue(Files.readString(Path.of(out+".score-manifest.json")).contains("\"control_value\": \"0\""));
        for(String conditions:List.of("missing","g1,1:100:A:G","g1,g1")) {
            Path bad=temporary.resolve("bad"+conditions.hashCode()+".tsv");List<String> invalid=arguments("binomial",bad,true);
            invalid.addAll(List.of("--condition-on",conditions));
            assertNotEquals(0,run(invalid).code);assertFalse(Files.exists(Path.of(bad+".score-manifest.json")));
        }
        Path alias=temporary.resolve("alias.tsv");List<String> valid=arguments("binomial",alias,true);
        valid.addAll(List.of("--condition-on","1:100:A:G"));
        result=run(valid);assertEquals(0,result.code,result.console);
        assertEquals(3,DelimitedData.read(alias).rows().size());
    }

    @Test void invalidExportControlsDoNotPublishSuccess() throws Exception {
        for(List<String> additions:List.of(List.of("--resume"),List.of("--score-block-size","513"),List.of("--score-genome-build",""))) {
            Path out=temporary.resolve("reject"+additions.hashCode()+".tsv");List<String> args=arguments("binomial",out,true);args.addAll(additions);
            assertNotEquals(0,run(args).code);assertFalse(Files.exists(Path.of(out+".score-manifest.json")));
        }
        List<String> args=arguments("binomial",temporary.resolve("noswitch.tsv"),false);
        args.addAll(List.of("--condition-on","g1"));assertNotEquals(0,run(args).code);
    }

    @Test void invalidOrderingFailsWithoutLeavingAStaleSuccessManifest() throws Exception {
        Path out=temporary.resolve("overwrite.tsv");
        Result first=run(arguments("binomial",out,true));assertEquals(0,first.code,first.console);
        List<String> rows=new ArrayList<>(Files.readAllLines(FIXTURE.resolve("variants.vcf")));
        int firstVariant=0;while(rows.get(firstVariant).startsWith("#"))firstVariant++;
        rows.add(firstVariant+1,rows.get(firstVariant).replace("\tg1\t","\tduplicate\t"));
        Path duplicates=temporary.resolve("duplicates.vcf");Files.write(duplicates,rows);
        List<String> args=arguments("binomial",out,true);
        args.set(args.indexOf("--omics")+1,duplicates.toString());
        args.addAll(List.of("--score-block-size","1","--overwrite"));
        Result result=run(args);assertNotEquals(0,result.code);
        assertTrue(result.console.contains("duplicate score variant key"),result.console);
        assertFalse(Files.exists(Path.of(out+".score-manifest.json")));
        assertTrue(Files.readString(Path.of(out+".log")).contains("failed"));
    }

    @Test void filtersLeaveScoreFieldsBlankAndOnlyAcceptedRowsHaveCovariance() throws Exception {
        Path out=temporary.resolve("filtered.tsv");List<String> args=arguments("binomial",out,true);
        args.addAll(List.of("--max-marker-missing","0"));
        Result result=run(args);assertEquals(0,result.code,result.console);
        DelimitedData summary=DelimitedData.read(out);
        assertEquals(4,summary.rows().size());
        String[] missing=summary.rows().get(2);
        assertEquals("filtered",missing[summary.column("status")]);
        for(String column:ConditionalGwasExport.COLUMNS)assertEquals("",missing[summary.column(column)]);
        assertEquals(6,DelimitedData.read(Path.of(out+".score-cov.tsv")).rows().size());
    }

    @Test void fractionalBinaryAndUnsupportedNullsAreRejected() throws Exception {
        Path phenotype=temporary.resolve("fractional.tsv");
        DelimitedData original=DelimitedData.read(FIXTURE.resolve("phenotype.tsv"));
        List<String> rows=new ArrayList<>();rows.add(String.join("\t",original.header()));
        for(String[] row:original.rows()){row[original.column("binary")]=row[original.column("binary")].equals("1")?"0.8":"0.2";rows.add(String.join("\t",row));}
        Files.write(phenotype,rows);Path out=temporary.resolve("fractional-out.tsv");
        List<String> args=arguments("binomial",out,true);args.set(args.indexOf("--pheno")+1,phenotype.toString());
        Result result=run(args);assertNotEquals(0,result.code);assertTrue(result.console.contains("binomial response"),result.console);
        assertFalse(Files.exists(Path.of(out+".score-manifest.json")));
        args=arguments("binomial",temporary.resolve("quasi.tsv"),true);
        args.set(args.indexOf("--family")+1,"quasi-binomial");
        result=run(args);assertNotEquals(0,result.code);assertTrue(result.console.contains("separate score exporters"),result.console);
    }

    @Test void coxCategoricalCovariatesUseIdentifiableReferenceContrasts() throws Exception {
        List<String> original=Files.readAllLines(FIXTURE.resolve("phenotype.tsv"));List<String> rows=new ArrayList<>();
        rows.add(original.get(0)+"\tgroup\tindicator");
        for(int i=1;i<original.size();i++)rows.add(original.get(i)+(i%2==0?"\ta\t0":"\tb\t1"));
        Path phenotype=temporary.resolve("factors.tsv");Files.write(phenotype,rows);
        DelimitedData[] results=new DelimitedData[2];
        for(int j=0;j<2;j++) {
            Path out=temporary.resolve("factor"+j+".tsv");List<String> args=arguments("cox-efron",out,true);
            args.set(args.indexOf("--pheno")+1,phenotype.toString());
            args.set(args.indexOf("--formula")+1,"Surv(start,stop,event)~x+"+(j==0?"group":"indicator")+"+offset(o)+<omics>");
            Result result=run(args);assertEquals(0,result.code,result.console);results[j]=DelimitedData.read(out);
        }
        for(int j=0;j<4;j++)for(String column:List.of("score_u","score_variance","beta_score"))
            assertEquals(value(results[0],results[0].rows().get(j),column),value(results[1],results[1].rows().get(j),column),1e-10);
        List<String> repeated=arguments("cox-efron",temporary.resolve("repeated.tsv"),true);
        repeated.set(repeated.indexOf("--pheno")+1,phenotype.toString());
        repeated.addAll(List.of("--individual-id","group"));
        Result rejected=run(repeated);assertNotEquals(0,rejected.code);assertTrue(rejected.console.contains("distinct individuals"),rejected.console);
    }

    private static List<String> arguments(String model,Path out,boolean export) {
        boolean cox=model.startsWith("cox");
        List<String> args=new ArrayList<>(List.of("--omics",FIXTURE.resolve("variants.vcf").toString(),
            "--pheno",FIXTURE.resolve("phenotype.tsv").toString(),"--id","IID","--model",cox?"cox":"glm",
            "--formula",(cox?(model.contains("right")?"Surv(stop,event)":"Surv(start,stop,event)"):model.equals("binomial")?"binary":"count")+" ~ x + offset(o) + <omics>",
            "--threads","1","--out",out.toString()));
        if(cox)args.addAll(List.of("--ties",model.substring(model.lastIndexOf('-')+1)));else args.addAll(List.of("--family",model));
        if(export)args.addAll(List.of("--conditional-gwas-summary","--score-genome-build","GRCh38"));
        return args;
    }
    private static double value(DelimitedData table,String[] row,String column){return Double.parseDouble(row[table.column(column)]);}
    private static double[] numbers(String value){return Arrays.stream(value.split(",")).mapToDouble(Double::parseDouble).toArray();}
    private static Result run(List<String> args){ByteArrayOutputStream bytes=new ByteArrayOutputStream();PrintStream stream=new PrintStream(bytes);return new Result(JLinAlgCli.run(args.toArray(String[]::new),stream,stream),bytes.toString());}
    private record Result(int code,String console) { }
}
