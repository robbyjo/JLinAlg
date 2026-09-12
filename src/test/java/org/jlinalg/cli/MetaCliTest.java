/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.zip.GZIPOutputStream;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class MetaCliTest {
    @Test void fixedAndAllRandomEstimatorsMatchIndependentMetaforFixtures() throws Exception {
        var reference=new Properties();
        try(var in=getClass().getResourceAsStream("/meta/cli-reference.properties")) { reference.load(in); }
        double[] y={.2,.5,.1,.7},se={.1,.2,.15,.25};
        List<Path> inputs=new ArrayList<>();
        for(int i=0;i<4;i++)inputs.add(file("ref"+i+".tsv","feature_id\tbeta\tse\ngene\t"+y[i]+"\t"+se[i]+"\n"));
        Path mods=file("ref-mods.tsv","cohort\tdose\nc1\t-1\nc2\t0\nc3\t1\nc4\t2\n");
        for(String kind:List.of("pool","regression"))for(String method:List.of("FE","REML","DL","PM")) {
            Path output=directory.resolve(kind+method+".tsv");
            List<String> args=command(kind.equals("pool")?"meta-analysis":"meta-regression",output,inputs,
                "--model",method.equals("FE")?"fixed":"random","--tau-estimator",method.equals("FE")?"reml":method.toLowerCase(Locale.ROOT));
            if(kind.equals("regression"))args.addAll(List.of("--moderator-file",mods.toString(),"--moderators","dose"));
            assertEquals(0,run(args.toArray(String[]::new)),error);
            int terms=kind.equals("pool")?1:2;
            for(int term=0;term<terms;term++) {
                var values=row(output,"gene",kind.equals("pool")?"pooled":term==0?"(Intercept)":"dose");
                for(String field:List.of("beta","se","p","tau2","q")) {
                    String column=field.equals("p")?"p_value":field.equals("tau2")?"tau_squared":field;
                    int referenceTerm=field.equals("tau2")||field.equals("q")?0:term;
                    double expected=Double.parseDouble(reference.getProperty(kind+"."+method+"."+field+"."+referenceTerm));
                    assertEquals(expected,Double.parseDouble(values.get(column)),3e-7,kind+" "+method+" "+field);
                }
            }
        }
    }
    @TempDir Path directory;
    private String error;
    private int run(String... args) {
        var errors=new ByteArrayOutputStream();
        int code=JLinAlgCli.run(args,new PrintStream(new ByteArrayOutputStream()),new PrintStream(errors));
        error=errors.toString(); return code;
    }
    private Path file(String name,String content) throws IOException {
        Path path=directory.resolve(name); Files.writeString(path,content); return path;
    }
    private List<String> command(String command,Path output,List<Path> inputs,String... options) {
        var args=new ArrayList<String>(List.of(command,"--out",output.toString(),"--sort-chunk-rows","1","--block-rows","1"));
        for(int i=0;i<inputs.size();i++)args.addAll(List.of("--cohort","c"+(i+1)+"="+inputs.get(i)));
        args.addAll(List.of(options));return args;
    }
    private Map<String,String> row(Path output,String feature,String term) throws IOException {
        var table=DelimitedData.read(output);
        for(String[] row:table.rows())if(row[0].equals(feature)&&row[1].equals(term)) {
            Map<String,String> result=new HashMap<>();
            for(int i=0;i<row.length;i++)result.put(table.header().get(i),row[i]);
            return result;
        }
        fail("missing result: "+feature+" "+term);return null;
    }
    @Test void joinsUnsortedCohortsAndPreservesDirectionAndSingleCohortSemantics() throws Exception {
        List<Path> inputs=new ArrayList<>();
        for(int i=0;i<6;i++) {
            String text="feature_id\tbeta\tse\n";
            if(i==0)text+="z\t0\t.2\nsingle\t.4\t.1\n";
            if(i!=4)text+="a\t"+(i<2?1:-1)+"\t1\n";
            else text+="a\tNA\t1\n";
            inputs.add(file("c"+i+".tsv",text));
        }
        Path output=directory.resolve("result.tsv");
        assertEquals(0,run(command("meta-analysis",output,inputs,"--model","fixed").toArray(String[]::new)),error);
        var fit=row(output,"a","pooled");
        assertEquals("++--?-",fit.get("direction"));assertEquals("5",fit.get("n_cohorts"));
        assertEquals(-.2,Double.parseDouble(fit.get("beta")),1e-14);
        assertEquals(1/Math.sqrt(5),Double.parseDouble(fit.get("se")),1e-14);
        assertEquals("single_cohort",row(output,"single","pooled").get("status"));
        assertEquals("NA",row(output,"single","pooled").get("tau_squared"));
        assertEquals("0?????",row(output,"z","pooled").get("direction"));
        assertTrue(Files.readString(Path.of(output+".log")).contains("cohort_order=c1,c2,c3,c4,c5,c6"));
        assertEquals("a",DelimitedData.read(output).rows().get(0)[0]);
        Path filtered=directory.resolve("filtered.tsv");
        assertEquals(0,run(command("meta-analysis",filtered,inputs,"--min-cohorts","2").toArray(String[]::new)),error);
        assertEquals("below_min_cohorts",row(filtered,"single","pooled").get("status"));
        assertEquals("NA",row(filtered,"single","pooled").get("beta"));
    }
    @Test void randomEffectsRecoverAnalyticEqualVarianceHeterogeneity() throws Exception {
        List<Path> inputs=new ArrayList<>();
        for(int i=0;i<5;i++)inputs.add(file("r"+i+".tsv","feature_id\tbeta\tse\ngene\t"+(i-2)+"\t.1\n"));
        for(String estimator:List.of("reml","pm","dl")) {
            Path output=directory.resolve(estimator+".tsv");
            assertEquals(0,run(command("meta-analysis",output,inputs,"--tau-estimator",estimator).toArray(String[]::new)),error);
            var fit=row(output,"gene","pooled");
            assertEquals(2.49,Double.parseDouble(fit.get("tau_squared")),2e-6);
            assertEquals(Math.sqrt(.5),Double.parseDouble(fit.get("se")),3e-7);
            assertEquals("--0++",fit.get("direction"));
        }
    }
    @Test void regressionAlignsModeratorNamesAndReportsUnestimableRows() throws Exception {
        List<Path> inputs=new ArrayList<>();
        double[] effects={.2,.5,.1,.7},ses={.1,.2,.15,.25};
        for(int i=0;i<4;i++)inputs.add(file("g"+i+".tsv","feature_id\tbeta\tse\na\t"+effects[i]+"\t"+ses[i]+"\n"
            +(i<2?"few\t1\t.1\n":"")));
        Path mods=file("mods.tsv","cohort\tdose\nc4\t2\nc2\t0\nc1\t-1\nc3\t1\n");
        Path output=directory.resolve("reg.tsv");
        assertEquals(0,run(command("meta-regression",output,inputs,"--moderator-file",mods.toString(),
            "--moderators","dose","--tau-estimator","dl").toArray(String[]::new)),error);
        assertEquals(.3004978,Double.parseDouble(row(output,"a","(Intercept)").get("beta")),1e-7);
        assertEquals(.09178625,Double.parseDouble(row(output,"a","dose").get("beta")),1e-7);
        assertEquals("insufficient_cohorts",row(output,"few","NA").get("status"));
        Path constant=file("constant.tsv","cohort\tdose\nc1\t1\nc2\t1\nc3\t1\nc4\t1\n");
        Path deficient=directory.resolve("rank.tsv");
        assertEquals(0,run(command("meta-regression",deficient,inputs,"--moderator-file",constant.toString(),
            "--moderators","dose").toArray(String[]::new)),error);
        assertEquals("rank_deficient",row(deficient,"a","NA").get("status"));
    }
    @Test void readsGzipCsvAndEscapesFeatureIdsInCsvOutput() throws Exception {
        Path input=directory.resolve("cohort.csv.gz");
        try(var out=new GZIPOutputStream(Files.newOutputStream(input))) {
            out.write("id,b,se\n\"gene,one\",.2,.1\n".getBytes(java.nio.charset.StandardCharsets.UTF_8));
        }
        Path output=directory.resolve("quoted.csv");
        assertEquals(0,run(command("meta-analysis",output,List.of(input),"--id-column","id","--effect-column","b").toArray(String[]::new)),error);
        assertEquals("single_cohort",row(output,"gene,one","pooled").get("status"));
    }
    @Test void mergesMoreThanOneFanInWithoutKeepingAllFeaturesInMemory() throws Exception {
        StringBuilder text=new StringBuilder("feature_id\tbeta\tse\n");
        for(int i=69;i>=0;i--)text.append(String.format(Locale.ROOT,"f%03d\t.1\t.2\n",i));
        Path input=file("many.tsv",text.toString()),output=directory.resolve("many-out.tsv");
        assertEquals(0,run(command("meta-analysis",output,List.of(input)).toArray(String[]::new)),error);
        var rows=DelimitedData.read(output).rows();assertEquals(70,rows.size());
        assertEquals("f000",rows.get(0)[0]);assertEquals("f069",rows.get(69)[0]);
    }
    @Test void rejectsDuplicateIdsBadNumbersAndExistingOutputsWithoutPublishingPartialResults() throws Exception {
        for(String rows:List.of("a\t.1\t.2\nb\t.2\t.2\na\t.3\t.2\n","a\t.1\t0\n","a\toops\t.1\n")) {
            Path input=file("bad.tsv","feature_id\tbeta\tse\n"+rows),output=directory.resolve("bad-out.tsv");
            assertEquals(2,run(command("meta-analysis",output,List.of(input)).toArray(String[]::new)));
            assertFalse(Files.exists(output));assertFalse(Files.exists(Path.of(output+".log")));
            try(var files=Files.list(directory)) { assertFalse(files.anyMatch(p->p.getFileName().toString().startsWith(".jlinalg-meta-"))); }
        }
        Path input=file("valid.tsv","feature_id\tbeta\tse\na\t.1\t.2\n"),output=file("existing.tsv","preserve me");
        assertEquals(2,run(command("meta-analysis",output,List.of(input)).toArray(String[]::new)));
        assertEquals("preserve me",Files.readString(output));
        assertEquals(2,run("meta-analysis","--cohort","a="+input,"--out",directory.resolve("x").toString(),"--min-cohorts","0"));
        assertEquals(0,run("meta-analysis","--help"));assertEquals(0,run("meta-regression","--help"));
    }
}
