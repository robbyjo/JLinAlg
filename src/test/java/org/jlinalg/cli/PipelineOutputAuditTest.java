package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import java.util.*;
import org.jlinalg.pipeline.OmicsAssociationEstimate;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class PipelineOutputAuditTest {
    @TempDir Path directory;
    @Test void commandsWithoutOverwriteOptionPreserveExistingOutputs() throws Exception {
        Path input=directory.resolve("beta.csv");
        Files.writeString(input,"y,x\n0.20,0\n0.30,1\n0.42,2\n0.55,3\n0.68,4\n0.78,5\n");
        Path output=directory.resolve("existing.tsv");Files.writeString(output,"retain me");
        var stream=new java.io.PrintStream(java.io.OutputStream.nullOutputStream());
        assertEquals(2,JLinAlgCli.run(new String[]{"beta-regression","--input",input.toString(),
            "--response","y","--mean","x","--out",output.toString()},stream,stream));
        assertEquals(2,JLinAlgCli.run(new String[]{"mr-estimate","--input","src/test/resources/r-reference/mr-conditional-input.tsv",
            "--method","ivw","--out",output.toString()},stream,stream));
        assertEquals("retain me",Files.readString(output));
    }
    @Test void cliSmokeAlignsShuffledPhenotypesAndEmitsOneRowPerFeature() throws Exception {
        Path output=directory.resolve("smoke.tsv");
        var stream=new java.io.PrintStream(java.io.OutputStream.nullOutputStream());
        assertEquals(0,JLinAlgCli.run(new String[]{"--pheno","src/test/resources/pipeline-audit/smoke-pheno.tsv",
            "--omics","src/test/resources/pipeline-audit/smoke-omics.tsv","--id","id","--formula","y ~ <omics>",
            "--model","ols","--backend","cpu","--omics-type","generic","--annot","src/test/resources/pipeline-audit/smoke-annotation.tsv",
            "--annot-cols","all","--out",output.toString(),"--threads","1","--no-log"},stream,stream));
        var rows=Files.readAllLines(output);assertEquals(4,rows.size());
        String[] header=rows.get(0).split("\t"),first=rows.get(1).split("\t",-1);
        assertEquals("ok",first[0]);assertEquals("signal",first[1]);
        assertFalse(Arrays.asList(header).contains("omics_type"));
        assertFalse(Arrays.asList(header).contains("position"));
        assertFalse(Arrays.asList(header).contains("statistic_type"));
        assertFalse(Arrays.asList(header).contains("df_method"));
        assertFalse(Arrays.asList(header).contains("partial_r2_method"));
        assertFalse(Arrays.asList(header).contains("filter_reason"));
        assertFalse(Arrays.asList(header).contains("error_type"));
        assertFalse(Arrays.asList(header).contains("message"));
        assertEquals(1.0,Double.parseDouble(first[Arrays.asList(header).indexOf("beta")]),.02);
        assertTrue(rows.get(2).startsWith("failed\t"));assertTrue(rows.get(3).startsWith("ok\t"));
        assertTrue(Files.readString(Path.of(output+".manifest.json")).contains("\"tested_features\": \"2\""));
    }
    @Test void cliRejectsInputOutputAndLogAliasesBeforeWriting() throws Exception {
        Path pheno=directory.resolve("pheno.tsv");String original="id\ty\na\t1\nb\t2\nc\t3\n";
        Files.writeString(pheno,original);
        for(String[] extra:new String[][]{{"--out",pheno.toString(),"--overwrite"},
                {"--out",directory.resolve("out.tsv").toString(),"--log",pheno.toString()}}) {
            List<String> args=new ArrayList<>(List.of("--pheno",pheno.toString(),"--id","id","--formula","y ~ 1","--backend","cpu"));
            args.addAll(List.of(extra));
            var stream=new java.io.PrintStream(java.io.OutputStream.nullOutputStream());
            assertNotEquals(0,JLinAlgCli.run(args.toArray(String[]::new),stream,stream));
            assertEquals(original,Files.readString(pheno));
        }
    }
    @Test void resumeMustNotCertifyAnUnverifiedExistingFile() throws Exception {
        Path output=directory.resolve("unrelated.tsv");Files.writeString(output,"not a result");
        var stream=new java.io.PrintStream(java.io.OutputStream.nullOutputStream());
        assertNotEquals(0,JLinAlgCli.run(new String[]{"--pheno",directory.resolve("absent.tsv").toString(),
            "--id","id","--formula","y ~ 1","--out",output.toString(),"--resume","--no-log"},stream,stream));
        assertEquals("not a result",Files.readString(output));
    }
    @Test void bhMatchesFrozenRWithExplicitTestedFamily() throws Exception {
        List<String> fixture;
        try(var reader=new java.io.BufferedReader(new java.io.InputStreamReader(
                getClass().getResourceAsStream("/pipeline-audit/bh-reference.tsv"),java.nio.charset.StandardCharsets.UTF_8))) {
            fixture=reader.lines().skip(1).toList();
        }
        Path target=directory.resolve("r.tsv");
        try(var bh=new ExternalBh(target,false,2)) {
            bh.writeHeader(List.of("id"));
            for(int i=0;i<fixture.size();i++)bh.write(List.of(Integer.toString(i)),Double.parseDouble(fixture.get(i).split("\t")[0]));
            bh.finish();
        }
        List<String> actual=Files.readAllLines(target);
        for(int i=0;i<fixture.size();i++) {
            double expected=Double.parseDouble(fixture.get(i).split("\t")[1]);
            String found=actual.get(i+1).split("\t",-1)[1];
            if(Double.isNaN(expected))assertEquals("",found);
            else assertEquals(expected,Double.parseDouble(found),2e-15);
        }
    }

    @Test void multilineFieldsAndHeadersPreserveRecordToQAlignment() throws Exception {
        Path target=directory.resolve("out.tsv");
        try(var bh=new ExternalBh(target,false)) {
            bh.writeHeader(List.of("id","a\tb"));
            bh.write(List.of("one","first\r\nsecond\n\"quoted\""),.01);
            bh.write(List.of("two","ok"),.04);
            bh.write(List.of("failed",""),Double.NaN);
            assertEquals(2,bh.tests()); bh.finish(); bh.finish();
        }
        String text=Files.readString(target);
        assertTrue(text.startsWith("id\t\"a\tb\""));
        assertTrue(text.contains("\"first\r\nsecond\n\"\"quoted\"\"\"\t0.02"));
        assertTrue(text.contains("two\tok\t0.04"));
        assertTrue(text.contains("failed\t\t"));
        try(var files=Files.list(directory)){assertEquals(1,files.count());}
    }

    @Test void csvExtensionUsesCommaDelimiterAndCsvQuotingThroughBh() throws Exception {
        Path target=directory.resolve("out.csv");
        try(var bh=new ExternalBh(target,false,2)) {
            bh.writeHeader(List.of("id","label"));
            bh.write(List.of("one","a,b"),.01);
            bh.write(List.of("two","quoted \"value\""),.04);
            bh.finish();
        }
        String text=Files.readString(target);
        assertTrue(text.startsWith("id,label,fdr_bh"));
        assertTrue(text.contains("one,\"a,b\",0.02"));
        assertTrue(text.contains("two,\"quoted \"\"value\"\"\",0.04"));
        assertFalse(text.startsWith("id\t"));
    }

    @Test void missingAnnotationIsBlankRatherThanFatal() throws Exception {
        Path annotation=directory.resolve("annotation.tsv");
        Files.writeString(annotation,"id\tgene\na\tA\n");
        Path target=directory.resolve("out.tsv");
        try(var sink=new CliResultSink(target,false,false,"t",
                AnnotationLookup.read(annotation,null,List.of("all")),null)) {
            sink.acceptEstimate(new OmicsAssociationEstimate("missing",1,.5,2,10,.07,-1,1));
            sink.finish();
        }
        String[] fields=Files.readAllLines(target).get(1).split("\t",-1);
        assertEquals("",fields[fields.length-2]);
        assertEquals(.07,Double.parseDouble(fields[fields.length-1]),0);
    }

    @Test void closeCleansScratchButRetainsRecoverablePartial() throws Exception {
        Path target=directory.resolve("out.tsv");
        try(var bh=new ExternalBh(target,false)) {
            bh.writeHeader(List.of("id")); bh.write(List.of("one"),.1);
        }
        try(var files=Files.list(directory)) {
            assertEquals(List.of("out.tsv.partial"),files.map(p->p.getFileName().toString()).toList());
        }
    }

    @Test void finishWithoutOverwriteDoesNotClobberConcurrentOutput() throws Exception {
        Path target=directory.resolve("out.tsv");
        try(var bh=new ExternalBh(target,false)) {
            bh.writeHeader(List.of("id"));bh.write(List.of("one"),.1);
            Files.writeString(target,"other job");
            assertThrows(FileAlreadyExistsException.class,bh::finish);
        }
        assertEquals("other job",Files.readString(target));
    }

    @Test void externalChunksMatchAnalyticBhWithTiesInvalidAndMissingRows() throws Exception {
        Path target=directory.resolve("large.tsv");
        int n=10031;
        double[] p=new double[n];
        for(int i=0;i<n;i++)p[i]=i%101==0?Double.NaN:(i%97)/100.0;
        double[] finite=Arrays.stream(p).filter(Double::isFinite).sorted().toArray();
        Map<Double,Double> q=new HashMap<>();double running=1;
        for(int i=finite.length-1;i>=0;i--){running=Math.min(running,finite[i]*finite.length/(i+1));q.put(finite[i],running);}
        try(var bh=new ExternalBh(target,false,101)) {
            bh.writeHeader(List.of("id"));
            for(int i=0;i<n;i++)bh.write(List.of(Integer.toString(i)),p[i]);
            bh.finish();assertEquals(finite.length,bh.tests());
        }
        List<String> lines=Files.readAllLines(target);assertEquals(n+1,lines.size());
        for(int i=0;i<n;i++) {
            String[] fields=lines.get(i+1).split("\t",-1);
            assertEquals(Integer.toString(i),fields[0]);
            if(Double.isNaN(p[i]))assertEquals("",fields[1]);
            else assertEquals(q.get(p[i]),Double.parseDouble(fields[1]),1e-15);
        }
    }
}
