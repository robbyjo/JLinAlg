/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;
import static org.junit.jupiter.api.Assertions.*;
import java.io.*;
import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class RareMetaAdvancedCliTest {
    @TempDir Path dir;
    private void fixture()throws Exception {
        String h="##AnalyzedSamples=100\n#CHROM POS REF ALT N_INFORMATIVE ALL_AF CALL_RATE HWE_PVALUE U_STAT SQRT_V_STAT\n";
        Files.writeString(dir.resolve("a.score"),h+"1 10 A G 100 .01 1 1 2 2\n1 20 C T 100 .02 1 1 -3 3\n1 30 A C 100 .1 1 1 4 2\n");
        Files.writeString(dir.resolve("b.score"),h+"1 10 G A 100 .99 1 1 -4 2\n1 20 C T 100 .02 1 1 1 3\n1 30 C A 100 .9 1 1 -2 2\n");
        Files.writeString(dir.resolve("a.cov"),"#CHROM CURRENT_POS MARKERS_IN_WINDOW COV_MATRICES\n1 10 10,20,30 .04,.01,.01\n1 20 20,30 .09,.02\n1 30 30 .04\n");
        Files.writeString(dir.resolve("b.cov"),"#CHROM CURRENT_POS MARKERS_IN_WINDOW COV_MATRICES\n1 10 10,20,30 .04,-.01,.01\n1 20 20,30 .09,-.02\n1 30 30 .04\n");
        Files.writeString(dir.resolve("cohorts"),"cohort\tscores\tcovariance\na\ta.score\ta.cov\nb\tb.score\tb.cov\n");
        Files.writeString(dir.resolve("groups"),"gene 1:10:A:G 1:20:C:T\none 1:10:A:G\n");
    }
    private int run(String prefix,String tests,String... extra) {
        List<String> args=new ArrayList<>(List.of("rare-meta","--cohorts",dir.resolve("cohorts").toString(),"--groups",dir.resolve("groups").toString(),"--genome-build","test","--test",tests,"--weights","equal","--out",dir.resolve(prefix).toString()));
        args.addAll(List.of(extra));var error=new ByteArrayOutputStream();
        int code=JLinAlgCli.run(args.toArray(String[]::new),new PrintStream(OutputStream.nullOutputStream()),new PrintStream(error));
        if(code!=0)System.err.println(error);return code;
    }
    private List<Map<String,String>> rows(String file)throws Exception {
        List<String> lines=Files.readAllLines(dir.resolve(file));String[] header=lines.get(0).split("\t");List<Map<String,String>> result=new ArrayList<>();
        for(String line:lines.subList(1,lines.size())) {String[] f=line.split("\t",-1);assertEquals(header.length,f.length,line);Map<String,String> row=new HashMap<>();for(int i=0;i<f.length;i++)row.put(header[i],f[i]);result.add(row);}return result;
    }
    @Test void allTestsDiagnosticsAndParallelOutputAreConsistent()throws Exception {
        fixture();String tests="burden,skat,skat-o,vt,het-skat,het-skat-o,burden-fixed,burden-random";
        for(String prefix:List.of("serial","parallel"))assertEquals(0,run(prefix,tests,"--leave-variant-out","--leave-cohort-out","--cohort-results","--simulations","2000","--threads",prefix.equals("serial")?"1":"2","--cache-mb",prefix.equals("serial")?"0":"16"));
        for(String test:(tests+",qc").split(","))assertEquals(Files.readString(dir.resolve("serial."+test+".tsv")),Files.readString(dir.resolve("parallel."+test+".tsv")),test);
        var burden=rows("serial.burden.tsv");
        var omitted=burden.stream().filter(r->r.get("scope").equals("leave_variant:1:10:A:G")&&r.get("feature_id").equals("gene")).findFirst().orElseThrow();
        assertEquals(-2.0/18,Double.parseDouble(omitted.get("beta")),1e-14);
        assertTrue(burden.stream().anyMatch(r->r.get("feature_id").equals("one")&&r.get("status").equals("no_informative_variants")));
        var het=rows("serial.het-skat.tsv").stream().filter(r->r.get("scope").equals("meta")&&r.get("feature_id").equals("gene")).findFirst().orElseThrow();
        assertEquals("2",het.get("n_variants"));assertEquals("4",het.get("n_effect_dimensions"));assertEquals(30,Double.parseDouble(het.get("q")));
        for(String t:tests.split(","))rows("serial."+t+".tsv");
        assertTrue(rows("serial.qc.tsv").stream().allMatch(r->r.get("n_samples").equals("100.0")));
    }
    @Test void conditionalAnalysisOrientsBothAxesAndHandlesMissingConditions()throws Exception {
        fixture();Files.writeString(dir.resolve("conditions"),"gene 1:30:A:C\none 1:30:A:C\n");
        assertEquals(0,run("conditional","burden,skat-o","--condition",dir.resolve("conditions").toString(),"--skato-calibration","deterministic"));
        var gene=rows("conditional.burden.tsv").stream().filter(r->r.get("feature_id").equals("gene")).findFirst().orElseThrow();
        // Cohort Schur score sum: 4 - (1+2)*(4+2)/4 = -.5.
        // Information sum: 2*(4+1+1+9) - 2*(1+2)^2/4 = 25.5.
        assertEquals(-.5/25.5,Double.parseDouble(gene.get("beta")),1e-14);
        Files.writeString(dir.resolve("b.score"),Files.readString(dir.resolve("b.score")).replace("1 30 C A 100 .9 1 1 -2 2\n",""));
        assertEquals(2,run("missing","burden","--condition",dir.resolve("conditions").toString()));
        assertFalse(Files.exists(dir.resolve("missing.burden.tsv")));
        assertEquals(0,run("exclude","burden","--condition",dir.resolve("conditions").toString(),"--condition-missing","exclude"));
        assertTrue(rows("exclude.burden.tsv").stream().allMatch(r->r.get("n_cohorts").equals("1")));
        assertTrue(rows("exclude.qc.tsv").stream().anyMatch(r->r.get("reason").equals("missing_condition")));
    }
    @Test void unknownConditionGroupsAndAbsentCrossCovarianceFail()throws Exception {
        fixture();Files.writeString(dir.resolve("conditions"),"typo 1:30:A:C\n");
        assertEquals(2,run("typo","burden","--condition",dir.resolve("conditions").toString()));
        Files.writeString(dir.resolve("conditions"),"gene 1:30:A:C\n");
        Files.writeString(dir.resolve("a.cov"),Files.readString(dir.resolve("a.cov")).replace("10,20,30 .04,.01,.01","10,20 .04,.01"));
        assertEquals(2,run("cross","burden","--condition",dir.resolve("conditions").toString()));
    }
}
