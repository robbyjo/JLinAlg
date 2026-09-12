/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.*;
import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class RareMetaCliTest {
    @TempDir Path directory;
    private int run(String... args)throws Exception {
        ByteArrayOutputStream error=new ByteArrayOutputStream();
        int status=JLinAlgCli.run(args,new PrintStream(OutputStream.nullOutputStream()),new PrintStream(error));
        if(status!=0)System.err.println(error);
        return status;
    }
    @Test void upstreamRaremMetalFixtures()throws Exception {
        Path resources=Path.of("src/test/resources/raremetal").toAbsolutePath();
        Path manifest=directory.resolve("cohorts.tsv");
        Files.writeString(manifest,"cohort\tscores\tcovariance\nA\t"+resources.resolve("STUDY1.QT1.singlevar.score.txt.gz")+"\t"+resources.resolve("STUDY1.QT1.singlevar.cov.txt.gz")
            +"\nB\t"+resources.resolve("STUDY2.QT1.singlevar.score.txt.gz")+"\t"+resources.resolve("STUDY2.QT1.singlevar.cov.txt.gz")+"\n");
        Path prefix=directory.resolve("meta");
        assertEquals(0,run("rare-meta","--cohorts",manifest.toString(),"--out",prefix.toString(),"--test","single,burden,skat","--genome-build","tutorial",
            "--groups",resources.resolve("group.file").toString(),"--af-policy","raremetal","--hwe","0.00001","--call-rate","0.95"));
        var observed=table(Path.of(prefix+".burden.tsv"));
        for(String[] expected:reference(resources.resolve("COMBINED.QT1.meta.burden.results"))) {
            String[] actual=observed.get(expected[0]);assertNotNull(actual,expected[0]);
            assertEquals(Double.parseDouble(expected[9]),Double.parseDouble(actual[7]),1e-5,expected[0]+" beta");
            assertEquals(Double.parseDouble(expected[10]),Double.parseDouble(actual[10]),2e-5,expected[0]+" p");
        }
        observed=table(Path.of(prefix+".skat.tsv"));
        for(String[] expected:reference(resources.resolve("COMBINED.QT1.meta.SKAT_.results"))) {
            String[] actual=observed.get(expected[0]);assertNotNull(actual);
            assertEquals(Double.parseDouble(expected[9]),Double.parseDouble(actual[6]),.1,expected[0]+" Q");
            assertEquals(Double.parseDouble(expected[10]),Double.parseDouble(actual[7]),3e-5,expected[0]+" Davies p");
        }
        assertTrue(Files.readString(Path.of(prefix+".log")).contains("elapsed="));
        observed=table(Path.of(prefix+".single.tsv"));int matched=0;
        for(String[] expected:reference(resources.resolve("COMBINED.QT1.meta.singlevar.results"))) {
            // Upstream contains malformed unknown-ALT rows with shifted columns;
            // these are deliberately excluded rather than copied as associations.
            if(expected[7].equals("NA")||!expected[2].matches("[ACGT]+")||!expected[3].matches("[ACGT]+"))continue;
            String[] actual=observed.get(expected[0]+":"+expected[1]+":"+expected[2]+":"+expected[3]);
            int sign=1;
            if(actual==null){actual=observed.get(expected[0]+":"+expected[1]+":"+expected[3]+":"+expected[2]);sign=-1;}
            assertNotNull(actual,String.join(":",Arrays.copyOf(expected,4)));matched++;
            assertEquals(Double.parseDouble(expected[7]),sign*Double.parseDouble(actual[10]),Math.max(1e-5,Math.abs(Double.parseDouble(expected[7]))*2e-5));
            assertEquals(Double.parseDouble(expected[10]),Double.parseDouble(actual[13]),1e-5);
        }
        assertTrue(matched>500);
        Path mb=directory.resolve("mb");
        assertEquals(0,run("rare-meta","--cohorts",manifest.toString(),"--out",mb.toString(),"--test","burden","--weights","mb","--genome-build","tutorial",
            "--groups",resources.resolve("group.file").toString(),"--af-policy","raremetal","--hwe","0.00001","--call-rate","0.95"));
        observed=table(Path.of(mb+".burden.tsv"));
        for(String[] expected:reference(resources.resolve("COMBINED.QT1.meta.MB.results"))) {
            String[] actual=observed.get(expected[0]);assertNotNull(actual);
            assertEquals(Double.parseDouble(expected[9]),Double.parseDouble(actual[7]),1e-6);
            assertEquals(Double.parseDouble(expected[10]),Double.parseDouble(actual[10]),2e-5);
        }
    }
    private static List<String[]> reference(Path p)throws IOException{return Files.readAllLines(p).stream().filter(s->!s.startsWith("#")&&!s.isBlank()).map(s->s.split("\t")).toList();}
    private static Map<String,String[]> table(Path p)throws IOException {
        Map<String,String[]> result=new HashMap<>();for(String s:Files.readAllLines(p).subList(1,Files.readAllLines(p).size())){String[] f=s.split("\t");result.put(f[1],f);}return result;
    }
    @Test void selectedTestsMissingCohortsAndAlleleFlip()throws Exception {
        String header="##AnalyzedSamples=100\n#CHROM POS REF ALT N_INFORMATIVE ALL_AF CALL_RATE HWE_PVALUE U_STAT SQRT_V_STAT\n";
        Files.writeString(directory.resolve("a.score"),header+"1 10 A G 100 .01 1 1 2 2\n1 20 C T 100 .02 1 1 -3 3\n");
        Files.writeString(directory.resolve("b.score"),header+"1 10 G A 100 .99 1 1 -4 4\n");
        Files.writeString(directory.resolve("a.cov"),"#CHROM CURRENT_POS MARKERS_IN_WINDOW COV_MATRICES\n1 10 10,20 .04,.01\n1 20 20 .09\n");
        Files.writeString(directory.resolve("b.cov"),"#CHROM CURRENT_POS MARKERS_IN_WINDOW COV_MATRICES\n1 10 10 .16\n");
        Files.writeString(directory.resolve("cohorts.tsv"),"cohort\tscores\tcovariance\na\ta.score\ta.cov\nb\tb.score\tb.cov\n");
        Files.writeString(directory.resolve("groups"),"gene 1:10:A:G 1:20:C:T\n");
        Path out=directory.resolve("out");
        assertEquals(0,run("rare-meta","--cohorts",directory.resolve("cohorts.tsv").toString(),"--genome-build","GRCh38","--test","single,burden","--groups",directory.resolve("groups").toString(),"--out",out.toString()));
        var rows=table(Path.of(out+".single.tsv"));assertEquals("++",rows.get("1:10:A:G")[7]);assertEquals("-?",rows.get("1:20:C:T")[7]);
        assertEquals(.3,Double.parseDouble(rows.get("1:10:A:G")[10]),1e-15);
        String[] burden=table(Path.of(out+".burden.tsv")).get("gene");assertEquals(3.0/31,Double.parseDouble(burden[7]),1e-15);
        assertFalse(Files.exists(Path.of(out+".skat.tsv")));
        assertEquals(0,run("rare-meta","--cohorts",directory.resolve("cohorts.tsv").toString(),"--genome-build","GRCh38","--test","single","--min-cohorts","2","--out",directory.resolve("min").toString()));
        assertEquals("below_min_cohorts",table(directory.resolve("min.single.tsv")).get("1:20:C:T")[4]);
        assertEquals(0,run("rare-meta","--cohorts",directory.resolve("cohorts.tsv").toString(),"--genome-build","GRCh38","--test","burden","--window-size","10","--window-step","10","--out",directory.resolve("windows").toString()));
        assertEquals(Set.of("1:1-10","1:11-20"),table(directory.resolve("windows.burden.tsv")).keySet());
        Files.writeString(directory.resolve("a.cov"),"#CHROM CURRENT_POS MARKERS_IN_WINDOW COV_MATRICES\n1 10 10 .04\n1 20 20 .09\n");
        assertEquals(2,run("rare-meta","--cohorts",directory.resolve("cohorts.tsv").toString(),"--genome-build","GRCh38","--test","burden","--groups",directory.resolve("groups").toString(),"--out",directory.resolve("bad").toString()));
        assertFalse(Files.exists(directory.resolve("bad.burden.tsv")));
    }

    @Test void participantScoreExportRoundTripsIndexedSummaries()throws Exception {
        Path vcf=directory.resolve("cohort.vcf");
        Files.writeString(vcf,"##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            +"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\td\te\tf\n"
            +"1\t10\tv1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t0/0\t0/1\t1/1\t0/0\n"
            +"1\t20\tv2\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/0\t0/1\t0/0\t0/0\t0/0\n");
        Path pheno=directory.resolve("pheno.tsv");Files.writeString(pheno,"sample\ty\nf\t2\ne\t6\nd\t5\nc\t1\nb\t4\na\t3\n");
        assertEquals(0,run("rare-score","--vcf",vcf.toString(),"--pheno",pheno.toString(),"--id","sample","--response","y","--genome-build","GRCh38","--out",directory.resolve("cohort").toString()));
        try(var study=new org.jlinalg.raremetal.RareMetalStudy(directory.resolve("cohort.score.txt.gz"),directory.resolve("cohort.cov.txt.gz"),0,0)) {
            assertTrue(study.indexed());var rows=study.region("1",10,20);assertEquals(2,rows.size());
            // Independent intercept-only OLS: RSS=17.5, sigma^2=3.5.
            assertEquals(7/3.5,rows.get(10L).score(),1e-12);
            double[] covariance=study.covariance(List.of(rows.get(10L),rows.get(20L)));
            assertEquals((-4.0/3)/3.5,covariance[1],1e-12);
        }
        Files.writeString(directory.resolve("roundtrip.tsv"),"cohort\tscores\tcovariance\nA\tcohort.score.txt.gz\tcohort.cov.txt.gz\n");
        assertEquals(0,run("rare-meta","--cohorts",directory.resolve("roundtrip.tsv").toString(),"--genome-build","GRCh38","--test","single,burden,skat,skat-o","--maf","0.5","--weights","equal","--window-size","100","--out",directory.resolve("roundtrip").toString()));
    }
    @Test void rvtestsIndexedCovarianceLayout()throws Exception {
        Path root=Path.of("src/test/resources/raremetal").toAbsolutePath();
        Path manifest=directory.resolve("rvtests.tsv");
        Files.writeString(manifest,"cohort\tscores\tcovariance\nA\t"+root.resolve("STUDY1.rvtests.MetaScore.assoc.gz")+"\t"+root.resolve("STUDY1.rvtests.MetaCov.assoc.gz")
            +"\nB\t"+root.resolve("STUDY2.rvtests.MetaScore.assoc.gz")+"\t"+root.resolve("STUDY2.rvtests.MetaCov.assoc.gz")+"\n");
        assertEquals(0,run("rare-meta","--cohorts",manifest.toString(),"--genome-build","tutorial","--test","burden,skat","--af-policy","raremetal",
            "--groups",root.resolve("group.file").toString(),"--out",directory.resolve("rvtests").toString()));
        assertTrue(table(directory.resolve("rvtests.skat.tsv")).containsKey("GENE3"));
    }
    @Test void unknownAllelesDoNotHideAnInformativeLaterCohort()throws Exception {
        String header="##AnalyzedSamples=100\n#CHROM POS REF ALT N_INFORMATIVE ALL_AF U_STAT SQRT_V_STAT\n";
        Files.writeString(directory.resolve("a.score"),header+"1 10 A . 100 0 NA NA\n1 20 C . 100 .01 2 2\n");
        Files.writeString(directory.resolve("b.score"),header+"1 10 A G 100 .02 4 2\n1 20 C T 100 .02 4 2\n");
        Files.writeString(directory.resolve("cohorts.tsv"),"cohort\tscores\na\ta.score\nb\tb.score\n");
        assertEquals(0,run("rare-meta","--cohorts",directory.resolve("cohorts.tsv").toString(),"--genome-build","test",
            "--test","single,burden","--window-size","10","--out",directory.resolve("unknown").toString()));
        var single=table(directory.resolve("unknown.single.tsv"));
        assertEquals("?+",single.get("1:10:A:G")[7]);assertEquals("?+",single.get("1:20:C:T")[7]);
        var groups=table(directory.resolve("unknown.burden.tsv"));
        assertEquals(Set.of("1:1-10","1:11-20"),groups.keySet());
        assertEquals(1,Double.parseDouble(groups.get("1:1-10")[7]),1e-15);
    }
}
