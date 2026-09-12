/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;
import java.io.*;
import java.nio.file.*;
import java.util.*;

/** Dense LD block CLI benchmark; identical outputs gate every timed comparison. */
public final class RareMetaEngineeringBenchmark {
    private RareMetaEngineeringBenchmark() { }
    public static void main(String[] args)throws Exception {
        Path root=Files.createTempDirectory(Path.of("build"),"rare-engineering-");
        List<String> report=new ArrayList<>(List.of("variants,cohorts,groups,threads,cache_mb,repeat,seconds"));
        for(int n:new int[]{128,256}) {
            Path dir=Files.createDirectory(root.resolve("n"+n));StringBuilder manifest=new StringBuilder("cohort\tscores\tcovariance\n");
            for(int c=0;c<4;c++) {
                StringBuilder score=new StringBuilder("##AnalyzedSamples=1000\n#CHROM POS REF ALT N_INFORMATIVE ALL_AF CALL_RATE HWE_PVALUE U_STAT SQRT_V_STAT\n"),cov=new StringBuilder("#CHROM CURRENT_POS MARKERS_IN_WINDOW COV_MATRICES\n");
                for(int i=0;i<n;i++) {
                    score.append("1 ").append(i+1).append(" A G 1000 .01 1 1 ").append(.25*Math.sin(i+c)).append(" 1\n");
                    cov.append("1 ").append(i+1).append(' ');
                    for(int j=i;j<n;j++){if(j>i)cov.append(',');cov.append(j+1);}cov.append(' ');
                    for(int j=i;j<n;j++){if(j>i)cov.append(',');cov.append(Math.pow(.3,j-i)/1000);}cov.append('\n');
                }
                Files.writeString(dir.resolve(c+".score"),score);Files.writeString(dir.resolve(c+".cov"),cov);
                manifest.append("c").append(c).append('\t').append(c).append(".score\t").append(c).append(".cov\n");
            }
            Files.writeString(dir.resolve("cohorts"),manifest);StringBuilder groups=new StringBuilder();
            for(int g=0;g<6;g++){groups.append("g").append(g);for(int i=0;i<n;i++)groups.append(" 1:").append(i+1).append(":A:G");groups.append('\n');}
            Files.writeString(dir.resolve("groups"),groups);Map<String,String> reference=new HashMap<>();
            for(int repeat=0;repeat<3;repeat++)for(int mode=0;mode<3;mode++) {
                int threads=mode==2?2:1,cache=mode==0?0:16;Path out=dir.resolve("run"+repeat+"-"+mode);
                long start=System.nanoTime();ByteArrayOutputStream error=new ByteArrayOutputStream();
                int status=JLinAlgCli.run(new String[]{"rare-meta","--cohorts",dir.resolve("cohorts").toString(),"--groups",dir.resolve("groups").toString(),"--genome-build","synthetic","--test","burden,skat","--threads",""+threads,"--cache-mb",""+cache,"--out",out.toString()},new PrintStream(OutputStream.nullOutputStream()),new PrintStream(error));
                double seconds=(System.nanoTime()-start)/1e9;
                if(status!=0)throw new IllegalStateException(error.toString());
                for(String test:List.of("burden","skat","qc")) {
                    String value=Files.readString(Path.of(out+"."+test+".tsv"));
                    if(reference.putIfAbsent(test,value)!=null&&!reference.get(test).equals(value))throw new AssertionError("output differs: "+test);
                }
                String row=n+",4,6,"+threads+","+cache+","+repeat+","+seconds;report.add(row);System.out.println(row);
            }
        }
        Files.write(root.resolve("timings.csv"),report);System.out.println("Report: "+root.resolve("timings.csv"));
    }
}
