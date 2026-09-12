/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;
import static org.junit.jupiter.api.Assertions.*;
import java.io.*;
import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.jlinalg.raremetal.RareMetalStudy;

class RareScoreRelatedTest {
    @TempDir Path dir;
    @Test void relatedExportMatchesBalancedRandomInterceptClosedForm()throws Exception {
        double[] y={-2.2,-1.8,-1.3,-.7,-.1,.1,.6,1.4,1.8,2.2,2.7,3.3};
        int[][] g={{0,1,0,0,1,0,0,2,0,1,0,0},{1,0,0,1,0,0,0,0,1,0,0,1}};
        StringBuilder pheno=new StringBuilder("sample\ty\n"),vcf=new StringBuilder("##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        for(int i=0;i<12;i++)vcf.append("\ts").append(i);vcf.append('\n');
        for(int i=11;i>=0;i--)pheno.append("s").append(i).append('\t').append(y[i]).append('\n');
        for(int j=0;j<2;j++){vcf.append("1\t").append(10+j*10).append("\tv").append(j).append("\tA\tG\t.\tPASS\t.\tGT");for(int call:g[j])vcf.append('\t').append(call==0?"0/0":call==1?"0/1":"1/1");vcf.append('\n');}
        StringBuilder grm=new StringBuilder("sample");for(int i=11;i>=0;i--)grm.append("\ts").append(i);grm.append('\n');
        for(int i=11;i>=0;i--){grm.append("s").append(i);for(int j=11;j>=0;j--)grm.append('\t').append(i/2==j/2?1:0);grm.append('\n');}
        Files.writeString(dir.resolve("cohort.vcf"),vcf);Files.writeString(dir.resolve("pheno"),pheno);Files.writeString(dir.resolve("grm"),grm);
        ByteArrayOutputStream error=new ByteArrayOutputStream();
        int status=JLinAlgCli.run(new String[]{"rare-score","--vcf",dir.resolve("cohort.vcf").toString(),"--pheno",dir.resolve("pheno").toString(),"--id","sample","--response","y","--grm",dir.resolve("grm").toString(),"--genome-build","test","--out",dir.resolve("result").toString()},new PrintStream(OutputStream.nullOutputStream()),new PrintStream(error));
        assertEquals(0,status,error.toString());
        // Balanced two-person random-intercept model: independent ANOVA REML.
        double residual=.86/6,genetic=3.5-residual/2,a=genetic+residual,b=genetic;
        double[][] inverse=new double[12][12];double denom=a*a-b*b;
        for(int i=0;i<12;i++){inverse[i][i]=a/denom;inverse[i][i^1]=-b/denom;}
        double[] py=new double[12];for(int i=0;i<12;i++)for(int j=0;j<12;j++)py[i]+=inverse[i][j]*(y[j]-.5);
        try(var study=new RareMetalStudy(dir.resolve("result.score.txt.gz"),dir.resolve("result.cov.txt.gz"),0,1e-10)) {
            var rows=study.region("1",10,20);List<RareMetalStudy.Score> records=List.of(rows.get(10L),rows.get(20L));
            double[] actual=study.covariance(records);
            for(int h=0;h<2;h++) {
                double expected=0;for(int i=0;i<12;i++)expected+=g[h][i]*py[i];
                assertEquals(expected,records.get(h).score(),Math.max(1,Math.abs(expected))*2e-6);
                assertTrue(Double.isFinite(records.get(h).hwePValue()));
                for(int j=0;j<2;j++) {
                    double expectedV=0;double meanH=Arrays.stream(g[h]).average().orElseThrow(),meanJ=Arrays.stream(g[j]).average().orElseThrow();
                    for(int l=0;l<12;l++)for(int r=0;r<12;r++)expectedV+=(g[h][l]-meanH)*inverse[l][r]*(g[j][r]-meanJ);
                    assertEquals(expectedV,actual[h*2+j],Math.max(1,Math.abs(expectedV))*2e-6);
                }
            }
        }
        assertTrue(Files.readString(dir.resolve("result.log")).contains("related_Gaussian_REML_variance_components"));
    }
    @Test void hweAvoidsRareTailAbsoluteToleranceAndRejectsDosageApproximation() {
        double[] calls=new double[100];Arrays.fill(calls,0,50,0);Arrays.fill(calls,50,100,2);
        var hwe=org.jlinalg.genetics.HardyWeinberg.calculate(calls,null,-1);
        assertTrue(hwe.pValue()<1e-20);
        calls[0]=.2;assertTrue(Double.isNaN(org.jlinalg.genetics.HardyWeinberg.calculate(calls,null,-1).pValue()));
    }
}
