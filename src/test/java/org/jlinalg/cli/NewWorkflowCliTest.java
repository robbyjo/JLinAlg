/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Locale;
import java.util.Random;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class NewWorkflowCliTest {
    @TempDir Path temporary;

    @Test
    void pcaAndCombatCommandsWriteAuditableArtifacts() throws Exception {
        Path matrix=temporary.resolve("omics.tsv"),pheno=temporary.resolve("pheno.tsv");
        StringBuilder omics=new StringBuilder("feature\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts8\n");
        for(int f=0;f<24;f++){omics.append("g").append(f);for(int s=0;s<8;s++)omics.append('\t').append((s<4?-2:3)+Math.sin(f+s*.2));omics.append('\n');}
        Files.writeString(matrix,omics);Files.writeString(pheno,"id\tbatch\tcase\ns1\tA\t0\ns2\tA\t1\ns3\tA\t0\ns4\tA\t1\ns5\tB\t0\ns6\tB\t1\ns7\tB\t0\ns8\tB\t1\n");
        Run pca=run("confounders","--method","pca","--omics",matrix.toString(),"--factors","2","--out",temporary.resolve("pca").toString(),"--write-adjusted");
        assertEquals(0,pca.status,pca.error);assertTrue(Files.exists(temporary.resolve("pca.factors.tsv")));assertTrue(Files.readString(temporary.resolve("pca.manifest.tsv")).contains("factors\t2"));
        Run combat=run("batch-adjust","--method","combat","--omics",matrix.toString(),"--pheno",pheno.toString(),"--id","id","--batch","batch","--preserve","case","--out",temporary.resolve("combat").toString());
        assertEquals(0,combat.status,combat.error);assertTrue(Files.exists(temporary.resolve("combat.adjusted.tsv")));assertTrue(combat.output.contains("aligned_samples=8"));
    }

    @Test
    void predictionAndIvCommandsExposeEarlierFeatures() throws Exception {
        Path prediction=temporary.resolve("prediction.tsv");StringBuilder p=new StringBuilder("y\tx\tz\n");for(int i=0;i<40;i++){double x=(i-20)/8.0;p.append(x+(i%3)*.1).append('\t').append(x).append('\t').append(i%2).append('\n');}Files.writeString(prediction,p);
        Run predicted=run("glm-predict","--input",prediction.toString(),"--response","y","--predictors","x,z","--family","gaussian","--estimand","difference","--scenario-column","z","--first","1","--second","0");
        assertEquals(0,predicted.status,predicted.error);assertTrue(predicted.output.contains("scenario_difference"));

        Path iv=temporary.resolve("iv.tsv");StringBuilder rows=new StringBuilder("y\tw\tx\tz1\tz2\n");Random random=new Random(4);for(int i=0;i<80;i++){double w=(i%5)-2,z1=Math.sin(i*.7),z2=Math.cos(i*.31),u=random.nextGaussian()*.2,x=.8*z1-.5*z2+.2*w+u,y=1+2*x+.3*w+u+random.nextGaussian()*.1;rows.append(y).append('\t').append(w).append('\t').append(x).append('\t').append(z1).append('\t').append(z2).append('\n');}Files.writeString(iv,rows);
        Run fitted=run("iv-regression","--input",iv.toString(),"--response","y","--exogenous","w","--endogenous","x","--instruments","z1,z2","--out",temporary.resolve("ivfit").toString(),"--backend","cpu");
        assertEquals(0,fitted.status,fitted.error);assertTrue(Files.exists(temporary.resolve("ivfit.coefficients.tsv")));assertTrue(Files.readString(temporary.resolve("ivfit.strength.tsv")).contains("partial_r2"));
    }

    @Test
    void arimaRegressionCommandWritesScalableSmoothing() throws Exception {
        Path input=temporary.resolve("series.tsv");StringBuilder text=new StringBuilder("y\tt\n");double state=0;for(int i=0;i<80;i++){state=.65*state+Math.sin(i*.8)*.2;text.append(1+.03*i+state).append('\t').append(i).append('\n');}Files.writeString(input,text);
        Run result=run("arima-regression","--input",input.toString(),"--response","y","--predictors","t","--order","1,0,0","--out",temporary.resolve("arima").toString(),"--smooth");
        assertEquals(0,result.status,result.error);assertTrue(Files.exists(temporary.resolve("arima.smoothed.tsv")));assertTrue(Files.readString(temporary.resolve("arima.model.tsv")).contains("effective_ar"));
    }

    @Test
    void allNewCommandsHaveHelp() {
        for(String command:new String[]{"confounders","batch-adjust","glm-predict","iv-regression","arima-regression"}){
            Run result=run(command,"--help");assertEquals(0,result.status,command+": "+result.error);assertTrue(result.output.contains("Usage:"));
        }
    }

    private static Run run(String...args){ByteArrayOutputStream out=new ByteArrayOutputStream(),err=new ByteArrayOutputStream();int status=JLinAlgCli.run(args,new PrintStream(out,true,StandardCharsets.UTF_8),new PrintStream(err,true,StandardCharsets.UTF_8));return new Run(status,out.toString(StandardCharsets.UTF_8),err.toString(StandardCharsets.UTF_8));}
    private record Run(int status,String output,String error){}
}
