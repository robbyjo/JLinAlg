/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.*;
import java.nio.file.*;
import java.util.*;
import org.jlinalg.xwas.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class XwasWorkflowTest {
    @TempDir Path directory;
    private static final List<String> TRAITS=List.of("A","B","C","D");
    private int run(String... args) {
        ByteArrayOutputStream err=new ByteArrayOutputStream();
        int status=JLinAlgCli.run(args,new PrintStream(new ByteArrayOutputStream()),new PrintStream(err));
        if(status!=0)System.err.println(err);return status;
    }
    private static double[] expected(String name)throws IOException {
        return DelimitedData.read(Path.of("src/test/resources/xwas/"+name)).rows().stream().mapToDouble(r->Double.parseDouble(r[0])).toArray();
    }
    private static double[] matrix(String name,List<String> labels)throws IOException {
        return XwasFiles.matrix(Path.of("examples/xwas/"+name),labels);
    }
    @Test void ldscMatchesIndependentRQrAndCommonDeleteCovariance()throws Exception {
        Path out=directory.resolve("ldsc.tsv");
        assertEquals(0,run("ldsc","--input","examples/xwas/ldsc.tsv","--traits","A,B,C,D","--reference-variants","100000","--blocks","30","--out",out.toString()));
        assertArrayEquals(expected("ldsc-S.tsv"),XwasFiles.matrix(Path.of(out+".S.tsv"),TRAITS),1e-11);
        assertArrayEquals(expected("ldsc-V.tsv"),XwasFiles.matrix(Path.of(out+".V.tsv"),XwasFiles.pairs(TRAITS)),1e-11);
        assertArrayEquals(expected("ldsc-intercepts.tsv"),XwasFiles.matrix(Path.of(out+".intercepts.tsv"),TRAITS),1e-11);
        assertTrue(Files.readString(Path.of(out+".log")).contains("status=complete"));
        String prior=Files.readString(out);
        assertEquals(2,run("ldsc","--input","examples/xwas/ldsc.tsv","--traits","A,B,C,D","--reference-variants","100000","--blocks","30","--out",out.toString()));
        assertEquals(prior,Files.readString(out));
        // Actual LDSC outputs are admitted by the factor workflow without a format conversion.
        assertEquals(0,run("genomic-factor","--s",out+".S.tsv","--v",out+".V.tsv","--traits","A,B,C,D","--out",directory.resolve("from-ldsc.tsv").toString()));
    }
    @Test void molecularCliAlignsReversedAllelesAndMatchesRJointTest()throws Exception {
        for(String command:List.of("twas","pwas")) {
            Path out=directory.resolve(command+".tsv");
            assertEquals(0,run(command,"--gwas","examples/xwas/gwas.tsv","--weights","examples/xwas/weights.tsv","--reference","examples/xwas/reference.tsv","--ld","examples/xwas/ld.tsv","--joint","true","--out",out.toString()));
            var t=DelimitedData.read(out);String[] r=t.rows().get(0);double[] e=expected("twas.tsv");
            assertEquals(e[0],XwasFiles.number(t,r,"z"),1e-12);assertEquals(e[1],XwasFiles.number(t,r,"p"),1e-12);assertEquals(e[2],XwasFiles.number(t,r,"predicted_variance"),1e-12);
            var j=DelimitedData.read(Path.of(out+".joint.tsv"));assertEquals(expected("twas-joint.tsv")[0],XwasFiles.number(j,j.rows().get(0),"chi_square"),1e-12);
        }
    }
    @Test void factorFitAndSamplingOverlapMatchIndependentROptim()throws Exception {
        var fit=GenomicFactor.fit(matrix("factor-S.tsv",TRAITS),4,matrix("factor-V.tsv",XwasFiles.pairs(TRAITS)));
        double[] expected=expected("factor-fit.tsv");
        assertArrayEquals(Arrays.copyOfRange(expected,0,4),fit.loadings(),2e-6);
        assertArrayEquals(Arrays.copyOfRange(expected,4,8),fit.residualVariances(),2e-6);
        assertEquals(expected[8],fit.chiSquare(),1e-8);
        assertArrayEquals(expected("factor-covariance.tsv"),fit.parameterCovariance(),1e-8);
        Path out=directory.resolve("factor.tsv");
        assertEquals(0,run("genomic-factor","--s","examples/xwas/factor-S.tsv","--v","examples/xwas/factor-V.tsv","--traits","A,B,C,D","--gwas","examples/xwas/factor-gwas.tsv","--sampling-correlation","examples/xwas/sampling-correlation.tsv","--out",out.toString()));
        var t=DelimitedData.read(Path.of(out+".gwas.tsv"));String[] row=t.rows().get(0);double[] e=expected("factor-association.tsv");
        assertEquals(e[0],XwasFiles.number(t,row,"beta"),1e-6);assertEquals(e[1],XwasFiles.number(t,row,"se"),1e-6);assertEquals(e[2],XwasFiles.number(t,row,"q_snp"),1e-5);
    }
    @Test void trainedScoreRoundTripsAndEvaluatesIndependentCohort()throws Exception {
        Path weights=directory.resolve("score.tsv"),out=directory.resolve("predictions.tsv");
        assertEquals(0,run("score-train","--input","examples/xwas/train.tsv","--features","x1,x2","--outcome","y","--lambdas","0.1","--out",weights.toString()));
        double[] e=expected("score.tsv");var w=DelimitedData.read(weights);
        assertEquals(e[0],XwasFiles.number(w,w.rows().get(0),"intercept"),1e-8);
        assertEquals(e[1],XwasFiles.number(w,w.rows().get(0),"weight"),1e-8);assertEquals(e[2],XwasFiles.number(w,w.rows().get(1),"weight"),1e-8);
        assertEquals(0,run("score-apply","--input","examples/xwas/test.tsv","--weights",weights.toString(),"--outcome","y","--out",out.toString()));
        var t=DelimitedData.read(Path.of(out+".evaluation.tsv"));String[] r=t.rows().get(0);
        assertEquals(e[4],XwasFiles.number(t,r,"rmse"),1e-8);assertEquals(e[5],XwasFiles.number(t,r,"predictive_r2"),1e-8);
        assertEquals(e[6],XwasFiles.number(t,r,"calibration_intercept"),1e-8);assertEquals(e[7],XwasFiles.number(t,r,"calibration_slope"),1e-8);
        Path bad=directory.resolve("leakage.tsv");
        assertEquals(2,run("score-apply","--input","examples/xwas/train.tsv","--weights",weights.toString(),"--outcome","y","--out",bad.toString()));
        assertFalse(Files.exists(bad));
    }
    @Test void trainingCvIsReproducibleAndNeverReadsHeldOutOutcomes()throws Exception {
        Path a=directory.resolve("a.tsv"),b=directory.resolve("b.tsv");
        for(Path out:List.of(a,b))assertEquals(0,run("score-train","--input","examples/xwas/train.tsv","--features","x1,x2","--outcome","y","--lambdas","1,0.1,0.01","--folds","4","--out",out.toString()));
        assertEquals(Files.readString(a),Files.readString(b));
        assertEquals(2,run("score-train","--input","examples/xwas/train.tsv","--features","x1,y","--outcome","y","--lambdas","1","--out",directory.resolve("bad.tsv").toString()));
    }
    @Test void externalPolygenicScoreUsesTwoMinusDosageOnSwap()throws Exception {
        Path weights=directory.resolve("weights.tsv"),input=directory.resolve("dosage.tsv"),alleles=directory.resolve("alleles.tsv"),out=directory.resolve("scores.tsv");
        Files.writeString(weights,"variant\tea\toa\tweight\nv1\tA\tC\t2\n");
        Files.writeString(input,"sample\tv1\ty\na\t0\t1\nb\t1\t0\nc\t2\t0\nd\t0\t1\n");
        Files.writeString(alleles,"variant\tea\toa\nv1\tC\tA\n");
        assertEquals(0,run("score-apply","--input",input.toString(),"--weights",weights.toString(),"--alleles",alleles.toString(),"--outcome","y","--family","binary","--out",out.toString()));
        var result=DelimitedData.read(out);assertEquals(4,XwasFiles.number(result,result.rows().get(0),"score"));assertEquals(0,XwasFiles.number(result,result.rows().get(2),"score"));
        var eval=DelimitedData.read(Path.of(out+".evaluation.tsv"));assertEquals(1,XwasFiles.number(eval,eval.rows().get(0),"auc"));
    }
    @Test void degenerateMatricesAndInvalidAllelesFailExplicitly() {
        assertThrows(IllegalArgumentException.class,()->PredictedOmics.test(new double[]{1,1},new double[]{1,-1},new double[]{1,1},new double[]{1,1,1,1}));
        assertThrows(IllegalArgumentException.class,()->PredictedOmics.test(new double[]{1,1},new double[]{1,1},new double[]{1,1},new double[]{1,2,2,1}));
        assertThrows(IllegalArgumentException.class,()->PredictedOmics.joint(new double[]{1,1},new double[][]{{1,0},{2,0}},new double[]{1,1},new double[]{1,0,0,1}));
        assertThrows(IllegalArgumentException.class,()->PredictedOmics.alleleSign("A","C","A","G"));
        assertThrows(IllegalArgumentException.class,()->GenomicFactor.associate(new double[]{.1,.2},new double[]{1,1,1,1},new double[]{1,1}));
        assertThrows(IllegalArgumentException.class,()->GenomicFactor.fit(new double[]{1,0,0,1},2,new double[9]));
    }
    @Test void failedMolecularCoveragePublishesNoResult()throws Exception {
        Path bad=directory.resolve("gwas.tsv"),out=directory.resolve("bad.tsv");Files.writeString(bad,"variant\tea\toa\tz\nv1\tA\tC\t2\n");
        assertEquals(2,run("twas","--gwas",bad.toString(),"--weights","examples/xwas/weights.tsv","--reference","examples/xwas/reference.tsv","--ld","examples/xwas/ld.tsv","--out",out.toString()));assertFalse(Files.exists(out));
        Files.writeString(bad,"variant\tea\toa\tz\nv1\tA\tC\t2\nv1\tA\tC\t3\n");
        assertEquals(2,run("twas","--no-log","--gwas",bad.toString(),"--weights","examples/xwas/weights.tsv","--reference","examples/xwas/reference.tsv","--ld","examples/xwas/ld.tsv","--out",out.toString()));
    }
    @Test void scoreAucHandlesTiesAndR2CanBeNegative() {
        var e=PredictionScores.evaluate(new double[]{0,1,0,1},new double[]{1,1,1,1},new double[]{.5,.5,.5,.5},true);
        assertEquals(.5,e.auc());assertTrue(Double.isNaN(e.calibrationSlope()));
        var poor=PredictionScores.evaluate(new double[]{1,2,3,4},new double[]{9,8,7,6},new double[]{2.5,2.5,2.5,2.5},false);
        assertTrue(poor.predictiveR2()<0);
    }
    @Test void matrixLabelsAreAlignedAndSidecarCollisionsPreserveExistingFiles()throws Exception {
        Path permuted=directory.resolve("permuted.tsv");Files.writeString(permuted,"row\tb\ta\na\t.2\t1\nb\t1\t.2\n");
        assertArrayEquals(new double[]{1,.2,.2,1},XwasFiles.matrix(permuted,List.of("a","b")));
        Path out=directory.resolve("m.tsv"),sidecar=Path.of(out+".metadata.tsv");Files.writeString(sidecar,"preserve");
        assertEquals(2,run("twas","--gwas","examples/xwas/gwas.tsv","--weights","examples/xwas/weights.tsv","--reference","examples/xwas/reference.tsv","--ld","examples/xwas/ld.tsv","--out",out.toString()));
        assertFalse(Files.exists(out));assertEquals("preserve",Files.readString(sidecar));
    }
    @Test void helpIsWiredForAllCommands() {
        for(String name:List.of("ldsc","twas","pwas","genomic-factor","score-train","score-apply"))assertEquals(0,run(name,"--help"));
    }
    @Test void ldscVariableSampleSizesRecoverAnalyticSlopeAndNegativeH2IsUnavailable() {
        int n=24;double[] ld=new double[n],wld=new double[n];double[][] z=new double[n][1],ns=new double[n][1];int[] block=new int[n];
        for(int i=0;i<n;i++){ld[i]=1+i%7;wld[i]=2+i%4;ns[i][0]=1000+31*i;block[i]=i/4;z[i][0]=Math.sqrt(1.1+.3*ns[i][0]*ld[i]/10000);}
        var fit=LdScoreRegression.fit(ld,wld,z,ns,block,10000);
        assertEquals(.3,fit.geneticCovariance()[0],1e-12);assertEquals(1.1,fit.intercepts()[0],1e-12);
        for(int i=0;i<n;i++)z[i][0]=Math.sqrt(1-.1*ns[i][0]*ld[i]/10000);
        var negative=LdScoreRegression.fit(ld,wld,z,ns,block,10000);
        assertEquals(-.1,negative.geneticCovariance()[0],1e-12);assertTrue(Double.isNaN(negative.geneticCorrelations()[0]));
        ns[0][0]=0;assertThrows(IllegalArgumentException.class,()->LdScoreRegression.fit(ld,wld,z,ns,block,10000));
    }
    @Test void molecularTestIsInvariantToPositiveWeightScaleAndTracksOrientation() {
        double[] z={2,-1,3},sd={.6,.7,.5},r={1,.2,-.1,.2,1,.3,-.1,.3,1};
        var a=PredictedOmics.test(z,new double[]{.4,-.2,.1},sd,r);
        var b=PredictedOmics.test(z,new double[]{4,-2,1},sd,r);
        var c=PredictedOmics.test(z,new double[]{-.4,.2,-.1},sd,r);
        assertEquals(a.z(),b.z(),1e-14);assertEquals(a.z(),-c.z(),1e-14);assertEquals(a.pValue(),c.pValue(),1e-14);
        assertThrows(IllegalArgumentException.class,()->PredictedOmics.alleleSign("a","c","a","c"));
        assertThrows(IllegalArgumentException.class,()->XwasFiles.id("a\"b"));
    }
    @Test void trainedGenotypeToMolecularWeightsFeedTwasWithoutConversion()throws Exception {
        Path train=directory.resolve("genetic.tsv"),alleles=directory.resolve("alleles.tsv"),weights=directory.resolve("molecular.tsv"),out=directory.resolve("twas.tsv");
        StringBuilder table=new StringBuilder("sample\texpression\tv1\tv2\tv3\n");
        for(int i=0;i<30;i++){double a=i%3,b=(i/3)%3,c=(i/9)%3;table.append("s").append(i).append('\t').append(1+.3*a-.2*b+.1*c+.01*Math.sin(i)).append('\t').append(a).append('\t').append(b).append('\t').append(c).append('\n');}
        Files.writeString(train,table);Files.writeString(alleles,"variant\tea\toa\nv1\tA\tC\nv2\tA\tC\nv3\tA\tC\n");
        assertEquals(0,run("score-train","--input",train.toString(),"--features","v1,v2,v3","--outcome","expression","--lambdas","0.01","--alleles",alleles.toString(),"--out",weights.toString()));
        assertEquals(0,run("twas","--gwas","examples/xwas/gwas.tsv","--weights",weights.toString(),"--reference","examples/xwas/reference.tsv","--ld","examples/xwas/ld.tsv","--out",out.toString()));
        assertEquals(1,DelimitedData.read(out).rows().size());
        assertTrue(XwasFiles.number(DelimitedData.read(out),DelimitedData.read(out).rows().get(0),"predicted_variance")>0);
    }
}
