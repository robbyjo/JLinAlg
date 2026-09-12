/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.*;
import java.nio.file.*;
import java.util.*;
import org.jlinalg.raremetal.RareMetalStudy;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class RareMetaMetadataTest {
    @TempDir Path dir;
    private static final String SCORE="##AnalyzedSamples=100\n#CHROM POS REF ALT N_INFORMATIVE ALL_AF U_STAT SQRT_V_STAT\n1 10 A G 100 .01 2 2\n1 20 C T 100 .02 -3 3\n";
    private static final String COV="#CHROM CURRENT_POS MARKERS_IN_WINDOW COV_MATRICES\n1 10 10,20 .04,.01\n1 20 20 .09\n";
    private static final String MODEL="##SummaryMetadataVersion=1\n##TraitType=quantitative\n##TraitId=height\n##TraitUnits=cm\n##TraitTransformation=none\n##NullModel=unrelated-Gaussian\n##GenotypeModel=diploid-additive\n##EffectScale=trait-per-alt-allele\n##ScoreCalibration=asymptotic-normal\n##CovarianceModel=model-based\n##CovarianceScale=per-sample\n";
    private RareMetalStudy study(String score,String cov,boolean strict)throws IOException {
        Path s=dir.resolve("score"),c=dir.resolve("cov");
        Files.writeString(s,score);if(cov!=null)Files.writeString(c,cov);
        return new RareMetalStudy(s,cov==null?null:c,0,0,strict);
    }
    @Test void strictScaleRoundTripHasIndependentCovarianceAndBurdenReference()throws Exception {
        for(boolean scoreScale:new boolean[]{false,true}) {
            String model=scoreScale?MODEL.replace("=per-sample","=score"):MODEL;
            String cov=scoreScale?COV.replace(".04,.01","4,1").replace(".09","9"):COV;
            try(var s=study(model+SCORE,model+cov,true)) {
                var rows=s.region("1",10,20);
                assertTrue(s.metadata().complete());
                assertArrayEquals(new double[]{4,1,1,9},s.covariance(List.of(rows.get(10L),rows.get(20L))),1e-14);
                var fit=org.jlinalg.settest.SummarySetTests.burden("gene",
                    new org.jlinalg.settest.SetTestScoreState(new double[]{2,-3},s.covariance(List.of(rows.get(10L),rows.get(20L))),2),new double[]{1,1});
                assertEquals(-1.0/15,fit.beta(),1e-14);
                assertEquals(1/Math.sqrt(15),fit.standardError(),1e-14);
            }
        }
    }
    @Test void legacyAssumptionsAreExplicitAndStrictRejectsMissingDeclarations()throws Exception {
        try(var s=study(SCORE,COV,false)){assertEquals("legacy-assumed-quantitative",s.metadata().status());}
        assertTrue(assertThrows(IOException.class,()->study(SCORE,COV,true)).getMessage().contains("strict model metadata"));
        assertThrows(IOException.class,()->study(MODEL+SCORE,COV,true));
    }
    @Test void declaredUnsupportedModelsAndCalibrationNeverFallBackToLegacy() {
        for(String declaration:List.of("TraitType=binary","NullModel=unrelated-logistic","NullModel=Cox",
                "ScoreCalibration=saddlepoint","CovarianceScale=unknown","SummaryMetadataVersion=2",
                "BinaryTrait=True","GenotypeModel=dominant","CovarianceModel=cluster-robust")) {
            IOException error=assertThrows(IOException.class,()->study("##"+declaration+"\n"+SCORE,null,false));
            assertTrue(error.getMessage().contains("unsupported"),error.getMessage());
        }
    }
    @Test void conflictingPairDuplicateAndCompressedHeadersFailBeforeReadingData() {
        for(String declaration:List.of("TraitId=weight","TraitUnits=kg","AnalyzedSamples=99","CovarianceScale=score",
                "NullModel=related-Gaussian-REML","TraitTransformation=inverse-normal"))
            assertTrue(assertThrows(IOException.class,()->study(MODEL+SCORE,"##"+declaration+"\n"+COV,false)).getMessage().contains("incompatible metadata"));
        assertThrows(IOException.class,()->study("##AnalyzedSamples=100\n"+SCORE,null,false));
        assertThrows(IOException.class,()->study("##TraitId= \n"+SCORE,null,false));
        for(String header:List.of("#CHROM CURRENT_POS EXP COV_MATRICES","#CHROM CURRENT_POS REF ALT EXP COV_MATRICES"))
            assertTrue(assertThrows(IOException.class,()->study(SCORE,header+"\n1 10 10 4\n",false)).getMessage().contains("dedicated decoder"));
    }
    @Test void missingFirstCohortCannotHideIncompatibleLaterCohorts()throws Exception {
        Files.writeString(dir.resolve("a"),SCORE);
        Files.writeString(dir.resolve("b"),MODEL+SCORE);
        Files.writeString(dir.resolve("c"),MODEL.replace("TraitId=height","TraitId=weight")+SCORE);
        Files.writeString(dir.resolve("manifest"),"cohort\tscores\na\ta\nb\tb\nc\tc\n");
        ByteArrayOutputStream errors=new ByteArrayOutputStream();
        int result=JLinAlgCli.run(new String[]{"rare-meta","--cohorts",dir.resolve("manifest").toString(),"--genome-build","test","--test","single","--out",dir.resolve("out").toString()},new PrintStream(OutputStream.nullOutputStream()),new PrintStream(errors));
        assertEquals(2,result);assertTrue(errors.toString().contains("TraitId"));
        assertFalse(Files.exists(dir.resolve("out.single.tsv")));
    }
    @Test void metadataFromCovarianceAlsoConstrainsCohortCompatibility()throws Exception {
        try(var a=study(SCORE,"##TraitId=height\n"+COV,false)) {
            var metadata=a.metadata();assertEquals("height",metadata.value("TraitId"));
            try(var b=study("##TraitId=weight\n"+SCORE,null,false)) {
                assertThrows(IOException.class,()->metadata.requireCompatible(b.metadata()));
            }
        }
    }
    @Test void gaussianRemlCanPoolWithGaussianOlsButTransformationsMustAgree()throws Exception {
        try(var a=study(MODEL+SCORE,null,false)) {
            var metadata=a.metadata();
            try(var b=study(MODEL.replace("unrelated-Gaussian","related-Gaussian-REML")+SCORE,null,false)) {
                assertDoesNotThrow(()->metadata.requireCompatible(b.metadata()));
            }
            try(var b=study("##InverseNormal=ON\n"+SCORE,null,false)) {
                assertThrows(IOException.class,()->metadata.requireCompatible(b.metadata()));
            }
        }
        assertThrows(IOException.class,()->study("##InverseNormal=ON\n"+MODEL+SCORE,null,false));
    }
    @Test void duplicatePositionsRemainAmbiguousEvenForOneRequestedAllele()throws Exception {
        try(var s=study(SCORE+"1 20 C G 100 .01 2 2\n",COV,false)) {
            assertThrows(IOException.class,()->s.region("1",20,20));
            try(var cursor=s.cursor()){cursor.next();cursor.next();assertThrows(IOException.class,cursor::next);}
        }
    }
    @Test void rareCaseRFixtureIsNotMisrepresentedAsCalibratedGaussianInference()throws Exception {
        String[] f=Files.readAllLines(Path.of("src/test/resources/raremetal/binary-calibration-boundary.tsv")).get(1).split("\t");
        double u=Double.parseDouble(f[4]),v=Double.parseDouble(f[5]),normal=Double.parseDouble(f[6]),exact=Double.parseDouble(f[7]);
        assertEquals(.01,exact,1e-15);
        assertEquals(normal,org.jlinalg.settest.SummarySetTests.singleVariant("rare-case",u,v).pValue(),normal*1e-12);
        String row="##AnalyzedSamples=1000\n##TraitType=binary\n#CHROM POS REF ALT N_INFORMATIVE ALL_AF U_STAT SQRT_V_STAT\n"
            +"1 10 A G 1000 .005 "+u+" "+Math.sqrt(v)+"\n";
        assertThrows(IOException.class,()->study(row,null,false));
    }
    @Test void cliExportDeclaresTraitIdentityAndStrictImportPreservesNumbers()throws Exception {
        Path examples=Path.of("examples/rare-meta").toAbsolutePath();
        ByteArrayOutputStream errors=new ByteArrayOutputStream();
        PrintStream sink=new PrintStream(OutputStream.nullOutputStream());
        String[] args={"rare-score","--vcf",examples.resolve("cohort-a.vcf").toString(),"--pheno",examples.resolve("cohort-a.tsv").toString(),
            "--id","sample","--response","trait","--trait-id","height","--trait-units","cm","--genome-build","GRCh38","--out",dir.resolve("export").toString()};
        assertEquals(0,JLinAlgCli.run(args,sink,new PrintStream(errors)),errors.toString());
        try(var study=new RareMetalStudy(dir.resolve("export.score.txt.gz"),dir.resolve("export.cov.txt.gz"),0,0,true)) {
            assertEquals("height",study.metadata().value("TraitId"));assertEquals("cm",study.metadata().value("TraitUnits"));
        }
        Files.writeString(dir.resolve("manifest"),"cohort\tscores\tcovariance\na\texport.score.txt.gz\texport.cov.txt.gz\n");
        for(String policy:List.of("legacy","strict")) {
            String[] meta={"rare-meta","--cohorts",dir.resolve("manifest").toString(),"--genome-build","GRCh38","--test","single,burden","--window-size","100","--maf","0.5","--model-metadata",policy,"--out",dir.resolve(policy).toString()};
            assertEquals(0,JLinAlgCli.run(meta,sink,new PrintStream(errors)),errors.toString());
            assertTrue(Files.readString(dir.resolve(policy+".log")).contains("model_metadata=declared-quantitative"));
        }
        for(String test:List.of("single","burden","qc"))assertEquals(Files.readString(dir.resolve("legacy."+test+".tsv")),Files.readString(dir.resolve("strict."+test+".tsv")));
    }
}
