package org.jlinalg.association;
import static org.junit.jupiter.api.Assertions.*;
import java.util.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.ols.OlsOptions;
import org.junit.jupiter.api.Test;
class AssociationScaleAuditTest {
    @Test void extremeFiniteUnitsPreserveOlsAndP3dInference() {
        int n=40;double[] y=new double[n];double[][] fixed=new double[n][1],g=new double[n][3];
        double[] scales={1,1e-160,1e160};
        for(int i=0;i<n;i++) {fixed[i][0]=1;y[i]=1+.7*Math.sin(.7*i)+.2*Math.cos(1.3*i);
            for(int j=0;j<3;j++)g[i][j]=scales[j]*Math.sin(.7*i);}
        var options=AssociationEngineOptions.defaults().withBackendPolicy(BackendPolicy.CPU)
            .withFailurePolicy(AssociationFailurePolicy.RECORD_NAN);
        var ols=FastOlsAssociation.scanPredictors(y,fixed,g,null,null,null,OlsOptions.defaults(),options);
        assertTrue(ols.failures().isEmpty());
        var scanner=org.jlinalg.gwas.RemlAssociationScanner.prepare(y,fixed,
            List.of(org.jlinalg.reml.VarianceComponent.identity("e",n)),org.jlinalg.reml.RemlOptions.defaults(),BackendPolicy.CPU);
        var mixed=scanner.scan(g,null);
        for(int j=0;j<3;j++) {
            assertEquals(.696556541560302,ols.estimate(j).beta()*scales[j],2e-14);
            assertEquals(.03282048034188775,ols.estimate(j).standardError()*scales[j],2e-14);
            assertEquals(1,ols.estimate(j).pValue()/1.1339376191093256e-22,2e-12);
            assertEquals(mixed.beta()[0],mixed.beta()[j]*scales[j],2e-14);
            assertEquals(mixed.standardErrors()[0],mixed.standardErrors()[j]*scales[j],2e-14);
            assertEquals(1,mixed.pValues()[j]/mixed.pValues()[0],2e-12);
        }
        assertEquals(scales[2]*Math.sin(.7),g[1][2]); // no caller mutation
    }
    @Test void unrepresentableExportIsRecordedAsFailure() {
        double[] y={1,2,2,4,5,7,7,9};double[][] fixed=new double[8][1],g=new double[8][1];
        for(int i=0;i<8;i++){fixed[i][0]=1;g[i][0]=i*1e-320;}
        var options=AssociationEngineOptions.defaults().withBackendPolicy(BackendPolicy.CPU)
            .withFailurePolicy(AssociationFailurePolicy.RECORD_NAN);
        var fit=FastOlsAssociation.scanPredictors(y,fixed,g,null,null,null,OlsOptions.defaults(),options);
        assertEquals(1,fit.failures().size());assertTrue(Double.isNaN(fit.estimate(0).beta()));
        assertThrows(IllegalArgumentException.class,()->FastOlsAssociation.scanPredictors(y,fixed,g,null,null,null,
            OlsOptions.defaults(),options.withFailurePolicy(AssociationFailurePolicy.FAIL_FAST)));
    }
    @Test void fullScanMatchesFrozenRIncludingImputationAndInference() throws Exception {
        int n=400,m=128;double[] y=new double[n];double[][] x=new double[n][2],g=new double[n][m];
        for(int i=0;i<n;i++) {
            y[i]=Math.sin(.13*i)+.3*Math.cos(.07*i);x[i][0]=1;x[i][1]=Math.cos(.07*i);
            for(int j=0;j<m;j++)g[i][j]=(i+j)%97==0?Double.NaN:Math.sin(.023*(i+1)*(j+1))+.2*Math.cos(.031*(i+1)*(j+1));
        }
        var result=FastOlsAssociation.scanPredictors(y,x,g,null,null,null,OlsOptions.defaults(),
            AssociationEngineOptions.defaults().withBackendPolicy(BackendPolicy.CPU));
        double maxBeta=0,maxSe=0,maxP=0;
        try(var reader=new java.io.BufferedReader(new java.io.InputStreamReader(
                getClass().getResourceAsStream("/pipeline-audit/ols-scan-reference.tsv"),java.nio.charset.StandardCharsets.UTF_8))) {
            reader.readLine();
            for(int j=0;j<m;j++) {
                String[] fields=reader.readLine().split("\t");var e=result.estimate(j);
                maxBeta=Math.max(maxBeta,Math.abs(e.beta()-Double.parseDouble(fields[1])));
                maxSe=Math.max(maxSe,Math.abs(e.standardError()-Double.parseDouble(fields[2])));
                maxP=Math.max(maxP,Math.abs(e.pValue()-Double.parseDouble(fields[3])));
            }
            assertNull(reader.readLine());
        }
        assertTrue(result.failures().isEmpty());
        assertTrue(maxBeta<1e-12);assertTrue(maxSe<1e-12);assertTrue(maxP<1e-12);
        System.out.println("R OLS scan max errors beta="+maxBeta+" SE="+maxSe+" p="+maxP);
    }
    @Test void repeatedRefitAdapterDoesNotDiscardRemlNonconvergence() {
        var fitter=AssociationModels.reml(List.of(org.jlinalg.reml.VarianceComponent.identity("e",8)),
            org.jlinalg.reml.RemlOptions.builder().initialVariances(.01).maximumIterations(1).build());
        assertThrows(IllegalArgumentException.class,()->fitter.fit(new double[]{1,2,2,4,5,7,7,9},
            new double[]{1,1,1,1,1,1,1,1},8,1,BackendPolicy.CPU));
    }
    @Test void glmScoreIsInvariantToPredictorUnits() {
        double[] y={0,0,1,0,1,0,1,1};double[][] x=new double[8][1],g=new double[8][2];
        for(int i=0;i<8;i++){x[i][0]=1;g[i][0]=i;g[i][1]=i*1e-8;}
        var options=AssociationEngineOptions.defaults().withBackendPolicy(BackendPolicy.CPU)
            .withFailurePolicy(AssociationFailurePolicy.RECORD_NAN);
        var fit=FastGlmAssociation.prepare(y,x,org.jlinalg.glm.GlmFamilies.binomial(),null,null,
            org.jlinalg.glm.GlmOptions.defaults(),options).scan(g,null,options);
        assertEquals(fit.estimate(0).pValue(),fit.estimate(1).pValue(),1e-14);
    }
    @Test void predictorUnitsDoNotChangeOLSInferenceAndCollinearityStillFails() {
        double[] y={1,2,2,4,5,7,7,9};double[][] x=new double[8][1],g=new double[8][3];
        for(int i=0;i<8;i++){x[i][0]=1;g[i][0]=i/2;g[i][1]=(i/2)*1e-8;g[i][2]=1e-8;}
        var options=AssociationEngineOptions.defaults().withBackendPolicy(BackendPolicy.CPU)
            .withFailurePolicy(AssociationFailurePolicy.RECORD_NAN);
        var result=FastOlsAssociation.scanPredictors(y,x,g,null,null,null,OlsOptions.defaults(),options);
        assertEquals(2.25,result.estimate(0).beta(),1e-12);
        assertEquals(result.estimate(0).beta(),result.estimate(1).beta()*1e-8,1e-12);
        assertEquals(result.estimate(0).standardError(),result.estimate(1).standardError()*1e-8,1e-12);
        assertEquals(result.estimate(0).pValue(),result.estimate(1).pValue(),1e-14);
        assertTrue(Double.isNaN(result.estimate(2).beta()));assertEquals(1,result.failures().size());
    }
    @Test void nearlyPerfectSignalRetainsPositiveResidualUncertainty() {
        double[] y=new double[20];double[][] x=new double[20][1],g=new double[20][1];
        for(int i=0;i<20;i++){x[i][0]=1;g[i][0]=i;y[i]=i+(i%2==0?1:-1)*1e-7;}
        var r=FastOlsAssociation.scanPredictors(y,x,g,null,null,null,OlsOptions.defaults(),
            AssociationEngineOptions.defaults().withBackendPolicy(BackendPolicy.CPU));
        assertTrue(r.estimate(0).standardError()>1e-10);
        assertEquals(4.0721997167553101e-9,r.estimate(0).standardError(),2e-16);
    }
}
