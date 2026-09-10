package org.jlinalg.settest;
import static org.junit.jupiter.api.Assertions.*;
import java.util.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.ols.OlsOptions;
import org.jlinalg.pipeline.*;
import org.junit.jupiter.api.Test;
class SetScaleAuditTest {
    @Test void higherRankTailsMatchPositiveRSeriesAndRespectStochasticBounds() throws Exception {
        try(var reader=new java.io.BufferedReader(new java.io.InputStreamReader(
                getClass().getResourceAsStream("/pipeline-audit/rank-three-reference.tsv"),java.nio.charset.StandardCharsets.UTF_8))) {
            reader.readLine();String line;double previous=1;
            while((line=reader.readLine())!=null) {
                String[] f=line.split("\t");double q=Double.parseDouble(f[0]), expected=Double.parseDouble(f[1]);
                var tail=QuadraticFormDistribution.survival(q,new double[]{1,.5,.25});
                assertEquals("positive-gamma-series",tail.method());
                assertEquals(1,tail.pValue()/expected,2e-11);
                assertTrue(tail.pValue()<=Double.parseDouble(f[2]));
                assertTrue(tail.pValue()<previous);previous=tail.pValue();
                assertEquals(1,QuadraticFormDistribution.survival(q*1e100,new double[]{1e100,5e99,2.5e99}).pValue()/expected,2e-11);
            }
        }
    }
    @Test void unresolvableHigherRankSpectrumUsesLabeledSaddlepointFallback() {
        var tail=QuadraticFormDistribution.survival(100,
            new double[]{1,.5,1e-20});
        assertEquals("lugannani-rice-saddlepoint",tail.method());
        assertTrue(tail.pValue()>0&&tail.pValue()<1e-15);
        double rankTwo=QuadraticFormDistribution.survival(100,
            new double[]{1,.5}).pValue();
        assertEquals(1,tail.pValue()/rankTwo,.12);
    }
    @Test void paddedExactRankQuantilesRemainAccurateNearOne() {
        double p=Math.nextDown(1.0);
        assertEquals(jdistlib.ChiSquare.quantile(p,1,false,false),
            QuadraticFormDistribution.critical(new double[]{1,0,0,0,0,0,0,0,0,0},p));
        assertEquals(3*jdistlib.ChiSquare.quantile(p,2,false,false),
            QuadraticFormDistribution.critical(new double[]{0,3,3,0},p));
        // At small q: P(Q<=q) ~ q/(2*sqrt(lambda1*lambda2)).
        double actual=QuadraticFormDistribution.critical(new double[]{1,.5,0},p);
        assertEquals(1,actual/(2*Math.sqrt(.5)*(1-p)),2e-12);
        assertThrows(IllegalArgumentException.class,()->QuadraticFormDistribution.critical(new double[]{-1},1));
    }
    @Test void nearlyRankOneExtremeTailsHaveRelativeAccuracy() {
        for(double q:new double[]{50,200,700}) {
            double expected=jdistlib.ChiSquare.cumulative(q,1,false,false);
            for(double unit:new double[]{1,1e-100,1e100}) {
                double actual=QuadraticFormDistribution.survival(q*unit,new double[]{unit,1e-12*unit}).pValue();
                assertEquals(1,actual/expected,2e-9,"relative tail error at q="+q);
            }
        }
    }
    @Test void rankTwoTailAndCriticalMatchIndependentRPolarIntegral() {
        assertEquals(.16908111409965279,
            QuadraticFormDistribution.survival(20,new double[]{1,10}).pValue(),2e-12);
        assertEquals(39.484170286858287,QuadraticFormDistribution.critical(new double[]{1,10},.05),2e-9);
        assertEquals(.05,QuadraticFormDistribution.survival(
            QuadraticFormDistribution.critical(new double[]{1,10},.05),new double[]{1,10}).pValue(),2e-12);
    }
    @Test void mixtureScaleDoesNotChangeTailAndInvalidSpectrumIsRejected() {
        for(double scale:new double[]{1e-100,1e-16,1,1e100})
            assertEquals(.18587673236587593,
                QuadraticFormDistribution.survival(3.5*scale,new double[]{2*scale}).pValue(),2e-15);
        assertThrows(IllegalArgumentException.class,()->QuadraticFormDistribution.survival(1,new double[]{1,Double.NaN}));
        assertThrows(IllegalArgumentException.class,()->QuadraticFormDistribution.survival(1,new double[]{1,-1}));
    }
    @Test void kernelWeightsInDifferentUnitsDoNotChangeSingleVariantSkat() {
        double[] y={1,2,2,4,5,7,7,9},g={0,0,1,1,1,2,2,2};double[][] x=new double[8][1];
        for(double[] row:x)row[0]=1;
        var model=LinearSetTestNullModel.prepare(y,x,OlsOptions.defaults(),BackendPolicy.CPU);
        double reference=Double.NaN, burdenReference=Double.NaN, scoreReference=Double.NaN;
        for(double weight:new double[]{1,1e-8,1e8}) {
            var set=new VariantSet("s",List.of(new WeightedVariant(
                new VariantRecord("v","",0,"","",g,Double.NaN),EffectAllele.ALTERNATE,weight)));
            var result=SetTests.skat(set,model,SetTestOptions.defaults());
            if(Double.isNaN(reference))reference=result.pValue();
            assertEquals(reference,result.pValue(),2e-14);
            var burden=SetTests.burden(set,model,SetTestOptions.defaults());
            var score=SetTests.burden(SetTests.prepare(set,8,SetTestOptions.defaults()),model);
            if(Double.isNaN(burdenReference)){burdenReference=burden.pValue();scoreReference=score.pValue();}
            assertEquals(burdenReference,burden.pValue(),2e-14);
            assertEquals(scoreReference,score.pValue(),2e-14);
        }
    }
}
