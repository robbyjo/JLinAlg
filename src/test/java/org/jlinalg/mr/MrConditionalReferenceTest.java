/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import static org.junit.jupiter.api.Assertions.*;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;
import javax.xml.parsers.DocumentBuilderFactory;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.genetics.ConditionalAssociation;
import org.jlinalg.genetics.ConditionalAssociationModel;
import org.jlinalg.genetics.SecondarySignalClumper;
import org.junit.jupiter.api.Test;
import org.w3c.dom.Element;

class MrConditionalReferenceTest {
    private static final double[] BX = {.22,-.31,.45,.19,-.38};
    private static final double[] SX = {.03,.04,.035,.025,.045};
    private static final double[] SY = {.04,.05,.045,.035,.055};
    private static List<HarmonizedInstrument> instruments(boolean exact) {
        double[] noise = {.01,-.02,.015,.025,-.01};
        List<HarmonizedInstrument> result = new ArrayList<>();
        for (int i=0;i<BX.length;i++) result.add(new HarmonizedInstrument("rs"+(i+1),"A","C",BX[i],SX[i],
            Math.signum(BX[i])*(.08+.7*Math.abs(BX[i])+(exact?0:noise[i])),SY[i],.3,.3,false,false));
        return result;
    }
    private static double[][] ld() {
        double[][] r = new double[5][5];
        for(int i=0;i<5;i++) for(int j=0;j<5;j++) r[i][j]=Math.pow(.35,Math.abs(i-j));
        return r;
    }
    private static Properties reference() throws Exception {
        Properties properties = new Properties();
        try(InputStream in=MrConditionalReferenceTest.class.getResourceAsStream("/r-reference/mr-conditional.properties")) {
            assertNotNull(in); properties.load(in);
        }
        return properties;
    }
    private static void compare(Properties ref,String key,double b,double se,double p,double low,double high) {
        double[] actual={b,se,p,low,high}; String[] names={"beta","se","p","lower","upper"};
        for(int i=0;i<actual.length;i++) assertEquals(Double.parseDouble(ref.getProperty(key+"."+names[i])),actual[i],2e-11,key+"."+names[i]);
    }
    private static void compare(Properties ref,String key,MrEstimate fit) {
        compare(ref,key,fit.estimate(),fit.standardError(),fit.pValue(),fit.confidenceLower(),fit.confidenceUpper());
    }

    @Test void conditionalAndJointInferenceMatchRAndKnownVarianceRegression() throws Exception {
        Properties ref=reference();
        var values=instruments(false).stream().map(x->new ConditionalAssociation(x.variantId(),x.exposureEffect(),x.exposureStandardError(),"locus")).toList();
        var model=new ConditionalAssociationModel(values,ld());
        for(boolean joint : new boolean[]{false,true}) {
            var fits=model.condition(joint?new int[]{0,1,2,3,4}:new int[]{0,2},.9);
            for(int i=0;i<fits.size();i++) {
                var fit=fits.get(i);
                compare(ref,(joint?"joint.":"conditional.")+(i+1),fit.effect(),fit.standardError(),fit.pValue(),fit.confidenceLower(),fit.confidenceUpper());
            }
        }
    }

    @Test void generalizedAndOverlapInferenceMatchRIncludingConfidence() throws Exception {
        Properties ref=reference(); var values=instruments(false);
        for(boolean random:new boolean[]{false,true}) compare(ref,"ivw."+(random?"random":"fixed"),
            CorrelatedMendelianRandomization.ivw(values,ld(),random,.9,BackendPolicy.CPU).estimate());
        var egger=CorrelatedMendelianRandomization.egger(values,ld(),.9,BackendPolicy.CPU).estimate();
        compare(ref,"egger.slope",egger.slope());
        compare(ref,"egger.intercept",egger.intercept(),egger.interceptStandardError(),egger.interceptPValue(),egger.interceptConfidenceLower(),egger.interceptConfidenceUpper());
        double[] covariance=new double[5]; for(int i=0;i<5;i++) covariance[i]=.2*SX[i]*SY[i];
        var overlap=OverlapAwareMendelianRandomization.ivw(values,covariance,.9);
        assertTrue(overlap.converged()); compare(ref,"overlap",overlap.estimate());
        var wide=OverlapAwareMendelianRandomization.ivw(values,covariance,.99).estimate();
        assertEquals(overlap.estimate().estimate(),wide.estimate());
        assertTrue(wide.confidenceUpper()>overlap.estimate().confidenceUpper());
        assertThrows(IllegalArgumentException.class,()->OverlapAwareMendelianRandomization.ivw(values,covariance,1));
        covariance[0]=1; assertThrows(IllegalArgumentException.class,()->OverlapAwareMendelianRandomization.ivw(values,covariance,.95));
    }

    @Test void changingLeadSelectionMatchesIndependentRRefits() throws Exception {
        Properties ref=reference();
        var values=List.of(new ConditionalAssociation("proxy",1.3,.1,"locus"),
            new ConditionalAssociation("causal1",1,.1,"locus"),new ConditionalAssociation("causal2",1,.1,"locus"));
        var fit=SecondarySignalClumper.select(values,new double[][]{{1,.65,.65},{.65,1,0},{.65,0,1}},.05,3,.95);
        assertEquals(ref.getProperty("selection.history"),String.join(",",fit.changes()));
        assertEquals(ref.getProperty("selection.ids"),String.join(",",fit.selectedVariantIds()));
        for(int i=0;i<3;i++) {
            var row=fit.associations().get(i);
            compare(ref,"selection."+(i+1),row.effect(),row.standardError(),row.pValue(),row.confidenceLower(),row.confidenceUpper());
        }
    }

    @Test void identityLdPreservesIndependentEstimatorsAndConfidenceChangesOnlyIntervals() {
        var values=instruments(false); double[][] identity=new double[5][5];
        for(int i=0;i<5;i++) identity[i][i]=1;
        var independent=MendelianRandomization.egger(values,.9);
        var correlated=CorrelatedMendelianRandomization.egger(values,identity,.9,BackendPolicy.CPU).estimate();
        assertEquals(independent.intercept(),correlated.intercept(),1e-12);
        assertEquals(independent.slope().estimate(),correlated.slope().estimate(),1e-12);
        for(boolean random:new boolean[]{false,true}) {
            var plain=MendelianRandomization.ivw(values,random,.9);
            var gls=CorrelatedMendelianRandomization.ivw(values,identity,random,.9,BackendPolicy.CPU).estimate();
            assertEquals(plain.estimate(),gls.estimate(),1e-12);
            assertEquals(plain.standardError(),gls.standardError(),1e-12);
            var wide=CorrelatedMendelianRandomization.ivw(values,identity,random,.99,BackendPolicy.CPU).estimate();
            assertEquals(gls.pValue(),wide.pValue()); assertTrue(wide.confidenceUpper()>gls.confidenceUpper());
        }
    }

    @Test void eggerSvgLinePassesThroughOrientedPointsWithNonzeroIntercept() throws Exception {
        var values=instruments(true);
        for(boolean generalized:new boolean[]{false,true}) {
            var fit=generalized?CorrelatedMendelianRandomization.egger(values,ld(),.95,BackendPolicy.CPU).estimate():MendelianRandomization.egger(values,.95);
            assertEquals(.08,fit.intercept(),1e-12); assertEquals(.7,fit.slope().estimate(),1e-12);
            Path path=Files.createTempFile("mr-egger-", ".svg");
            try {
                MrPlot.writeSvg(path,values,fit);
                var xml=DocumentBuilderFactory.newInstance().newDocumentBuilder().parse(path.toFile());
                Element line=(Element)xml.getElementsByTagName("line").item(0);
                double y0=Double.parseDouble(line.getAttribute("y1")),y1=Double.parseDouble(line.getAttribute("y2"));
                assertTrue(y0>=0 && y0<=430 && y1>=0 && y1<=430);
                var circles=xml.getElementsByTagName("circle");
                double firstX=Double.parseDouble(((Element)circles.item(0)).getAttribute("cx"));
                for(int i=0;i<circles.getLength();i++) {
                    Element circle=(Element)circles.item(i);
                    double x=Double.parseDouble(circle.getAttribute("cx")),y=Double.parseDouble(circle.getAttribute("cy"));
                    assertEquals(y0+(y1-y0)*x/610,y,1e-9);
                    if(i==1) assertTrue(x>firstX,"negative exposure must be reoriented before plotting");
                }
                assertThrows(IllegalArgumentException.class,()->MrPlot.writeSvg(path,values,fit.slope()));
            } finally {Files.deleteIfExists(path);}
        }
    }
}
