/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.loess;
import static org.junit.jupiter.api.Assertions.*;
import java.util.Arrays;
import org.junit.jupiter.api.Test;

class LoessScaleAuditTest {
    @Test void commonWeightScalePreservesFitsAndLeverage() {
        double[]x={0,1,2,3,4,5,6,7,8,9},y={1,2,1,3,4,6,5,8,7,9};
        var expected=Loess.fit(x,y);
        for(double unit:new double[]{1e-200,1e-20,1e200,1e308}){
            double[]w=new double[x.length];Arrays.fill(w,unit);
            var actual=Loess.fit(x,y,w,LoessOptions.defaults());
            assertArrayEquals(expected.fittedValues(),actual.fittedValues(),1e-12);
            assertArrayEquals(expected.leverage(),actual.leverage(),1e-12);
            assertEquals(expected.traceHat(),actual.traceHat(),1e-12);
        }
    }
    @Test void largeResponseValuesDoNotOverflowTheLocalRightHandSide() {
        double[]x=new double[100],y=new double[100];
        for(int i=0;i<x.length;i++){x[i]=i;y[i]=.5+.001*i;}
        var expected=Loess.fit(x,y);
        double[]scaled=y.clone();for(int i=0;i<scaled.length;i++)scaled[i]*=1e308;
        var actual=Loess.fit(x,scaled);
        double[]values=actual.fittedValues();
        for(int i=0;i<values.length;i++)assertEquals(expected.fittedValues()[i],values[i]/1e308,1e-12);
        assertEquals(expected.predict(new double[]{40.5})[0],actual.predict(new double[]{40.5})[0]/1e308,1e-12);
    }
    @Test void extremePredictorRangePreservesNormalizedGeometry() {
        double[]x=new double[20],large=new double[20],y=new double[20];
        for(int i=0;i<x.length;i++){x[i]=-1+2*i/19.0;large[i]=x[i]*1e308;y[i]=x[i]*x[i]+.1*Math.sin(i);}
        var options=LoessOptions.defaults().withSpan(1);
        assertArrayEquals(Loess.fit(x,y,options).fittedValues(),Loess.fit(large,y,options).fittedValues(),1e-12);
    }
}
