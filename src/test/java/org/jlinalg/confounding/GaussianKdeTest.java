/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.confounding;
import static org.junit.jupiter.api.Assertions.*;
import java.util.Random;
import org.junit.jupiter.api.Test;

class GaussianKdeTest {
    @Test void fftDensityMatchesDirectKernelIncludingSparseTails() {
        Random random = new Random(42); double[] x = new double[2500];
        for(int i=0;i<x.length;i++)x[i]=random.nextGaussian()+(i%3==0?3:0);
        x[0]=-8; x[1]=10;
        for(double h:new double[]{.08,.3,1.2}) {
            double[] density=GaussianKde.atObservations(x,h);
            for(int i=0;i<x.length;i+=11) {
                double direct=0;
                for(double v:x)direct+=Math.exp(-.5*Math.pow((x[i]-v)/h,2));
                direct/=x.length*h*Math.sqrt(2*Math.PI);
                assertEquals(1,density[i]/direct,.0005,"bandwidth="+h+" index="+i);
            }
        }
    }
    @Test void identicalFeaturesDoNotRequirePairwiseWork() {
        double[] x=new double[850000];
        double[] density=GaussianKde.atObservations(x,.5);
        assertEquals(1/(.5*Math.sqrt(2*Math.PI)),density[0],1e-15);
        assertEquals(density[0],density[x.length-1]);
    }
}
