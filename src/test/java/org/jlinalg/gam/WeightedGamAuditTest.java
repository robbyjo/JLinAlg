/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.gam;
import static org.junit.jupiter.api.Assertions.*;
import java.util.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

class WeightedGamAuditTest {
    @Test void offsetsAreRestoredWithoutChangingResidualsOrCovarianceEstimates() {
        int n=30;double[] y=new double[n],adjusted=new double[n],w=new double[n],o=new double[n],x=new double[n];double[][] fixed=new double[n][1];
        for(int r=0;r<n;r++){fixed[r][0]=1;x[r]=r/(double)(n-1);w[r]=.5+(r%4)*.2;o[r]=10+r*.2;
            adjusted[r]=2+Math.sin(6*x[r])+.1*Math.cos(3*r);y[r]=adjusted[r]+o[r];adjusted[r]=y[r]-o[r];}
        var smooth=List.of(PSplineTerm.of("s",x,6));
        var fit=WeightedGam.fitGaussian(y,fixed,smooth,w,o,RemlOptions.defaults(),BackendPolicy.CPU);
        var control=WeightedGam.fitGaussian(adjusted,fixed,smooth,w,null,RemlOptions.defaults(),BackendPolicy.CPU);
        assertArrayEquals(control.reml().varianceComponents(),fit.reml().varianceComponents(),1e-6);
        assertArrayEquals(control.residuals(),fit.residuals(),1e-6);
        for(int r=0;r<n;r++) {
            assertEquals(y[r],fit.fittedValues()[r]+fit.residuals()[r],1e-12);
            assertEquals(control.fittedValues()[r]+o[r],fit.fittedValues()[r],1e-6);
        }
    }
}
