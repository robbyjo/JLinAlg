/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.xwas;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;

class PredictedOmicsTest {
    @Test void preparedModelsMatchAnalyticCorrelatedJointTestAndOwnTheirInputs() {
        double[] z={2,-1},sd={2,3},ld={1,.4,.4,1};double[][] w={{1,0},{0,1}};
        var fit=PredictedOmics.prepare(z,w,sd,ld);
        var a=fit.associations();
        assertEquals(2,a.get(0).z(),1e-14);assertEquals(-1,a.get(1).z(),1e-14);
        assertEquals(4,a.get(0).predictedVariance());assertEquals(9,a.get(1).predictedVariance());
        double expected=(4+1+1.6)/(1-.16);
        assertEquals(expected,fit.joint().chiSquare(),1e-13);
        assertEquals(Math.exp(-expected/2),fit.joint().pValue(),1e-14);
        assertEquals(PredictedOmics.test(z,w[0],sd,ld),a.get(0));
        z[0]=100;sd[0]=100;ld[0]=100;w[0][0]=100;
        assertEquals(expected,fit.joint().chiSquare(),1e-13);
        assertThrows(UnsupportedOperationException.class,()->a.clear());
    }
    @Test void dependentModelsRetainMarginalsButRejectJointInference() {
        var fit=PredictedOmics.prepare(new double[]{1,2},new double[][]{{1,0},{2,0}},new double[]{1,1},new double[]{1,.3,.3,1});
        assertEquals(2,fit.associations().size());assertEquals(fit.associations().get(0).z(),fit.associations().get(1).z());
        assertThrows(IllegalArgumentException.class,fit::joint);
        assertThrows(IllegalArgumentException.class,()->PredictedOmics.prepare(new double[]{1,2},new double[][]{{1,0}},new double[]{1,1},new double[]{1,2,2,1}));
        assertThrows(IllegalArgumentException.class,()->PredictedOmics.prepare(new double[]{1,2},new double[][]{{1,-1}},new double[]{1,1},new double[]{1,1,1,1}));
    }
}
