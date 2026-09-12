/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

import static org.junit.jupiter.api.Assertions.*;
import java.util.function.DoubleUnaryOperator;
import org.junit.jupiter.api.Test;

class MixedModelProfileTest {
    private static final DoubleUnaryOperator PROFILE=t->-.5*t*t;

    @Test void boundaryProfileUsesTheMixtureCriticalValue() {
        var result=MixedModelProfile.boundaryInterval(PROFILE,0,0,.95,0,4,12);
        // 50:50 mixture: chi-square_1 90th percentile = normal_0.95 squared.
        assertEquals(1.6448536269514722,result.upper(),1e-8);
        assertEquals(0,result.lower());
        assertTrue(result.upperFound());assertFalse(result.lowerFound());
    }

    @Test void boundaryProfileRejectsInvalidConfidenceBeforeCallingTheRefit() {
        for(double confidence:new double[]{-.1,0,.5,1,1.2,Double.NaN,Double.POSITIVE_INFINITY})
            assertThrows(IllegalArgumentException.class,()->MixedModelProfile.boundaryInterval(
                t->{fail("invalid input must not call refit");return 0;},0,0,confidence,0,4,12));
    }

    @Test void bothProfileApisValidateAllCommonInputs() {
        for(boolean boundary:new boolean[]{false,true}) {
            assertThrows(IllegalArgumentException.class,()->call(boundary,null,0,0,0,4,12));
            for(double invalid:new double[]{Double.NaN,Double.NEGATIVE_INFINITY,Double.POSITIVE_INFINITY}) {
                assertThrows(IllegalArgumentException.class,()->call(boundary,PROFILE,invalid,0,0,4,12));
                assertThrows(IllegalArgumentException.class,()->call(boundary,PROFILE,0,invalid,0,4,12));
                assertThrows(IllegalArgumentException.class,()->call(boundary,PROFILE,0,0,invalid,4,12));
                assertThrows(IllegalArgumentException.class,()->call(boundary,PROFILE,0,0,0,invalid,12));
            }
            assertThrows(IllegalArgumentException.class,()->call(boundary,PROFILE,-1,0,0,4,12));
            assertThrows(IllegalArgumentException.class,()->call(boundary,PROFILE,5,0,0,4,12));
            assertThrows(IllegalArgumentException.class,()->call(boundary,PROFILE,0,0,0,0,12));
            assertThrows(IllegalArgumentException.class,()->call(boundary,PROFILE,0,0,1,-1,12));
            assertThrows(IllegalArgumentException.class,()->call(boundary,PROFILE,0,0,0,4,7));
            assertThrows(IllegalStateException.class,()->call(boundary,t->Double.NaN,0,0,0,4,12));
        }
    }

    private static ProfileLikelihoodInterval call(boolean boundary,DoubleUnaryOperator f,
            double estimate,double likelihood,double lower,double upper,int grid) {
        return boundary?MixedModelProfile.boundaryInterval(f,estimate,likelihood,.95,lower,upper,grid)
            :MixedModelProfile.interval(f,estimate,likelihood,.95,lower,upper,grid);
    }
}
