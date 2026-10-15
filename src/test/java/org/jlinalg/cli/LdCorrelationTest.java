package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.IOException;
import org.junit.jupiter.api.Test;

class LdCorrelationTest {
    @Test void normalizesOnlyRoundingBeyondPerfectCorrelation() throws Exception {
        assertEquals(1.0,LdCorrelation.bounded(Math.nextUp(1.0)));
        assertEquals(-1.0,LdCorrelation.bounded(Math.nextDown(-1.0)));
        for(double value:new double[]{-1,-.9,0,.001,.9,1})
            assertEquals(value,LdCorrelation.bounded(value));
    }
    @Test void rejectsMaterialBoundViolationsAndNonfiniteValues() {
        for(double value:new double[]{1+1e-11,-1-1e-11,Double.NaN,
                Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY})
            assertThrows(IOException.class,()->LdCorrelation.bounded(value));
    }
}
