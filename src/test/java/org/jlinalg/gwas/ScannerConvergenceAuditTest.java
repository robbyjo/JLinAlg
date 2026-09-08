package org.jlinalg.gwas;
import static org.junit.jupiter.api.Assertions.*;
import java.util.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.reml.*;
import org.junit.jupiter.api.Test;
class ScannerConvergenceAuditTest {
    @Test void markerUnitsDoNotChangeP3dInference() {
        double[] y={1,2,2,4,5,7,7,9};double[][] x=new double[8][1],g=new double[8][3];
        for(int i=0;i<8;i++){x[i][0]=1;g[i][0]=i/2;g[i][1]=(i/2)*1e-8;g[i][2]=1;}
        var scanner=RemlAssociationScanner.prepare(y,x,List.of(VarianceComponent.identity("e",8)),
            RemlOptions.defaults(),BackendPolicy.CPU);
        var result=scanner.scan(g,null);
        assertEquals(result.pValues()[0],result.pValues()[1],1e-14);
        assertTrue(Double.isNaN(result.pValues()[2]));
    }
    @Test void failedNullFitCannotProduceAssociationPValues() {
        double[] y={1,2,2,4,5,7,7,9};double[][] x=new double[8][1];
        for(double[] row:x)row[0]=1;
        var ex=assertThrows(IllegalArgumentException.class,()->RemlAssociationScanner.prepare(y,x,
            List.of(VarianceComponent.identity("e",8)),
            RemlOptions.builder().initialVariances(.01).maximumIterations(1).build(),BackendPolicy.CPU));
        assertTrue(ex.getMessage().contains("did not converge"));
    }
}
