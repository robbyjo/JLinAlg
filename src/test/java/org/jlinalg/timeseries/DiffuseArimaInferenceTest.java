/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Properties;
import org.junit.jupiter.api.Test;

class DiffuseArimaInferenceTest {
    @Test void covarianceMatchesStatsArimaAcrossIntegratedSeasonalMissingAndDrift() throws Exception {
        for(String name:new String[]{"integrated-ar","integrated-ma","integrated-missing","seasonal",
                "seasonal-only","twice-integrated","drift","seasonal-drift","monthly-missing"}) {
            var result=fixture(name);assertTrue(result.coefficientInferenceAvailable(),name);
            Properties ref=new Properties();try(var in=Files.newInputStream(Path.of("src/test/resources/feasible-extensions/"+name+".properties"))){ref.load(in);}
            double[] expected=Arrays.stream(ref.getProperty("covariance").split(",")).mapToDouble(Double::parseDouble).toArray();
            double[] actual=result.coefficientCovariance();int k=result.standardErrors().length;
            for(int i=0;i<k;i++)for(int j=0;j<k;j++) {
                double scale=Math.sqrt(expected[i*k+i]*expected[j*k+j]);
                assertEquals(expected[i*k+j],actual[i*k+j],2e-3*scale,name+" covariance "+i+","+j);
            }
        }
    }
    @Test void irregularRandomWalkDriftHasExactGlsVariance() {
        double[] y={0,1.2,Double.NaN,4.5,3.8,Double.NaN,Double.NaN,7.9,8.1};
        var fit=DiffuseArima.fit(y,new ArimaOrder(0,1,0),ArimaOptions.builder().includeDrift(true).build());
        assertTrue(fit.coefficientInferenceAvailable());
        assertEquals((y[8]-y[0])/8,fit.fit().location(),1e-12);
        assertEquals(fit.fit().innovationVariance()/8,fit.coefficientCovariance()[0],1e-8);
        double[] scaled=y.clone();for(int i=0;i<scaled.length;i++)scaled[i]*=1e5;
        var large=DiffuseArima.fit(scaled,new ArimaOrder(0,1,0),ArimaOptions.builder().includeDrift(true).build());
        assertTrue(large.coefficientInferenceAvailable());
        assertEquals(fit.coefficientCovariance()[0],large.coefficientCovariance()[0]/1e10,2e-7);
        double[] copy=fit.coefficientCovariance();copy[0]=-1;assertTrue(fit.coefficientCovariance()[0]>0);
        var compatibility=new DiffuseArima.Result(fit.fit(),1,fit.likelihoodObservations(),true);
        assertFalse(compatibility.coefficientInferenceAvailable());assertTrue(Double.isNaN(compatibility.standardErrors()[0]));
    }
    @Test void noCoefficientsHasEmptyCovarianceAndNonconvergenceSuppressesInference() throws Exception {
        double[] y={0,1,1.5,3,2,4};
        var empty=DiffuseArima.fit(y,new ArimaOrder(0,1,0),ArimaOptions.defaults());
        assertEquals(0,empty.coefficientCovariance().length);assertTrue(empty.coefficientInferenceAvailable());
        double[] series=series("integrated-missing");
        var stopped=DiffuseArima.fit(series,new ArimaOrder(1,1,1),ArimaOptions.builder().maximumFunctionEvaluations(20).optimizationTolerance(1e-14).build());
        assertFalse(stopped.fit().converged());assertFalse(stopped.coefficientInferenceAvailable());
        assertTrue(Arrays.stream(stopped.coefficientCovariance()).allMatch(Double::isNaN));
    }
    private static DiffuseArima.Result fixture(String name) throws Exception {
        Map<String,Integer> spec=new HashMap<>();
        for(String line:Files.readAllLines(Path.of("src/test/resources/timeseries/"+name+"-reference.csv")).subList(1,9)) {
            String[] cell=line.split(",");spec.put(cell[0],Integer.parseInt(cell[1]));
        }
        return DiffuseArima.fit(series(name),new ArimaOrder(spec.get("p"),spec.get("d"),spec.get("q")),ArimaOptions.builder()
            .seasonalOrder(new SeasonalArimaOrder(spec.get("P"),spec.get("D"),spec.get("Q"),spec.get("period")))
            .includeDrift(spec.get("drift")==1).build());
    }
    private static double[] series(String name) throws Exception {
        return Files.readAllLines(Path.of("src/test/resources/timeseries/"+name+"-series.csv")).stream().skip(1).mapToDouble(Double::parseDouble).toArray();
    }
}
