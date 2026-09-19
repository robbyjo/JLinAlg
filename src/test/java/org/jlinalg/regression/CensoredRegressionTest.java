/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

import static org.junit.jupiter.api.Assertions.*;
import java.io.*;
import java.util.*;
import org.junit.jupiter.api.Test;

class CensoredRegressionTest {
    @Test void allDistributionsMatchSurvregIncludingObservedCovariance() throws IOException {
        double[][] data;
        try(var reader=new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("/r-reference/regression-inference/censored.tsv")))) {
            data=reader.lines().skip(1).map(s->Arrays.stream(s.split("\t")).mapToDouble(Double::parseDouble).toArray()).toArray(double[][]::new);
        }
        Properties ref=new Properties();try(var in=getClass().getResourceAsStream("/r-reference/regression-inference/reference.properties")){ref.load(in);}
        for(var distribution:CensoredRegression.Distribution.values()) {
            boolean gaussian=distribution==CensoredRegression.Distribution.GAUSSIAN;
            double[] y=Arrays.stream(data).mapToDouble(r->r[gaussian?0:4]).toArray();
            int[] censor=Arrays.stream(data).mapToInt(r->(int)r[gaussian?3:5]).toArray();
            double[][] x=Arrays.stream(data).map(r->new double[]{1,r[1],r[2]}).toArray(double[][]::new);
            var fit=CensoredRegression.fit(y,x,censor,distribution);String prefix="censored."+distribution.name().toLowerCase(Locale.ROOT);
            assertArrayEquals(values(ref,prefix+".beta"),fit.coefficients(),2e-6,prefix);
            assertArrayEquals(values(ref,prefix+".cov"),fit.covariance(),2e-7,prefix);
            assertEquals(values(ref,prefix+".scale")[0],fit.scale(),2e-6,prefix);
            assertEquals(values(ref,prefix+".logLik")[0],fit.logLikelihood(),2e-7,prefix);
            if(gaussian) {
                assertTrue(fit.observedMean(new double[]{1,0,0},-.4,2.5)>-.4);
                assertTrue(fit.observedMean(new double[]{1,0,0},-.4,2.5)<2.5);
                assertEquals(fit.location(new double[]{1,0,0}),fit.observedMean(new double[]{1,0,0},Double.NEGATIVE_INFINITY,Double.POSITIVE_INFINITY),1e-12);
            } else {
                double[] row={1,.2,-.1};double median=fit.timeQuantile(row,.5);
                assertEquals(.5,fit.survival(row,median),1e-12);
                assertArrayEquals(new double[]{1,1,1},fit.timeRatio(row,row,.95),1e-12);
                assertEquals(Math.exp(fit.coefficients()[1]),fit.timeRatio(new double[]{1,0,0},new double[]{1,1,0},.95)[0],1e-12);
                double[] ci=fit.timeQuantileInterval(row,.25,.95);assertTrue(ci[1]<ci[0]&&ci[0]<ci[2]);
                double[] rescaled=y.clone();for(int i=0;i<rescaled.length;i++)rescaled[i]*=10;
                var again=CensoredRegression.fit(rescaled,x,censor,distribution);
                assertEquals(fit.coefficients()[0]+Math.log(10),again.coefficients()[0],2e-6);
                assertEquals(10*median,again.timeQuantile(row,.5),2e-5);
            }
        }
    }
    @Test void failsForUnidentifiedOrInvalidData() {
        double[] y={1,2,3,4};double[][] x={{1},{1},{1},{1}};
        assertThrows(IllegalArgumentException.class,()->CensoredRegression.fit(y,x,new int[]{1,1,1,1},CensoredRegression.Distribution.GAUSSIAN));
        assertThrows(IllegalArgumentException.class,()->CensoredRegression.fit(new double[]{0,2,3,4},x,new int[4],CensoredRegression.Distribution.LOGNORMAL));
        assertThrows(IllegalArgumentException.class,()->CensoredRegression.fit(y,new double[][]{{1,1},{1,1},{1,1},{1,1}},new int[4],CensoredRegression.Distribution.GAUSSIAN));
    }
    private static double[] values(Properties p,String key){return Arrays.stream(p.getProperty(key).split(",")).mapToDouble(Double::parseDouble).toArray();}
}
