/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

import static org.junit.jupiter.api.Assertions.*;
import java.io.*;
import java.util.*;
import org.junit.jupiter.api.Test;

class OrdinalRegressionTest {
    @Test void orderedLogitAndProbitMatchMassPolr() throws IOException {
        double[][] data;
        try(var reader=new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("/r-reference/regression-inference/ordinal.tsv")))){
            data=reader.lines().skip(1).map(s->Arrays.stream(s.split("\t")).mapToDouble(Double::parseDouble).toArray()).toArray(double[][]::new);
        }
        Properties ref=new Properties();try(var in=getClass().getResourceAsStream("/r-reference/regression-inference/reference.properties")){ref.load(in);}
        int[] y=Arrays.stream(data).mapToInt(r->(int)r[0]).toArray();double[][] x=Arrays.stream(data).map(r->new double[]{r[1],r[2]}).toArray(double[][]::new);
        for(var link:OrdinalRegression.Link.values()) {
            var fit=OrdinalRegression.fit(y,x,4,link);String key="ordinal."+link.name().toLowerCase(Locale.ROOT);
            assertArrayEquals(values(ref,key+".beta"),fit.coefficients(),3e-6);
            assertArrayEquals(values(ref,key+".thresholds"),fit.thresholds(),3e-6);
            assertArrayEquals(values(ref,key+".cov"),fit.covariance(),3e-6);
            assertEquals(values(ref,key+".logLik")[0],fit.logLikelihood(),1e-7);
            for(double value:new double[]{-30,-5,0,5,30}){
                double[] probs=fit.probabilities(new double[]{value,0});assertEquals(1,Arrays.stream(probs).sum(),1e-12);for(double prob:probs)assertTrue(prob>=0&&prob<=1);
            }
        }
        assertThrows(IllegalArgumentException.class,()->OrdinalRegression.fit(y,Arrays.stream(x).map(r->new double[]{1,r[0]}).toArray(double[][]::new),4,OrdinalRegression.Link.LOGIT));
        assertThrows(IllegalArgumentException.class,()->OrdinalRegression.fit(y,x,5,OrdinalRegression.Link.LOGIT));
    }
    private static double[] values(Properties p,String key){return Arrays.stream(p.getProperty(key).split(",")).mapToDouble(Double::parseDouble).toArray();}
    @Test void rejectsSeparatedBinaryOrdinalFit() {
        double[][] x={{-4},{-3},{-2},{-1},{1},{2},{3},{4}};int[] y={0,0,0,0,1,1,1,1};
        for(var link:OrdinalRegression.Link.values())assertThrows(IllegalArgumentException.class,()->OrdinalRegression.fit(y,x,2,link));
    }
}
