/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.differential;
import static org.junit.jupiter.api.Assertions.*;
import java.util.Arrays;
import java.util.Random;
import org.junit.jupiter.api.Test;

class VoomTrendTest {
    @Test void interpolatedTrendTracksIndependentBruteForceWeights() {
        Random random = new Random(481);
        double[] x = new double[1500], y = new double[x.length], query = new double[15000];
        for (int i = 0; i < x.length; i++) { x[i] = random.nextGaussian() * 2; y[i] = .4 + random.nextDouble(); }
        for (int i = 0; i < query.length; i++) query[i] = random.nextGaussian() * 2.5;
        VoomTrend trend = new VoomTrend(x, y, query);
        for (int i = 0; i < query.length; i += 37) {
            double expected = direct(query[i], x, y);
            // Precision weights use the fourth power; check that estimand too.
            assertEquals(1, Math.pow(expected / trend.at(query[i]), 4), .002);
        }
    }
    @Test void tiedLocationsAndBoundaryQueriesAreStable() {
        double[] x = new double[300], y = new double[300]; Arrays.fill(y, 2);
        VoomTrend trend = new VoomTrend(x, y, new double[] {-100,0,100});
        for (double q : new double[] {-100,-.1,0,.1,100}) assertEquals(2, trend.at(q), 1e-12);
    }
    private static double direct(double q, double[] x, double[] y) {
        double[] distance = Arrays.stream(x).map(v -> Math.abs(v-q)).sorted().toArray();
        double radius = Math.max(1e-12, distance[(x.length+1)/2-1]);
        double numerator = 0, denominator = 0;
        for(int i=0;i<x.length;i++) {
            double u=Math.abs(x[i]-q)/radius;
            if(u>1)continue;
            double weight=Math.pow(1-Math.pow(u,3),3);
            numerator+=weight*y[i];denominator+=weight;
        }
        return numerator/denominator;
    }
}
