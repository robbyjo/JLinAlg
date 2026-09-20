/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.spatial;

import java.util.*;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class SpatialStatisticsTest {
    @Test void pathStatisticsMatchHandCalculatedAndExhaustiveNull() {
        var edges=List.of(new SpatialStatistics.Edge(0,1,1),new SpatialStatistics.Edge(1,2,1),new SpatialStatistics.Edge(2,3,1));
        double[] values={1,2,3,4};double[] actual=SpatialStatistics.statistics(values,edges);
        assertEquals(1.0/3,actual[0],1e-15);assertEquals(.3,actual[1],1e-15);
        // Enumerate all 24 assignments independently, directly using centered products.
        int count=0,extremeI=0,extremeC=0;
        for(int a=1;a<=4;a++)for(int b=1;b<=4;b++)for(int c=1;c<=4;c++)for(int d=1;d<=4;d++) {
            if(Set.of(a).contains(b)||a==c||a==d||b==c||b==d||c==d)continue;
            double i=4*((a-2.5)*(b-2.5)+(b-2.5)*(c-2.5)+(c-2.5)*(d-2.5))/15;
            double g=((a-b)*(a-b)+(b-c)*(b-c)+(c-d)*(c-d))/10.0;
            count++;if(Math.abs(i+1.0/3)>=2.0/3-1e-12)extremeI++;if(Math.abs(g-1)>=.7-1e-12)extremeC++;
        }
        assertEquals(24,count);
        var fit=SpatialStatistics.test(values,edges,new String[]{"a","a","a","a"},19999,51);
        assertEquals(extremeI/24.0,fit.moranP(),.015);assertEquals(extremeC/24.0,fit.gearyP(),.015);
        assertEquals(fit,SpatialStatistics.test(values,edges,new String[]{"a","a","a","a"},19999,51));
    }
    @Test void graphSeparatesSectionsGapsAndCompartmentWithThreeDimensions() {
        double[][] xyz={{0,0,0},{1,0,0},{0,0,2},{0,0,0},{100,0,0}};
        String[] group={"s1","s1","s1","s2","s1"};
        var radius=SpatialStatistics.graph(xyz,group,"radius",1.5,1);
        assertEquals(List.of(new SpatialStatistics.Edge(0,1,1)),radius);
        var knn=SpatialStatistics.graph(xyz,group,"knn",2.1,1);
        assertEquals(2,knn.size());assertTrue(knn.stream().noneMatch(e->e.source()==3||e.target()>=3));
    }
    @Test void singletonPermutationStrataGiveNoEvidenceAndBadGraphsReject() {
        var edges=List.of(new SpatialStatistics.Edge(0,1,1),new SpatialStatistics.Edge(1,2,1));
        var result=SpatialStatistics.test(new double[]{1,2,9},edges,new String[]{"a","b","c"},99,3);
        assertEquals(1,result.moranP());assertEquals(1,result.gearyP());
        assertThrows(IllegalArgumentException.class,()->SpatialStatistics.statistics(new double[]{1,1,1},edges));
        assertThrows(IllegalArgumentException.class,()->SpatialStatistics.statistics(new double[]{1,2,3},List.of(edges.get(0),edges.get(0))));
    }
}
