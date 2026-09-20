/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.spatial;

import java.util.*;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class SpatialStatisticsTest {
    @Test void constantDecimalMeasurementsRejectDespiteMeanRounding() {
        double[] constant=new double[7];Arrays.fill(constant,.1);
        var edges=List.of(new SpatialStatistics.Edge(0,1,1),new SpatialStatistics.Edge(1,2,1));
        assertThrows(IllegalArgumentException.class,()->SpatialStatistics.statistics(constant,edges));
    }
    @Test void statisticsAreInvariantToExtremeFiniteScales() {
        var edges=List.of(new SpatialStatistics.Edge(0,1,1),new SpatialStatistics.Edge(1,2,1),new SpatialStatistics.Edge(2,3,1));
        for(double scale:new double[]{1e-200,1e200}) {
            double[] actual=SpatialStatistics.statistics(new double[]{scale,2*scale,3*scale,4*scale},edges);
            assertArrayEquals(new double[]{1.0/3,.3},actual,1e-14);
        }
        assertArrayEquals(SpatialStatistics.statistics(new double[]{-1,-.5,.5,1},edges),
            SpatialStatistics.statistics(new double[]{-1e308,-5e307,5e307,1e308},edges),1e-14);
    }
    @Test void slopesHandleExtremeUnitsAndRejectConstantDecimalDistances() {
        for(double scale:new double[]{1e-200,1,1e200}) {
            double actual=SpatialStatistics.linearSlope(new double[]{scale,2*scale,3*scale,4*scale},new double[]{2,4,6,8});
            assertEquals(2,actual*scale,1e-14);
        }
        double[] x=new double[7];Arrays.fill(x,.1);
        assertThrows(IllegalArgumentException.class,()->SpatialStatistics.linearSlope(x,new double[]{1,2,3,4,5,6,7}));
        assertEquals(0,SpatialStatistics.linearSlope(new double[]{0,1e-320,2e-320},new double[]{.1,.1,.1}));
        assertEquals(Math.scalb(1.0,1000),SpatialStatistics.linearSlope(
            new double[]{0,Math.scalb(1.0,-1000),Math.scalb(1.0,-999)},
            new double[]{Math.scalb(1.0,40),0,Math.scalb(1.0,40)+2}),Math.scalb(1.0,960));
    }
    @Test void optimizedGraphsMatchIndependentAllPairReferenceIncludingTies() {
        Random random=new Random(32);int n=60;double[][] xyz=new double[n][3];String[] groups=new String[n];
        for(int i=0;i<n;i++){groups[i]="section"+(i%3);for(int j=0;j<3;j++)xyz[i][j]=random.nextInt(6);}
        for(String method:List.of("radius","knn"))for(int k:List.of(1,4,20)) {
            Set<String> expected=new TreeSet<>();
            for(int i=0;i<n;i++) {
                final int source=i;List<Integer> near=new ArrayList<>();
                for(int j=0;j<n;j++)if(i!=j&&groups[i].equals(groups[j])&&squared(xyz[i],xyz[j])<=9)near.add(j);
                near.sort(Comparator.<Integer>comparingDouble(j->squared(xyz[source],xyz[j])).thenComparingInt(j->j));
                for(int j:near.subList(0,method.equals("radius")?near.size():Math.min(k,near.size())))expected.add(Math.min(i,j)+":"+Math.max(i,j));
            }
            Set<String> actual=new TreeSet<>();for(var e:SpatialStatistics.graph(xyz,groups,method,3,k))actual.add(e.source()+":"+e.target());
            assertEquals(expected,actual);
        }
    }
    private static double squared(double[] a,double[] b) {
        double sum=0;for(int i=0;i<a.length;i++)sum+=(a[i]-b[i])*(a[i]-b[i]);return sum;
    }
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
        var tiny=SpatialStatistics.graph(new double[][]{{0,0},{1e-160,1e-160}},new String[]{"s","s"},"radius",1e-159,1);
        assertEquals(Math.sqrt(2)*1e-160,tiny.get(0).distance(),1e-174);
    }
    @Test void singletonPermutationStrataGiveNoEvidenceAndBadGraphsReject() {
        var edges=List.of(new SpatialStatistics.Edge(0,1,1),new SpatialStatistics.Edge(1,2,1));
        var result=SpatialStatistics.test(new double[]{1,2,9},edges,new String[]{"a","b","c"},99,3);
        assertEquals(1,result.moranP());assertEquals(1,result.gearyP());
        assertThrows(IllegalArgumentException.class,()->SpatialStatistics.statistics(new double[]{1,1,1},edges));
        assertThrows(IllegalArgumentException.class,()->SpatialStatistics.statistics(new double[]{1,2,3},List.of(edges.get(0),edges.get(0))));
    }
}
