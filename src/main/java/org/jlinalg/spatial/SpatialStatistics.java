/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.spatial;

import java.util.*;

/** Bounded physical-coordinate graphs and explicitly stratified random-label spatial tests. */
public final class SpatialStatistics {
    private SpatialStatistics() { }
    /** An undirected, unweighted edge; each unordered pair appears once. */
    public record Edge(int source, int target, double distance) { }
    /** Within-section autocorrelation and Monte Carlo random-label p values. */
    public record Autocorrelation(double moran, double geary, double moranP, double gearyP) { }

    /** Radius or union-kNN graph, never crossing a supplied section/compartment key. */
    public static List<Edge> graph(double[][] coordinates, String[] compartments,
            String method, double radius, int k) {
        int n=coordinates.length;
        if(n<2 || compartments.length!=n || !Set.of("radius","knn").contains(method)
                || !(radius>0) || !Double.isFinite(radius) || k<1)
            throw new IllegalArgumentException("Invalid graph inputs/method/radius/k");
        int dim=coordinates[0].length;
        if(dim!=2 && dim!=3)throw new IllegalArgumentException("Coordinates must have two or three dimensions");
        Map<String,List<Integer>> groups=new LinkedHashMap<>();
        for(int i=0;i<n;i++) {
            if(coordinates[i].length!=dim || compartments[i]==null || compartments[i].isBlank())throw new IllegalArgumentException("Invalid coordinates/compartments");
            for(double v:coordinates[i])if(!Double.isFinite(v))throw new IllegalArgumentException("Nonfinite coordinate");
            groups.computeIfAbsent(compartments[i],key->new ArrayList<>()).add(i);
        }
        long work=0;for(var group:groups.values())work+=(long)group.size()*group.size();
        if(work>25_000_000)throw new IllegalArgumentException("Graph exceeds 25 million within-compartment distance evaluations; partition the data");
        Map<Long,Edge> edges=new TreeMap<>();
        for(var group:groups.values())for(int i:group) {
            List<Edge> near=new ArrayList<>();
            for(int j:group)if(i!=j) {
                double d=0;for(int axis=0;axis<dim;axis++)d=Math.hypot(d,coordinates[i][axis]-coordinates[j][axis]);
                if(d<=radius)near.add(new Edge(i,j,d));
            }
            near.sort(Comparator.comparingDouble(Edge::distance).thenComparingInt(Edge::target));
            int take=method.equals("knn")?Math.min(k,near.size()):near.size();
            for(int a=0;a<take;a++) {
                Edge e=near.get(a);int s=Math.min(i,e.target),t=Math.max(i,e.target);
                edges.putIfAbsent(((long)s<<32)|t,new Edge(s,t,e.distance));
                if(edges.size()>1_000_000)throw new IllegalArgumentException("Graph exceeds one million edges");
            }
        }
        return List.copyOf(edges.values());
    }

    /** Unweighted symmetric I/C; includes isolated observations in the declared section universe. */
    public static double[] statistics(double[] values, List<Edge> edges) {
        if(values.length<3 || edges.isEmpty())throw new IllegalArgumentException("Need >=3 observations and at least one edge");
        double mean=0;
        for(double v:values){if(!Double.isFinite(v))throw new IllegalArgumentException("Nonfinite measurement");mean+=v/values.length;}
        double variance=0,cross=0,difference=0;
        for(double v:values)variance+=(v-mean)*(v-mean);
        if(!(variance>0) || !Double.isFinite(variance))throw new IllegalArgumentException("Constant or overflowing spatial feature");
        Set<Long> seen=new HashSet<>();
        for(Edge e:edges) {
            if(e.source<0 || e.target>=values.length || e.source>=e.target
                    || !seen.add(((long)e.source<<32)|e.target))throw new IllegalArgumentException("Invalid/duplicate undirected edge");
            cross+=(values[e.source]-mean)*(values[e.target]-mean);
            difference+=Math.pow(values[e.source]-values[e.target],2);
        }
        return new double[]{values.length*cross/(edges.size()*variance),
            (values.length-1)*difference/(2*edges.size()*variance)};
    }

    /** Two-sided tests around the unrestricted random-label expectations -1/(n-1) and 1.
     * Restricted permutations use the same fixed test statistics, with stratum labels held fixed. */
    public static Autocorrelation test(double[] values,List<Edge> edges,String[] strata,int permutations,long seed) {
        if(strata.length!=values.length || permutations<1 || permutations>100_000
                || (long)(values.length+edges.size())*permutations>100_000_000)
            throw new IllegalArgumentException("Invalid or excessive permutation work");
        List<int[]> groups=groups(strata);
        double[] observed=statistics(values,edges);int mi=0,gc=0;
        double center=-1.0/(values.length-1);Random random=new Random(seed);
        for(int b=0;b<permutations;b++) {
            double[] shuffled=values.clone();shuffle(shuffled,groups,random);
            double[] stat=statisticsUnchecked(shuffled,edges);
            if(Math.abs(stat[0]-center)>=Math.abs(observed[0]-center)-1e-12)mi++;
            if(Math.abs(stat[1]-1)>=Math.abs(observed[1]-1)-1e-12)gc++;
        }
        return new Autocorrelation(observed[0],observed[1],(mi+1.0)/(permutations+1),(gc+1.0)/(permutations+1));
    }
    private static double[] statisticsUnchecked(double[] values,List<Edge> edges) {
        double mean=Arrays.stream(values).average().orElseThrow(),v=0,c=0,d=0;
        for(double x:values)v+=(x-mean)*(x-mean);
        for(Edge e:edges){c+=(values[e.source]-mean)*(values[e.target]-mean);d+=Math.pow(values[e.source]-values[e.target],2);}
        return new double[]{values.length*c/(edges.size()*v),(values.length-1)*d/(2*edges.size()*v)};
    }
    /** Exchangeability strata for random-label tests. */
    public static List<int[]> groups(String[] strata) {
        Map<String,List<Integer>> groups=new LinkedHashMap<>();
        for(int i=0;i<strata.length;i++) {
            if(strata[i]==null || strata[i].isBlank())throw new IllegalArgumentException("Blank permutation stratum");
            groups.computeIfAbsent(strata[i],s->new ArrayList<>()).add(i);
        }
        return groups.values().stream().map(g->g.stream().mapToInt(Integer::intValue).toArray()).toList();
    }
    /** Fisher-Yates permutation within each stratum. */
    public static void shuffle(double[] values,List<int[]> groups,Random random) {
        for(int[] g:groups)for(int i=g.length-1;i>0;i--){int j=random.nextInt(i+1);double v=values[g[i]];values[g[i]]=values[g[j]];values[g[j]]=v;}
    }
}
