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
        if(method.equals("radius")) {
            List<Edge> edges=new ArrayList<>();
            for(var group:groups.values())for(int a=0;a<group.size();a++)for(int b=a+1;b<group.size();b++) {
                int i=group.get(a),j=group.get(b);double d=distance(coordinates[i],coordinates[j]);
                if(d<=radius){edges.add(new Edge(i,j,d));if(edges.size()>1_000_000)throw new IllegalArgumentException("Graph exceeds one million edges");}
            }
            edges.sort(Comparator.comparingInt(Edge::source).thenComparingInt(Edge::target));
            return List.copyOf(edges);
        }
        Map<Long,Edge> edges=new TreeMap<>();
        Comparator<Edge> order=Comparator.comparingDouble(Edge::distance).thenComparingInt(Edge::target);
        for(var group:groups.values())for(int i:group) {
            PriorityQueue<Edge> near=new PriorityQueue<>(order.reversed());
            for(int j:group)if(i!=j) {
                double d=distance(coordinates[i],coordinates[j]);
                if(d<=radius) {
                    Edge candidate=new Edge(i,j,d);
                    if(near.size()<k)near.add(candidate);
                    else if(order.compare(candidate,near.element())<0){near.remove();near.add(candidate);}
                }
            }
            for(Edge e:near) {
                int s=Math.min(i,e.target),t=Math.max(i,e.target);
                edges.putIfAbsent(((long)s<<32)|t,new Edge(s,t,e.distance));
                if(edges.size()>1_000_000)throw new IllegalArgumentException("Graph exceeds one million edges");
            }
        }
        return List.copyOf(edges.values());
    }
    private static double distance(double[] a,double[] b) {
        double sum=0,max=0;
        for(int j=0;j<a.length;j++){double delta=a[j]-b[j];sum+=delta*delta;max=Math.max(max,Math.abs(delta));}
        if(sum>=Double.MIN_NORMAL && Double.isFinite(sum))return Math.sqrt(sum);
        if(max==0 || !Double.isFinite(max))return max;
        sum=0;for(int j=0;j<a.length;j++){double delta=(a[j]-b[j])/max;sum+=delta*delta;}
        return max*Math.sqrt(sum);
    }

    private record Centered(double[] values,double scale,double sumSquares) { }
    private static Centered center(double[] values) {
        if(values.length==0)throw new IllegalArgumentException("Empty measurements");
        double max=0;for(double v:values){if(!Double.isFinite(v))throw new IllegalArgumentException("Nonfinite measurement");max=Math.max(max,Math.abs(v));}
        double[] z=new double[values.length];double scale=0;
        for(int i=0;i<z.length;i++){z[i]=values[i]-values[0];scale=Math.max(scale,Math.abs(z[i]));}
        if(scale==0)return new Centered(z,1,0);
        if(Double.isInfinite(scale)) {
            // Opposite-sign finite extremes can overflow their difference.
            scale=max;for(int i=0;i<z.length;i++)z[i]=values[i]/scale-values[0]/scale;
        } else for(int i=0;i<z.length;i++)z[i]/=scale;
        double mean=0,correction=0;
        for(double v:z){double term=v/z.length-correction,total=mean+term;correction=(total-mean)-term;mean=total;}
        double variance=0;for(int i=0;i<z.length;i++){z[i]-=mean;variance+=z[i]*z[i];}
        return new Centered(z,scale,variance);
    }
    private static void validateEdges(int n,List<Edge> edges) {
        if(n<3 || edges.isEmpty())throw new IllegalArgumentException("Need >=3 observations and at least one edge");
        Set<Long> seen=new HashSet<>();
        for(Edge e:edges)if(e.source<0 || e.target>=n || e.source>=e.target
                || !Double.isFinite(e.distance) || e.distance<0 || !seen.add(((long)e.source<<32)|e.target))
            throw new IllegalArgumentException("Invalid/duplicate undirected edge");
    }

    /** Unweighted symmetric I/C; includes isolated observations in the declared section universe. */
    public static double[] statistics(double[] values, List<Edge> edges) {
        validateEdges(values.length,edges);Centered centered=center(values);
        if(!(centered.sumSquares>0))throw new IllegalArgumentException("Constant spatial feature");
        return statisticsUnchecked(centered.values,edges,centered.sumSquares);
    }

    /** Two-sided tests around the unrestricted random-label expectations -1/(n-1) and 1.
     * Restricted permutations use the same fixed test statistics, with stratum labels held fixed. */
    public static Autocorrelation test(double[] values,List<Edge> edges,String[] strata,int permutations,long seed) {
        if(strata.length!=values.length || permutations<1 || permutations>100_000
                || (long)(values.length+edges.size())*permutations>100_000_000)
            throw new IllegalArgumentException("Invalid or excessive permutation work");
        List<int[]> groups=groups(strata);validateEdges(values.length,edges);Centered centered=center(values);
        if(!(centered.sumSquares>0))throw new IllegalArgumentException("Constant spatial feature");
        double[] observed=statisticsUnchecked(centered.values,edges,centered.sumSquares);int mi=0,gc=0;
        double center=-1.0/(values.length-1);Random random=new Random(seed);
        double[] shuffled=new double[values.length];
        for(int b=0;b<permutations;b++) {
            System.arraycopy(centered.values,0,shuffled,0,shuffled.length);shuffle(shuffled,groups,random);
            double[] stat=statisticsUnchecked(shuffled,edges,centered.sumSquares);
            if(Math.abs(stat[0]-center)>=Math.abs(observed[0]-center)-1e-12)mi++;
            if(Math.abs(stat[1]-1)>=Math.abs(observed[1]-1)-1e-12)gc++;
        }
        return new Autocorrelation(observed[0],observed[1],(mi+1.0)/(permutations+1),(gc+1.0)/(permutations+1));
    }
    private static double[] statisticsUnchecked(double[] values,List<Edge> edges,double v) {
        double c=0,d=0;
        for(Edge e:edges){c+=values[e.source]*values[e.target];double delta=values[e.source]-values[e.target];d+=delta*delta;}
        return new double[]{values.length*c/(edges.size()*v),(values.length-1)*d/(2*edges.size()*v)};
    }

    /** Prespecified descriptive slope, centered and scaled before cross-products.
     * No independent-observation inference is attached to this within-sample summary. */
    public static double linearSlope(double[] distance,double[] response) {
        if(distance.length<3 || distance.length!=response.length)throw new IllegalArgumentException("Slope requires >=3 aligned observations");
        Centered x=center(distance),y=center(response);
        if(!(x.sumSquares>0))throw new IllegalArgumentException("Distance has no within-sample variation");
        double cross=0;for(int i=0;i<distance.length;i++)cross+=x.values[i]*y.values[i];
        double coefficient=cross/x.sumSquares;
        if(coefficient==0)return 0;
        // Combine binary exponents separately so an intermediate scale ratio cannot
        // overflow or underflow when the final slope is representable.
        int a=Math.getExponent(coefficient),b=Math.getExponent(y.scale),c=Math.getExponent(x.scale);
        double slope=Math.scalb(Math.scalb(coefficient,-a)*Math.scalb(y.scale,-b)/Math.scalb(x.scale,-c),a+b-c);
        if(!Double.isFinite(slope))throw new IllegalArgumentException("Distance slope is not representable; rescale physical units");
        return slope;
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
