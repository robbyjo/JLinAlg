/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.List;

/** Graph-Laplacian smooth over discrete regions or network nodes. */
public final class MarkovRandomFieldSmoothTerm {
    private MarkovRandomFieldSmoothTerm() { }

    public record Edge(int first, int second, double weight) {
        public Edge {
            if (first < 0 || second < 0 || first == second
                    || !(weight > 0.0) || !Double.isFinite(weight)) {
                throw new IllegalArgumentException("invalid MRF edge");
            }
        }
    }

    /** Creates a one-hot node basis and weighted graph-Laplacian penalty. */
    public static QuadraticSmoothTerm of(
            String name, int[] observationNodes, int nodeCount, List<Edge> edges) {
        if (observationNodes == null || observationNodes.length < 2
                || nodeCount < 2 || edges == null || edges.isEmpty()) {
            throw new IllegalArgumentException("nodes and edges are required");
        }
        double[] design = new double[observationNodes.length * nodeCount];
        for (int row = 0; row < observationNodes.length; row++) {
            int node = observationNodes[row];
            if (node < 0 || node >= nodeCount) {
                throw new IllegalArgumentException("observation node is out of range");
            }
            design[row * nodeCount + node] = 1.0;
        }
        double[] penalty = new double[nodeCount * nodeCount];
        for (Edge edge : edges) {
            if (edge.first() >= nodeCount || edge.second() >= nodeCount) {
                throw new IllegalArgumentException("edge node is out of range");
            }
            penalty[edge.first() * nodeCount + edge.first()] += edge.weight();
            penalty[edge.second() * nodeCount + edge.second()] += edge.weight();
            penalty[edge.first() * nodeCount + edge.second()] -= edge.weight();
            penalty[edge.second() * nodeCount + edge.first()] -= edge.weight();
        }
        return new QuadraticSmoothTerm(name, observationNodes.length, nodeCount,
            design, List.of(penalty));
    }
}
