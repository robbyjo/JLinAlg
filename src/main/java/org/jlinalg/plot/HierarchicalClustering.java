/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.plot;

import java.util.ArrayList;
import java.util.List;

/** Deterministic Euclidean average-linkage clustering for plot ordering. */
final class HierarchicalClustering {
    private HierarchicalClustering() { }

    static Node cluster(double[][] observations) {
        validate(observations);
        int count = observations.length;
        if (count == 1) return Node.leaf(0);
        int capacity = 2 * count - 1;
        double[][] distance = new double[capacity][capacity];
        Node[] nodes = new Node[capacity];
        List<Integer> active = new ArrayList<>();
        for (int first = 0; first < count; first++) {
            nodes[first] = Node.leaf(first);
            active.add(first);
            for (int second = 0; second < first; second++) {
                distance[first][second] = euclidean(observations[first],
                    observations[second]);
                distance[second][first] = distance[first][second];
            }
        }
        int next = count;
        while (active.size() > 1) {
            int bestFirst = 0, bestSecond = 1;
            double best = Double.POSITIVE_INFINITY;
            for (int first = 0; first < active.size(); first++)
                for (int second = first + 1; second < active.size(); second++) {
                    double candidate = distance[active.get(first)][active.get(second)];
                    if (candidate < best) {
                        best = candidate;
                        bestFirst = first;
                        bestSecond = second;
                    }
                }
            int first = active.get(bestFirst), second = active.get(bestSecond);
            Node left = nodes[first], right = nodes[second];
            if (left.minimumLeaf() > right.minimumLeaf()) {
                Node swap = left;
                left = right;
                right = swap;
            }
            nodes[next] = Node.branch(left, right, best);
            for (int candidate : active) {
                if (candidate == first || candidate == second) continue;
                double value = (nodes[first].size() * distance[first][candidate]
                    + nodes[second].size() * distance[second][candidate])
                    / (nodes[first].size() + nodes[second].size());
                distance[next][candidate] = value;
                distance[candidate][next] = value;
            }
            active.remove(bestSecond);
            active.remove(bestFirst);
            active.add(next++);
        }
        return nodes[active.get(0)];
    }

    private static void validate(double[][] observations) {
        if (observations == null || observations.length == 0
                || observations[0] == null || observations[0].length == 0)
            throw new IllegalArgumentException("clustering observations are required");
        int dimensions = observations[0].length;
        for (double[] row : observations) {
            if (row == null || row.length != dimensions)
                throw new IllegalArgumentException("clustering matrix must be rectangular");
            for (double value : row)
                if (!Double.isFinite(value))
                    throw new IllegalArgumentException("clustering values must be finite");
        }
    }

    private static double euclidean(double[] first, double[] second) {
        double sum = 0.0;
        for (int index = 0; index < first.length; index++) {
            double difference = first[index] - second[index];
            sum += difference * difference;
        }
        return Math.sqrt(sum);
    }

    record Node(int leaf, Node left, Node right, double height, int size,
            int minimumLeaf) {
        static Node leaf(int index) {
            return new Node(index, null, null, 0.0, 1, index);
        }

        static Node branch(Node left, Node right, double height) {
            return new Node(-1, left, right, height, left.size + right.size,
                Math.min(left.minimumLeaf, right.minimumLeaf));
        }

        boolean isLeaf() { return leaf >= 0; }

        int[] order() {
            int[] result = new int[size];
            append(result, 0);
            return result;
        }

        private int append(int[] result, int offset) {
            if (isLeaf()) {
                result[offset] = leaf;
                return offset + 1;
            }
            return right.append(result, left.append(result, offset));
        }
    }
}
