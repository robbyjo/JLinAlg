/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.Arrays;

/** Independent grouped random coefficients sharing an estimated covariance matrix.
 * For nesting, use composite group labels, e.g. study + '/' + outcome.
 * Random-design rows explicitly include a column of ones when an intercept is wanted. */
public final class MetaRandomEffect {
    public enum Structure { DIAGONAL, UNSTRUCTURED }
    private final String name;
    private final String[] groups;
    private final double[][] design;
    private final Structure structure;

    public MetaRandomEffect(String name, String[] groups, double[][] design, Structure structure) {
        if (name == null || name.isBlank() || groups == null || design == null || structure == null
                || groups.length != design.length || design.length < 2 || design[0] == null || design[0].length == 0)
            throw new IllegalArgumentException("random-effect name, groups, design and structure required");
        this.name = name; this.groups = groups.clone(); this.structure = structure;
        this.design = new double[design.length][];
        for (int i = 0; i < design.length; i++) {
            if (groups[i] == null || groups[i].isBlank() || design[i] == null || design[i].length != design[0].length)
                throw new IllegalArgumentException("random-effect groups and design rows are invalid");
            this.design[i] = design[i].clone();
            for (double v : design[i]) if (!Double.isFinite(v)) throw new IllegalArgumentException("random design must be finite");
        }
        if (Arrays.stream(groups).distinct().count() < 2) throw new IllegalArgumentException("at least two random-effect groups required");
    }
    public static MetaRandomEffect intercept(String name, String[] groups) {
        double[][] z = new double[groups.length][1];
        for (double[] row : z) row[0] = 1;
        return new MetaRandomEffect(name, groups, z, Structure.DIAGONAL);
    }
    public String name() { return name; }
    public String[] groups() { return groups.clone(); }
    public double[][] design() { return Arrays.stream(design).map(double[]::clone).toArray(double[][]::new); }
    public Structure structure() { return structure; }
    int rows() { return groups.length; }
    int columns() { return design[0].length; }
    int parameters() { int q = columns(); return structure == Structure.DIAGONAL ? q : q * (q + 1) / 2; }
    boolean sameGroup(int i, int j) { return groups[i].equals(groups[j]); }
    double z(int i, int j) { return design[i][j]; }

    double[] covariance(double[] theta, int offset) {
        int q = columns(); double[] l = new double[q * q], g = new double[q * q];
        for (int i = 0; i < q; i++) for (int j = 0; j <= i; j++)
            if (structure == Structure.UNSTRUCTURED || i == j) l[i * q + j] = theta[offset++];
        for (int i = 0; i < q; i++) for (int j = 0; j < q; j++)
            for (int k = 0; k < q; k++) g[i * q + j] += l[i * q + k] * l[j * q + k];
        return g;
    }
}
