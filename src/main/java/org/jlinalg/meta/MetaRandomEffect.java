/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.Arrays;
import java.util.List;

/** Independent grouped random coefficients sharing an estimated covariance matrix.
 * For nesting, use composite group labels, e.g. study + '/' + outcome.
 * Random-design rows explicitly include a column of ones when an intercept is wanted. */
public final class MetaRandomEffect {
    public enum Structure {
        DIAGONAL,
        UNSTRUCTURED,
        /** One common variance and one exchangeable correlation. */
        COMPOUND_SYMMETRY,
        /** One common variance and AR(1) correlation over design columns. */
        AR1
    }
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
        if ((structure == Structure.COMPOUND_SYMMETRY
                || structure == Structure.AR1) && design[0].length < 2)
            throw new IllegalArgumentException("structured correlation requires at least two random coefficients");
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
    int parameters() {
        int q = columns();
        return switch (structure) {
            case DIAGONAL -> q;
            case UNSTRUCTURED -> q * (q + 1) / 2;
            case COMPOUND_SYMMETRY, AR1 -> 2;
        };
    }
    boolean sameGroup(int i, int j) { return groups[i].equals(groups[j]); }
    double z(int i, int j) { return design[i][j]; }

    double[] covariance(double[] theta, int offset) {
        int q = columns();
        if (structure == Structure.COMPOUND_SYMMETRY
                || structure == Structure.AR1) {
            double variance = theta[offset] * theta[offset];
            double raw = theta[offset + 1];
            double rho;
            if (structure == Structure.AR1) rho = Math.tanh(raw);
            else {
                double lower = -1.0 / (q - 1.0);
                rho = lower + (1.0 - lower) / (1.0 + Math.exp(-raw));
            }
            double[] g = new double[q * q];
            for (int row = 0; row < q; row++) for (int column = 0; column < q; column++)
                g[row * q + column] = variance * (structure == Structure.AR1
                    ? Math.pow(rho, Math.abs(row - column))
                    : row == column ? 1.0 : rho);
            return g;
        }
        double[] l = new double[q * q], g = new double[q * q];
        for (int i = 0; i < q; i++) for (int j = 0; j <= i; j++)
            if (structure == Structure.UNSTRUCTURED || i == j) l[i * q + j] = theta[offset++];
        for (int i = 0; i < q; i++) for (int j = 0; j < q; j++)
            for (int k = 0; k < q; k++) g[i * q + j] += l[i * q + k] * l[j * q + k];
        return g;
    }

    void initialize(double scale, double[] scales, double[] start, int offset) {
        int q = columns();
        if (structure == Structure.COMPOUND_SYMMETRY || structure == Structure.AR1) {
            double norm = 0;
            for (int a = 0; a < q; a++) for (int i = 0; i < rows(); i++) norm += z(i,a)*z(i,a)/(rows()*q);
            scales[offset] = scale / Math.sqrt(Math.max(1e-12, norm));
            scales[offset+1] = 1;
            start[offset] = .5;
            start[offset+1] = structure == Structure.COMPOUND_SYMMETRY
                ? -Math.log(q - 1.0) : 0;
            return;
        }
        int k = offset;
        for (int a = 0; a < q; a++) {
            double norm = 0;
            for (int i = 0; i < rows(); i++) norm += z(i,a)*z(i,a)/rows();
            for (int b = 0; b <= a; b++) if (structure == Structure.UNSTRUCTURED || a == b) {
                scales[k] = scale / Math.sqrt(Math.max(1e-12, norm));
                start[k++] = a == b ? .5 : 0;
            }
        }
    }

    List<double[]> covarianceBases() {
        int q = columns();
        java.util.ArrayList<double[]> result = new java.util.ArrayList<>();
        if (structure == Structure.COMPOUND_SYMMETRY || structure == Structure.AR1) {
            double[] identity = new double[q*q], correlation = new double[q*q];
            for (int i=0;i<q;i++) identity[i*q+i]=1;
            if (structure == Structure.COMPOUND_SYMMETRY) {
                for(int i=0;i<q;i++)for(int j=0;j<q;j++)if(i!=j)correlation[i*q+j]=1;
            } else for(int i=0;i<q-1;i++) correlation[i*q+i+1]=correlation[(i+1)*q+i]=1;
            result.add(identity); result.add(correlation); return result;
        }
        for (int a=0;a<q;a++) for(int b=0;b<=a;b++)
            if(structure==Structure.UNSTRUCTURED||a==b){
                double[] basis=new double[q*q];
                basis[a*q+b]=1; basis[b*q+a]=1;
                result.add(basis);
            }
        return result;
    }
}
