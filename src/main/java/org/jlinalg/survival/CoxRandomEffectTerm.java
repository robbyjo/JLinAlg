/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import java.util.List;
import java.util.Objects;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.mixed.RandomEffectTerm;

/** Random-effect design and unit-variance Gaussian precision for Cox frailty. */
public final class CoxRandomEffectTerm {
    private final String name;
    private final int observations;
    private final int coefficients;
    private final double[] design;
    private final double[] precision;
    private final List<String> coefficientNames;

    public CoxRandomEffectTerm(
            String name, double[][] design, double[] precision,
            List<String> coefficientNames) {
        if (design == null || design.length == 0 || design[0] == null)
            throw new IllegalArgumentException("frailty design is required");
        this.name = Objects.requireNonNull(name, "name");
        observations = design.length;
        coefficients = design[0].length;
        this.design = MatrixOps.rowMajor(design, observations);
        if (precision == null || precision.length != coefficients * coefficients
                || coefficientNames == null
                || coefficientNames.size() != coefficients)
            throw new IllegalArgumentException(
                "frailty precision and coefficient names have invalid dimensions");
        this.precision = MatrixOps.finiteCopy(precision, "frailty precision");
        this.coefficientNames = List.copyOf(coefficientNames);
        validatePrecision();
    }

    private CoxRandomEffectTerm(
            String name, int observations, int coefficients,
            double[] design, double[] precision,
            List<String> coefficientNames) {
        this.name = name;
        this.observations = observations;
        this.coefficients = coefficients;
        this.design = design;
        this.precision = precision;
        this.coefficientNames = coefficientNames;
        validatePrecision();
    }

    public static CoxRandomEffectTerm independent(RandomEffectTerm term) {
        if (term == null)
            throw new IllegalArgumentException("random-effect term is required");
        return new CoxRandomEffectTerm(term.name(), term.observations(),
            term.coefficients(), term.design(),
            MatrixOps.identity(term.coefficients()), term.coefficientNames());
    }

    public String name() { return name; }
    public int observations() { return observations; }
    public int coefficients() { return coefficients; }
    public double[] design() { return design.clone(); }
    public double[] precision() { return precision.clone(); }
    public List<String> coefficientNames() { return coefficientNames; }
    double[] designView() { return design; }
    double[] precisionView() { return precision; }

    private void validatePrecision() {
        if (name == null || name.isBlank() || coefficients < 1)
            throw new IllegalArgumentException("frailty term name/dimensions are invalid");
        for (int row = 0; row < coefficients; row++) {
            if (!(precision[row * coefficients + row] > 0))
                throw new IllegalArgumentException(
                    "frailty precision diagonal must be positive");
            for (int column = 0; column < row; column++) {
                double left = precision[row * coefficients + column];
                double right = precision[column * coefficients + row];
                if (Math.abs(left - right) > 1e-10
                        * Math.max(1, Math.max(Math.abs(left), Math.abs(right))))
                    throw new IllegalArgumentException(
                        "frailty precision must be symmetric");
            }
        }
    }
}
