/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.List;
import java.util.Objects;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.RemlOptions;

/**
 * Reusable Gaussian GAM structure for repeated-response fits.
 *
 * <p>Spline knots, basis matrices, and penalty factors are constructed once.
 * Refits additionally warm-start variance components from the previous fit.</p>
 */
public final class PreparedGam {
    private final double[] parametricDesign;
    private final int observations;
    private final int parametricColumns;
    private final List<PSplineTerm> smoothTerms;
    private final RemlOptions options;
    private final BackendPolicy backendPolicy;

    /** Prepares a model from conventional fixed-effect arrays. */
    public PreparedGam(
            double[][] parametricDesign,
            List<PSplineTerm> smoothTerms,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        if (parametricDesign == null || parametricDesign.length == 0
                || parametricDesign[0] == null) {
            throw new IllegalArgumentException("parametric design is required");
        }
        observations = parametricDesign.length;
        parametricColumns = parametricDesign[0].length;
        this.parametricDesign = MatrixOps.rowMajor(
            parametricDesign, observations);
        this.smoothTerms = validateTerms(smoothTerms, observations);
        this.options = Objects.requireNonNull(options, "options");
        this.backendPolicy = Objects.requireNonNull(
            backendPolicy, "backendPolicy");
    }

    /** Prepares with default REML controls and preferred acceleration. */
    public static PreparedGam of(
            double[][] parametricDesign, List<PSplineTerm> smoothTerms) {
        return new PreparedGam(parametricDesign, smoothTerms,
            RemlOptions.defaults(), BackendPolicy.PREFERRED);
    }

    /** Fits a response against the retained design and smooth bases. */
    public GamResult fit(double[] response) {
        validateResponse(response);
        return Gam.fitGaussian(response, parametricDesign,
            observations, parametricColumns, smoothTerms,
            options, backendPolicy);
    }

    /** Refits using the previous variance components as REML starting values. */
    public GamResult refit(GamResult previous, double[] response) {
        Objects.requireNonNull(previous, "previous");
        validateResponse(response);
        RemlOptions warm = options.toBuilder().initialVariances(
            previous.mixedModel().reml().varianceComponents()).build();
        return Gam.fitGaussian(response, parametricDesign,
            observations, parametricColumns, smoothTerms,
            warm, backendPolicy);
    }

    public int observations() { return observations; }
    public int parametricColumns() { return parametricColumns; }
    public List<PSplineTerm> smoothTerms() { return smoothTerms; }

    private void validateResponse(double[] response) {
        if (response == null || response.length != observations) {
            throw new IllegalArgumentException(
                "response length must equal prepared observations");
        }
        MatrixOps.requireFinite(response, "response");
    }

    private static List<PSplineTerm> validateTerms(
            List<PSplineTerm> terms, int observations) {
        if (terms == null || terms.isEmpty()) {
            throw new IllegalArgumentException("at least one smooth term is required");
        }
        List<PSplineTerm> result = List.copyOf(terms);
        for (PSplineTerm term : result) {
            if (term == null || term.observations() != observations) {
                throw new IllegalArgumentException(
                    "smooth terms must match prepared observations");
            }
        }
        return result;
    }
}
