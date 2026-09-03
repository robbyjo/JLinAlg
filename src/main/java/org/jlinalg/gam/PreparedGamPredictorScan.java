/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.gam;

import java.util.List;
import java.util.Objects;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.RemlOptions;

/**
 * Retains response, parametric design, and compute context while the smooth
 * covariate changes across a high-throughput predictor scan.
 */
public final class PreparedGamPredictorScan implements AutoCloseable {
    private final double[] response;
    private final double[] parametricDesign;
    private final int observations;
    private final int parametricColumns;
    private final RemlOptions options;
    private final BackendContext context;
    private volatile boolean closed;

    /** Creates a scan from a conventional row-by-column parametric design. */
    public PreparedGamPredictorScan(
            double[] response,
            double[][] parametricDesign,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        this(response, MatrixOps.rowMajor(parametricDesign, response.length),
            response.length, parametricDesign[0].length,
            options, backendPolicy);
    }

    /** Creates a scan from a contiguous row-major parametric design. */
    public PreparedGamPredictorScan(
            double[] response,
            double[] parametricDesign,
            int observations,
            int parametricColumns,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        MatrixOps.validateModelData(
            response, parametricDesign, observations, parametricColumns);
        this.response = response.clone();
        this.parametricDesign = parametricDesign.clone();
        this.observations = observations;
        this.parametricColumns = parametricColumns;
        this.options = Objects.requireNonNull(options, "options");
        this.context = BackendContext.select(
            Objects.requireNonNull(backendPolicy, "backendPolicy"));
    }

    /** Fits a cubic second-difference P-spline for one changing predictor. */
    public GamResult fit(
            String termName, double[] covariate, int basisDimension) {
        requireOpen();
        return fit(List.of(
            PSplineTerm.of(termName, covariate, basisDimension)));
    }

    /** Fits one or more changing smooth terms with the retained null design. */
    public GamResult fit(List<PSplineTerm> smoothTerms) {
        requireOpen();
        return Gam.fitGaussian(response, parametricDesign, observations,
            parametricColumns, smoothTerms, options, context);
    }

    public int observations() { return observations; }
    public int parametricColumns() { return parametricColumns; }

    private void requireOpen() {
        if (closed)
            throw new IllegalStateException("prepared GAM predictor scan is closed");
    }

    @Override
    public void close() {
        if (closed) return;
        closed = true;
        context.close();
    }
}
