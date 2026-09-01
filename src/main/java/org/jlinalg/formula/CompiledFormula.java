/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.formula;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.Glm;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.glm.GlmOptions;
import org.jlinalg.glm.GlmResult;
import org.jlinalg.ols.Ols;
import org.jlinalg.ols.OlsOptions;
import org.jlinalg.ols.OlsResult;

/** A formula compiled into reusable contiguous arrays. */
public final class CompiledFormula {
    private final double[] response;
    private final double[] design;
    private final int rows;
    private final int columns;
    private final List<String> coefficientNames;
    private final double[] weights;
    private final double[] offset;

    CompiledFormula(
            double[] response, double[] design, int rows, int columns,
            List<String> coefficientNames,
            double[] weights, double[] offset) {
        this.response = response;
        this.design = design;
        this.rows = rows;
        this.columns = columns;
        this.coefficientNames = List.copyOf(coefficientNames);
        this.weights = weights;
        this.offset = offset;
    }

    public int rows() { return rows; }
    public int columns() { return columns; }
    public List<String> coefficientNames() { return coefficientNames; }
    public double[] response() { return response.clone(); }
    public double[] design() { return design.clone(); }
    public double[] weights() { return weights == null ? null : weights.clone(); }
    public double[] offset() { return offset == null ? null : offset.clone(); }

    double[] responseView() { return response; }
    double[] designView() { return design; }
    double[] weightsView() { return weights; }
    double[] offsetView() { return offset; }

    /** Fits OLS without reparsing or rebuilding the model matrix. */
    public OlsResult fitOls(OlsOptions options, BackendPolicy backendPolicy) {
        return Ols.fit(response, design, rows, columns,
            weights, offset, options, backendPolicy);
    }

    /** Fits a GLM without reparsing or rebuilding the model matrix. */
    public GlmResult fitGlm(
            GlmFamily family, GlmOptions options, BackendPolicy backendPolicy) {
        return Glm.fit(response, design, rows, columns,
            family, weights, offset, options, backendPolicy);
    }
}
