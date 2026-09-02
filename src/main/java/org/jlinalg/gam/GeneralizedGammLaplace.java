/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.glmm.GlmmLaplace;
import org.jlinalg.glmm.GlmmLaplaceOptions;
import org.jlinalg.glmm.GlmmLaplaceResult;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.VarianceComponent;

/** Accuracy-oriented generalized additive mixed models using Laplace integration. */
public final class GeneralizedGammLaplace {
    private GeneralizedGammLaplace() { }

    /** Fits smooth, grouped, pedigree, GRM, or cryptic-relatedness components jointly. */
    public static GeneralizedGammLaplaceResult fit(
            double[] response,
            double[][] parametricDesign,
            List<PSplineTerm> smoothTerms,
            GlmFamily family,
            List<VarianceComponent> randomComponents,
            double[] priorWeights,
            double[] offset,
            GlmmLaplaceOptions options,
            BackendPolicy backendPolicy) {
        if (response == null) throw new IllegalArgumentException("response is required");
        double[] fixed = MatrixOps.rowMajor(parametricDesign, response.length);
        int parametricColumns = parametricDesign[0].length;
        List<VarianceComponent> additional = randomComponents == null
            ? List.of() : List.copyOf(randomComponents);
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            PSplineMixedModelCompiler.Compiled compiled =
                PSplineMixedModelCompiler.compile(fixed, response.length,
                    parametricColumns, smoothTerms, backend);
            List<VarianceComponent> components = new ArrayList<>();
            HashSet<String> names = new HashSet<>();
            for (PSplineMixedModelCompiler.Term term : compiled.terms()) {
                if (!names.add(term.term().name())) {
                    throw new IllegalArgumentException("duplicate component name: "
                        + term.term().name());
                }
                components.add(new VarianceComponent(term.term().name(),
                    response.length,
                    PSplineMixedModelCompiler.covarianceBasis(term)));
            }
            for (VarianceComponent component : additional) {
                if (component == null || component.dimension() != response.length
                        || !names.add(component.name())) {
                    throw new IllegalArgumentException(
                        "random components need unique names and matching observations");
                }
                components.add(component);
            }
            GlmmLaplaceResult fit = GlmmLaplace.fit(response,
                compiled.fixedDesign(), response.length, compiled.fixedColumns(),
                family, components, priorWeights, offset, options, backendPolicy);
            Map<String, double[]> smooth = new LinkedHashMap<>();
            Map<String, Double> smoothing = new LinkedHashMap<>();
            Map<String, double[]> random = new LinkedHashMap<>();
            double[] variances = fit.varianceComponents();
            int index = 0;
            for (PSplineMixedModelCompiler.Term term : compiled.terms()) {
                smooth.put(term.term().name(),
                    fit.componentPredictor(term.term().name()));
                smoothing.put(term.term().name(), 1.0 / variances[index++]);
            }
            for (VarianceComponent component : additional) {
                random.put(component.name(), fit.componentPredictor(component.name()));
            }
            return new GeneralizedGammLaplaceResult(fit,
                parametricColumns, smooth, smoothing, random);
        }
    }
}
