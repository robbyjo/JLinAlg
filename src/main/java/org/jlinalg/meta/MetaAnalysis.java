/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.List;
import jdistlib.ChiSquare;
import jdistlib.Normal;
import jdistlib.T;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.DegreesOfFreedomMethod;

/** Fixed- and random-effects inverse-variance meta-analysis. */
public final class MetaAnalysis {
    private MetaAnalysis() { }

    public static MetaAnalysisResult fit(List<MetaStudy> studies) {
        return fit(studies, MetaAnalysisOptions.randomEffects(),
            BackendPolicy.PREFERRED);
    }

    public static MetaAnalysisResult fit(
            List<MetaStudy> studies, MetaAnalysisOptions options,
            BackendPolicy backendPolicy) {
        if (options == null || backendPolicy == null)
            throw new IllegalArgumentException("options and backend policy are required");
        MetaMath.Data data = MetaMath.data(studies);
        int count = studies.size();
        double[] intercept = new double[count];
        java.util.Arrays.fill(intercept, 1.0);
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            double tauSquared = MetaMath.estimateTauSquared(
                data, intercept, 1, options, context.backend());
            MetaMath.Fit fitted = MetaMath.fit(
                data, intercept, 1, tauSquared, context.backend());
            MetaMath.Fit fixed = MetaMath.fit(
                data, intercept, 1, 0.0, context.backend());
            double residualDf = count - 1.0;
            double scale = inferenceScale(options, fitted.qe(), residualDf);
            double standardError = Math.sqrt(fitted.covariance()[0] * scale);
            AssociationStatistics statistics = statistics(options,
                fitted.beta(), new double[] {standardError}, residualDf);
            double critical = critical(options, residualDf);
            double lower = fitted.beta()[0] - critical * standardError;
            double upper = fitted.beta()[0] + critical * standardError;
            double predictionLower = Double.NaN;
            double predictionUpper = Double.NaN;
            if (options.method() == MetaAnalysisMethod.RANDOM_EFFECT) {
                double predictionStandardError = Math.sqrt(
                    tauSquared + standardError * standardError);
                predictionLower = fitted.beta()[0] - critical * predictionStandardError;
                predictionUpper = fitted.beta()[0] + critical * predictionStandardError;
            }
            double q = fixed.qe();
            double qP = ChiSquare.cumulative(q, residualDf, false, false);
            double iSquared = q > 0.0
                ? 100.0 * Math.max(0.0, (q - residualDf) / q) : 0.0;
            double hSquared = Math.max(1.0, q / residualDf);
            double[] normalized = normalize(fitted.weights());
            List<String> names = studies.stream().map(MetaStudy::name).toList();
            return new MetaAnalysisResult(names, options.method(),
                options.tauSquaredEstimator(), options.inferenceMethod(), statistics,
                lower, upper, predictionLower, predictionUpper, q, residualDf,
                qP, tauSquared, iSquared, hSquared, normalized,
                context.provenance());
        }
    }

    static double inferenceScale(
            MetaAnalysisOptions options, double q, double degreesOfFreedom) {
        return switch (options.inferenceMethod()) {
            case NORMAL, STUDENT_T -> 1.0;
            case HARTUNG_KNAPP -> q / degreesOfFreedom;
            case MODIFIED_HARTUNG_KNAPP -> Math.max(1.0, q / degreesOfFreedom);
        };
    }

    static AssociationStatistics statistics(
            MetaAnalysisOptions options, double[] beta, double[] errors,
            double degreesOfFreedom) {
        return options.inferenceMethod() == MetaInferenceMethod.NORMAL
            ? AssociationStatistics.normal(beta, errors)
            : AssociationStatistics.studentT(beta, errors, degreesOfFreedom,
                DegreesOfFreedomMethod.RESIDUAL);
    }

    static double critical(MetaAnalysisOptions options, double degreesOfFreedom) {
        double probability = 0.5 + options.confidenceLevel() / 2.0;
        return options.inferenceMethod() == MetaInferenceMethod.NORMAL
            ? Normal.quantile(probability, 0.0, 1.0, true, false)
            : T.quantile(probability, degreesOfFreedom, true, false);
    }

    static double[] normalize(double[] weights) {
        double sum = 0.0;
        for (double value : weights) sum += value;
        double[] result = weights.clone();
        for (int index = 0; index < result.length; index++) result[index] /= sum;
        return result;
    }
}
