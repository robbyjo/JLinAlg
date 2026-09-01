/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.ArrayList;
import java.util.List;
import jdistlib.ChiSquare;
import jdistlib.accelerator.CholeskyFactor;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.internal.MatrixOps;

/** Weighted fixed- and random-effects meta-regression with multiple moderators. */
public final class MetaRegression {
    private MetaRegression() { }

    /** Fits a model with an intercept and the supplied moderator columns. */
    public static MetaRegressionResult fit(
            List<MetaStudy> studies,
            double[][] moderators,
            List<String> moderatorNames) {
        return fit(studies, moderators, moderatorNames, true,
            MetaAnalysisOptions.randomEffects(), BackendPolicy.PREFERRED);
    }

    /** Fits the requested meta-regression design and heterogeneity model. */
    public static MetaRegressionResult fit(
            List<MetaStudy> studies,
            double[][] moderators,
            List<String> moderatorNames,
            boolean includeIntercept,
            MetaAnalysisOptions options,
            BackendPolicy backendPolicy) {
        MetaMath.Data data = MetaMath.data(studies);
        if (moderators == null || moderatorNames == null || options == null
                || backendPolicy == null || moderators.length != studies.size())
            throw new IllegalArgumentException(
                "moderators, names, options, and backend policy are required");
        int moderatorCount = moderators.length == 0 || moderators[0] == null
            ? 0 : moderators[0].length;
        if (moderatorCount < 1 || moderatorNames.size() != moderatorCount)
            throw new IllegalArgumentException(
                "one name is required for each moderator column");
        double[] moderatorData = MatrixOps.rowMajor(moderators, studies.size());
        int columns = moderatorCount + (includeIntercept ? 1 : 0);
        if (studies.size() <= columns)
            throw new IllegalArgumentException(
                "meta-regression needs more studies than coefficients");
        double[] design = new double[studies.size() * columns];
        for (int row = 0; row < studies.size(); row++) {
            int offset = 0;
            if (includeIntercept) design[row * columns + offset++] = 1.0;
            System.arraycopy(moderatorData, row * moderatorCount,
                design, row * columns + offset, moderatorCount);
        }
        List<String> coefficientNames = new ArrayList<>(columns);
        if (includeIntercept) coefficientNames.add("(Intercept)");
        for (String name : moderatorNames) {
            if (name == null || name.isBlank())
                throw new IllegalArgumentException("moderator names must not be blank");
            coefficientNames.add(name);
        }

        try (BackendContext context = BackendContext.select(backendPolicy)) {
            double tauSquared = MetaMath.estimateTauSquared(
                data, design, columns, options, context.backend());
            MetaMath.Fit fitted = MetaMath.fit(
                data, design, columns, tauSquared, context.backend());
            MetaMath.Fit fixed = MetaMath.fit(
                data, design, columns, 0.0, context.backend());
            double degrees = studies.size() - columns;
            double scale = MetaAnalysis.inferenceScale(
                options, fitted.qe(), degrees);
            double[] covariance = fitted.covariance().clone();
            for (int index = 0; index < covariance.length; index++)
                covariance[index] *= scale;
            double[] errors = new double[columns];
            for (int index = 0; index < columns; index++)
                errors[index] = Math.sqrt(covariance[index * columns + index]);
            AssociationStatistics statistics = MetaAnalysis.statistics(
                options, fitted.beta(), errors, degrees);

            int moderatorStart = includeIntercept ? 1 : 0;
            // Q_M remains an asymptotic Wald chi-square diagnostic; optional
            // Knapp-Hartung scaling applies to coefficient t inference only.
            ModeratorTest moderatorTest = moderatorTest(fitted.beta(),
                fitted.covariance(), columns, moderatorStart, context.backend());
            double residualP = ChiSquare.cumulative(
                fixed.qe(), degrees, false, false);
            double iSquared = fixed.qe() > 0.0
                ? 100.0 * Math.max(0.0, (fixed.qe() - degrees) / fixed.qe())
                : 0.0;
            double hSquared = Math.max(1.0, fixed.qe() / degrees);

            double[] intercept = new double[studies.size()];
            java.util.Arrays.fill(intercept, 1.0);
            double nullTau = MetaMath.estimateTauSquared(
                data, intercept, 1, options, context.backend());
            double rSquared = nullTau > 0.0
                ? 100.0 * Math.max(0.0, 1.0 - tauSquared / nullTau) : 0.0;

            return new MetaRegressionResult(coefficientNames, options.method(),
                options.tauSquaredEstimator(), options.inferenceMethod(), statistics,
                covariance, fixed.qe(), degrees, residualP,
                moderatorTest.statistic(), moderatorTest.degreesOfFreedom(),
                moderatorTest.pValue(), tauSquared, iSquared, hSquared, rSquared,
                MetaAnalysis.normalize(fitted.weights()), context.provenance());
        }
    }

    private static ModeratorTest moderatorTest(
            double[] beta, double[] covariance, int columns, int start,
            jdistlib.accelerator.ComputeBackend backend) {
        int count = columns - start;
        double[] selected = new double[count];
        double[] selectedCovariance = new double[count * count];
        for (int row = 0; row < count; row++) {
            selected[row] = beta[start + row];
            for (int column = 0; column < count; column++)
                selectedCovariance[row * count + column] = covariance[
                    (start + row) * columns + start + column];
        }
        CholeskyFactor factor = backend.dpotrf(selectedCovariance, count);
        double[] solved = factor.solve(selected);
        double statistic = 0.0;
        for (int index = 0; index < count; index++)
            statistic += selected[index] * solved[index];
        double pValue = ChiSquare.cumulative(statistic, count, false, false);
        return new ModeratorTest(statistic, count, pValue);
    }

    private record ModeratorTest(
        double statistic, double degreesOfFreedom, double pValue) { }
}
