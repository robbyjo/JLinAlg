/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.prediction;

import java.util.Objects;
import jdistlib.Normal;
import jdistlib.T;
import org.jlinalg.gee.GeeInference;
import org.jlinalg.gee.GeeResult;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.glm.GlmResult;
import org.jlinalg.internal.MatrixOps;

/**
 * Common response-scale prediction, scenario, and marginal-effect inference
 * for fitted GLMs and marginal GEE models.
 *
 * <p>All intervals quantify coefficient uncertainty in an expected response or
 * derived estimand. They deliberately exclude residual, process, random-effect,
 * and future-outcome variation.</p>
 */
public final class ModelPredictions {
    private ModelPredictions() { }

    /** Expected responses for GLM design rows, using zero offsets and 95% CIs. */
    public static ExpectedResponse[] expectedResponses(
            GlmResult fit, GlmFamily family, double[][] design) {
        return expectedResponses(model(fit), family, design, null, 0.95);
    }

    /** Expected responses for GLM design rows and supplied link-scale offsets. */
    public static ExpectedResponse[] expectedResponses(
            GlmResult fit, GlmFamily family, double[][] design,
            double[] offset, double confidenceLevel) {
        return expectedResponses(model(fit), family, design, offset, confidenceLevel);
    }

    /** Expected responses for GEE design rows, using zero offsets and 95% CIs. */
    public static ExpectedResponse[] expectedResponses(
            GeeResult fit, GlmFamily family, double[][] design) {
        return expectedResponses(model(fit), family, design, null, 0.95);
    }

    /** Expected responses for GEE design rows and supplied link-scale offsets. */
    public static ExpectedResponse[] expectedResponses(
            GeeResult fit, GlmFamily family, double[][] design,
            double[] offset, double confidenceLevel) {
        return expectedResponses(model(fit), family, design, offset, confidenceLevel);
    }

    /** Population-average response difference between two GLM scenarios. */
    public static ScenarioDifference scenarioDifference(
            GlmResult fit, GlmFamily family,
            double[][] first, double[][] second) {
        return scenarioDifference(model(fit), family, first, second,
            null, null, Averaging.POPULATION_AVERAGE, 0.95);
    }

    /** Response difference between two GLM scenarios under explicit semantics. */
    public static ScenarioDifference scenarioDifference(
            GlmResult fit, GlmFamily family,
            double[][] first, double[][] second,
            double[] firstOffset, double[] secondOffset,
            Averaging averaging, double confidenceLevel) {
        return scenarioDifference(model(fit), family, first, second,
            firstOffset, secondOffset, averaging, confidenceLevel);
    }

    /** Population-average response difference between two GEE scenarios. */
    public static ScenarioDifference scenarioDifference(
            GeeResult fit, GlmFamily family,
            double[][] first, double[][] second) {
        return scenarioDifference(model(fit), family, first, second,
            null, null, Averaging.POPULATION_AVERAGE, 0.95);
    }

    /** Response difference between two GEE scenarios under explicit semantics. */
    public static ScenarioDifference scenarioDifference(
            GeeResult fit, GlmFamily family,
            double[][] first, double[][] second,
            double[] firstOffset, double[] secondOffset,
            Averaging averaging, double confidenceLevel) {
        return scenarioDifference(model(fit), family, first, second,
            firstOffset, secondOffset, averaging, confidenceLevel);
    }

    /** Population-average ratio of positive GLM response means. */
    public static RiskRatio riskRatio(
            GlmResult fit, GlmFamily family,
            double[][] first, double[][] second) {
        return riskRatio(model(fit), family, first, second,
            null, null, Averaging.POPULATION_AVERAGE, 0.95);
    }

    /** Ratio of positive GLM response means under explicit semantics. */
    public static RiskRatio riskRatio(
            GlmResult fit, GlmFamily family,
            double[][] first, double[][] second,
            double[] firstOffset, double[] secondOffset,
            Averaging averaging, double confidenceLevel) {
        return riskRatio(model(fit), family, first, second,
            firstOffset, secondOffset, averaging, confidenceLevel);
    }

    /** Population-average ratio of positive GEE response means. */
    public static RiskRatio riskRatio(
            GeeResult fit, GlmFamily family,
            double[][] first, double[][] second) {
        return riskRatio(model(fit), family, first, second,
            null, null, Averaging.POPULATION_AVERAGE, 0.95);
    }

    /** Ratio of positive GEE response means under explicit semantics. */
    public static RiskRatio riskRatio(
            GeeResult fit, GlmFamily family,
            double[][] first, double[][] second,
            double[] firstOffset, double[] secondOffset,
            Averaging averaging, double confidenceLevel) {
        return riskRatio(model(fit), family, first, second,
            firstOffset, secondOffset, averaging, confidenceLevel);
    }

    /** Population-average GLM marginal effect for one linear design column. */
    public static AverageMarginalEffect averageMarginalEffect(
            GlmResult fit, GlmFamily family,
            double[][] design, int predictorColumn) {
        return averageMarginalEffect(model(fit), family, design, null,
            predictorColumn, Averaging.POPULATION_AVERAGE, 0.95);
    }

    /** GLM marginal effect with explicit averaging and link-scale offsets. */
    public static AverageMarginalEffect averageMarginalEffect(
            GlmResult fit, GlmFamily family, double[][] design, double[] offset,
            int predictorColumn, Averaging averaging, double confidenceLevel) {
        return averageMarginalEffect(model(fit), family, design, offset,
            predictorColumn, averaging, confidenceLevel);
    }

    /** Population-average GEE marginal effect for one linear design column. */
    public static AverageMarginalEffect averageMarginalEffect(
            GeeResult fit, GlmFamily family,
            double[][] design, int predictorColumn) {
        return averageMarginalEffect(model(fit), family, design, null,
            predictorColumn, Averaging.POPULATION_AVERAGE, 0.95);
    }

    /** GEE marginal effect with explicit averaging and link-scale offsets. */
    public static AverageMarginalEffect averageMarginalEffect(
            GeeResult fit, GlmFamily family, double[][] design, double[] offset,
            int predictorColumn, Averaging averaging, double confidenceLevel) {
        return averageMarginalEffect(model(fit), family, design, offset,
            predictorColumn, averaging, confidenceLevel);
    }

    private static ExpectedResponse[] expectedResponses(
            Model model, GlmFamily family, double[][] design,
            double[] offset, double confidenceLevel) {
        Design data = validateDesign(design, offset, model.beta().length);
        requireFamily(model, family);
        requireEstimable(model, data.rows());
        double critical = critical(model, confidenceLevel);
        ExpectedResponse[] result = new ExpectedResponse[data.rows().length];
        for (int row = 0; row < result.length; row++) {
            double[] x = data.rows()[row];
            double eta = dot(x, model.beta()) + data.offset()[row];
            double linkSe = standardError(x, model.covariance());
            double mean = family.inverseLink(eta);
            double[] gradient = scale(x, family.meanDerivative(eta));
            double meanSe = standardError(gradient, model.covariance());
            double lower = family.inverseLink(eta - critical * linkSe);
            double upper = family.inverseLink(eta + critical * linkSe);
            result[row] = new ExpectedResponse(eta, linkSe,
                new MeanEstimate(mean, meanSe,
                    Math.min(lower, upper), Math.max(lower, upper)));
        }
        return result;
    }

    private static ScenarioDifference scenarioDifference(
            Model model, GlmFamily family,
            double[][] first, double[][] second,
            double[] firstOffset, double[] secondOffset,
            Averaging averaging, double confidenceLevel) {
        requireFamily(model, family);
        Objects.requireNonNull(averaging, "averaging");
        double critical = critical(model, confidenceLevel);
        Scenario left = scenario(model, family, first, firstOffset, averaging);
        Scenario right = scenario(model, family, second, secondOffset, averaging);
        double covariance = bilinear(left.gradient(), model.covariance(),
            right.gradient());
        double[] gradient = subtract(left.gradient(), right.gradient());
        double estimate = left.mean() - right.mean();
        double standardError = standardError(gradient, model.covariance());
        return new ScenarioDifference(averaging,
            meanEstimate(left, model, critical),
            meanEstimate(right, model, critical), covariance,
            estimate, standardError,
            estimate - critical * standardError,
            estimate + critical * standardError);
    }

    private static RiskRatio riskRatio(
            Model model, GlmFamily family,
            double[][] first, double[][] second,
            double[] firstOffset, double[] secondOffset,
            Averaging averaging, double confidenceLevel) {
        requireFamily(model, family);
        Objects.requireNonNull(averaging, "averaging");
        double critical = critical(model, confidenceLevel);
        Scenario left = scenario(model, family, first, firstOffset, averaging);
        Scenario right = scenario(model, family, second, secondOffset, averaging);
        if (!(left.mean() > 0.0) || !(right.mean() > 0.0)) {
            throw new IllegalArgumentException(
                "risk ratios require positive scenario response means");
        }
        double covariance = bilinear(left.gradient(), model.covariance(),
            right.gradient());
        double[] logGradient = subtract(scale(left.gradient(), 1.0 / left.mean()),
            scale(right.gradient(), 1.0 / right.mean()));
        double logStandardError = standardError(logGradient, model.covariance());
        double logRatio = Math.log(left.mean()) - Math.log(right.mean());
        double ratio = Math.exp(logRatio);
        return new RiskRatio(averaging,
            meanEstimate(left, model, critical),
            meanEstimate(right, model, critical), covariance,
            ratio, ratio * logStandardError, logStandardError,
            Math.exp(logRatio - critical * logStandardError),
            Math.exp(logRatio + critical * logStandardError));
    }

    private static AverageMarginalEffect averageMarginalEffect(
            Model model, GlmFamily family, double[][] design, double[] offset,
            int predictorColumn, Averaging averaging, double confidenceLevel) {
        requireFamily(model, family);
        Objects.requireNonNull(averaging, "averaging");
        Design data = validateDesign(design, offset, model.beta().length);
        if (predictorColumn < 0 || predictorColumn >= model.beta().length) {
            throw new IllegalArgumentException(
                "predictor column must identify a fitted coefficient");
        }
        requireEstimableCoordinate(model, predictorColumn);
        double critical = critical(model, confidenceLevel);
        int parameters = model.beta().length;
        double estimate = 0.0;
        double[] gradient = new double[parameters];
        double coefficient = model.beta()[predictorColumn];
        double denominator;
        double[][] rows;
        double[] offsets;
        if (averaging == Averaging.AT_AVERAGE_COVARIATES) {
            rows = new double[][] {columnMeans(data.rows())};
            offsets = new double[] {mean(data.offset())};
            denominator = 1.0;
        } else {
            rows = data.rows();
            offsets = data.offset();
            denominator = rows.length;
        }
        requireEstimable(model, rows);
        for (int row = 0; row < rows.length; row++) {
            double eta = dot(rows[row], model.beta()) + offsets[row];
            double first = family.meanDerivative(eta);
            double second = family.meanSecondDerivative(eta);
            estimate += first * coefficient / denominator;
            for (int column = 0; column < parameters; column++) {
                double direct = column == predictorColumn ? first : 0.0;
                gradient[column] += (direct
                    + second * rows[row][column] * coefficient) / denominator;
            }
        }
        if (!Double.isFinite(estimate) || !finite(gradient)) {
            throw new IllegalArgumentException(
                "family derivatives produced a non-finite marginal effect");
        }
        double standardError = standardError(gradient, model.covariance());
        return new AverageMarginalEffect(averaging, predictorColumn,
            estimate, standardError,
            estimate - critical * standardError,
            estimate + critical * standardError);
    }

    private static Scenario scenario(
            Model model, GlmFamily family, double[][] design, double[] offset,
            Averaging averaging) {
        Design data = validateDesign(design, offset, model.beta().length);
        int parameters = model.beta().length;
        double[] gradient = new double[parameters];
        if (averaging == Averaging.AT_AVERAGE_COVARIATES) {
            double[] x = columnMeans(data.rows());
            requireEstimable(model, new double[][] {x});
            double eta = dot(x, model.beta()) + mean(data.offset());
            double value = family.inverseLink(eta);
            return new Scenario(value, scale(x, family.meanDerivative(eta)));
        }
        requireEstimable(model, data.rows());
        double value = 0.0;
        for (int row = 0; row < data.rows().length; row++) {
            double eta = dot(data.rows()[row], model.beta()) + data.offset()[row];
            value += family.inverseLink(eta) / data.rows().length;
            double derivative = family.meanDerivative(eta) / data.rows().length;
            for (int column = 0; column < parameters; column++) {
                gradient[column] += derivative * data.rows()[row][column];
            }
        }
        return new Scenario(value, gradient);
    }

    private static MeanEstimate meanEstimate(
            Scenario scenario, Model model, double critical) {
        double standardError = standardError(
            scenario.gradient(), model.covariance());
        return new MeanEstimate(scenario.mean(), standardError,
            scenario.mean() - critical * standardError,
            scenario.mean() + critical * standardError);
    }

    private static Model model(GlmResult fit) {
        Objects.requireNonNull(fit, "fit");
        if (!fit.converged()) {
            throw new IllegalArgumentException(
                "GLM predictions require a converged fit: "
                    + fit.convergenceMessage());
        }
        return new Model(fit.coefficients(), fit.covariance(), fit.family(),
            fit.estimatedDispersion(), fit.residualDegreesOfFreedom(), fit);
    }

    private static Model model(GeeResult fit) {
        Objects.requireNonNull(fit, "fit");
        if (!fit.converged()) {
            throw new IllegalArgumentException(
                "GEE predictions require a converged fit: "
                    + fit.convergenceMessage());
        }
        return new Model(fit.coefficients(), fit.covariance(), fit.family(),
            fit.inference() == GeeInference.CLUSTER_T,
            fit.degreesOfFreedom(), null);
    }

    private static void requireFamily(Model model, GlmFamily family) {
        Objects.requireNonNull(family, "family");
        if (!model.familyName().equals(family.name())) {
            throw new IllegalArgumentException(
                "supplied family " + family.name()
                    + " does not match fitted family " + model.familyName());
        }
    }

    private static void requireEstimable(Model model, double[][] combinations) {
        if (model.glmFit() != null) model.glmFit().requireEstimable(combinations);
    }

    private static void requireEstimableCoordinate(Model model, int column) {
        if (model.glmFit() != null) {
            model.glmFit().requireEstimableCoordinate(column);
        }
    }

    private static Design validateDesign(
            double[][] design, double[] offset, int columns) {
        if (design == null || design.length == 0
                || offset != null && offset.length != design.length) {
            throw new IllegalArgumentException("prediction design is invalid");
        }
        double[][] copy = new double[design.length][];
        for (int row = 0; row < design.length; row++) {
            if (design[row] == null || design[row].length != columns) {
                throw new IllegalArgumentException(
                    "prediction design columns must equal coefficient count");
            }
            copy[row] = MatrixOps.finiteCopy(
                design[row], "prediction design row");
        }
        double[] offsets = offset == null ? new double[design.length]
            : MatrixOps.finiteCopy(offset, "prediction offset");
        return new Design(copy, offsets);
    }

    private static double critical(Model model, double confidenceLevel) {
        if (!(confidenceLevel > 0.0 && confidenceLevel < 1.0)) {
            throw new IllegalArgumentException(
                "confidence level must lie strictly between zero and one");
        }
        double probability = 0.5 + confidenceLevel / 2.0;
        return model.studentT()
            ? T.quantile(probability, model.degreesOfFreedom(), true, false)
            : Normal.quantile(probability, 0.0, 1.0, true, false);
    }

    private static double standardError(double[] gradient, double[] covariance) {
        double variance = bilinear(gradient, covariance, gradient);
        double tolerance = 64.0 * Math.ulp(1.0)
            * Math.max(1.0, maximumAbsolute(covariance)
                * squaredNorm(gradient));
        if (variance < -tolerance || !Double.isFinite(variance)) {
            throw new IllegalArgumentException(
                "coefficient covariance produced invalid estimand variance");
        }
        return Math.sqrt(Math.max(0.0, variance));
    }

    private static double bilinear(
            double[] left, double[] matrix, double[] right) {
        int dimension = left.length;
        double result = 0.0;
        for (int row = 0; row < dimension; row++) {
            for (int column = 0; column < dimension; column++) {
                result += left[row] * matrix[row * dimension + column]
                    * right[column];
            }
        }
        return result;
    }

    private static double dot(double[] left, double[] right) {
        double result = 0.0;
        for (int index = 0; index < left.length; index++) {
            result += left[index] * right[index];
        }
        return result;
    }

    private static double[] scale(double[] values, double scale) {
        double[] result = new double[values.length];
        for (int index = 0; index < values.length; index++) {
            result[index] = scale * values[index];
        }
        return result;
    }

    private static double[] subtract(double[] left, double[] right) {
        double[] result = new double[left.length];
        for (int index = 0; index < left.length; index++) {
            result[index] = left[index] - right[index];
        }
        return result;
    }

    private static double[] columnMeans(double[][] matrix) {
        double[] result = new double[matrix[0].length];
        for (double[] row : matrix) {
            for (int column = 0; column < result.length; column++) {
                result[column] += row[column] / matrix.length;
            }
        }
        return result;
    }

    private static double mean(double[] values) {
        double result = 0.0;
        for (double value : values) result += value / values.length;
        return result;
    }

    private static double maximumAbsolute(double[] values) {
        double result = 0.0;
        for (double value : values) result = Math.max(result, Math.abs(value));
        return result;
    }

    private static double squaredNorm(double[] values) {
        double result = 0.0;
        for (double value : values) result += value * value;
        return result;
    }

    private static boolean finite(double[] values) {
        for (double value : values) if (!Double.isFinite(value)) return false;
        return true;
    }

    private record Model(
        double[] beta, double[] covariance, String familyName,
        boolean studentT, double degreesOfFreedom, GlmResult glmFit) { }
    private record Design(double[][] rows, double[] offset) { }
    private record Scenario(double mean, double[] gradient) { }
}
