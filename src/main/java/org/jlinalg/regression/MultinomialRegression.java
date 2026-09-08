/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

/** Baseline-category multinomial logistic regression by monotone gradient ascent. */
public final class MultinomialRegression {
    private MultinomialRegression() { }

    public static MultinomialRegressionResult fit(int[] response,
            double[][] predictors, int classCount, MultinomialOptions options) {
        if (response == null || predictors == null || response.length < 2
                || predictors.length != response.length || classCount < 2
                || options == null || predictors[0] == null)
            throw new IllegalArgumentException("multinomial inputs are invalid");
        int rows = response.length, columns = predictors[0].length;
        if(columns < 1) throw new IllegalArgumentException("at least one predictor is required");
        double[] x = new double[rows * columns];
        for (int row = 0; row < rows; row++) {
            if (response[row] < 0 || response[row] >= classCount
                    || predictors[row] == null || predictors[row].length != columns)
                throw new IllegalArgumentException("response classes or predictor rows are invalid");
            for (int column = 0; column < columns; column++) {
                x[row * columns + column] = predictors[row][column];
                if (!Double.isFinite(x[row * columns + column])) throw new IllegalArgumentException("predictors must be finite");
            }
        }
        double[] scales=RegressionOptimizer.scaleColumns(x,rows,columns);
        var optimized=RegressionOptimizer.minimize(point -> {
            double[] gradient=gradient(response,x,point,rows,columns,classCount);
            for(int j=0;j<gradient.length;j++)gradient[j]/=-rows;
            return new RegressionOptimizer.Evaluation(-logLikelihood(response,x,point,rows,columns,classCount)/rows,gradient);
        },(classCount-1)*columns,options.maximumIterations(),options.relativeTolerance(),options.initialStep());
        double[] beta=optimized.point();
        double likelihood=logLikelihood(response,x,beta,rows,columns,classCount);
        double[] probabilities=probabilities(x,beta,rows,columns,classCount);
        for(int j=0;j<beta.length;j++) beta[j]/=scales[j%columns];
        return new MultinomialRegressionResult(beta, probabilities,
            likelihood, rows, columns, classCount, optimized.iterations(), optimized.converged());
    }

    public static MultinomialRegressionResult fit(int[] response, double[][] predictors, int classCount) {
        return fit(response, predictors, classCount, MultinomialOptions.defaults());
    }

    private static double[] gradient(int[] y, double[] x, double[] beta, int rows, int columns, int classes) {
        double[] result = new double[beta.length], probabilities = probabilities(x, beta, rows, columns, classes);
        for (int row = 0; row < rows; row++) for (int cls = 1; cls < classes; cls++) {
            double residual = (y[row] == cls ? 1.0 : 0.0) - probabilities[row * classes + cls];
            for (int column = 0; column < columns; column++) result[(cls - 1) * columns + column] += x[row * columns + column] * residual;
        }
        return result;
    }
    private static double logLikelihood(int[] y, double[] x, double[] beta, int rows, int columns, int classes) {
        double result=0; double[] logits=new double[classes];
        for(int row=0;row<rows;row++) {
            double maximum=0;
            for(int cls=1;cls<classes;cls++) {
                logits[cls]=0;
                for(int column=0;column<columns;column++)logits[cls]+=x[row*columns+column]*beta[(cls-1)*columns+column];
                maximum=Math.max(maximum,logits[cls]);
            }
            double sum=0;for(double logit:logits)sum+=Math.exp(logit-maximum);
            result+=logits[y[row]]-maximum-Math.log(sum);
        }
        return result;
    }
    private static double[] probabilities(double[] x, double[] beta, int rows, int columns, int classes) {
        double[] result = new double[rows * classes];
        for (int row = 0; row < rows; row++) {
            double maximum = 0.0;
            for (int cls = 1; cls < classes; cls++) { double value = 0; for (int column = 0; column < columns; column++) value += x[row * columns + column] * beta[(cls - 1) * columns + column]; maximum = Math.max(maximum, value); }
            double denominator = Math.exp(-maximum);
            for (int cls = 1; cls < classes; cls++) { double value = 0; for (int column = 0; column < columns; column++) value += x[row * columns + column] * beta[(cls - 1) * columns + column]; result[row * classes + cls] = Math.exp(value - maximum); denominator += result[row * classes + cls]; }
            result[row * classes] = Math.exp(-maximum) / denominator;
            for (int cls = 1; cls < classes; cls++) result[row * classes + cls] /= denominator;
        }
        return result;
    }
}
