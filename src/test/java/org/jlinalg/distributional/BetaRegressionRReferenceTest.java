/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

/** Comparisons with the GasolineYield and FoodExpenditure betareg examples. */
final class BetaRegressionRReferenceTest {
    private static final double[] GASOLINE_COEFFICIENTS = {
        -6.159571047010208, 1.727728875069054, 1.322596915621213,
        1.572309886594047, 1.059714112762823, 1.133751781080070,
        1.040161812360569, 0.543692226082664, 0.495900661509114,
        0.385792958032884, 0.010966874176797
    };
    private static final double[] GASOLINE_STANDARD_ERRORS = {
        0.182324675698894, 0.101229390444847, 0.117902041897632,
        0.116104500626117, 0.102359826066170, 0.103523238546190,
        0.106036474164913, 0.109127466666666, 0.108925669305345,
        0.118593267790192, 0.000412647504395
    };

    @Test
    void gasolineBatchAndTemperatureMatchesBetaregTableOne()
            throws IOException {
        double[][] rows = read("/r-reference/betareg-gasoline.tsv", 6);
        double[] response = column(rows, 0);
        double[][] meanDesign = gasolineDesign(rows);

        BetaRegressionResult fit = BetaRegression.fit(response, meanDesign,
            intercept(rows.length),
            BetaRegressionOptions.constantPrecisionDefaults(),
            BackendPolicy.CPU);

        assertTrue(fit.converged(), fit.convergenceMessage());
        assertArrayClose(GASOLINE_COEFFICIENTS,
            fit.meanCoefficients(), 2e-8);
        assertArrayClose(GASOLINE_STANDARD_ERRORS,
            meanStandardErrors(fit, meanDesign[0].length), 2e-8);
        assertEquals(440.2783885611994,
            fit.precisionCoefficients()[0], 2e-6);
        assertEquals(110.0256249613124,
            fit.standardErrors()[meanDesign[0].length], 2e-6);
        assertEquals(84.79755796201518, fit.logLikelihood(), 2e-9);
        assertTrue(fit.iterations() < 20,
            "specialized scoring should need fewer than betareg's 54 updates");
    }

    @Test
    void variablePrecisionGasolineMatchesBetaregExample()
            throws IOException {
        double[][] rows = read("/r-reference/betareg-gasoline.tsv", 6);
        double[] response = column(rows, 0);
        double[][] precisionDesign = new double[rows.length][2];
        for (int row = 0; row < rows.length; row++) {
            precisionDesign[row][0] = 1.0;
            precisionDesign[row][1] = rows[row][4];
        }

        BetaRegressionResult fit = BetaRegression.fit(response,
            gasolineDesign(rows), precisionDesign,
            BetaRegressionOptions.variablePrecisionDefaults(),
            BackendPolicy.CPU);

        assertTrue(fit.converged(), fit.convergenceMessage());
        assertArrayClose(new double[] {
            -5.923236136239194, 1.601987749622551, 1.297266254630967,
            1.565338274749211, 1.030071970491647, 1.154163041862477,
            1.019444647902441, 0.622259050245802, 0.564582997139931,
            0.359438984247131, 0.010359481692118
        }, fit.meanCoefficients(), 2e-8);
        assertEquals(1.364088821346420,
            fit.precisionCoefficients()[0], 1e-7);
        assertEquals(0.014570318309691,
            fit.precisionCoefficients()[1], 2e-10);
        int start = fit.meanCoefficients().length;
        assertEquals(1.225781237264785,
            fit.standardErrors()[start], 1e-7);
        assertEquals(0.003618284536182,
            fit.standardErrors()[start + 1], 2e-10);
        assertEquals(86.97706518350222, fit.logLikelihood(), 2e-9);
    }

    @Test
    void householdFoodExpenditureMatchesBetaregTableTwo()
            throws IOException {
        double[][] rows = read("/r-reference/betareg-food.tsv", 3);
        double[] response = new double[rows.length];
        double[][] design = new double[rows.length][3];
        for (int row = 0; row < rows.length; row++) {
            response[row] = rows[row][0] / rows[row][1];
            design[row][0] = 1.0;
            design[row][1] = rows[row][1];
            design[row][2] = rows[row][2];
        }

        BetaRegressionResult fit = BetaRegression.fit(response, design,
            intercept(rows.length),
            BetaRegressionOptions.constantPrecisionDefaults(),
            BackendPolicy.CPU);

        assertTrue(fit.converged(), fit.convergenceMessage());
        assertArrayClose(new double[] {
            -0.622548056189129, -0.012298840534926, 0.118462097692880
        }, fit.meanCoefficients(), 2e-8);
        assertArrayClose(new double[] {
            0.223853539305629, 0.003035584648663, 0.035340667012308
        }, meanStandardErrors(fit, 3), 2e-8);
        assertEquals(35.60975032891577,
            fit.precisionCoefficients()[0], 2e-7);
        assertEquals(8.079598247959055,
            fit.standardErrors()[3], 2e-7);
        assertEquals(45.333509321238, fit.logLikelihood(), 2e-9);
    }

    @Test
    void validatesOpenIntervalAndSupportsRowMajorInput() {
        double[] response = {0.2, 0.3, 0.4, 0.5};
        double[] mean = {1.0, 0.0, 1.0, 1.0, 1.0, 2.0, 1.0, 3.0};
        double[] precision = {1.0, 1.0, 1.0, 1.0};
        BetaRegressionResult fit = BetaRegression.fit(response,
            mean, 2, precision, 1,
            BetaRegressionOptions.constantPrecisionDefaults(),
            BackendPolicy.CPU);
        assertTrue(fit.converged(), fit.convergenceMessage());
        assertEquals(4, fit.observations());
        assertThrows(IllegalArgumentException.class,
            () -> BetaRegression.fit(new double[] {0.0, 0.5},
                new double[][] {{1.0}, {1.0}}));
    }

    private static double[][] gasolineDesign(double[][] rows) {
        double[][] design = new double[rows.length][11];
        for (int row = 0; row < rows.length; row++) {
            design[row][0] = 1.0;
            int batch = (int) rows[row][5];
            if (batch < 10) design[row][batch] = 1.0;
            design[row][10] = rows[row][4];
        }
        return design;
    }

    private static double[][] intercept(int rows) {
        double[][] result = new double[rows][1];
        for (double[] row : result) row[0] = 1.0;
        return result;
    }

    private static double[] column(double[][] rows, int column) {
        double[] result = new double[rows.length];
        for (int row = 0; row < rows.length; row++) {
            result[row] = rows[row][column];
        }
        return result;
    }

    private static double[] meanStandardErrors(
            BetaRegressionResult fit, int columns) {
        double[] result = new double[columns];
        System.arraycopy(fit.standardErrors(), 0, result, 0, columns);
        return result;
    }

    private static void assertArrayClose(
            double[] expected, double[] actual, double tolerance) {
        assertEquals(expected.length, actual.length);
        for (int index = 0; index < expected.length; index++) {
            assertEquals(expected[index], actual[index], tolerance,
                "coefficient " + index);
        }
    }

    private static double[][] read(String resource, int columns)
            throws IOException {
        InputStream input = BetaRegressionRReferenceTest.class
            .getResourceAsStream(resource);
        if (input == null) throw new IOException("missing " + resource);
        List<double[]> rows = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(
                input, StandardCharsets.UTF_8))) {
            String line = reader.readLine();
            while ((line = reader.readLine()) != null) {
                String[] fields = line.split("\\t");
                double[] values = new double[columns];
                for (int column = 0; column < columns; column++) {
                    values[column] = Double.parseDouble(fields[column]);
                }
                rows.add(values);
            }
        }
        return rows.toArray(double[][]::new);
    }
}
