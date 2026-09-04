/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.sem.Sem;
import org.jlinalg.sem.SemFitResult;
import org.jlinalg.sem.SemModel;
import org.jlinalg.sem.SemOptions;

/** TOPMed cardiometabolic observed-variable SEM benchmark shared with lavaan. */
public final class TopmedSemBenchmark {
    private static final String[] VARIABLES = {
        "Sex", "Age", "BMI", "Waist", "Systolic_BP", "Diastolic_BP",
        "Glucose", "HDL", "LnTG", "LogInsulin", "CRP", "eGFR"
    };
    private static volatile double checksum;

    private TopmedSemBenchmark() { }

    public static void main(String[] arguments) throws IOException {
        Path input = Path.of(stringProperty("jlinalg.benchmark.input",
            "D:/Research/topmed/splicing-bmi/new/mastermat-batch1234-wcbc-ext.csv"));
        int measurements = integerProperty("jlinalg.benchmark.measurements", 10);
        double[][] data = standardizedCompleteRows(input);
        SemModel model = model();
        for (int iteration = 0; iteration < 3; iteration++)
            checksum += Sem.fit(data, model).logLikelihood();
        double[] seconds = new double[measurements];
        SemFitResult fit = null;
        for (int iteration = 0; iteration < measurements; iteration++) {
            long started = System.nanoTime();
            fit = Sem.fit(data, model, SemOptions.defaults(), BackendPolicy.CPU);
            seconds[iteration] = (System.nanoTime() - started) / 1e9;
            checksum += fit.logLikelihood();
        }
        java.util.Arrays.sort(seconds);
        double median = seconds[seconds.length / 2];
        System.out.println("runtime,observations,variables,parameters,evaluations,"
            + "converged,median_seconds,log_likelihood,chi_square,df,cfi,tli,rmsea,srmr");
        System.out.printf(Locale.ROOT,
            "JLinAlg,%d,%d,%d,%d,%s,%.9f,%.12g,%.12g,%d,%.12g,%.12g,%.12g,%.12g%n",
            fit.observations(), VARIABLES.length, fit.parameters().size(),
            fit.functionEvaluations(), fit.converged(), median, fit.logLikelihood(),
            fit.chiSquare(), fit.degreesOfFreedom(), fit.cfi(), fit.tli(),
            fit.rmsea(), fit.srmr());
        System.out.println("label,estimate,se");
        fit.parameters().forEach(parameter -> System.out.printf(Locale.ROOT,
            "%s,%.12g,%.12g%n", parameter.label(), parameter.estimate(),
            parameter.standardError()));
    }

    static SemModel model() {
        SemModel.Builder model = SemModel.builder(VARIABLES)
            .regression("BMI", "Age", 0.1)
            .regression("BMI", "Sex", 0.1)
            .regression("Waist", "BMI", 0.5)
            .regression("Waist", "Sex", 0.1)
            .regression("Systolic_BP", "Age", 0.2)
            .regression("Systolic_BP", "BMI", 0.1)
            .regression("Diastolic_BP", "Age", 0.1)
            .regression("Diastolic_BP", "BMI", 0.1)
            .regression("Glucose", "Age", 0.1)
            .regression("Glucose", "BMI", 0.1)
            .regression("HDL", "BMI", -0.1)
            .regression("HDL", "LnTG", -0.2)
            .regression("HDL", "Sex", 0.1)
            .regression("LnTG", "BMI", 0.2)
            .regression("LnTG", "Glucose", 0.1)
            .regression("LnTG", "Sex", -0.1)
            .regression("LogInsulin", "BMI", 0.2)
            .regression("LogInsulin", "Glucose", 0.2)
            .regression("CRP", "BMI", 0.2)
            .regression("CRP", "Sex", 0.1)
            .regression("eGFR", "Age", -0.4)
            .regression("eGFR", "Sex", 0.1)
            .regression("eGFR", "BMI", 0.1)
            .covariance("Age", "Sex", 0.0);
        for (String variable : VARIABLES) model.variance(variable, 1.0);
        return model.build();
    }

    private static double[][] standardizedCompleteRows(Path input) throws IOException {
        List<double[]> rows = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(input)) {
            String headerLine = reader.readLine();
            if (headerLine == null) throw new IllegalArgumentException("empty CSV");
            String[] header = headerLine.split(",", -1);
            Map<String, Integer> positions = new HashMap<>();
            for (int index = 0; index < header.length; index++)
                positions.put(header[index], index);
            int[] selected = new int[VARIABLES.length];
            for (int index = 0; index < VARIABLES.length; index++) {
                Integer position = positions.get(VARIABLES[index]);
                if (position == null)
                    throw new IllegalArgumentException("missing CSV column: " + VARIABLES[index]);
                selected[index] = position;
            }
            String line;
            while ((line = reader.readLine()) != null) {
                String[] fields = line.split(",", -1);
                double[] row = new double[VARIABLES.length];
                boolean complete = true;
                for (int index = 0; index < selected.length; index++) {
                    String field = selected[index] < fields.length ? fields[selected[index]] : "";
                    if (field.isBlank() || field.equalsIgnoreCase("NA")) {
                        complete = false;
                        break;
                    }
                    try {
                        row[index] = Double.parseDouble(field);
                    } catch (NumberFormatException exception) {
                        complete = false;
                        break;
                    }
                    complete &= Double.isFinite(row[index]);
                }
                if (complete) rows.add(row);
            }
        }
        if (rows.size() <= VARIABLES.length)
            throw new IllegalArgumentException("too few complete rows");
        double[][] result = rows.toArray(double[][]::new);
        for (int variable = 0; variable < VARIABLES.length; variable++) {
            double mean = 0.0;
            for (double[] row : result) mean += row[variable];
            mean /= result.length;
            double sumSquares = 0.0;
            for (double[] row : result) {
                double centered = row[variable] - mean;
                sumSquares += centered * centered;
            }
            double standardDeviation = Math.sqrt(sumSquares / (result.length - 1.0));
            if (!(standardDeviation > 0.0))
                throw new IllegalArgumentException("constant CSV column: " + VARIABLES[variable]);
            for (double[] row : result)
                row[variable] = (row[variable] - mean) / standardDeviation;
        }
        return result;
    }

    private static String stringProperty(String name, String defaultValue) {
        String value = System.getProperty(name, defaultValue);
        if (value.isBlank()) throw new IllegalArgumentException(name + " must be nonblank");
        return value;
    }

    private static int integerProperty(String name, int defaultValue) {
        int value = Integer.getInteger(name, defaultValue);
        if (value < 1) throw new IllegalArgumentException(name + " must be positive");
        return value;
    }
}
