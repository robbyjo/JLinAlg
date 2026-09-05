/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.distributional;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.pedigree.PedigreeIndividual;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;
import org.junit.jupiter.api.Test;

/** Independent TMB sparse-GMRF gate for correlated pedigree effects. */
final class ZeroInflatedPedigreeTmbReferenceTest {
    @Test
    void correlatedPedigreeMatchesIndependentTmbTemplate() throws Exception {
        Data data = data();
        Pedigree pedigree = pedigree();
        PedigreeRandomEffectTerm term = PedigreeRandomEffectTerm.of(
            "pedigree", data.ids(), pedigree);
        CorrelatedZeroInflatedRandomEffect correlated =
            CorrelatedZeroInflatedRandomEffect.pedigree("additive", term);
        ZeroInflatedMixedResult fit = SparseZeroInflatedMixedModel.fitPoisson(
            data.response(), data.design(), 2, data.design(), 2,
            List.of(), List.of(), List.of(), List.of(), List.of(correlated),
            null, new ZeroInflatedMixedOptions(1200, 140, 1e-7, 0.3,
                1e-8, 1e4, 1e-5, 1e5, 25.0,
                new double[] {0.4, 0.8}), BackendPolicy.CPU);
        Properties reference = properties();
        assertEquals(value(reference, "count_intercept"),
            fit.countCoefficients()[0], 0.045);
        assertEquals(value(reference, "count_x"),
            fit.countCoefficients()[1], 0.035);
        assertEquals(value(reference, "zero_intercept"),
            fit.zeroCoefficients()[0], 0.08);
        assertEquals(value(reference, "zero_x"),
            fit.zeroCoefficients()[1], 0.07);
        assertEquals(value(reference, "count_sd"),
            Math.sqrt(fit.varianceComponents()[0]), 0.06);
        assertEquals(value(reference, "zero_sd"),
            Math.sqrt(fit.varianceComponents()[1]), 0.09);
        assertEquals(value(reference, "correlation"),
            fit.correlations().get("additive"), 0.10);
        assertEquals(value(reference, "log_likelihood"),
            fit.marginalLogLikelihood(), 0.55);
        assertEquals(2 * pedigree.size(), fit.randomCoefficientCount());
        assertTrue(fit.converged(), fit.convergenceMessage());
    }

    private static Pedigree pedigree() {
        List<PedigreeIndividual> members = new ArrayList<>();
        for (int index = 1; index <= 20; index++)
            members.add(PedigreeIndividual.founder("p" + index));
        for (int index = 21; index <= 60; index++) {
            int family = (index - 21) / 4;
            members.add(new PedigreeIndividual("p" + index,
                "p" + (family + 1), "p" + (family + 11)));
        }
        return Pedigree.of(members);
    }

    private static Properties properties() throws Exception {
        Properties result = new Properties();
        try (InputStream input = ZeroInflatedPedigreeTmbReferenceTest.class
                .getResourceAsStream(
                    "/r-reference/zero-inflated-pedigree-tmb.properties")) {
            result.load(input);
        }
        return result;
    }

    private static double value(Properties properties, String name) {
        return Double.parseDouble(properties.getProperty(name));
    }

    private static Data data() throws Exception {
        List<Double> response = new ArrayList<>();
        List<Double> x = new ArrayList<>();
        List<String> ids = new ArrayList<>();
        try (InputStream input = ZeroInflatedPedigreeTmbReferenceTest.class
                .getResourceAsStream(
                    "/r-reference/zero-inflated-pedigree-tmb-data.tsv");
             BufferedReader reader = new BufferedReader(new InputStreamReader(
                 input, StandardCharsets.UTF_8))) {
            reader.readLine();
            for (String line = reader.readLine(); line != null;
                    line = reader.readLine()) {
                String[] fields = line.split("\t", -1);
                response.add(Double.valueOf(fields[0]));
                x.add(Double.valueOf(fields[1]));
                ids.add(fields[2]);
            }
        }
        double[] design = new double[2 * response.size()];
        for (int row = 0; row < response.size(); row++) {
            design[2 * row] = 1.0;
            design[2 * row + 1] = x.get(row);
        }
        return new Data(toArray(response), design, List.copyOf(ids));
    }

    private static double[] toArray(List<Double> values) {
        double[] result = new double[values.size()];
        for (int index = 0; index < result.length; index++)
            result[index] = values.get(index);
        return result;
    }

    private record Data(
            double[] response, double[] design, List<String> ids) { }
}
