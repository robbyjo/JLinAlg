/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.distributional;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
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
import java.util.Properties;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glmm.GlmmLaplaceOptions;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.pedigree.PedigreeIndividual;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;
import org.junit.jupiter.api.Test;

/** Independent comparisons with glmmTMB 1.1.14 beta-family Laplace fits. */
final class BetaMixedModelRReferenceTest {
    @Test
    void groupedRandomInterceptTracksGlmmTmb() throws IOException {
        Fixture fixture = fixture();
        BetaMixedModelResult fit = BetaMixedModel.fit(fixture.response,
            fixture.fixed, fixture.response.length, 2,
            List.of(RandomEffectTerm.randomIntercept("group", fixture.groups)),
            null, null, controls(), BackendPolicy.CPU);

        assertTrue(fit.converged(), fit.convergenceMessage());
        assertEquals(fixture.reference("beta_0"),
            fit.meanCoefficients()[0], 0.003);
        assertEquals(fixture.reference("beta_1"),
            fit.meanCoefficients()[1], 0.006);
        assertEquals(fixture.reference("precision"), fit.precision(), 0.01);
        assertEquals(fixture.reference("variance"),
            fit.varianceComponents()[0], 0.002);
        assertEquals(fixture.reference("log_likelihood"),
            fit.marginalLogLikelihood(), 0.13);
        assertEquals(fixture.groups.stream().distinct().count(),
            fit.randomEffects("group").length);
    }

    @Test
    void unrelatedFounderPedigreeEqualsGroupedRandomIntercept()
            throws IOException {
        Fixture fixture = fixture();
        List<PedigreeIndividual> founders = fixture.groups.stream().distinct()
            .map(PedigreeIndividual::founder).toList();
        PedigreeRandomEffectTerm pedigree = PedigreeRandomEffectTerm.of(
            "animal", fixture.groups, Pedigree.of(founders));

        BetaMixedModelResult grouped = BetaMixedModel.fit(fixture.response,
            fixture.fixed, fixture.response.length, 2,
            List.of(RandomEffectTerm.randomIntercept("animal", fixture.groups)),
            null, null, controls(), BackendPolicy.CPU);
        BetaMixedModelResult genetic = BetaMixedModel.fitPedigree(
            fixture.response, fixture.fixed, fixture.response.length, 2,
            pedigree, null, null, controls(), BackendPolicy.CPU);

        assertArrayEquals(grouped.meanCoefficients(),
            genetic.meanCoefficients(), 2e-8);
        assertArrayEquals(grouped.varianceComponents(),
            genetic.varianceComponents(), 2e-8);
        assertArrayEquals(grouped.randomEffects("animal"),
            genetic.randomEffects("animal"), 2e-8);
        assertEquals(grouped.precision(), genetic.precision(), 2e-8);
    }

    @Test
    void pedigreeRetainsUnobservedAncestorsAndRejectsBoundaryData() {
        Pedigree pedigree = Pedigree.of(List.of(
            PedigreeIndividual.founder("sire"),
            PedigreeIndividual.founder("dam"),
            new PedigreeIndividual("offspring", "sire", "dam")));
        List<String> animals = List.of("offspring", "offspring", "offspring",
            "offspring", "offspring", "offspring");
        PedigreeRandomEffectTerm term = PedigreeRandomEffectTerm.of(
            "animal", animals, pedigree);
        double[] fixed = {1, -1, 1, -.6, 1, -.2, 1, .2, 1, .6, 1, 1};
        BetaMixedModelResult fit = BetaMixedModel.fitPedigree(
            new double[] {.15, .24, .31, .43, .57, .70}, fixed, 6, 2,
            term, null, null,
            new BetaMixedModelOptions(new GlmmLaplaceOptions(
                8, 60, 1e-6, 1.0, 1e-8, 100.0, null),
            10.0, 1e-6, 1e8), BackendPolicy.CPU);
        assertEquals(3, fit.randomEffects("animal").length);
        assertThrows(IllegalArgumentException.class,
            () -> BetaMixedModel.fit(new double[] {0.0, 0.5},
                new double[][] {{1.0}, {1.0}},
                List.of(RandomEffectTerm.randomIntercept("g",
                    List.of("a", "b")))));
    }

    private static BetaMixedModelOptions controls() {
        return new BetaMixedModelOptions(new GlmmLaplaceOptions(
            35, 100, 1e-7, 0.75, 1e-8, 100.0, null),
            10.0, 1e-6, 1e8);
    }

    private static Fixture fixture() throws IOException {
        InputStream input = BetaMixedModelRReferenceTest.class
            .getResourceAsStream("/r-reference/beta-mixed-glmmtmb.tsv");
        if (input == null) throw new IOException("missing beta mixed fixture");
        List<Double> response = new ArrayList<>();
        List<Double> fixed = new ArrayList<>();
        List<String> groups = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(
                input, StandardCharsets.UTF_8))) {
            String line = reader.readLine();
            while ((line = reader.readLine()) != null) {
                String[] fields = line.split("\t");
                response.add(Double.valueOf(fields[0]));
                fixed.add(1.0);
                fixed.add(Double.valueOf(fields[1]));
                groups.add(fields[2]);
            }
        }
        Properties reference = new Properties();
        try (InputStream properties = BetaMixedModelRReferenceTest.class
                .getResourceAsStream(
                    "/r-reference/beta-mixed-glmmtmb.properties")) {
            if (properties == null) throw new IOException("missing properties");
            reference.load(properties);
        }
        return new Fixture(response.stream().mapToDouble(Double::doubleValue)
            .toArray(), fixed.stream().mapToDouble(Double::doubleValue).toArray(),
            groups, reference);
    }

    private record Fixture(double[] response, double[] fixed,
            List<String> groups, Properties properties) {
        double reference(String name) {
            return Double.parseDouble(properties.getProperty(name));
        }
    }
}
