/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import static org.junit.jupiter.api.Assertions.*;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

/** Frozen metafor 5.0-1 / clubSandwich 0.7.0 comparisons; regenerate with resources/meta/generate-metafor.R. */
public class MetaMetaforParityTest {
    private static final BackendPolicy CPU = BackendPolicy.CPU;
    private static final Path ROOT = Path.of("src/test/resources/meta");
    private final Map<String, Double> reference = reference();

    @Test void trimAndFillMatchesMetaforCountsEffectsAndStandardErrors() throws IOException {
        List<MetaStudy> studies = bias();
        for (String method : List.of("FE", "REML", "DL", "PM")) {
            MetaAnalysisOptions options = options(method, MetaInferenceMethod.NORMAL);
            for (MetaBiasCorrections.TrimEstimator estimator : MetaBiasCorrections.TrimEstimator.values())
                for (MetaBiasCorrections.Side side : List.of(MetaBiasCorrections.Side.LEFT, MetaBiasCorrections.Side.RIGHT)) {
                    String key = "trim." + method + "." + estimator + "." + side.name().toLowerCase();
                    if (reference.containsKey(key + ".undefined.0")) {
                        assertThrows(ArithmeticException.class, () -> MetaBiasCorrections.trimAndFill(studies, options, side, estimator, CPU));
                        continue;
                    }
                    var f = MetaBiasCorrections.trimAndFill(studies, options, side, estimator, CPU);
                    close(key + ".k0", new double[] {f.imputedStudyCount()}, 0);
                    close(key + ".beta", new double[] {f.adjustedEffect()}, 2e-7);
                    close(key + ".se", new double[] {f.adjustedFit().standardError()}, 2e-7);
                    assertEquals(studies.size() + f.imputedStudyCount(), f.augmentedStudies().size());
                }
            var auto = MetaBiasCorrections.trimAndFill(studies, options, MetaBiasCorrections.Side.AUTO, MetaBiasCorrections.TrimEstimator.L0, CPU);
            close("trim." + method + ".auto.beta", new double[] {auto.adjustedEffect()}, 2e-7);
            close("trim." + method + ".auto.k0", new double[] {auto.imputedStudyCount()}, 0);
        }
        assertTrue(MetaBiasCorrections.trimAndFill(studies).imputedStudyCount() > 1);
        List<MetaStudy> reflected = studies.stream().map(s -> new MetaStudy(s.name(), -s.effectSize(), s.standardError())).toList();
        assertEquals(-MetaBiasCorrections.trimAndFill(studies).adjustedEffect(),
            MetaBiasCorrections.trimAndFill(reflected).adjustedEffect(), 1e-10);
    }

    @Test void petPeeseAndEggerUseEstimatedDispersionAndResidualT() throws IOException {
        for (String key : List.of("pet", "peese")) {
            var f = key.equals("pet") ? MetaBiasCorrections.pet(bias(), .9) : MetaBiasCorrections.peese(bias(), .9);
            close(key + ".beta", new double[] {f.intercept(), f.slope()}, 1e-11);
            close(key + ".cov", f.covariance(), 1e-11); close(key + ".p", f.pValues(), 1e-11);
            close(key + ".ci", new double[] {f.confidenceLower()[0], f.confidenceUpper()[0], f.confidenceLower()[1], f.confidenceUpper()[1]}, 1e-11);
            assertArrayEquals(new double[] {10, 10}, f.associationStatistics().degreesOfFreedom());
        }
        var egger = MetaPublicationBias.diagnose(bias());
        close("egger", new double[] {egger.eggerStatistic(), egger.eggerPValue(), egger.eggerIntercept()}, 1e-10);
        assertThrows(IllegalArgumentException.class, () -> MetaPublicationBias.diagnose(bias().subList(0, 2)));
    }

    @Test void trimfillWithPositiveHeterogeneityIncludesFiniteQ0Case() throws IOException {
        List<MetaStudy> studies = bias("trim-moderate.csv");
        for (var estimator : MetaBiasCorrections.TrimEstimator.values()) {
            var f = MetaBiasCorrections.trimAndFill(studies,options("REML",MetaInferenceMethod.NORMAL),
                MetaBiasCorrections.Side.LEFT,estimator,CPU);
            close("moderate." + estimator,new double[] {f.imputedStudyCount(),f.adjustedEffect(),
                f.adjustedFit().standardError(),f.adjustedFit().tauSquared()},2e-7);
        }
    }

    @Test void diagonalRandomSlopesMatchMetafor() throws IOException {
        Data d = data();
        var f = MetaMultilevelRegression.fit(d.studies,d.mods,List.of("x"),d.v,true,
            List.of(new MetaRandomEffect("study",d.groups,d.random,MetaRandomEffect.Structure.DIAGONAL)),
            MetaMultilevelRegression.Estimation.REML,options("REML",MetaInferenceMethod.NORMAL),CPU);
        assertTrue(f.converged());
        close("diagonal.REML.G",f.randomCovariances().get(0),2e-5);
        close("diagonal.REML.beta",f.beta(),2e-6);
        close("diagonal.REML.ll",new double[] {f.logLikelihood()},1e-8);
    }

    @Test void structuredCovariancesAndVarianceProfilesAreEstimable() throws IOException {
        Data d=data();
        for(var structure:List.of(MetaRandomEffect.Structure.COMPOUND_SYMMETRY,
                MetaRandomEffect.Structure.AR1)){
            var fit=MetaMultilevelRegression.fit(d.studies,d.mods,List.of("x"),
                d.v,true,List.of(new MetaRandomEffect("study",d.groups,
                    d.random,structure)),MetaMultilevelRegression.Estimation.REML,
                options("REML",MetaInferenceMethod.NORMAL),CPU);
            assertTrue(fit.converged());
            double[] covariance=fit.randomCovariances().get(0);
            assertEquals(covariance[0],covariance[3],1e-12);
            assertTrue(Math.abs(covariance[1])<covariance[0]);
        }
        var interval=MetaMultilevelRegression.profileRandomStandardDeviation(
            d.studies,d.mods,List.of("x"),d.v,true,
            List.of(new MetaRandomEffect("study",d.groups,d.random,
                MetaRandomEffect.Structure.DIAGONAL)),0,0,
            MetaMultilevelRegression.Estimation.REML,.90,2,
            options("REML",MetaInferenceMethod.NORMAL),CPU);
        assertTrue(interval.estimate()>0);
        assertTrue(interval.lowerFound()||interval.lower()==0);
        assertTrue(interval.upperFound());
    }

    @Test void knownCovarianceLikelihoodIncludesGaussianConstant() throws IOException {
        Data d = data();
        var fit = MetaMultilevelRegression.fit(d.studies, d.mods, List.of("x"), d.v, true, CPU);
        close("known.beta", fit.beta(), 1e-11); close("known.cov", fit.covariance(), 1e-11);
        close("known.ll", new double[] {fit.logLikelihood()}, 1e-10);
    }

    @Test void jointRandomCoefficientsAndNestedVariancesMatchMlAndReml() throws IOException {
        Data d = data();
        for (var estimation : MetaMultilevelRegression.Estimation.values()) {
            var fit = hierarchical(d, estimation, false);
            assertTrue(fit.converged(), "hierarchical optimizer failed: " + estimation);
            String key = "hier." + estimation;
            close(key + ".beta", fit.beta(), 2e-6);
            close(key + ".cov", fit.model().covariance(), 2e-6);
            close(key + ".G", fit.randomCovariances().get(0), 2e-5);
            close(key + ".ll", new double[] {fit.logLikelihood()}, 1e-8);
            close(key + ".p", fit.associationStatistics().pValues(), 2e-6);
            var nested = hierarchical(d, estimation, true);
            assertTrue(nested.converged()); key = "nested." + estimation;
            close(key + ".beta", nested.beta(), 2e-6);
            close(key + ".variance", new double[] {nested.randomCovariances().get(0)[0], nested.randomCovariances().get(1)[0]}, 2e-5);
            close(key + ".ll", new double[] {nested.logLikelihood()}, 1e-8);
            if (estimation == MetaMultilevelRegression.Estimation.REML) {
                var cr = MetaClusterRobust.fit(d.studies, d.groups, d.mods, List.of("x"), fit,
                    options("REML", MetaInferenceMethod.STUDENT_T), MetaRobustCorrection.CR2, CPU);
                close("hier.CR2.cov", cr.covariance(), 2e-7);
                close("hier.CR2.df", cr.associationStatistics().degreesOfFreedom(), 2e-5);
                close("hier.CR2.p", cr.associationStatistics().pValues(), 2e-6);
            }
        }
    }

    @Test void allClusterCorrectionsMatchClubSandwich() throws IOException {
        Data d = data();
        for (var correction : MetaRobustCorrection.values()) {
            var f = MetaClusterRobust.fit(d.studies, d.groups, d.mods, List.of("x"), d.v,
                options("FE", MetaInferenceMethod.STUDENT_T), correction, CPU);
            String key = correction == MetaRobustCorrection.CR1P ? "CR1p" : correction.name();
            close(key + ".cov", f.covariance(), 1e-10);
            close(key + ".p", f.associationStatistics().pValues(), 1e-10);
            if (correction == MetaRobustCorrection.CR2) close("CR2.df", f.associationStatistics().degreesOfFreedom(), 1e-9);
            assertTrue(f.confidenceLower()[0] < f.beta()[0]);
        }
    }

    @Test void correlatedEstimatorsAndInferenceAreNotIgnored() throws IOException {
        Data d = data();
        for (String method : List.of("REML", "DL", "PM")) for (var inference : MetaInferenceMethod.values()) {
            var f = MetaCorrelatedAnalysis.fit(d.studies, d.v, options(method, inference), CPU);
            String key = "correlated." + method;
            assertEquals(reference.get(key + ".0"), f.tauSquared(), 2e-6);
            assertEquals(reference.get(key + ".1"), f.pooledEffectSize(), 2e-8);
            double scale = switch (inference) {
                case NORMAL, STUDENT_T -> 1;
                case HARTUNG_KNAPP -> reference.get(key + ".3") / 35;
                case MODIFIED_HARTUNG_KNAPP -> Math.max(1, reference.get(key + ".3") / 35);
            };
            assertEquals(reference.get(key + ".2") * scale, f.covariance()[0], 1e-7);
            assertEquals(reference.get(key + ".4"), f.generalizedQ(), 1e-9);
        }
    }

    @Test void invalidClustersAndCovarianceStructuresAreRejected() throws IOException {
        Data d = data(); String[] single = new String[d.groups.length]; Arrays.fill(single, "one");
        assertThrows(IllegalArgumentException.class, () -> MetaClusterRobust.fit(d.studies, single, d.mods,
            List.of("x"), MetaAnalysisOptions.fixedEffect(), CPU));
        assertThrows(IllegalArgumentException.class, () -> MetaClusterRobust.fit(d.studies, d.groups, d.mods,
            List.of("x"), options("FE", MetaInferenceMethod.HARTUNG_KNAPP), CPU));
        double[][] malformed = new double[d.v.length][];
        assertThrows(IllegalArgumentException.class, () -> MetaCorrelatedAnalysis.fit(d.studies, malformed, MetaAnalysisOptions.fixedEffect(), CPU));
        List<MetaRandomEffect> duplicate = List.of(MetaRandomEffect.intercept("a", d.groups), MetaRandomEffect.intercept("b", d.groups));
        assertThrows(IllegalArgumentException.class, () -> MetaMultilevelRegression.fit(d.studies, d.mods, List.of("x"), d.v, true,
            duplicate, MetaMultilevelRegression.Estimation.REML, MetaAnalysisOptions.randomEffects(), CPU));
        d.v[0][3] = d.v[3][0] = .001;
        assertThrows(IllegalArgumentException.class, () -> MetaClusterRobust.fit(d.studies, d.groups, d.mods,
            List.of("x"), d.v, MetaAnalysisOptions.fixedEffect(), MetaRobustCorrection.CR2, CPU));
    }

    @Test void hedgesSampleProductsAndSumsDoNotOverflow() {
        for (int n : new int[] {50000, Integer.MAX_VALUE}) {
            double df = 2.0 * n - 2, j = 1 - 3 / (4 * df - 1);
            var effect = MetaEffectSizes.standardizedMeanDifference(2, 1, 1, 1, n, n);
            assertEquals(j, effect.effectSize(), 1e-15);
            assertEquals(j * Math.sqrt(2.0 / n + 1 / (2 * df)), effect.standardError(), 1e-15);
        }
    }

    @Test void zeroVarianceBoundaryAndDefaultBudgetDoNotFalselyConverge() throws IOException {
        Data d = data();
        List<MetaStudy> noiseless = new ArrayList<>();
        for (int i = 0; i < d.studies.size(); i++) noiseless.add(new MetaStudy("e" + i, .4 + .3 * d.mods[i][0], d.studies.get(i).standardError()));
        var fit = MetaMultilevelRegression.fit(noiseless, d.mods, List.of("x"), d.v, true,
            List.of(MetaRandomEffect.intercept("study", d.groups)), MetaMultilevelRegression.Estimation.REML,
            MetaAnalysisOptions.randomEffects(), CPU);
        assertTrue(fit.converged()); assertTrue(fit.randomCovariances().get(0)[0] < 1e-10);
        assertArrayEquals(new double[] {.4,.3}, fit.beta(), 1e-10);
        var regular = MetaMultilevelRegression.fit(d.studies, d.mods, List.of("x"), d.v, true,
            List.of(new MetaRandomEffect("study", d.groups, d.random, MetaRandomEffect.Structure.UNSTRUCTURED)),
            MetaMultilevelRegression.Estimation.REML, MetaAnalysisOptions.randomEffects(), CPU);
        assertTrue(regular.converged()); close("hier.REML.G",regular.randomCovariances().get(0),2e-5);
        var exhausted = MetaMultilevelRegression.fit(d.studies, d.mods, List.of("x"), d.v, true,
            List.of(new MetaRandomEffect("study", d.groups, d.random, MetaRandomEffect.Structure.UNSTRUCTURED)),
            MetaMultilevelRegression.Estimation.REML, MetaAnalysisOptions.builder().maximumIterations(10).build(), CPU);
        assertFalse(exhausted.converged(), "an exhausted covariance search must not report success");
    }

    @Test void clusterPermutationAndConfidenceLevelAreRespected() throws IOException {
        Data d = data(); int n = d.studies.size();
        var options = MetaAnalysisOptions.builder().method(MetaAnalysisMethod.FIXED_EFFECT).inferenceMethod(MetaInferenceMethod.STUDENT_T).confidenceLevel(.9).build();
        var original = MetaClusterRobust.fit(d.studies,d.groups,d.mods,List.of("x"),d.v,options,MetaRobustCorrection.CR2,CPU);
        List<MetaStudy> studies = new ArrayList<>(); String[] groups = new String[n]; double[][] mods = new double[n][1], v = new double[n][n];
        for (int i = 0; i < n; i++) {
            int a = n - 1 - i; studies.add(d.studies.get(a)); groups[i] = d.groups[a]; mods[i] = d.mods[a];
            for (int j = 0; j < n; j++) v[i][j] = d.v[a][n-1-j];
        }
        var reordered = MetaClusterRobust.fit(studies,groups,mods,List.of("x"),v,options,MetaRobustCorrection.CR2,CPU);
        assertArrayEquals(original.covariance(),reordered.covariance(),1e-12);
        assertArrayEquals(original.associationStatistics().degreesOfFreedom(),reordered.associationStatistics().degreesOfFreedom(),1e-10);
        var narrow = MetaCorrelatedAnalysis.fit(d.studies,d.v,options,CPU);
        var wide = MetaCorrelatedAnalysis.fit(d.studies,d.v,options("FE",MetaInferenceMethod.STUDENT_T),CPU);
        assertTrue(narrow.confidenceLower()>wide.confidenceLower());
        assertTrue(narrow.confidenceUpper()<wide.confidenceUpper());
    }

    private void close(String key, double[] actual, double tolerance) {
        for (int i = 0; i < actual.length; i++) {
            assertNotNull(reference.get(key + "." + i), "missing fixture " + key);
            assertEquals(reference.get(key + "." + i), actual[i], tolerance, key + "." + i);
        }
    }
    static MetaAnalysisOptions options(String method, MetaInferenceMethod inference) {
        return MetaAnalysisOptions.builder().method(method.equals("FE") ? MetaAnalysisMethod.FIXED_EFFECT : MetaAnalysisMethod.RANDOM_EFFECT)
            .tauSquaredEstimator(switch (method) { case "DL" -> TauSquaredEstimator.DERSIMONIAN_LAIRD; case "PM" -> TauSquaredEstimator.PAULE_MANDEL; default -> TauSquaredEstimator.REML; })
            .inferenceMethod(inference).maximumIterations(1000).tolerance(1e-11).build();
    }
    static MetaHierarchicalResult hierarchical(Data d, MetaMultilevelRegression.Estimation estimation, boolean nested) {
        List<MetaRandomEffect> terms = nested ? List.of(MetaRandomEffect.intercept("study", d.groups), MetaRandomEffect.intercept("effect", d.ids))
            : List.of(new MetaRandomEffect("study", d.groups, d.random, MetaRandomEffect.Structure.UNSTRUCTURED));
        return MetaMultilevelRegression.fit(d.studies, d.mods, List.of("x"), d.v, true, terms,
            estimation, options("REML", MetaInferenceMethod.NORMAL), CPU);
    }
    private static List<MetaStudy> bias() throws IOException {
        return bias("bias.csv");
    }
    private static List<MetaStudy> bias(String filename) throws IOException {
        List<String> lines = Files.readAllLines(ROOT.resolve(filename)); List<MetaStudy> studies = new ArrayList<>();
        for (int i = 1; i < lines.size(); i++) { String[] a = lines.get(i).split(","); studies.add(new MetaStudy("b" + i, Double.parseDouble(a[0]), Double.parseDouble(a[1]))); }
        return studies;
    }
    static Data data() throws IOException {
        List<String> lines = Files.readAllLines(ROOT.resolve("hierarchical.csv")), cov = Files.readAllLines(ROOT.resolve("sampling.csv"));
        int n = lines.size() - 1; List<MetaStudy> studies = new ArrayList<>();
        double[][] mods = new double[n][1], random = new double[n][2], v = new double[n][n]; String[] groups = new String[n], ids = new String[n];
        for (int i = 0; i < n; i++) {
            String[] a = lines.get(i + 1).split(","), b = cov.get(i).split(",");
            studies.add(new MetaStudy("e" + i, Double.parseDouble(a[0]), Math.sqrt(Double.parseDouble(a[1]))));
            mods[i][0] = Double.parseDouble(a[2]); random[i][0] = 1; random[i][1] = mods[i][0]; groups[i] = a[3]; ids[i] = a[4];
            for (int j = 0; j < n; j++) v[i][j] = Double.parseDouble(b[j]);
        }
        return new Data(studies, mods, random, v, groups, ids);
    }
    private static Map<String, Double> reference() {
        try {
            Map<String, Double> map = new HashMap<>(); List<String> lines = Files.readAllLines(ROOT.resolve("metafor-reference.csv"));
            for (String line : lines.subList(1, lines.size())) { String[] a = line.split(","); map.put(a[0], Double.parseDouble(a[1])); }
            return map;
        } catch (IOException e) { throw new java.io.UncheckedIOException(e); }
    }
    record Data(List<MetaStudy> studies, double[][] mods, double[][] random, double[][] v, String[] groups, String[] ids) { }

    /** Isolated runner for shared-workspace development; the same methods run through JUnit in Gradle. */
    public static void main(String[] args) throws Exception {
        int count = 0;
        for (Class<?> type : List.of(MetaMetaforParityTest.class, MetaAnalysisTest.class, MetaAdvancedTest.class, MetaExtensionsTest.class)) {
            Object test = type.getDeclaredConstructor().newInstance();
            for (var method : type.getDeclaredMethods()) if (method.isAnnotationPresent(Test.class)) {
                System.out.println("Checking " + method.getName()); method.invoke(test); count++;
            }
        }
        System.out.println(count + " meta tests passed.");
    }
}
