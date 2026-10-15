package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.jlinalg.study.StudyOperationExecutor.Context;
import org.jlinalg.study.StudyOperationExecutor.State;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class StudyGeneticEvidenceTest {
    @TempDir Path project;

    @Test void runsMrHarmonizationAndSensitivityFromReviewedRegionalSources() throws Exception {
        Context context = context("mendelian_randomization", Map.of());
        var result = new StudyGeneticEvidence().execute(context);
        assertEquals(State.COMPLETE, result.state(), result.details().toString());
        assertTrue(Files.readString(context.outputDirectory().resolve("pair/mr-ivw-fixed.tsv")).contains("0.5"));
        assertTrue(Files.isRegularFile(context.outputDirectory().resolve("pair/mr-egger.tsv")));
        assertTrue(Files.isRegularFile(context.outputDirectory().resolve("pair/mr-loo-0-estimate.tsv")));
        assertTrue(Files.readString(context.outputDirectory().resolve("pair/harmonization.tsv")).contains("rs2\t1\t200\tA\tG\t0.2\t0.1\t-1\t1"));
        assertTrue(Files.readString(context.outputDirectory().resolve("pair/exposure-source.tsv.source.yaml")).contains("sha256"));
    }

    @Test void runsFineMappingAndColocalizationWithPriorSensitivity() throws Exception {
        Context context = context("colocalization", Map.of());
        var result = new StudyGeneticEvidence().execute(context);
        assertEquals(State.COMPLETE, result.state(), result.details().toString());
        assertTrue(Files.isRegularFile(context.outputDirectory().resolve("pair/exposure-fine-mapping.tsv")));
        assertTrue(Files.isRegularFile(context.outputDirectory().resolve("pair/colocalization.tsv")));
        assertTrue(Files.isRegularFile(context.outputDirectory().resolve("pair/colocalization-prior-sensitivity.tsv")));
    }

    @Test void overlapUnapprovedAndWrongFeatureCannotProduceInference() throws Exception {
        for (Map<String, Object> change : List.of(Map.<String, Object>of("sample_overlap", "unknown"),
                Map.<String, Object>of("approved", false), Map.<String, Object>of("feature_id", "wrong-feature"))) {
            Context context = context("mendelian_randomization", change);
            var result = new StudyGeneticEvidence().execute(context);
            assertEquals(State.WAITING_FOR_INPUT, result.state());
            assertFalse(Files.exists(context.outputDirectory().resolve("pair/mr-ivw-fixed.tsv")));
            assertFalse(((List<?>) result.details().get("needs")).isEmpty());
        }
    }

    @Test void mismatchedBuildAndIndefiniteLdStopBeforeInference() throws Exception {
        Context context = context("mendelian_randomization", Map.of("outcome", source("outcome", "GRCh37")));
        var result = new StudyGeneticEvidence().execute(context);
        assertEquals(State.WAITING_FOR_INPUT, result.state());
        context = context("mendelian_randomization", Map.of());
        Files.writeString(project.resolve("ld.tsv"), "variant_id\trs1\trs2\trs3\nrs1\t1\t.9\t.9\nrs2\t.9\t1\t-.9\nrs3\t.9\t-.9\t1\n");
        result = new StudyGeneticEvidence().execute(context);
        assertEquals(State.WAITING_FOR_INPUT, result.state());
        assertTrue(result.details().toString().contains("not positive semidefinite"));
    }

    @Test void genotypeDerivedMrLdCanSkipOnlyTheCubicPsdRecheck() throws Exception {
        Path matrix = project.resolve("indefinite-ld.tsv");
        Files.writeString(matrix, "variant_id\trs1\trs2\trs3\nrs1\t1\t.9\t.9\nrs2\t.9\t1\t-.9\nrs3\t.9\t-.9\t1\n");
        List<StudyGeneticEvidence.Aligned> aligned = List.of(
            aligned("rs1", 100), aligned("rs2", 200), aligned("rs3", 300));
        assertDoesNotThrow(() -> StudyGeneticEvidence.matrix(matrix, aligned, false));
        var failure = assertThrows(java.io.IOException.class,
            () -> StudyGeneticEvidence.matrix(matrix, aligned, true));
        assertTrue(failure.getMessage().contains("not positive semidefinite"));
    }

    @Test void streamedLdReordersRowsAndColumnsAndOrientsTheSelectedSubmatrix()
            throws Exception {
        Path matrix = project.resolve("reordered.csv");
        Files.writeString(matrix, "variant_id,\"rs2\",rs3,rs1\n"
            + "rs3,.2,1,.3\nrs1,.1,.3,1\nrs2,1,.2,.1\n");
        var rs3 = aligned("rs3", 300);
        var flipped = new StudyGeneticEvidence.Aligned(rs3.exposure(),
            rs3.outcome(), 1, -1);
        double[][] actual = StudyGeneticEvidence.matrix(matrix,
            List.of(flipped, aligned("rs1", 100)), true);
        assertArrayEquals(new double[] {1, -.3}, actual[0], 1e-15);
        assertArrayEquals(new double[] {-.3, 1}, actual[1], 1e-15);
    }

    @Test void streamedLdStillRejectsMalformedOrMismatchedRows() throws Exception {
        Path matrix = project.resolve("malformed.tsv");
        for (String rows : List.of("rs1\t1\t0\nrs1\t0\t1\n",
                "rs1\t1\t0\n", "rs1\t1\t0\nrsX\t0\t1\n",
                "rs1\t1\nrs2\t0\t1\n")) {
            Files.writeString(matrix, "variant_id\trs1\trs2\n" + rows);
            assertThrows(java.io.IOException.class, () -> StudyGeneticEvidence.matrix(
                matrix, List.of(aligned("rs1", 100)), false));
        }
    }

    @Test void denseClumpingLdCanExceedFineMappingBudgetWithoutKeepingStrings()
            throws Exception {
        int n = 5001;
        Path matrix = project.resolve("dense-clumping.tsv");
        try (var writer = Files.newBufferedWriter(matrix)) {
            writer.write("variant_id");
            for (int i = 1; i <= n; i++) writer.write("\trs" + i);
            writer.newLine();
            String zeros = "\t0".repeat(n);
            for (int i = 1; i <= n; i++) {
                StringBuilder row = new StringBuilder("rs" + i).append(zeros);
                row.setCharAt(("rs" + i).length() + 2 * i - 1, '1');
                writer.append(row).append('\n');
            }
        }
        double[][] actual = StudyGeneticEvidence.matrix(matrix,
            List.of(aligned("rs5001", 5001)), false);
        assertEquals(1, actual[0][0]);
        var failure = assertThrows(java.io.IOException.class,
            () -> StudyGeneticEvidence.matrix(matrix,
                List.of(aligned("rs5001", 5001)), true));
        assertTrue(failure.getMessage().contains("at most 5000"));
    }

    @Test void palindromeAndWeakInstrumentsAreNotSilentlyUsed() throws Exception {
        Context context = context("mendelian_randomization", Map.of());
        Files.writeString(project.resolve("exposure.tsv"), Files.readString(project.resolve("exposure.tsv")).replace("rs2\t1\t200\tA\tG", "rs2\t1\t200\tA\tT"));
        var result = new StudyGeneticEvidence().execute(context);
        assertEquals(State.COMPLETE, result.state(), result.details().toString());
        assertTrue(result.details().toString().contains("palindromic"));
        assertFalse(Files.readString(context.outputDirectory().resolve("pair/mr-instruments.tsv")).contains("rs2"));
        context = context("mendelian_randomization", Map.of());
        Files.writeString(project.resolve("exposure.tsv"), Files.readString(project.resolve("exposure.tsv")).replace("1e-12", ".1"));
        result = new StudyGeneticEvidence().execute(context);
        assertEquals(State.COMPLETE, result.state());
        assertTrue(result.details().toString().contains("not_estimable"));
        assertTrue(result.details().toString().contains("No independent strong instrument"));
        assertTrue(((List<?>) result.details().get("needs")).isEmpty());
    }

    @Test void oneIndependentStrongInstrumentProducesLimitedWaldRatio() throws Exception {
        Context context = context("mendelian_randomization", Map.of());
        Files.writeString(project.resolve("ld.tsv"), "variant_id\trs1\trs2\trs3\nrs1\t1\t.9\t.9\nrs2\t.9\t1\t.9\nrs3\t.9\t.9\t1\n");
        var result = new StudyGeneticEvidence().execute(context);
        assertEquals(State.COMPLETE, result.state(), result.details().toString());
        Path output = context.outputDirectory().resolve("pair/mr-wald-ratio.tsv");
        assertTrue(Files.isRegularFile(output));
        assertTrue(Files.readString(output).contains("WALD_RATIO\t0.5"));
        assertFalse(Files.exists(context.outputDirectory().resolve("pair/mr-ivw-fixed.tsv")));
        assertTrue(result.details().toString().contains("Only one independent strong instrument remained"));
        assertTrue(result.details().toString().contains("not estimable"));
    }

    @Test void oneAlignedVariantCanProduceWaldRatio() throws Exception {
        Context context = context("mendelian_randomization", Map.of());
        Files.writeString(project.resolve("outcome.tsv"),
            "variant_id\tchromosome\tposition\teffect_allele\tother_allele\tbeta\tse\tp_value\n"
                + "rs1\t1\t100\tA\tG\t0.05\t0.01\t1e-8\n");
        var result = new StudyGeneticEvidence().execute(context);
        assertEquals(State.COMPLETE, result.state(), result.details().toString());
        assertTrue(Files.isRegularFile(context.outputDirectory().resolve("pair/mr-wald-ratio.tsv")));
        assertTrue(result.details().toString().contains("aligned_variants=1"));
    }

    @Test void partialColocalizationFiltersLargeOutcomeToExposureCandidates() throws Exception {
        Map<String, Object> exposure = source("exposure", "GRCh38");
        exposure.put("complete_region", false); exposure.put("export_p_threshold", 1e-7);
        StringBuilder outcome = new StringBuilder("variant_id\tchromosome\tposition\teffect_allele\tother_allele\tbeta\tse\tp_value\n");
        outcome.append("rs1\t1\t100\tA\tG\t0.05\t0.01\t1e-8\nrs2\t1\t200\tG\tA\t-0.1\t0.01\t1e-8\nrs3\t1\t300\tA\tG\t0.15\t0.01\t1e-8\n");
        for (int i = 0; i < 5001; i++) outcome.append("extra").append(i)
            .append("\t1\t400\tA\tG\t0.001\t0.01\t0.5\n");
        Context context = context("colocalization", Map.of("exposure", exposure,
            "partial_colocalization", true, "minimum_exposure_overlap", 0.5));
        Files.writeString(project.resolve("outcome.tsv"), outcome);
        var result = new StudyGeneticEvidence().execute(context);
        assertEquals(State.COMPLETE, result.state(), result.details().toString());
        assertFalse(result.details().toString().contains("5000-variant"));
    }

    @Test void convergedFineMappingWithoutCredibleSignalIsCompletedNotEstimable() throws Exception {
        Context context = context("colocalization", Map.of());
        Files.writeString(project.resolve("outcome.tsv"),
            "variant_id\tchromosome\tposition\teffect_allele\tother_allele\tbeta\tse\tp_value\n"
                + "rs1\t1\t100\tA\tG\t1e-8\t1\t0.99\n"
                + "rs2\t1\t200\tG\tA\t-1e-8\t1\t0.99\n"
                + "rs3\t1\t300\tA\tG\t1e-8\t1\t0.99\n");
        var result = new StudyGeneticEvidence().execute(context);
        assertEquals(State.COMPLETE, result.state(), result.details().toString());
        assertTrue(result.details().toString().contains("not_estimable"));
        assertTrue(result.details().toString().contains("no credible signal"));
        assertTrue(((List<?>) result.details().get("needs")).isEmpty());
    }

    @Test void exhaustedPartialRegionSourcesAreCompletedNotEstimable() throws Exception {
        Map<String, Object> exposure = source("exposure", "GRCh38");
        exposure.put("complete_region", false); exposure.put("export_p_threshold", 1e-7);
        Context context = context("colocalization", Map.of("exposure", exposure,
            "partial_colocalization", true, "minimum_exposure_overlap", 0.5));
        Files.writeString(project.resolve("outcome.tsv"),
            "variant_id\tchromosome\tposition\teffect_allele\tother_allele\tbeta\tse\tp_value\n"
                + "rs1\t1\t100\tA\tG\t0.05\t0.01\t1e-8\n");
        var result = new StudyGeneticEvidence().execute(context);
        assertEquals(State.COMPLETE, result.state(), result.details().toString());
        assertTrue(result.details().toString().contains("not_estimable"));
        assertTrue(result.details().toString().contains("fewer than two variants could be aligned"));
        assertTrue(((List<?>) result.details().get("needs")).isEmpty());
    }

    @Test void studyFineMappingRecordsExtendedIterationCeiling() throws Exception {
        Context context = context("colocalization", Map.of());
        var result = new StudyGeneticEvidence().execute(context);
        assertEquals(State.COMPLETE, result.state(), result.details().toString());
        assertTrue(Files.readString(context.outputDirectory().resolve("pair/exposure-command.log"))
            .contains("--max-iterations, 2000"));
        assertTrue(result.details().toString().contains("maximum_iterations_per_attempt=2000"));
        assertTrue(result.details().toString().contains("requested_maximum_effects=3"));
        assertTrue(result.details().toString().contains("selected_effects=3"));
        assertTrue(result.details().toString().contains("exploratory_largest_convergent_effect_count_fallback=false"));
    }

    @Test void automaticLdNeverGuessesAncestryAndDoesNotLiftBuild() throws Exception {
        var result = new StudyGeneticEvidence().execute(context("fine_mapping", Map.of("ld", Map.of())));
        assertEquals(State.WAITING_FOR_INPUT, result.state());
        assertTrue(result.details().toString().contains("ancestry is not guessed"));
        result = new StudyGeneticEvidence().execute(context("fine_mapping", Map.of("ld", Map.of("population", "EUR"))));
        assertEquals(State.WAITING_FOR_INPUT, result.state());
        assertTrue(result.details().toString().contains("no implicit liftover"));
    }

    @Test void filteredExposureCanSupportReviewedMrButNotColocalization() throws Exception {
        Map<String, Object> exposure = source("exposure", "GRCh38");
        exposure.put("complete_region", false); exposure.put("export_p_threshold", 1e-7);
        var mr = new StudyGeneticEvidence().execute(context("mendelian_randomization", Map.of("exposure", exposure)));
        assertEquals(State.COMPLETE, mr.state(), mr.details().toString());
        var coloc = new StudyGeneticEvidence().execute(context("colocalization", Map.of("exposure", exposure)));
        assertEquals(State.WAITING_FOR_INPUT, coloc.state());
        assertTrue(coloc.details().toString().contains("full regional statistics"));
        var partial = new StudyGeneticEvidence().execute(context("colocalization", Map.of(
            "exposure", exposure, "partial_colocalization", true, "minimum_exposure_overlap", 0.5)));
        assertEquals(State.COMPLETE, partial.state(), partial.details().toString());
        assertTrue(partial.details().toString().contains("exploratory_partial_region_major_limitation"));
        assertTrue(partial.details().toString().contains("absence of colocalization is not evidence of absence"));
        assertTrue(partial.details().toString().contains("exploratory_largest_convergent_effect_count_fallback=true"));
        exposure.put("export_p_threshold", 1e-10);
        mr = new StudyGeneticEvidence().execute(context("mendelian_randomization", Map.of("exposure", exposure)));
        assertEquals(State.WAITING_FOR_INPUT, mr.state());
    }

    @Test void mrReadsIrqtlLogProbabilityAndCompositeRsidWithoutUnderflow() throws Exception {
        Context context = context("mendelian_randomization", Map.of());
        String exposure = Files.readString(project.resolve("exposure.tsv"))
            .replace("p_value", "log10_p")
            .replace("rs1\t1\t100\tA\tG\t0.1\t0.01\t1e-12", "chr1:100:G:A:rs1\t1\t100\tA\tG\t0.1\t0.01\t-12")
            .replace("rs2\t1\t200\tA\tG\t0.2\t0.01\t1e-12", "chr1:200:G:A:rs2\t1\t200\tA\tG\t0.2\t0.01\t-12")
            .replace("rs3\t1\t300\tA\tG\t0.3\t0.01\t1e-12", "chr1:300:G:A:rs3\t1\t300\tA\tG\t0.3\t0.01\t-400");
        Files.writeString(project.resolve("exposure.tsv"), exposure);
        var result = new StudyGeneticEvidence().execute(context);
        assertEquals(State.COMPLETE, result.state(), result.details().toString());
        String instruments = Files.readString(context.outputDirectory().resolve("pair/mr-instruments.tsv"));
        assertTrue(instruments.contains("rs1"));
        assertFalse(instruments.contains("chr1:100"));
    }

    private Map<String, Object> source(String role, String build) {
        Map<String, Object> source = new LinkedHashMap<>();
        source.put("dataset_id", "d1"); source.put("file_role", role); source.put("source_id", "synthetic-" + role);
        source.put("citation", "Synthetic test fixture, not biological evidence"); source.put("genome_build", build);
        source.put("ancestry", "European"); source.put("sample_size", 10000); source.put("complete_region", true);
        source.put("trait", role.equals("outcome") ? "BMI" : "splicing"); source.put("tissue", "blood"); source.put("assay", "sQTL");
        source.put("feature_id", "probe1");
        return source;
    }

    private static StudyGeneticEvidence.Aligned aligned(String id, long position) {
        var variant = new StudyGeneticEvidence.Variant(id, "1", position,
            "A", "G", 0.1, 0.01, 1e-12);
        return new StudyGeneticEvidence.Aligned(variant, variant, 1, 1);
    }
    private Context context(String operation, Map<String, Object> changes) throws Exception {
        String header = "variant_id\tchromosome\tposition\teffect_allele\tother_allele\tbeta\tse\tp_value\n";
        Files.writeString(project.resolve("exposure.tsv"), header + "rs1\t1\t100\tA\tG\t0.1\t0.01\t1e-12\nrs2\t1\t200\tA\tG\t0.2\t0.01\t1e-12\nrs3\t1\t300\tA\tG\t0.3\t0.01\t1e-12\n");
        Files.writeString(project.resolve("outcome.tsv"), header + "rs1\t1\t100\tA\tG\t0.05\t0.01\t1e-8\nrs2\t1\t200\tG\tA\t-0.1\t0.01\t1e-8\nrs3\t1\t300\tA\tG\t0.15\t0.01\t1e-8\n");
        Files.writeString(project.resolve("ld.tsv"), "variant_id\trs1\trs2\trs3\nrs1\t1\t0\t0\nrs2\t0\t1\t0\nrs3\t0\t0\t1\n");
        Files.writeString(project.resolve("alleles.tsv"), "variant_id\teffect_allele\tother_allele\nrs1\tA\tG\nrs2\tA\tG\nrs3\tA\tG\n");
        Map<String, Object> ld = source("ld", "GRCh38"); ld.put("alleles_file_role", "alleles");
        Map<String, Object> pair = new LinkedHashMap<>();
        pair.put("id", "pair"); pair.put("feature_id", "probe1"); pair.put("approved", true); pair.put("rationale", "Synthetic test");
        pair.put("sample_overlap", "none"); pair.put("sample_overlap_basis", "Separate simulated samples"); pair.put("instrument_validity_rationale", "Simulation only");
        pair.put("region", Map.of("chromosome", "1", "start", 1, "end", 500));
        pair.put("exposure", source("exposure", "GRCh38")); pair.put("outcome", source("outcome", "GRCh38")); pair.put("ld", ld); pair.putAll(changes);
        Map<String, Object> study = Map.of("resources", Map.of("network", "offline"),
            "datasets", List.of(Map.of("id", "d1", "files", Map.of("exposure", "exposure.tsv", "outcome", "outcome.tsv", "ld", "ld.tsv", "alleles", "alleles.tsv"))),
            "hypotheses", List.of(Map.of("id", "h1", "outcome", "BMI", "analysis", Map.of("bioinformatics", Map.of("genetic_pairs", List.of(pair))))));
        Path upstream = Files.createTempDirectory(project, "upstream-");
        Files.createDirectories(upstream.resolve("standard-results"));
        Files.writeString(upstream.resolve("manifest.yaml"), "operation_type: validate_harmonize\nstate: complete\n");
        Files.writeString(upstream.resolve("standard-results/d1.tsv"), "feature_id\tp_value\nprobe1\t1e-9\n");
        Files.writeString(upstream.resolve("feature-annotation.tsv"), "feature_id\tgene_symbol\tgene_id\tchromosome\tstart\tend\tgenome_build\tmapping_status\nprobe1\tTEST\tENSG1\t1\t100\t300\tGRCh38\tlocated\n");
        return new Context(project, Files.createTempDirectory(project, "output-"), Map.of("hypothesis", "h1", "type", operation, "datasets", List.of("d1")), study, Map.of(), Map.of("harmonize", upstream));
    }
}
