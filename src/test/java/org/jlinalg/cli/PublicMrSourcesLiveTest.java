package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import org.jlinalg.study.StudyOperationExecutor.Context;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.condition.EnabledIfEnvironmentVariable;
import org.junit.jupiter.api.io.TempDir;

/** Opt-in smoke test for public indexed sources; excluded from ordinary offline tests. */
@EnabledIfEnvironmentVariable(named = "JLINALG_LIVE_PUBLIC_TEST", matches = "true")
class PublicMrSourcesLiveTest {
    @TempDir Path directory;

    @Test void retrievesPanUkbEffectsAndGrch38LdForRealIrqtlVariants() throws Exception {
        Map<String, Object> study = Map.of(
            "project", Map.of("data_policy", Map.of("remote_transmission", "aggregate_results")),
            "resources", Map.of("downloads", "automatic", "network", "online"));
        Context context = new Context(directory, directory.resolve("output"),
            Map.of("hypothesis", "live", "type", "mendelian_randomization"),
            study, Map.of(), Map.of());
        Files.createDirectories(context.outputDirectory());
        List<PanUkbbGwasOutcome.Instrument> instruments = List.of(
            new PanUkbbGwasOutcome.Instrument("rs2749735", "1", 93158582, "T", "G"),
            new PanUkbbGwasOutcome.Instrument("rs2783506", "1", 93157232, "G", "A"));
        Path outcome = context.outputDirectory().resolve("outcome.tsv");
        Map<String, Object> retrieved = PanUkbbGwasOutcome.extract(context, instruments, outcome);
        assertEquals(2, ((Number) retrieved.get("retained_instruments")).intValue());
        assertEquals(3, Files.readAllLines(outcome).size());

        List<StudyGeneticEvidence.Variant> variants = List.of(
            new StudyGeneticEvidence.Variant("rs2749735", "1", 93158582, "T", "G", -0.019, 0.0018, 1e-25),
            new StudyGeneticEvidence.Variant("rs2783506", "1", 93157232, "G", "A", -0.019, 0.0018, 1e-25));
        var ld = ThousandGenomesGrch38Ld.extract(context, "EUR", variants,
            context.outputDirectory().resolve("ld"));
        assertTrue(Files.isRegularFile(ld.matrix()));
        assertEquals("GRCh38", ld.provenance().get("genome_build"));
        assertTrue(Files.readString(context.outputDirectory().resolve("ld/ld-reference.yaml")).contains("GRCh38"));
    }
}
