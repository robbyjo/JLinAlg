package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;
import org.jlinalg.study.StudyOperationExecutor.Context;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.condition.EnabledIfEnvironmentVariable;
import org.junit.jupiter.api.io.TempDir;

/** Opt-in smoke test for indexed harmonized GWAS Catalog regional retrieval. */
@EnabledIfEnvironmentVariable(named = "JLINALG_LIVE_PUBLIC_TEST", matches = "true")
class PublicColocalizationSourcesLiveTest {
    @TempDir Path directory;

    @Test void retrievesHarmonizedGrch38BmiRegion() throws Exception {
        Map<String, Object> study = Map.of("project", Map.of("data_policy", Map.of("remote_transmission", "aggregate_results")),
            "resources", Map.of("downloads", "automatic", "network", "online"));
        Context context = new Context(directory, directory.resolve("output"), Map.of("hypothesis", "live", "type", "colocalization"),
            study, Map.of(), Map.of());
        Path output = context.outputDirectory().resolve("bmi-region.tsv");
        var result = GwasCatalogRegionalOutcome.extract(context, "GCST90479523", "European", "22", 40250000, 40265000, output);
        assertTrue(result.retainedVariants() > 2);
        assertEquals("GRCh38", result.metadata().get("genome_assembly"));
        assertTrue(Files.readString(output).startsWith("variant_id\tchromosome\tposition"));
        assertTrue(Files.readString(output.resolveSibling("bmi-region.tsv.source.yaml")).contains("a2760408739b078b1be0cba92cc9f3dc"));
    }
}
