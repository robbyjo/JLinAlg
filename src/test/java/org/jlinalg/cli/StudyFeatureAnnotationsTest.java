package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.jlinalg.study.StudyOperationExecutor.Context;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class StudyFeatureAnnotationsTest {
    @TempDir Path project;

    @Test void resolvesLocalProbesConvertsBedIntervalsAndRetainsAmbiguity() throws Exception {
        Files.writeString(project.resolve("annotation.tsv"), """
            feature_id\tgene_id\tgene_symbol\tchr\tstart\tend
            cg1\tENSG1\tGENE1\t1\t99\t100
            cg2\tENSG2\tGENE2\t2\t199\t205
            cg2\tENSG3\tGENE3\t3\t299\t305
            """);
        Path out = assess("cg1\ncg2\ncg3\n", Map.of("files", Map.of("annotation", "annotation.tsv"),
            "genome_build", "hg38", "coordinate_system", "0-based-half-open"));
        String table = Files.readString(out.resolve("feature-annotation.tsv"));
        assertTrue(table.contains("cg1\tENSG1\tGENE1\t1\t100\t100\tGRCh38"));
        assertEquals(2, table.lines().filter(line -> line.startsWith("cg2\t") && line.endsWith("ambiguous")).count());
        assertTrue(table.lines().anyMatch(line -> line.startsWith("cg3\t\t\t") && line.endsWith("unmapped")));
        assertTrue(Files.readString(out.resolve("annotation-summary.yaml")).contains("input_requests"));
    }

    @Test void cannotUseLocationsWithoutBuildOrCoordinateConvention() throws Exception {
        Files.writeString(project.resolve("annotation.tsv"), "feature_id\tgene_symbol\tchr\tposition\np1\tG1\t1\t100\n");
        Path out = assess("p1\n", Map.of("files", Map.of("annotation", "annotation.tsv")));
        assertTrue(Files.readString(out.resolve("feature-annotation.tsv")).contains("missing_location"));
        assertFalse(Files.readString(out.resolve("feature-annotation.tsv")).contains("located"));
    }

    @Test void gtfMatchesExactTranscriptVersionsAndRejectsMalformedCoordinates() throws Exception {
        Files.writeString(project.resolve("annotation.gtf"), """
            chr1\tGENCODE\ttranscript\t10\t20\t.\t+\t.\tgene_id "ENSG00001.1"; transcript_id "ENST00001.2"; gene_name "G1";
            """);
        Path out = assess("ENST00001.2\nENST00001.1\nENST00001\n", Map.of(
            "files", Map.of("annotation", "annotation.gtf"), "genome_build", "GRCh38"));
        String table = Files.readString(out.resolve("feature-annotation.tsv"));
        assertTrue(table.lines().anyMatch(line -> line.startsWith("ENST00001.2\t") && line.endsWith("located")));
        assertTrue(table.lines().anyMatch(line -> line.startsWith("ENST00001\t") && line.endsWith("located")));
        assertTrue(table.lines().anyMatch(line -> line.startsWith("ENST00001.1\t") && line.endsWith("unmapped")));
        Files.writeString(project.resolve("invalid.tsv"), "feature_id\tchr\tstart\tend\np1\t1\t20\t10\n");
        assertThrows(java.io.IOException.class, () -> assess("p1\n", Map.of("files", Map.of("annotation", "invalid.tsv"),
            "genome_build", "GRCh38", "coordinate_system", "1-based-closed")));
    }

    @Test void unknownProbesAreNeverAssumedToBeSymbols() throws Exception {
        Path out = assess("cg123\n123_at\n", Map.of());
        String table = Files.readString(out.resolve("feature-annotation.tsv"));
        assertEquals(2, table.lines().filter(line -> line.endsWith("unmapped")).count());
        out = assess("MYGENE\n", Map.of("feature_id_type", "gene_symbol"));
        assertTrue(Files.readString(out.resolve("feature-annotation.tsv")).contains("MYGENE\t\tMYGENE"));
    }

    @Test void conflictingBuildInTableFailsClearly() throws Exception {
        Files.writeString(project.resolve("annotation.tsv"), "feature_id\tchr\tposition\tgenome_build\np1\t1\t100\tGRCh37\n");
        var error = assertThrows(java.io.IOException.class, () -> assess("p1\n", Map.of(
            "files", Map.of("annotation", "annotation.tsv"), "genome_build", "GRCh38", "coordinate_system", "1-based-closed")));
        assertTrue(error.getMessage().contains("Conflicting genome builds"));
    }

    @Test void geneIdentifierIsNotDisplayedAsAnHgncSymbol() throws Exception {
        Files.writeString(project.resolve("annotation.tsv"), "feature_id\tgene_id\tgene_symbol\np1\tENSG00001\tENSG00001\n");
        Path out = assess("p1\n", Map.of("files", Map.of("annotation", "annotation.tsv")));
        assertTrue(Files.readString(out.resolve("feature-annotation.tsv")).contains("p1\tENSG00001\t\t"));
    }

    private Path assess(String features, Map<String, Object> settings) throws Exception {
        Path out = Files.createTempDirectory(project, "output-");
        Path rows = out.resolve("standard.tsv");
        Files.writeString(rows, "feature_id\n" + features);
        Map<String, Object> dataset = new LinkedHashMap<>(settings);
        dataset.put("id", "d1");
        Map<String, Object> study = Map.of("datasets", List.of(dataset), "resources", Map.of("downloads", "deny"));
        StudyFeatureAnnotations.write(new Context(project, out, Map.of(), study, Map.of(), Map.of()), Map.of("d1", rows));
        return out;
    }
}
