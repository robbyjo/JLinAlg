/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class ConditionalScoreCliTest {
    private static final Path FIXTURE = Path.of(
        "src/test/resources/conditional-score-export");
    @TempDir Path temporary;

    @Test void importsExportSchemaAndMatchesHandSchurComplement() throws Exception {
        Path summary = export("complete", 64);
        Path cohorts = cohortManifest("one", summary);
        Path output = temporary.resolve("conditional.tsv");
        Result run = run(List.of("conditional-score", "--cohorts",
            cohorts.toString(), "--targets", "1:200:A:G,1:300:A:G",
            "--condition-on", "1:100:A:G", "--out", output.toString()));
        assertEquals(0, run.code, run.console);

        Properties expected = new Properties();
        try (var reader = Files.newBufferedReader(
                FIXTURE.resolve("binomial.properties"))) {
            expected.load(reader);
        }
        double[] u = numbers(expected.getProperty("u"));
        double[] v = numbers(expected.getProperty("v"));
        DelimitedData table = DelimitedData.read(output);
        assertEquals(2, table.rows().size());
        int[] target = {1,2};
        for (int i = 0; i < target.length; i++) {
            String[] row = table.rows().get(i); int a = target[i];
            double adjustedU = u[a] - v[a*4]*u[0]/v[0];
            double adjustedV = v[a*4+a] - v[a*4]*v[a]/v[0];
            assertEquals(adjustedU,
                Double.parseDouble(row[table.column("score_u")]), 2e-6);
            assertEquals(adjustedV,
                Double.parseDouble(row[table.column("score_variance")]), 2e-6);
            assertEquals("local_schur_one_step",
                row[table.column("inference_scope")]);
            assertEquals("1", row[table.column("cohorts")]);
        }
        DelimitedData covariance = DelimitedData.read(Path.of(
            output + ".score-cov.tsv"));
        assertEquals(3, covariance.rows().size());
        double expectedCross = v[1*4+2] - v[1*4]*v[2]/v[0];
        assertEquals(expectedCross, Double.parseDouble(
            covariance.rows().get(1)[covariance.column("score_covariance")]),
            2e-6);
        String metadata = Files.readString(Path.of(output + ".metadata.tsv"));
        assertTrue(metadata.contains("nonlinear_cohort_refit\tfalse"));
        assertTrue(metadata.contains("unknown and rejected; never filled with zero"));
    }

    @Test void crossBlockCovarianceRemainsUnknownAndIncompleteFilesReject()
            throws Exception {
        Path summary = export("blocks", 2);
        Path cohorts = cohortManifest("one", summary);
        Path output = temporary.resolve("cross-block.tsv");
        Result cross = run(List.of("conditional-score", "--cohorts",
            cohorts.toString(), "--targets", "1:300:A:G",
            "--condition-on", "1:100:A:G", "--out", output.toString()));
        assertEquals(2, cross.code);
        assertTrue(cross.console.contains("unknown, not zero"), cross.console);
        assertFalse(Files.exists(output));

        Path covariance = Path.of(summary + ".score-cov.tsv");
        List<String> rows = new ArrayList<>(Files.readAllLines(covariance));
        rows.remove(rows.size()-1);
        Path incomplete = temporary.resolve("incomplete.score-cov.tsv");
        Files.write(incomplete, rows);
        IOException failure = assertThrows(IOException.class,
            () -> ConditionalScoreImporter.read("one", summary, incomplete,
                Path.of(summary + ".score-manifest.json"),
                List.of("1:100:A:G","1:200:A:G")));
        assertTrue(failure.getMessage().contains("complete upper triangle"),
            failure.getMessage());
    }

    @Test void helpSeparatesScoreConditioningFromMrAndCohortRefits() {
        Result result = run(List.of("conditional-score", "--help"));
        assertEquals(0, result.code);
        assertTrue(result.console.contains("not a nonlinear cohort refit"));
        assertTrue(result.console.contains("not the separate mr-estimate"));
    }

    @Test void fractionalManifestSchemaVersionRejects() throws Exception {
        Path summary = export("fractional-schema", 64);
        Path original = Path.of(summary + ".score-manifest.json");
        String json = Files.readString(original);
        assertTrue(json.contains("\"schema_version\": 1,"));
        Path fractional = temporary.resolve("fractional-schema.json");
        Files.writeString(fractional, json.replace(
            "\"schema_version\": 1,", "\"schema_version\": 1.5,"));
        IOException failure = assertThrows(IOException.class,
            () -> ConditionalScoreImporter.read("one", summary,
                Path.of(summary + ".score-cov.tsv"), fractional,
                List.of("1:100:A:G","1:200:A:G")));
        assertTrue(failure.getMessage().contains("schema_version 1"),
            failure.getMessage());
    }

    private Path export(String name, int blockSize) {
        Path output = temporary.resolve(name + ".tsv");
        List<String> arguments = new ArrayList<>(List.of(
            "--omics", FIXTURE.resolve("variants.vcf").toString(),
            "--pheno", FIXTURE.resolve("phenotype.tsv").toString(),
            "--id", "IID", "--model", "glm",
            "--formula", "binary ~ x + offset(o) + <omics>",
            "--threads", "1", "--out", output.toString(),
            "--family", "binomial", "--conditional-gwas-summary",
            "--score-genome-build", "GRCh38", "--score-block-size",
            Integer.toString(blockSize)));
        Result result = run(arguments);
        assertEquals(0, result.code, result.console);
        return output;
    }

    private Path cohortManifest(String cohort, Path summary) throws Exception {
        Path path = temporary.resolve(summary.getFileName() + "-cohorts.tsv");
        Files.writeString(path, "cohort\tsummary\tcovariance\tmanifest\n"
            + cohort + "\t" + summary.getFileName() + "\t"
            + Path.of(summary + ".score-cov.tsv").getFileName() + "\t"
            + Path.of(summary + ".score-manifest.json").getFileName() + "\n");
        return path;
    }

    private static Result run(List<String> arguments) {
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        PrintStream stream = new PrintStream(bytes);
        int code = JLinAlgCli.run(arguments.toArray(String[]::new), stream, stream);
        return new Result(code, bytes.toString());
    }
    private static double[] numbers(String value) {
        return Arrays.stream(value.split(","))
            .mapToDouble(Double::parseDouble).toArray();
    }
    private record Result(int code, String console) { }
}
