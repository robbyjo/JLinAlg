/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

final class JointMrCliTest {
    @TempDir Path directory;
    private Path input, correlation;

    @BeforeEach void writeInputs() throws Exception {
        input = directory.resolve("joint.tsv");
        correlation = directory.resolve("outcome-correlation.tsv");
        StringBuilder text = new StringBuilder(
            "variant_id\tbeta_exposure_x1\tse_exposure_x1\t"
            + "beta_exposure_x2\tse_exposure_x2\t"
            + "beta_outcome_y1\tse_outcome_y1\t"
            + "beta_outcome_y2\tse_outcome_y2\n");
        double[][] x = {{.10,.02},{.02,.12},{.15,-.03},{-.04,.16},
            {.12,.10},{-.08,.07},{.18,-.06},{.04,-.11}};
        for (int index = 0; index < x.length; index++)
            text.append("v").append(index + 1).append('\t')
                .append(x[index][0]).append("\t.01\t")
                .append(x[index][1]).append("\t.012\t")
                .append(.7*x[index][0]-.4*x[index][1]).append("\t.02\t")
                .append(-.2*x[index][0]+.9*x[index][1]).append("\t.025\n");
        Files.writeString(input, text);
        Files.writeString(correlation, "1\t.4\n.4\t1\n");
    }

    @Test void mvmrAndMultivariateCommandsHaveDistinctContracts()
            throws Exception {
        Run mvmr = run("mr-mvmr", "--input", input.toString(),
            "--exposures", "x1,x2", "--outcome", "y1", "--backend", "cpu");
        assertEquals(0, mvmr.status(), mvmr.error());
        assertEquals(3, mvmr.output().lines().count());
        assertTrue(mvmr.output().startsWith(
            "method\texposure\toutcome\ttest\ttest_scope\ttest_df\ttest_p_value\t"));
        assertTrue(mvmr.output().contains("MVMR_IVW\tx1\ty1"));

        Run multivariate = run("mr-multivariate", "--input", input.toString(),
            "--exposures", "x1,x2", "--outcomes", "y1,y2",
            "--outcome-correlation", correlation.toString(),
            "--backend", "cpu", "--plot",
            directory.resolve("multivariate.svg").toString());
        assertEquals(0, multivariate.status(), multivariate.error());
        assertEquals(5, multivariate.output().lines().count());
        assertTrue(multivariate.output().startsWith(
            "method\texposure\toutcome\ttest\ttest_scope\ttest_df\ttest_p_value\t"));
        for (String row : multivariate.output().lines().toList())
            assertEquals(27, row.split("\t", -1).length, row);
        assertTrue(multivariate.output().contains("MULTIVARIATE_IVW_FIXED"));
        assertTrue(multivariate.output().contains("overall_p_value"));
        assertTrue(Files.readString(directory.resolve("multivariate.svg"))
            .contains("<svg"));

        assertEquals(2, run("mr-mvmr", "--input", input.toString(),
            "--exposures", "x1").status());
        assertEquals(2, run("mr-multivariate", "--input", input.toString(),
            "--exposures", "x1", "--outcomes", "y1",
            "--outcome-correlation", correlation.toString()).status());
    }

    @Test void pressoRequiresAndWritesAuditableDiagnostics() throws Exception {
        Path output = directory.resolve("presso.tsv");
        Path diagnostics = directory.resolve("presso-diagnostics.tsv");
        Run result = run("mr-multivariate", "--input", input.toString(),
            "--exposures", "x1,x2", "--outcomes", "y1,y2",
            "--outcome-correlation", correlation.toString(),
            "--method", "presso", "--simulations", "20", "--seed", "17",
            "--output", output.toString(), "--diagnostics", diagnostics.toString(),
            "--backend", "cpu");
        assertEquals(0, result.status(), result.error());
        assertTrue(Files.readString(output).contains("MULTIVARIATE_PRESSO_RAW"));
        String audit = Files.readString(diagnostics);
        assertTrue(audit.contains("mahalanobis_distance"));
        assertTrue(audit.contains("global_p_value"));
        assertEquals(9, audit.lines().count());
    }

    private Run run(String... arguments) {
        String[] logged = new String[arguments.length + 1];
        logged[0] = arguments[0];
        logged[1] = "--no-log";
        System.arraycopy(arguments, 1, logged, 2, arguments.length - 1);
        ByteArrayOutputStream output = new ByteArrayOutputStream();
        ByteArrayOutputStream error = new ByteArrayOutputStream();
        int status = JLinAlgCli.run(logged,
            new PrintStream(output, true, StandardCharsets.UTF_8),
            new PrintStream(error, true, StandardCharsets.UTF_8));
        return new Run(status, output.toString(StandardCharsets.UTF_8),
            error.toString(StandardCharsets.UTF_8));
    }

    private record Run(int status, String output, String error) { }
}
