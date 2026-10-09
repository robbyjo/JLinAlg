/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class PlotCliTest {
    @TempDir Path directory;

    @Test
    void commandRendersEveryDocumentedChartFamily() throws Exception {
        Path data = data();
        List<Invocation> invocations = List.of(
            new Invocation("scatter", "--x", "x", "--y", "y", "--group", "group"),
            new Invocation("line", "--x", "x", "--y", "y", "--group", "group"),
            new Invocation("connected-scatter", "--x", "x", "--y", "y"),
            new Invocation("area", "--x", "x", "--y", "y"),
            new Invocation("step", "--x", "x", "--y", "y"),
            new Invocation("bubble", "--x", "x", "--y", "y", "--size", "size",
                "--group", "group"),
            new Invocation("histogram", "--x", "y", "--bins", "4"),
            new Invocation("density", "--x", "y"),
            new Invocation("dot", "--x", "category", "--y", "value"),
            new Invocation("bar", "--x", "category", "--y", "value", "--group",
                "group", "--position", "stacked"),
            new Invocation("pie", "--x", "category", "--y", "value"),
            new Invocation("donut", "--x", "category", "--y", "value"),
            new Invocation("box", "--x", "group", "--y", "y"),
            new Invocation("violin", "--x", "group", "--y", "y"));
        for (Invocation invocation : invocations) {
            Path plot = directory.resolve(invocation.type + ".png");
            List<String> arguments = base(invocation.type, data, plot);
            arguments.addAll(invocation.arguments);
            Run result = run(arguments);
            assertEquals(0, result.status, invocation.type + ": " + result.messages);
            byte[] signature = Files.readAllBytes(plot);
            assertTrue(signature.length > 1_000, invocation.type);
            assertArrayEquals(new byte[] {(byte) 0x89, 0x50, 0x4e, 0x47},
                java.util.Arrays.copyOf(signature, 4), invocation.type);
        }
    }

    @Test
    void commandExportsRasterAndVectorFormats() throws Exception {
        Path data = data();
        for (String extension : List.of("png", "svg", "pdf")) {
            Path plot = directory.resolve("scatter." + extension);
            Run result = run(base("scatter", data, plot, "--x", "x", "--y", "y",
                "--title", "Format test", "--theme", "minimal"));
            assertEquals(0, result.status, result.messages);
            byte[] bytes = Files.readAllBytes(plot);
            assertTrue(bytes.length > 500);
            if (extension.equals("png"))
                assertArrayEquals(new byte[] {(byte) 0x89, 0x50, 0x4e, 0x47},
                    java.util.Arrays.copyOf(bytes, 4));
            if (extension.equals("svg"))
                assertTrue(new String(bytes, StandardCharsets.UTF_8).contains("<svg"));
            if (extension.equals("pdf"))
                assertEquals("%PDF", new String(bytes, 0, 4, StandardCharsets.US_ASCII));
        }
    }

    @Test
    void commandRendersClusteredHeatmapsAndExtendedGalleryFamilies() throws Exception {
        Path heatmap = directory.resolve("heatmap.csv");
        Files.writeString(heatmap, """
            column,row,value
            A,R1,1
            B,R1,4
            C,R1,2
            A,R2,5
            B,R2,2
            C,R2,6
            A,R3,2
            B,R3,7
            C,R3,3
            """);
        for (String mode : List.of("rows", "columns", "both"))
            assertPng(run(base("heatmap", heatmap,
                directory.resolve("heatmap-" + mode + ".png"),
                "--x", "column", "--y", "row", "--fill", "value",
                "--cluster", mode, "--colors", "#2166AC,#FFFFFF,#B2182B")),
                directory.resolve("heatmap-" + mode + ".png"), "heatmap-" + mode);

        Path spider = directory.resolve("spider.csv");
        Files.writeString(spider, """
            axis,value,series
            speed,0.8,A
            power,0.6,A
            range,0.9,A
            speed,0.5,B
            power,0.9,B
            range,0.4,B
            """);
        check("spider", spider, "--x", "axis", "--y", "value", "--group", "series");

        Path observations = directory.resolve("observations.csv");
        Files.writeString(observations, """
            id,group,a,b,c,d
            S1,G1,1,8,2,4
            S2,G1,2,7,3,5
            S3,G2,8,1,7,5
            S4,G2,7,2,8,4
            """);
        check("parallel", observations, "--id", "id", "--columns", "a,b,c,d",
            "--group", "group", "--line-pattern", "dash-dot", "--line-width", "2.5",
            "--point-shape", "diamond", "--point-size", "10");
        check("dendrogram", observations, "--id", "id", "--columns", "a,b,c,d");
        check("dendrogram-horizontal", "dendrogram", observations, "--id", "id",
            "--columns", "a,b,c,d", "--orientation", "horizontal");

        Path categories = data();
        check("circular-bar", categories, "--x", "category", "--y", "value",
            "--group", "group");
        check("treemap", categories, "--x", "category", "--y", "value",
            "--group", "group");

        Path edges = directory.resolve("edges.csv");
        Files.writeString(edges, """
            source,target,weight,kind
            A,B,5,first
            A,C,3,first
            B,D,4,second
            C,D,2,second
            C,E,1,third
            D,F,5,third
            E,F,1,third
            """);
        for (String type : List.of("network", "chord", "edge-bundling", "sankey"))
            check(type, edges, "--source", "source", "--target", "target",
                "--weight", "weight", "--group", "kind");

        Path memberships = directory.resolve("memberships.csv");
        Files.writeString(memberships, """
            feature,set
            A,Discovery
            B,Discovery
            B,Replication
            C,Replication
            """);
        check("venn", memberships, "--x", "feature", "--group", "set");
    }

    @Test
    void commandReportsMissingRowsAndProtectsExistingOutputs() throws Exception {
        Path data = data();
        Path plot = directory.resolve("missing.png");
        Run first = run(base("scatter", data, plot, "--x", "x", "--y", "missing"));
        assertEquals(0, first.status, first.messages);
        assertTrue(first.messages.contains("omitted 1 row"));
        byte[] original = Files.readAllBytes(plot);
        Run second = run(base("scatter", data, plot, "--x", "x", "--y", "y"));
        assertEquals(2, second.status);
        assertArrayEquals(original, Files.readAllBytes(plot));
        Path rejected = directory.resolve("rejected.png");
        Run strict = run(base("scatter", data, rejected, "--x", "x", "--y", "missing",
            "--na", "error"));
        assertEquals(2, strict.status);
        assertTrue(strict.messages.contains("missing value"));
    }

    @Test
    void helpListsMappingsAndFormats() {
        Run result = run(List.of("plot", "--help", "--no-config"));
        assertEquals(0, result.status);
        assertTrue(result.messages.contains("--group COLUMN"));
        assertTrue(result.messages.contains(".pdf"));
        assertTrue(result.messages.contains("violin"));
        assertTrue(result.messages.contains("venn"));
    }

    private Path data() throws Exception {
        Path path = directory.resolve("plot-data.csv");
        Files.writeString(path, """
            x,y,group,size,category,value,missing
            1,2,A,4,alpha,3,2
            2,3,A,9,beta,4,3
            3,2,B,16,alpha,2,NA
            4,5,B,25,beta,6,5
            5,4,B,36,gamma,5,4
            """);
        return path;
    }

    private static List<String> base(String type, Path input, Path output,
            String... additional) {
        List<String> result = new ArrayList<>(List.of("plot", "--type", type,
            "--input", input.toString(), "--out", output.toString(), "--no-config"));
        result.addAll(List.of(additional));
        return result;
    }

    private static Run run(List<String> arguments) {
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        PrintStream stream = new PrintStream(bytes, true, StandardCharsets.UTF_8);
        int status = JLinAlgCli.run(arguments.toArray(String[]::new), stream, stream);
        return new Run(status, bytes.toString(StandardCharsets.UTF_8));
    }

    private void check(String type, Path input, String... arguments) throws Exception {
        check(type, type, input, arguments);
    }

    private void check(String name, String type, Path input, String... arguments)
            throws Exception {
        Path plot = directory.resolve(name + ".png");
        assertPng(run(base(type, input, plot, arguments)), plot, name);
    }

    private static void assertPng(Run run, Path plot, String name) throws Exception {
        assertEquals(0, run.status, name + ": " + run.messages);
        byte[] bytes = Files.readAllBytes(plot);
        assertTrue(bytes.length > 1_000, name);
        assertArrayEquals(new byte[] {(byte) 0x89, 0x50, 0x4e, 0x47},
            java.util.Arrays.copyOf(bytes, 4), name);
    }

    private record Run(int status, String messages) { }

    private static final class Invocation {
        private final String type;
        private final List<String> arguments;
        private Invocation(String type, String... arguments) {
            this.type = type;
            this.arguments = List.of(arguments);
        }
    }
}
