/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

final class MrEstimatorExtensionTest {
    @TempDir Path directory;
    private record Run(int status,String output,String error) { }
    private Run run(String method,String...extra) throws Exception {
        Path input=Path.of("src/test/resources/r-reference/mr-conditional-input.tsv");
        Path ld=Path.of("src/test/resources/r-reference/mr-conditional-ld.tsv");
        List<String> args=new ArrayList<>(List.of("mr-estimate","--input",input.toString(),"--ld",ld.toString(),"--method",method,"--backend","cpu"));
        args.addAll(List.of(extra));
        ByteArrayOutputStream output=new ByteArrayOutputStream(),error=new ByteArrayOutputStream();
        int status=JLinAlgCli.run(args.toArray(String[]::new),new PrintStream(output,true,StandardCharsets.UTF_8),new PrintStream(error,true,StandardCharsets.UTF_8));
        return new Run(status,output.toString(StandardCharsets.UTF_8),error.toString(StandardCharsets.UTF_8));
    }
    @Test void conditionalRequiresExplicitConditioningAndReportsInstrumentTests() throws Exception {
        Run missing=run("conditional"); assertEquals(2,missing.status()); assertTrue(missing.error().contains("--condition-on"));
        Run conditional=run("conditional","--condition-on","rs1,rs3");
        assertEquals(0,conditional.status(),conditional.error());
        assertTrue(conditional.output().contains("conditional_instrument_p_value"));
        assertFalse(conditional.output().contains("IVW_GENERALIZED_FIXED"));
        assertEquals(6,conditional.output().lines().count());
        assertTrue(conditional.output().contains("rs1,rs3"));
        assertEquals(2,run("conditional","--condition-on","absent").status());
        assertEquals(0,run("conditional","--joint","--association","outcome").status());
        assertEquals(2,run("conditional","--joint","--plot",directory.resolve("bad.svg").toString()).status());
        assertEquals(2,run("ivw-generalized-fixed","--condition-on","rs1").status());
    }
    @Test void generalizedRetainsCausalEffectInference() throws Exception {
        Run result=run("ivw-generalized-fixed","--confidence","0.9");
        assertEquals(0,result.status(),result.error());
        assertTrue(result.output().contains("IVW_GENERALIZED_FIXED"));
        assertTrue(result.output().contains("causal_p_value"));
    }
    @Test void overlapConfidenceIsPropagatedIncludingAllPath() throws Exception {
        Path covariance=directory.resolve("cov.tsv");
        Files.writeString(covariance,"0.00024\n0.0004\n0.000315\n0.000175\n0.000495\n");
        for(String method:List.of("overlap-aware","all")) {
            Run narrow=run(method,"--sampling-covariance",covariance.toString(),"--confidence","0.9");
            Run wide=run(method,"--sampling-covariance",covariance.toString(),"--confidence","0.99");
            assertEquals(0,narrow.status(),narrow.error()); assertEquals(0,wide.status(),wide.error());
            String[] n=narrow.output().lines().filter(s->s.startsWith("OVERLAP_AWARE_IVW")).findFirst().orElseThrow().split("\t");
            String[] w=wide.output().lines().filter(s->s.startsWith("OVERLAP_AWARE_IVW")).findFirst().orElseThrow().split("\t");
            assertEquals(n[1],w[1]); assertEquals(n[3],w[3]); assertTrue(Double.parseDouble(w[5])>Double.parseDouble(n[5]));
        }
    }
    @Test void secondarySignalsRunWithoutOutcomeColumnsAndReturnSelection() throws Exception {
        Path input=directory.resolve("exposure.tsv"),ld=directory.resolve("ld.tsv");
        Files.writeString(input,"variant_id\tbeta\tse\nproxy\t1.3\t.1\ncausal1\t1\t.1\ncausal2\t1\t.1\n");
        Files.writeString(ld,"1\t.65\t.65\n.65\t1\t0\n.65\t0\t1\n");
        Run selection=run("secondary-signals","--input",input.toString(),"--ld",ld.toString(),"--p-threshold",".05");
        assertEquals(0,selection.status(),selection.error());
        var rows=selection.output().lines().skip(1).toList();
        assertTrue(rows.get(0).endsWith("causal1,causal2\tfalse"));
        assertTrue(rows.get(1).endsWith("causal2\ttrue"));
        assertTrue(rows.get(2).endsWith("causal1\ttrue"));
    }
    @Test void eggerCliPreservesFullFitForPlotting() throws Exception {
        for(String method:List.of("egger","egger-generalized")) {
            Path svg=directory.resolve(method+".svg");
            Run result=run(method,"--plot",svg.toString());
            assertEquals(0,result.status(),result.error());
            assertTrue(Files.readString(svg).contains("<line"));
        }
    }
}
