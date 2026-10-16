/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import org.jlinalg.study.StudyOperationExecutor.Context;
import org.jlinalg.study.StudyOperationExecutor.State;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class StudyEffectQcTest {
    @TempDir Path directory;

    private Path table(String name, double multiplier, double extreme) throws Exception {
        Path path = directory.resolve(name + ".tsv");
        StringBuilder rows = new StringBuilder("feature_id\tbeta\tse\tp_value\n");
        for (int i=0;i<30;i++) rows.append("feature-").append(i).append('\t').append((i-15)*.01*multiplier)
            .append('\t').append(.1*multiplier).append("\t.000001\n");
        rows.append("extreme\t").append(extreme*multiplier).append('\t').append(200000*multiplier).append("\t.000001\n");
        Files.writeString(path,rows); return path;
    }

    private StudyEffectQc.Result apply(Path path, Map<String,Object> settings) throws Exception {
        return StudyEffectQc.apply(Map.of("data",path),Map.of("effect_outliers",settings),directory);
    }

    @Test void excludesPositiveAndNegativeExtremesWithoutUsingSignificanceOrClipping() throws Exception {
        for (double extreme : List.of(600000.0,-600000.0)) {
            Path path=table("input",1,extreme);
            var qc=apply(path,Map.of());
            assertEquals(1,qc.document().get("excluded_rows"));
            assertFalse(Files.readString(path).contains("extreme"));
            assertTrue(Files.readString(path).contains("feature-29\t0.14"));
            String audit=Files.readString(directory.resolve("effect-outlier-exclusions.tsv"));
            assertTrue(audit.contains("data\textreme\t"+extreme));
            assertTrue(audit.contains("extreme_robust_effect"));
        }
    }

    @Test void selectionIsInvariantToUniformUnitsAndDoesNotRejectUniformlyLargeEffects() throws Exception {
        for (double multiplier : List.of(1e-6,1.0,1e6)) {
            Path path=table("scaled",multiplier,600000);
            assertEquals(1,apply(path,Map.of()).document().get("excluded_rows"));
        }
        Path path=table("large-units",1e8,.1);
        assertEquals(0,apply(path,Map.of()).document().get("excluded_rows"));
    }

    @Test void smallAndZeroScaleFamiliesAreAuditedRatherThanDividingByZero() throws Exception {
        Path path=directory.resolve("tiny.tsv");
        Files.writeString(path,"feature_id\tbeta\tse\tp_value\na\t1\t1\t.9\nb\t600000\t200000\t.01\n");
        var qc=apply(path,Map.of());
        assertEquals(0,qc.document().get("excluded_rows"));
        assertTrue(qc.document().toString().contains("insufficient_features"));
        StringBuilder rows=new StringBuilder("feature_id\tbeta\tse\tp_value\n");
        for(int i=0;i<30;i++)rows.append("f").append(i).append("\t0\t1\t.5\n");
        rows.append("extreme\t600000\t200000\t.01\n"); Files.writeString(path,rows);
        qc=apply(path,Map.of());
        assertEquals(0,qc.document().get("excluded_rows"));
        assertTrue(qc.document().toString().contains("zero_or_nonfinite_robust_scale"));
        assertEquals(1,apply(path,Map.of("max_abs_effect",1000)).document().get("excluded_rows"));
    }

    @Test void zeroMadUsesIqrFallbackAndAbsoluteBoundIsStrict() throws Exception {
        Path path=directory.resolve("fallback.tsv");
        StringBuilder rows=new StringBuilder("feature_id\tbeta\tse\tp_value\n");
        for(int i=0;i<30;i++)rows.append("f").append(i).append('\t').append(i<20 ? 0 : 1).append("\t1\t.5\n");
        rows.append("extreme\t600000\t200000\t.01\n"); Files.writeString(path,rows);
        var qc=apply(path,Map.of());
        assertEquals(1,qc.document().get("excluded_rows"));
        assertTrue(qc.document().toString().contains("normal_consistent_IQR"));
        Files.writeString(path,"feature_id\tbeta\tse\tp_value\na\t1\t1\t.01\nb\t-1\t1\t.01\nc\t1.001\t1\t.01\n");
        assertEquals(1,apply(path,Map.of("max_abs_effect",1)).document().get("excluded_rows"));
        assertTrue(Files.readString(path).contains("b\t-1"));
    }

    @Test void disabledPolicyPreservesRowsAndInvalidSettingsFailClearly() throws Exception {
        Path path=table("disabled",1,600000); String original=Files.readString(path);
        assertEquals(0,apply(path,Map.of("enabled",false,"max_abs_effect",1)).document().get("excluded_rows"));
        assertEquals(original,Files.readString(path));
        for(Map<String,Object> settings : List.of(Map.<String,Object>of("enabled","yes"),Map.<String,Object>of("robust_z",0),
                Map.<String,Object>of("robust_z",Double.NaN),Map.<String,Object>of("min_features",2),Map.<String,Object>of("min_features",3.5),Map.<String,Object>of("max_abs_effect",-1)))
            assertThrows(IllegalArgumentException.class,()->apply(path,settings));
        Files.writeString(path,"feature_id\tbeta\tse\tp_value\na\t1\t1\t.5\nb\t2\t1\t.5\n");
        assertThrows(IOException.class,()->apply(path,Map.of("max_abs_effect",1e-12)));
        assertTrue(Files.exists(directory.resolve("effect-outlier-qc.yaml")));
    }

    @Test void harmonizationPropagatesExclusionsThroughCommonUniverseComparisonsMetaAndMr() throws Exception {
        Path left=table("left",1,600000),right=table("right",1,.1);
        String original=Files.readString(left);
        Map<String,Object> analysis=Map.of("harmonization",Map.of("common_universe","intersection"),"bioinformatics",Map.of("fdr",.99));
        Map<String,Object> study=Map.of("resources",Map.of("downloads","deny"),"hypotheses",List.of(Map.of("id","h","analysis",analysis)),
            "datasets",List.of(Map.of("id","left","cohort","c1","stage","association_results","type","generic","files",Map.of("results",left.toString())),
                Map.of("id","right","cohort","c2","stage","association_results","type","generic","files",Map.of("results",right.toString()))));
        Path output=Files.createDirectory(directory.resolve("harmonized"));
        Context context=new Context(directory,output,Map.of("type","validate_harmonize","hypothesis","h","datasets",List.of("left","right")),study,Map.of(),Map.of());
        var result=new StudyCliExecutor().execute(context); assertEquals(State.COMPLETE,result.state());
        assertEquals(original,Files.readString(left),"The supplied source must remain unchanged");
        for(String dataset:List.of("left","right"))assertFalse(Files.readString(output.resolve("standard-results/"+dataset+".tsv")).contains("extreme"));
        assertFalse(Files.readString(output.resolve("feature-annotation.tsv")).contains("extreme"));
        assertTrue(result.details().containsKey("method_report"));
        Path comparisons=Files.createDirectory(directory.resolve("comparisons"));
        var compared=new StudyCliExecutor().execute(new Context(directory,comparisons,Map.of("type","effect_comparison","hypothesis","h"),study,Map.of(),Map.of("harmonized",output)));
        assertEquals(State.COMPLETE,compared.state());
        assertTrue(Files.readString(comparisons.resolve("effect-comparisons.tsv")).contains("\t30\t"));
        Path meta=Files.createDirectory(directory.resolve("meta"));
        assertEquals(State.COMPLETE,new StudyCliExecutor().execute(new Context(directory,meta,Map.of("type","meta_analysis","hypothesis","h"),study,Map.of(),Map.of("harmonized",output))).state());
        assertFalse(Files.readString(meta.resolve("meta-analysis.tsv")).contains("extreme"));
        Files.writeString(meta.resolve("manifest.yaml"),"operation_type: meta_analysis\nstate: complete\n");
        Files.writeString(output.resolve("manifest.yaml"),"operation_type: validate_harmonize\nstate: complete\n");
        Files.writeString(output.resolve("feature-annotation.tsv"),"feature_id\tgene_symbol\tgene_id\tchromosome\tstart\tend\tgenome_build\tmapping_status\nfeature-0\tTEST\tGENE\t1\t100\t300\tGRCh38\tlocated\n");
        Context mr=new Context(directory,meta,Map.of("type","mendelian_randomization","hypothesis","h"),study,Map.of(),Map.of("harmonized",output,"meta",meta));
        assertFalse(StudyBioinformatics.mrFeatures(mr).isEmpty());
        assertFalse(StudyBioinformatics.mrFeatures(mr).stream().anyMatch(feature->feature.id().equals("extreme")));
    }
}
