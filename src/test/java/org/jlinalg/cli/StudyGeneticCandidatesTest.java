package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;

import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class StudyGeneticCandidatesTest {
    @TempDir Path directory;

    @Test void filtersUnsupportedRowsAndKeepsTheStrongestUniqueRsid() throws Exception {
        Path input = directory.resolve("candidates.tsv");
        String header = "feature_id\tvariant_id\tchromosome\tposition\teffect_allele\tother_allele\tbeta\tse\tlog10_p\teaf\tn\treview_flags\n";
        Files.writeString(input, header
            + "tx1\tchr1:100:G:A:rs1\t1\t100\tA\tG\t.1\t.01\t-8\t.2\t4700\t.\n"
            + "tx1\tchr1:100:G:A:rs1\t1\t100\tA\tG\t.2\t.01\t-10\t.2\t4700\t.\n"
            + "tx1\tchr1:200:GA:G:rs2\t1\t200\tGA\tG\t.1\t.01\t-9\t.2\t4700\t.\n"
            + "tx1\tchr1:300:G:A:rs3\t1\t300\tA\tG\t.1\t0\t-9\t.2\t4700\t.\n"
            + "tx1\tchr1:400:G:A:rs4\t1\t400\tA\tG\t.1\t.01\t1\t.2\t4700\t.\n"
            + "tx1\tchr1:500:G:A:rs5\t1\t500\tA\tG\t.1\t.01\t-7\t.2\t4700\t.\n");
        var selection = StudyGeneticCandidates.read(input, 5e-8);
        assertEquals(1, selection.byFeature().get("tx1").size());
        assertEquals(-10, selection.byFeature().get("tx1").get(0).log10p());
        assertEquals(1, selection.exclusions().get("duplicate_variant_id"));
        assertEquals(1, selection.exclusions().get("unsupported_non_biallelic_snp_alleles"));
        assertEquals(1, selection.exclusions().get("nonpositive_standard_error"));
        assertEquals(1, selection.exclusions().get("invalid_log10_p"));
        assertEquals(1, selection.exclusions().get("above_instrument_p_threshold"));
        Path output = directory.resolve("exposure.tsv");
        StudyGeneticCandidates.writeExposure(output,
            selection.byFeature().get("tx1"));
        String written = Files.readString(output);
        assertEquals(1, written.lines().filter(line -> line.startsWith("tx1\t")).count());
        assertFalse(written.contains("GA\tG"));
    }
}
