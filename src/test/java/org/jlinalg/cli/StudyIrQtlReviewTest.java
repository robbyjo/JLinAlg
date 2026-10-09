package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;
import java.util.Set;
import org.apache.commons.compress.compressors.bzip2.BZip2CompressorOutputStream;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class StudyIrQtlReviewTest {
    @TempDir Path directory;
    @Test void compressedExactTranscriptReviewPreservesLogPAndBlocksFilteredColoc() throws Exception {
        Path source = directory.resolve("data.csv.bz2"), output = directory.resolve("candidates.tsv");
        String csv = "SNP,SNP_Chr,SNP_Pos,EAF,N,Tx,Is_Cis,RSq,Fx,T,log10P\n"
            + "chr1:123:G:A:rs1,1,123,0.002,612,ENST1.2,1,.99999,-.28,-7000000,-3186\n"
            + "chr1:456:G:C:rs2,1,456,0.2,1200,ENST1.1,1,.02,.1,7,-10\n";
        try (var out = new BZip2CompressorOutputStream(Files.newOutputStream(source))) { out.write(csv.getBytes(java.nio.charset.StandardCharsets.UTF_8)); }
        var result = StudyIrQtlReview.scan(source, Map.of("genome_build", "GRCh38", "effect_allele", "ALT", "complete_region", false,
            "export_p_threshold", 1e-7, "rationale", "Match measured isoform ratio to same transcript-version QTL"), Set.of("ENST1.2"), output);
        assertEquals(2L, result.get("rows_scanned")); assertEquals(1L, result.get("exact_feature_matches"));
        assertEquals("blocked_significance_filtered_irQTL", result.get("colocalization_status"));
        String table = Files.readString(output);
        assertTrue(table.contains("\tA\tG\t-0.28\t")); assertTrue(table.contains("-3186.0"));
        assertTrue(table.contains("near_perfect_fit")); assertFalse(table.contains("ENST1.1"));
        assertTrue(Files.readString(output.resolveSibling("candidates.tsv.source.yaml")).contains("scientific_rationale"));
    }
    @Test void undeclaredAllelesOrBuildCannotBeGuessed() {
        assertThrows(IllegalArgumentException.class, () -> StudyIrQtlReview.scan(directory.resolve("absent"), Map.of(), Set.of(), directory.resolve("out")));
    }
}
