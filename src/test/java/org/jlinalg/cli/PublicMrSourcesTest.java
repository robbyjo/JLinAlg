package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;

import htsjdk.samtools.liftover.LiftOver;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class PublicMrSourcesTest {
    @TempDir Path directory;
    @Test void bmiSelectsIndependentPanUkbEurOutcome() {
        assertTrue(PanUkbbGwasOutcome.supports("BMI", "EUR"));
        assertTrue(PanUkbbGwasOutcome.supports("body mass index", "EUR"));
        assertTrue(PanUkbbGwasOutcome.supports("obesity", "EUR"));
        assertEquals(419163, PanUkbbGwasOutcome.descriptor().get("sample_size"));
        assertEquals("GRCh37", PanUkbbGwasOutcome.descriptor().get("source_genome_build"));
        var obesity = PanUkbbGwasOutcome.profile("Obesity", "EUR");
        assertEquals("icd10-E66", obesity.trait());
        assertEquals(420531, obesity.sampleSize());
        assertEquals(15917, PanUkbbGwasOutcome.descriptor(obesity).get("cases"));
    }

    @Test void localChainProjectionPreservesOrComplementsAlleles()
            throws Exception {
        Path chain=directory.resolve("test.chain");
        Files.writeString(chain,"chain 1 chr1 100 + 0 100 chr1 100 + 0 100 1\n100\n\n"
            + "chain 1 chr2 100 + 0 100 chr2 100 - 0 100 2\n100\n");
        LiftOver lift=new LiftOver(chain.toFile());
        var forward=PanUkbbGwasOutcome.project(
            new PanUkbbGwasOutcome.Instrument("rs1","1",10,"T","G"),lift);
        assertEquals(10,forward.position());assertEquals("T",forward.effect());
        assertEquals("G",forward.other());
        var reverse=PanUkbbGwasOutcome.project(
            new PanUkbbGwasOutcome.Instrument("rs2","2",10,"A","C"),lift);
        assertEquals(91,reverse.position());assertEquals("T",reverse.effect());
        assertEquals("G",reverse.other());
    }

    @Test void grch38LdRequiresAnExplicitSupportedPopulation() {
        assertEquals("European", ThousandGenomesGrch38Ld.ancestry("EUR"));
        assertThrows(IllegalArgumentException.class, () -> ThousandGenomesGrch38Ld.ancestry("auto"));
    }
}
