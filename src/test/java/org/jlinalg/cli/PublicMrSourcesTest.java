package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

class PublicMrSourcesTest {
    @Test void bmiSelectsIndependentPanUkbEurOutcome() {
        assertTrue(PanUkbbGwasOutcome.supports("BMI", "EUR"));
        assertTrue(PanUkbbGwasOutcome.supports("body mass index", "EUR"));
        assertFalse(PanUkbbGwasOutcome.supports("obesity", "EUR"));
        assertEquals(419163, PanUkbbGwasOutcome.descriptor().get("sample_size"));
        assertEquals("GRCh37", PanUkbbGwasOutcome.descriptor().get("source_genome_build"));
    }

    @Test void ncbiSpdiProjectionPreservesTheMappedAllele() throws Exception {
        String json = "{\"data\":{\"spdis\":["
            + "{\"seq_id\":\"NC_000001.11\",\"position\":93158581,\"deleted_sequence\":\"G\",\"inserted_sequence\":\"T\"},"
            + "{\"seq_id\":\"NC_000001.10\",\"position\":93624138,\"deleted_sequence\":\"G\",\"inserted_sequence\":\"T\"}]}}";
        var input = new PanUkbbGwasOutcome.Instrument("rs2749735", "1", 93158582, "T", "G");
        var mapped = PanUkbbGwasOutcome.parseProjection(input, "NC_000001.10", "https://example.invalid", json);
        assertEquals(93624139, mapped.position());
        assertEquals("T", mapped.effect());
        assertEquals("G", mapped.other());
    }

    @Test void grch38LdRequiresAnExplicitSupportedPopulation() {
        assertEquals("European", ThousandGenomesGrch38Ld.ancestry("EUR"));
        assertThrows(IllegalArgumentException.class, () -> ThousandGenomesGrch38Ld.ancestry("auto"));
    }
}
