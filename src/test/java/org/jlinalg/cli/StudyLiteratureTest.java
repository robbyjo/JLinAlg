package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.util.Map;
import org.junit.jupiter.api.Test;

class StudyLiteratureTest {
    @Test void extractsLiteralBoundedTextAndNeverInventsHypothesisSupport() {
        String abstractText = "<p>We studied <b>PCBP2</b> in tissue and found no evidence of a replicated obesity association in this population.</p>";
        var record = StudyLiterature.abstractEvidence(Map.of("abstractText", abstractText, "doi", "10.test/example"), "PCBP2", "obesity", "queries/0001");
        String normalized = abstractText.replaceAll("<[^>]*>", " ").replaceAll("\\s+", " ").trim();
        assertEquals(record.get("excerpt"), normalized.substring((Integer) record.get("excerpt_start"), (Integer) record.get("excerpt_end")));
        assertEquals(12, record.get("excerpt_word_count"));
        assertTrue(record.get("excerpt").toString().contains("no evidence"));
        assertEquals("unreviewed_candidate_not_hypothesis_support", record.get("assessment"));
        assertFalse(record.containsKey("abstractText"));
    }
    @Test void absentAbstractIsNotReconstructedFromTitle() {
        var record = StudyLiterature.abstractEvidence(Map.of("title", "A title"), "GENE", "BMI", "queries/0001");
        assertEquals("abstract_unavailable", record.get("excerpt_status")); assertFalse(record.containsKey("excerpt"));
    }
}
