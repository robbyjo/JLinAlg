package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;

import java.io.StringReader;
import org.junit.jupiter.api.Test;

class StudyYamlTest {
    @Test void readsLargeGeneratedEvidenceWithinTheBoundedStudyLimit() {
        String payload = "x".repeat(3_300_000);
        Object parsed = StudyYaml.load(new StringReader("evidence: " + payload));
        assertEquals(payload, StudyPublicEvidence.map(parsed).get("evidence"));
    }

    @Test void stillRejectsDuplicateKeysAndAliases() {
        assertThrows(RuntimeException.class,
            () -> StudyYaml.load("value: 1\nvalue: 2\n"));
        assertThrows(RuntimeException.class,
            () -> StudyYaml.load("first: &x [1]\nsecond: *x\n"));
    }
}
