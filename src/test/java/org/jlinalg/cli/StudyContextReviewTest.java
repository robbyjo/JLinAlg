package org.jlinalg.cli;
import static org.junit.jupiter.api.Assertions.*;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;

class StudyContextReviewTest {
    @Test void tissueAndModelContrastsAreCandidatesNotCausesAndPrivateFieldsStayLocal() {
        var study = Map.<String, Object>of("datasets", List.of(
            Map.of("id", "a", "cohort", "fhs", "tissue", "whole_blood", "analysis_description", Map.of("model", "linear", "covariates", List.of("age", "sex")), "files", Map.of("phenotype", "PATIENT_SECRET.csv")),
            Map.of("id", "b", "cohort", "jhs", "tissue", "peripheral_blood_mononuclear_cells", "analysis_description", Map.of("model", "mixed", "relatedness", "pedigree"), "donor", "PRIVATE_PERSON", "metadata", Map.of("dictionary", "PRIVATE_TEXT")),
            Map.of("id", "c", "cohort", "other", "tissue", "brain")),
            "hypotheses", List.of(Map.of("id", "H1", "description", "PRIVATE_HYPOTHESIS", "inputs", Map.of("discovery", List.of("a"), "replication", List.of("b")), "outcome", "BMI")));
        var inventory = StudyContextReview.inventory(study, "H1"); String publicJson = StudyAssistance.json(StudyContextReview.publicContext(inventory));
        assertTrue(publicJson.contains("whole_blood")); assertTrue(publicJson.contains("peripheral_blood_mononuclear_cells"));
        assertTrue(publicJson.contains("analysis/model")); assertTrue(publicJson.contains("not a demonstrated cause"));
        for (String secret : List.of("PATIENT_SECRET", "PRIVATE_PERSON", "PRIVATE_TEXT", "PRIVATE_HYPOTHESIS", "brain")) assertFalse(publicJson.contains(secret));
        assertTrue(inventory.get("local_omitted_fields").toString().contains("description"));
        assertTrue(StudyContextReview.instructions().contains("ALL observed results"));
        assertTrue(StudyContextReview.instructions().contains("confounded with cohort"));
    }
    @Test void unknownMetadataDoesNotBecomeAnEqualOrFabricatedComparison() {
        var study = Map.<String, Object>of("datasets", List.of(Map.of("id", "a", "tissue", "whole_blood"), Map.of("id", "b", "files", Map.of("matrix", "secret"))));
        assertEquals(List.of(), StudyContextReview.publicContext(StudyContextReview.inventory(study, null)).get("declared_contrasts"));
    }
    @Test void unfamiliarRootContextRemainsVisibleForLocalReviewWithoutDisclosingValues() {
        var inventory = StudyContextReview.inventory(Map.of("project", Map.of("description", "PRIVATE_PROJECT"),
            "custom_assumptions", List.of("PRIVATE_ASSUMPTION"), "datasets", List.of()), null);
        String local = inventory.get("local_omitted_fields").toString();
        String disclosed = StudyAssistance.json(StudyContextReview.publicContext(inventory));
        assertTrue(local.contains("custom_assumptions")); assertTrue(local.contains("project"));
        assertFalse(local.contains("PRIVATE_")); assertFalse(disclosed.contains("PRIVATE_"));
        assertEquals(2, StudyContextReview.publicContext(inventory).get("unreviewed_field_count"));
    }
}
