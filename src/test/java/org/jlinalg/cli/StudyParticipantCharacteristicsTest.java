package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import org.jlinalg.study.StudyOperationExecutor.Context;
import org.jlinalg.study.StudyOperationExecutor.State;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class StudyParticipantCharacteristicsTest {
    @TempDir Path project;
    private Context context(Map<String, Object> settings) throws Exception {
        Files.writeString(project.resolve("phenotype.csv"), "id,Sex,Age,BMI\nPRIVATE_ID_A,1,20,22\nPRIVATE_ID_B,2,40,NA\nPRIVATE_ID_C,NA,NA,30\n");
        var dataset = Map.<String, Object>of("id", "d", "cohort", "c", "files", Map.of("phenotype", "phenotype.csv", "membership", "members.csv"),
            "columns", Map.of("characteristics", settings));
        return new Context(project, Files.createTempDirectory(project, "output-"), Map.of("type", "participant_characteristics", "hypothesis", "h", "datasets", List.of("d")),
            Map.of("datasets", List.of(dataset), "hypotheses", List.of(Map.of("id", "h", "outcome", "BMI", "description", "A scientific question"))), Map.of(), Map.of());
    }
    @Test void availableCasesUseSampleSdAndCorrectSexDenominatorWithoutLeakingIds() throws Exception {
        Context c = context(Map.of("id", "id", "male_value", 1, "female_value", 2));
        var result = new StudyParticipantCharacteristics().execute(c);
        assertEquals(State.WAITING_FOR_INPUT, result.state());
        String yaml = Files.readString(c.outputDirectory().resolve("participant-characteristics.yaml"));
        assertFalse(yaml.contains("PRIVATE_ID")); assertTrue(yaml.contains("male_percent: 50.0"));
        assertTrue(yaml.contains("14.142135623730951")); assertTrue(yaml.contains("mean: 26.0"));
        assertTrue(yaml.contains("analysis_population_verified: false"));
    }
    @Test void membershipAndDeclarationAreRequiredForPublicationReadyCounts() throws Exception {
        Files.writeString(project.resolve("members.csv"), "id\nPRIVATE_ID_A\nPRIVATE_ID_B\n");
        Context c = context(Map.of("id", "id", "male_value", 1, "female_value", 2, "membership_role", "membership", "population_rationale", "Exact analysis design-matrix sample IDs"));
        var result = new StudyParticipantCharacteristics().execute(c);
        assertEquals(State.COMPLETE, result.state(), result.details().toString());
        assertTrue(result.details().toString().contains("N=2"));
        assertTrue(result.details().toString().contains("analysis_population_verified=true"));
    }
    @Test void missingSdAndUnknownSexDoNotBecomeZeros() throws Exception {
        Context c = context(Map.of("id", "id", "filters", Map.of("id", List.of("PRIVATE_ID_C"))));
        var result = new StudyParticipantCharacteristics().execute(c);
        assertEquals(State.WAITING_FOR_INPUT, result.state());
        assertTrue(result.details().toString().contains("not_estimable"));
        assertFalse(result.details().toString().contains("PRIVATE_ID_C"));
    }
    @Test void duplicateParticipantsAndNonfiniteNumbersAreRejected() throws Exception {
        Context c = context(Map.of("id", "id"));
        Files.writeString(project.resolve("phenotype.csv"), "id,Age,BMI\nx,20,20\nx,30,30\n");
        var result = new StudyParticipantCharacteristics().execute(c);
        assertEquals(State.WAITING_FOR_INPUT, result.state());
        assertTrue(result.details().toString().contains("unique"));
        Files.writeString(project.resolve("phenotype.csv"), "id,Age,BMI\nx,Infinity,20\n");
        result = new StudyParticipantCharacteristics().execute(c);
        assertTrue(result.details().toString().contains("Nonfinite"));
    }
    @Test void explicitObesityThresholdIsStrictAndExcludesMissingBmi() throws Exception {
        Context c = context(Map.of("id", "id", "derived_binary", Map.of("Obesity", Map.of("column", "BMI", "operator", "gt", "threshold", 30, "definition", "BMI > 30"))));
        Files.writeString(project.resolve("phenotype.csv"), "id,Sex,Age,BMI\na,1,20,30\nb,2,30,30.01\nc,1,40,29\nd,2,50,NA\n");
        var result = new StudyParticipantCharacteristics().execute(c);
        var tables = StudyPublicEvidence.list(result.details().get("tables"));
        var summary = StudyPublicEvidence.map(tables.get(0));
        var obesity = StudyPublicEvidence.map(StudyPublicEvidence.map(summary.get("binary")).get("Obesity"));
        assertEquals(1, obesity.get("case_n")); assertEquals(3, obesity.get("denominator")); assertEquals(1, obesity.get("missing"));
        assertEquals(100.0 / 3, (Double) obesity.get("percent"));
        assertTrue(summary.containsKey("selection_sha256"));
    }
    @Test void invalidBinaryRuleDoesNotBecomeAFalseZero() throws Exception {
        Context c = context(Map.of("derived_binary", Map.of("Obesity", Map.of("column", "BMI", "operator", ">", "threshold", 30, "definition", "BMI > 30"))));
        var result = new StudyParticipantCharacteristics().execute(c);
        assertEquals(State.WAITING_FOR_INPUT, result.state());
        assertTrue(result.details().toString().contains("derived_binary requires"));
    }
}
