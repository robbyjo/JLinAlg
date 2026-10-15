package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;

import java.io.StringReader;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.io.TempDir;
import org.junit.jupiter.api.Test;

class StudyYamlTest {
    @TempDir Path directory;

    @Test void generatedEvidenceHasItsOwnBoundedBudgetAndOrdinaryContractsRemainStrict() throws Exception {
        Path file=directory.resolve("genetic-evidence.yaml");
        // Many short comment lines exercise the byte/code-point budget without
        // making this regression depend on scanning one enormous YAML scalar.
        Files.writeString(file,("#"+"x".repeat(1000)+"\n").repeat(16_000)+"pairs: []\n");
        assertThrows(RuntimeException.class,()->StudyYaml.loadMap(file));
        assertEquals(java.util.List.of(),StudyYaml.loadArtifact(file).get("pairs"));
        Files.writeString(file,"first: &x [1]\nsecond: *x\n");
        assertThrows(RuntimeException.class,()->StudyYaml.loadArtifact(file));
        Files.writeString(file,"value: 1\nvalue: 2\n");
        assertThrows(RuntimeException.class,()->StudyYaml.loadArtifact(file));
        try(var random=new java.io.RandomAccessFile(file.toFile(),"rw")){random.setLength(128L*1024*1024+1);}
        assertThrows(java.io.IOException.class,()->StudyYaml.loadArtifact(file));
    }
    @Test void readsLargeGeneratedEvidenceWithinTheBoundedStudyLimit() {
        String payload = "x".repeat(3_300_000);
        Object parsed = StudyYaml.load(new StringReader("evidence: " + payload));
        assertEquals(payload, StudyPublicEvidence.map(parsed).get("evidence"));
    }

    @Test void sharedReportingReaderAcceptsExpandedDependencyManifests() throws Exception {
        Path file=directory.resolve("manifest.yaml");
        Files.writeString(file,("#"+"x".repeat(1000)+"\n").repeat(16_000)
            +"operation_type: mendelian_randomization\nstate: complete\n");
        assertThrows(RuntimeException.class,()->StudyYaml.loadMap(file));
        var method=StudyCliExecutor.class.getDeclaredMethod("loadYaml",Path.class);
        method.setAccessible(true);
        var parsed=(java.util.Map<?,?>)method.invoke(null,file);
        assertEquals("mendelian_randomization",parsed.get("operation_type"));
        assertEquals("complete",parsed.get("state"));
        Files.writeString(file,"state: complete\nstate: failed\n");
        assertThrows(java.lang.reflect.InvocationTargetException.class,()->method.invoke(null,file));
    }

    @Test void stillRejectsDuplicateKeysAndAliases() {
        assertThrows(RuntimeException.class,
            () -> StudyYaml.load("value: 1\nvalue: 2\n"));
        assertThrows(RuntimeException.class,
            () -> StudyYaml.load("first: &x [1]\nsecond: *x\n"));
    }
}
