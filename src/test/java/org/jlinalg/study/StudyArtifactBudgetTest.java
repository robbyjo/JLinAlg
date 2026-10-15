package org.jlinalg.study;

import static org.junit.jupiter.api.Assertions.*;
import java.io.RandomAccessFile;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class StudyArtifactBudgetTest {
    @TempDir Path directory;

    @Test void generatedGeneticEvidenceCanExceedContractBudgetWithoutWeakeningContracts() throws Exception {
        Path file=directory.resolve("genetic-evidence.yaml");
        String large="a".repeat(1024);
        String padding=("#"+"x".repeat(1000)+"\n").repeat(16_000);
        Files.writeString(file,padding+"pairs: []\naudit: "+large+"\n");
        assertThrows(IllegalArgumentException.class,()->StudyIO.load(file));
        assertEquals(large,StudyIO.loadArtifact(file).get("audit"));
        assertEquals(large,StudyIO.map(StudyEvidenceIndex.keyFindings(directory,"mendelian_randomization"),"findings").get("audit"));
    }

    @Test void generatedArtifactsRetainAliasDuplicateAndDepthProtections() throws Exception {
        Path file=directory.resolve("manifest.yaml");
        for(String unsafe:new String[]{"x: &a [1]\ny: *a\n","x: 1\nx: 2\n","x: "+"[".repeat(100)+"1"+"]".repeat(100)}) {
            Files.writeString(file,unsafe);
            assertThrows(IllegalArgumentException.class,()->StudyIO.loadArtifact(file));
        }
    }

    @Test void generatedArtifactByteBudgetIsBoundedBeforeParsing() throws Exception {
        Path file=directory.resolve("oversized.yaml");
        try(var random=new RandomAccessFile(file.toFile(),"rw")){random.setLength(128L*1024*1024+1);}
        var error=assertThrows(java.io.IOException.class,()->StudyIO.loadArtifact(file));
        assertTrue(error.getMessage().contains("128 MiB"));
    }
}
