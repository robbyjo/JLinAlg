/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.study;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

/** Local evidence packaging, explicitly not a de-identification or submission service. */
public final class StudyEvidenceArchive {
    private StudyEvidenceArchive() { }

    public static String export(Path project, Path output, StudyProject.Actor actor,
            boolean acknowledgeSensitive) throws IOException {
        if (!acknowledgeSensitive) throw new IllegalArgumentException("Evidence, notebook text and intermediate outputs may be sensitive. Review locally and use --acknowledge-sensitive; export does not authorize transmission or journal/regulatory submission.");
        return StudyProjectLock.with(project, () -> exportLocked(project.toAbsolutePath().normalize(), output.toAbsolutePath().normalize(), actor));
    }

    private static String exportLocked(Path root, Path output, StudyProject.Actor actor) throws IOException {
        StudyProject.verify(root);
        if (Files.exists(output)) throw new IOException("Archive output already exists; choose a new filename");
        Set<String> directories = Set.of("ledger", "state", "plans", "amendments", "observations", "questions", "answers", "proposals", "runs", "index");
        Set<String> files = Set.of("study.yaml", "resolved-study.yaml", "SESSION_HISTORY.md", "QUESTIONS.md", "AGENTS.md");
        for (String directory : directories) if (output.startsWith(root.resolve(directory)))
            throw new IOException("Archive cannot be written inside an authoritative evidence directory");
        List<Path> selected = new ArrayList<>();
        for (String name : files) if (Files.isRegularFile(root.resolve(name), LinkOption.NOFOLLOW_LINKS)) selected.add(root.resolve(name));
        for (String name : directories) if (Files.isDirectory(root.resolve(name), LinkOption.NOFOLLOW_LINKS)) {
            try (var paths = Files.walk(root.resolve(name))) {
                for (Path path : paths.sorted().toList()) {
                    if (Files.isSymbolicLink(path)) throw new IOException("Evidence archive refuses symbolic links: " + path);
                    if (Files.isRegularFile(path, LinkOption.NOFOLLOW_LINKS) && !path.toString().contains(".operation-")) selected.add(path);
                }
            }
        }
        selected.sort(Path::compareTo);
        Map<String, Object> manifest = new LinkedHashMap<>();
        for (Path path : selected) manifest.put(portable(root.relativize(path)), StudyIO.hash(path));
        Map<String, Object> limitations = new LinkedHashMap<>();
        limitations.put("scope", "Local evidence snapshot; not automatically submitted anywhere");
        limitations.put("privacy", "May contain sensitive result files, intermediate analysis data, scientist answers, notebook entries and LLM interactions. Not certified anonymous. Human disclosure review remains required.");
        limitations.put("omitted", "Original dataset files outside the selected evidence directories, reference caches and software binaries are not bundled. Input hashes and source records identify reconstruction prerequisites. Include original scripts as declared operation artifacts before exporting.");
        limitations.put("regulatory_status", "No assertion of FDA compliance, validated computerized-system status, digital signature or tamper-proof storage. Hashes support integrity checking only.");
        limitations.put("files_sha256", manifest);
        Files.createDirectories(output.getParent());
        Path stage = Files.createTempFile(output.getParent(), ".evidence-archive-", ".zip");
        try {
            try (ZipOutputStream zip = new ZipOutputStream(Files.newOutputStream(stage), StandardCharsets.UTF_8)) {
                for (Path path : selected) {
                    String name = portable(root.relativize(path));
                    ZipEntry entry = new ZipEntry(name); entry.setTime(0); zip.putNextEntry(entry);
                    Files.copy(path, zip); zip.closeEntry();
                    if (!StudyIO.hash(path).equals(manifest.get(name))) throw new IOException("Evidence changed during export: " + name);
                }
                ZipEntry entry = new ZipEntry("EVIDENCE_BUNDLE.json"); entry.setTime(0); zip.putNextEntry(entry);
                zip.write(StudyIO.canonicalJson(limitations).getBytes(StandardCharsets.UTF_8)); zip.closeEntry();
            }
            StudyProject.verify(root);
            Files.move(stage, output);
        } finally { Files.deleteIfExists(stage); }
        String hash = StudyIO.hash(output);
        StudyProject.append(root, "evidence_archive_exported", actor, Map.of("summary", "Exported local evidence snapshot; sharing remains subject to review",
            "archive", output.toString(), "archive_sha256", hash, "files", selected.size()));
        StudyEvidenceIndex.rebuild(root);
        return hash;
    }
    private static String portable(Path path) { return path.toString().replace('\\', '/'); }
}
