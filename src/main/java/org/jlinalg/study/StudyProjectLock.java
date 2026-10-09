/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.study;

import java.io.IOException;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.nio.channels.OverlappingFileLockException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.HashSet;
import java.util.Set;

/** Cross-process project mutation lock, reentrant only on the owning thread. */
final class StudyProjectLock {
    private StudyProjectLock() { }
    @FunctionalInterface interface Action<T> { T run() throws IOException; }
    private static final ThreadLocal<Set<Path>> HELD = ThreadLocal.withInitial(HashSet::new);

    static <T> T with(Path project, Action<T> action) throws IOException {
        Files.createDirectories(project.toAbsolutePath().normalize());
        Path root = project.toRealPath();
        if (HELD.get().contains(root)) return action.run();
        try (FileChannel channel = FileChannel.open(root.resolve(".study-run.lock"),
                StandardOpenOption.CREATE, StandardOpenOption.WRITE)) {
            try (FileLock lock = channel.tryLock()) {
                if (lock == null) throw new IOException("Another study mutation or run is active");
                HELD.get().add(root);
                try { return action.run(); }
                finally { HELD.get().remove(root); }
            } catch (OverlappingFileLockException busy) {
                throw new IOException("Another study mutation or run is active", busy);
            }
        }
    }
}
