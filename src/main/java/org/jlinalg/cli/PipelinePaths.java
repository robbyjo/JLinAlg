package org.jlinalg.cli;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/** Preflight safeguards run before logs or result scratch files are opened. */
final class PipelinePaths {
    private PipelinePaths() { }

    /** Commands without an overwrite option must never truncate an existing file. */
    static void requireFreshOutputs(Path... paths) throws IOException {
        List<Path> seen = new ArrayList<>();
        for (Path path : paths) {
            if (path == null) continue;
            if (Files.exists(path)) throw new IOException("output exists; choose a new output path: " + path);
            for (Path previous : seen)
                if (same(path,previous)) throw new IOException("output paths must be distinct: " + path);
            seen.add(path);
        }
    }

    static void validate(CliOptions options) {
        List<Path> inputs = new ArrayList<>(Arrays.asList(options.phenotype, options.omics,
            options.annotation, options.pedigree, options.grm, options.bgenSamples));
        inputs.removeIf(java.util.Objects::isNull);
        inputs.addAll(options.transformPlugins);
        if (options.grm != null) {
            String value = options.grm.toString();
            String lower = value.toLowerCase(java.util.Locale.ROOT);
            if (lower.endsWith(".grm.bin")) value = value.substring(0,value.length()-8);
            else if (lower.endsWith(".grm.id")) value = value.substring(0,value.length()-7);
            inputs.add(Path.of(value + ".grm.bin"));
            inputs.add(Path.of(value + ".grm.id"));
        }
        List<Path> outputs = new ArrayList<>(List.of(options.output,
            Path.of(options.output + ".partial"), options.manifestPath()));
        if (!options.noLog) outputs.add(options.logPath());
        try {
            for (int index=0; index<outputs.size(); index++) {
                Path output = outputs.get(index);
                for (Path input : inputs)
                    if (same(output,input)) throw new IllegalArgumentException(
                        "output/log path aliases an input: " + output);
                for (int previous=0; previous<index; previous++)
                    if (same(output,outputs.get(previous))) throw new IllegalArgumentException(
                        "output, partial, manifest, and log paths must be distinct: " + output);
                if (!options.overwrite && !options.resume && Files.exists(output))
                    throw new IllegalArgumentException("output/log exists; use --overwrite: " + output);
            }
            if (options.resume && Files.exists(options.output))
                throw new IllegalArgumentException("cannot verify an existing result for --resume; "
                    + "use a new --out or --overwrite without --resume to rerun");
        } catch (IOException exception) {
            throw new IllegalArgumentException("cannot validate input/output paths",exception);
        }
    }

    private static boolean same(Path left, Path right) throws IOException {
        if (Files.exists(left) && Files.exists(right) && Files.isSameFile(left,right)) return true;
        return canonical(left).equals(canonical(right));
    }

    private static Path canonical(Path path) throws IOException {
        Path absolute = path.toAbsolutePath().normalize();
        if (Files.exists(absolute)) return absolute.toRealPath();
        Path parent = absolute.getParent();
        return parent == null ? absolute : canonical(parent).resolve(absolute.getFileName());
    }
}
