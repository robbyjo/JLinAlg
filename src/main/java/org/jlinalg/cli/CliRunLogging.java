/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.*;
import java.util.*;
import java.util.function.ToIntFunction;

/** Run lifecycle for subcommands that do not have a command-specific log. */
final class CliRunLogging {
    private CliRunLogging() { }
    private static final Set<String> COMMANDS = Set.of("ld-db", "mr-instruments", "clump",
        "mr-xwas", "mr-estimate", "beta-regression", "penalized-regression",
        "ldsc", "twas", "pwas", "genomic-factor", "score-train", "score-apply");

    static boolean accepts(String command) { return COMMANDS.contains(command); }

    static int run(String[] arguments, PrintStream error, ToIntFunction<String[]> action) {
        List<String> args = new ArrayList<>();
        Path explicitLog = null, output = null;
        boolean disabled = false, overwrite = false;
        try {
            if (Arrays.stream(arguments).anyMatch(a -> Set.of("--help", "-h", "--version").contains(a))
                    || arguments.length == 2 && arguments[1].equals("help"))
                return action.applyAsInt(arguments);
            for (int i = 0; i < arguments.length; i++) {
                String arg = arguments[i];
                if (arg.equals("--no-log")) { disabled = true; continue; }
                if (arg.equals("--log")) {
                    if (explicitLog != null || i+1 == arguments.length || arguments[i+1].startsWith("--"))
                        throw new IllegalArgumentException("--log requires one file path");
                    explicitLog = Path.of(arguments[++i]); continue;
                }
                args.add(arg);
                if (arg.equals("--overwrite")) overwrite = true;
                if ((arg.equals("--out") || arg.equals("--output")) && i+1 < arguments.length
                        && !arguments[i+1].startsWith("--")) output = Path.of(arguments[i+1]);
            }
            String[] forwarded = args.toArray(String[]::new);
            if (disabled) {
                if (explicitLog != null) throw new IllegalArgumentException("--log and --no-log cannot be combined");
                return action.applyAsInt(forwarded);
            }
            Path log = explicitLog != null ? explicitLog : output != null ? Path.of(output + ".log")
                : Path.of("logs", "jlinalg-" + arguments[0] + "-" + java.time.Instant.now().toEpochMilli()
                    + "-" + UUID.randomUUID() + ".log");
            // Never append timing data to an input or another requested output,
            // including aliases, even when the command allows --overwrite.
            for (String value : forwarded) {
                if (value.startsWith("--")) continue;
                Path candidate;
                try { candidate = Path.of(value); }
                catch (InvalidPathException ignored) { continue; }
                if (log.toAbsolutePath().normalize().equals(candidate.toAbsolutePath().normalize())
                        || Files.exists(log) && Files.exists(candidate) && Files.isSameFile(log,candidate))
                    throw new IllegalArgumentException("log path aliases a command input/output: " + log);
            }
            if (!overwrite) PipelinePaths.requireFreshOutputs(log);
            try (RunLog journal = RunLog.openPlain(log, overwrite)) {
                journal.info("command=" + String.join(" ", forwarded));
                journal.info("version=" + JLinAlgCli.version());
                // Keep log-location notices out of machine-readable stdout.
                error.println("Log: " + log.toAbsolutePath());
                int status = action.applyAsInt(forwarded);
                journal.info("exit_code=" + status);
                journal.complete(status == 0 ? "complete" : "failed");
                return status;
            }
        } catch (IOException | RuntimeException failure) {
            error.println("jlinalg: " + failure.getMessage());
            return 2;
        }
    }
}
