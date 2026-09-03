/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Path;
import java.util.Locale;

/** Commands for discovering and installing LD reference databases. */
final class LdDatabaseCli {
    private static final String DOWNLOAD_COMMAND =
        "java -jar jlinalg-<version>.jar ld-db download "
        + "--database 1000g-phase3 [--location DIRECTORY]";

    @FunctionalInterface
    interface Installer {
        void install(LdDatabaseSpec database, Path location,
            PrintStream output) throws IOException, InterruptedException;
    }

    private LdDatabaseCli() { }

    static int run(String[] arguments, PrintStream output,
            PrintStream errorOutput) {
        return run(arguments, output, errorOutput,
            Path.of("").toAbsolutePath().normalize(),
            new LdDatabaseInstaller()::install);
    }

    static int run(String[] arguments, PrintStream output,
            PrintStream errorOutput, Path currentDirectory,
            Installer installer) {
        if (arguments.length == 0) {
            error(errorOutput, "an ld-db command is required");
            guidance(errorOutput);
            return 2;
        }
        return switch (arguments[0].toLowerCase(Locale.ROOT)) {
            case "list" -> list(arguments, output, errorOutput);
            case "download" -> download(arguments, output, errorOutput,
                currentDirectory, installer);
            default -> {
                error(errorOutput, "unknown ld-db command: " + arguments[0]);
                guidance(errorOutput);
                yield 2;
            }
        };
    }

    private static int list(String[] arguments, PrintStream output,
            PrintStream errorOutput) {
        if (arguments.length != 1) {
            error(errorOutput, "ld-db list does not accept options");
            return 2;
        }
        output.println(LdDatabaseCatalog.choices());
        output.println("\nDownload with:");
        output.println("  " + DOWNLOAD_COMMAND);
        return 0;
    }

    private static int download(String[] arguments, PrintStream output,
            PrintStream errorOutput, Path currentDirectory,
            Installer installer) {
        String databaseId = null;
        String locationValue = null;
        try {
            for (int index = 1; index < arguments.length; index++) {
                String option = arguments[index];
                switch (option) {
                    case "--database" -> databaseId = value(
                        arguments, ++index, option);
                    case "--location" -> locationValue = valueAllowBlank(
                        arguments, ++index, option);
                    case "--help", "-h" -> {
                        output.println(downloadHelp());
                        return 0;
                    }
                    default -> throw new IllegalArgumentException(
                        "unknown ld-db download option: " + option);
                }
            }
        } catch (IllegalArgumentException exception) {
            error(errorOutput, exception.getMessage());
            guidance(errorOutput);
            return 2;
        }

        if (databaseId == null || databaseId.isBlank()) {
            error(errorOutput, "no LD database was specified; "
                + "--database is required");
            guidance(errorOutput);
            return 2;
        }
        LdDatabaseSpec database = LdDatabaseCatalog.find(databaseId);
        if (database == null) {
            error(errorOutput, "unknown LD database: " + databaseId);
            guidance(errorOutput);
            return 2;
        }

        Path location = locationValue == null || locationValue.isBlank()
            ? currentDirectory
            : currentDirectory.resolve(locationValue).normalize();
        try {
            installer.install(database, location.toAbsolutePath().normalize(),
                output);
            return 0;
        } catch (InterruptedException exception) {
            Thread.currentThread().interrupt();
            error(errorOutput, "LD database download was interrupted");
            return 1;
        } catch (IOException exception) {
            error(errorOutput, "failed to install " + database.id()
                + ": " + exception.getMessage());
            return 1;
        }
    }

    private static void guidance(PrintStream output) {
        output.println("Install a freely available LD database with:");
        output.println("  " + DOWNLOAD_COMMAND);
        output.println(LdDatabaseCatalog.choices());
        output.println("Omit --location, or pass an empty value, "
            + "to install in the current directory.");
    }

    private static String downloadHelp() {
        return "Usage:\n  " + DOWNLOAD_COMMAND + "\n\n"
            + LdDatabaseCatalog.choices() + "\n\n"
            + "--location defaults to the current directory.";
    }

    private static String value(String[] values, int index, String option) {
        String result = valueAllowBlank(values, index, option);
        if (result.isBlank())
            throw new IllegalArgumentException(option
                + " requires a non-blank value");
        return result;
    }

    private static String valueAllowBlank(String[] values, int index,
            String option) {
        if (index >= values.length)
            throw new IllegalArgumentException(option + " requires a value");
        return values[index];
    }

    private static void error(PrintStream output, String message) {
        output.println("jlinalg: " + message);
    }
}
