/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import org.jlinalg.compute.BackendPolicy;

/** Parsed command-line configuration with deliberately dependency-free syntax. */
final class CliOptions {
    Path omics;
    Path phenotype;
    Path output;
    Path log;
    Path annotation;
    Path pedigree;
    Path grm;
    Path bgenSamples;
    String idColumn;
    String formula;
    String model = "auto";
    String family = "gaussian";
    String link;
    String ties = "efron";
    String degreesOfFreedom = "auto";
    String omicsType = "auto";
    String varianceComponents = "null-model";
    BackendPolicy backend = BackendPolicy.PREFERRED;
    String annotationId;
    String individualId;
    String pedigreeId;
    String sireId;
    String damId;
    String pedigreeFamilyId;
    String caseValue;
    String controlValue;
    String hweSamples = "all";
    String hweFilterScope = "all";
    final List<String> annotationColumns = new ArrayList<>();
    final List<String> transforms = new ArrayList<>();
    final List<Path> transformPlugins = new ArrayList<>();
    double minimumMaf;
    double minimumMac;
    double maximumMissingRate = 1.0;
    double minimumInfo = Double.NaN;
    double minimumHweP;
    int threads = Math.max(1, Runtime.getRuntime().availableProcessors());
    int blockSize;
    int checkpointEvery = 1000;
    boolean resume;
    boolean dryRun;
    boolean explain;
    boolean overwrite;
    boolean noLog;
    boolean help;
    boolean version;

    static CliOptions parse(String[] arguments) {
        CliOptions result = new CliOptions();
        for (int index = 0; index < arguments.length; index++) {
            String option = arguments[index];
            switch (option) {
                case "--help", "-h" -> result.help = true;
                case "--version" -> result.version = true;
                case "--resume" -> result.resume = true;
                case "--dry-run" -> result.dryRun = true;
                case "--explain" -> result.explain = true;
                case "--overwrite" -> result.overwrite = true;
                case "--no-log" -> result.noLog = true;
                case "--omics" -> result.omics = path(value(arguments, ++index, option));
                case "--pheno" -> result.phenotype = path(value(arguments, ++index, option));
                case "--out" -> result.output = path(value(arguments, ++index, option));
                case "--log" -> result.log = path(value(arguments, ++index, option));
                case "--annot" -> result.annotation = path(value(arguments, ++index, option));
                case "--pedigree" ->
                    result.pedigree = path(value(arguments, ++index, option));
                case "--grm" ->
                    result.grm = path(value(arguments, ++index, option));
                case "--sample-file" ->
                    result.bgenSamples = path(value(arguments, ++index, option));
                case "--id" -> result.idColumn = value(arguments, ++index, option);
                case "--formula" -> result.formula = value(arguments, ++index, option);
                case "--model" -> result.model = lower(value(arguments, ++index, option));
                case "--family" -> result.family = lower(value(arguments, ++index, option));
                case "--link" -> result.link = lower(value(arguments, ++index, option));
                case "--ties" ->
                    result.ties = lower(value(arguments, ++index, option));
                case "--df" ->
                    result.degreesOfFreedom = lower(value(arguments, ++index, option));
                case "--omics-type" ->
                    result.omicsType = lower(value(arguments, ++index, option));
                case "--variance-components" ->
                    result.varianceComponents = lower(value(arguments, ++index, option));
                case "--backend" -> result.backend = backend(
                    value(arguments, ++index, option));
                case "--annot-id" ->
                    result.annotationId = value(arguments, ++index, option);
                case "--individual-id" ->
                    result.individualId = value(arguments, ++index, option);
                case "--pedigree-id" ->
                    result.pedigreeId = value(arguments, ++index, option);
                case "--sire-id" ->
                    result.sireId = value(arguments, ++index, option);
                case "--dam-id" ->
                    result.damId = value(arguments, ++index, option);
                case "--pedigree-family-id" ->
                    result.pedigreeFamilyId = value(arguments, ++index, option);
                case "--annot-cols" -> addList(result.annotationColumns,
                    value(arguments, ++index, option));
                case "--transform" ->
                    result.transforms.add(value(arguments, ++index, option));
                case "--transform-plugin" ->
                    result.transformPlugins.add(
                        path(value(arguments, ++index, option)));
                case "--case-value" ->
                    result.caseValue = value(arguments, ++index, option);
                case "--control-value" ->
                    result.controlValue = value(arguments, ++index, option);
                case "--hwe-samples" ->
                    result.hweSamples = lower(value(arguments, ++index, option));
                case "--hwe-filter-scope" ->
                    result.hweFilterScope = lower(value(arguments, ++index, option));
                case "--min-maf" ->
                    result.minimumMaf = number(arguments, ++index, option);
                case "--min-mac" ->
                    result.minimumMac = number(arguments, ++index, option);
                case "--max-marker-missing" ->
                    result.maximumMissingRate = number(arguments, ++index, option);
                case "--min-info" ->
                    result.minimumInfo = number(arguments, ++index, option);
                case "--min-hwe-p" ->
                    result.minimumHweP = number(arguments, ++index, option);
                case "--threads" ->
                    result.threads = integer(arguments, ++index, option);
                case "--checkpoint-every" ->
                    result.checkpointEvery = integer(arguments, ++index, option);
                case "--block-size" -> {
                    String requested = lower(value(arguments, ++index, option));
                    result.blockSize = requested.equals("auto") ? 0
                        : positiveInteger(requested, option);
                }
                default -> throw new IllegalArgumentException(
                    "unknown option: " + option);
            }
        }
        return result;
    }

    void validateForRun() {
        if (phenotype == null || idColumn == null || formula == null
                || output == null)
            throw new IllegalArgumentException(
                "--pheno, --id, --formula, and --out are required");
        if (threads < 1 || checkpointEvery < 1)
            throw new IllegalArgumentException(
                "--threads and --checkpoint-every must be positive");
        if (minimumMaf < 0 || minimumMaf > 0.5 || minimumMac < 0
                || maximumMissingRate < 0 || maximumMissingRate > 1
                || minimumHweP < 0 || minimumHweP > 1)
            throw new IllegalArgumentException("filter thresholds are invalid");
        if (!varianceComponents.equals("null-model")
                && !varianceComponents.equals("refit"))
            throw new IllegalArgumentException(
                "--variance-components must be null-model or refit");
        if (!hweSamples.equals("all"))
            throw new IllegalArgumentException(
                "--hwe-samples currently supports only its default, all");
        if (!hweFilterScope.equals("all")
                && !hweFilterScope.equals("controls")
                && !hweFilterScope.equals("cases"))
            throw new IllegalArgumentException(
                "--hwe-filter-scope must be all, controls, or cases");
        if (pedigree != null && (pedigreeId == null
                || sireId == null || damId == null))
            throw new IllegalArgumentException(
                "--pedigree requires --pedigree-id, --sire-id, and --dam-id");
    }

    Path logPath() {
        return log != null ? log
            : Path.of(output.toString() + ".log");
    }

    Path manifestPath() {
        return Path.of(output.toString() + ".manifest.json");
    }

    private static Path path(String value) { return Path.of(value); }
    private static String lower(String value) {
        return value.toLowerCase(Locale.ROOT);
    }
    private static BackendPolicy backend(String value) {
        try {
            return BackendPolicy.valueOf(value.toUpperCase(Locale.ROOT));
        } catch (IllegalArgumentException failure) {
            throw new IllegalArgumentException(
                "--backend must be preferred, cholmod, gpu, cuda, opencl, "
                    + "vulkan, onemkl, openblas, auto, or cpu: " + value,
                failure);
        }
    }
    private static void addList(List<String> target, String text) {
        for (String value : text.split(",", -1)) {
            if (!value.isBlank()) target.add(value.trim());
        }
    }
    private static double number(String[] values, int index, String option) {
        String token = value(values, index, option);
        try {
            double parsed = Double.parseDouble(token);
            if (!Double.isFinite(parsed))
                throw new NumberFormatException("not finite");
            return parsed;
        } catch (NumberFormatException exception) {
            throw new IllegalArgumentException(
                option + " requires a finite number: " + token, exception);
        }
    }
    private static int integer(String[] values, int index, String option) {
        return positiveInteger(value(values, index, option), option);
    }
    private static int positiveInteger(String token, String option) {
        try {
            int parsed = Integer.parseInt(token);
            if (parsed < 1) throw new NumberFormatException("not positive");
            return parsed;
        } catch (NumberFormatException exception) {
            throw new IllegalArgumentException(
                option + " requires a positive integer: " + token, exception);
        }
    }
    private static String value(String[] values, int index, String option) {
        if (index >= values.length)
            throw new IllegalArgumentException(option + " requires a value");
        return values[index];
    }
}
