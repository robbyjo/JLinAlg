/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.Arrays;

/** Executable, bounded-memory statistical association command line. */
public final class JLinAlgCli {
    private JLinAlgCli() { }

    public static void main(String[] arguments) {
        int status = run(arguments);
        if (status != 0) System.exit(status);
    }

    static int run(String[] arguments) {
        return run(arguments, System.out, System.err);
    }

    static int run(String[] arguments, PrintStream output,
            PrintStream errorOutput) {
        if (arguments.length > 0 && arguments[0].equals("ld-db"))
            return LdDatabaseCli.run(Arrays.copyOfRange(
                arguments, 1, arguments.length), output, errorOutput);
        if (arguments.length > 0 && arguments[0].equals("mr-instruments"))
            return MrInstrumentCli.run(Arrays.copyOfRange(
                arguments, 1, arguments.length), output, errorOutput);
        if (arguments.length > 0 && arguments[0].equals("clump"))
            return LdClumpCli.run(Arrays.copyOfRange(
                arguments, 1, arguments.length), output, errorOutput);
        if (arguments.length > 0 && arguments[0].equals("mr-xwas"))
            return MrXwasCli.run(Arrays.copyOfRange(
                arguments, 1, arguments.length), output, errorOutput);
        CliOptions options;
        try {
            options = CliOptions.parse(arguments);
            if (options.help) {
                output.println(help());
                return 0;
            }
            if (options.version) {
                output.println(version());
                return 0;
            }
            options.validateForRun();
        } catch (RuntimeException exception) {
            errorOutput.println("jlinalg: " + exception.getMessage());
            errorOutput.println("Use --help for usage.");
            return 2;
        }

        RunLog log = null;
        try {
            if (!options.noLog)
                log = RunLog.open(options.logPath(), options.resume);
            info(log, "version=" + version());
            info(log, "command=" + String.join(" ", arguments));
            info(log, "java=" + System.getProperty("java.version"));
            info(log, "max_heap_bytes=" + Runtime.getRuntime().maxMemory());
            FormulaPlan plan = FormulaPlan.parse(options.formula);
            if (plan.hasOmics() && options.omics == null)
                throw new IllegalArgumentException(
                    "formula contains <omics> but --omics is absent");
            if (!plan.hasOmics() && options.omics != null)
                throw new IllegalArgumentException(
                    "--omics requires a <omics> formula term");
            if (options.pedigree != null)
                throw new IllegalArgumentException(
                    "pedigree CLI execution is not yet available in this build");
            if (options.minimumHweP > 0)
                throw new IllegalArgumentException(
                    "--min-hwe-p filtering is not yet available; HWE is reported");
            if (options.resume && Files.exists(options.output)) {
                info(log, "resume_result=already-complete");
                complete(log, "complete");
                return 0;
            }
            AnalysisRunner runner = new AnalysisRunner(options, plan, log);
            int result = runner.execute();
            complete(log, "complete");
            return result;
        } catch (Exception exception) {
            try {
                error(log, exception.getClass().getSimpleName()
                    + ": " + exception.getMessage());
                complete(log, "failed");
            } catch (IOException ignored) {
                // Preserve the original failure.
            }
            errorOutput.println("jlinalg: " + exception.getMessage());
            if (log == null) exception.printStackTrace(errorOutput);
            return 1;
        } finally {
            if (log != null) {
                try {
                    log.close();
                } catch (IOException exception) {
                    errorOutput.println(
                        "jlinalg: failed to close log: " + exception.getMessage());
                }
            }
        }
    }

    static String version() {
        String value = JLinAlgCli.class.getPackage()
            .getImplementationVersion();
        return value == null ? "development" : value;
    }

    private static void info(RunLog log, String message) throws IOException {
        if (log != null) log.info(message);
    }
    private static void error(RunLog log, String message) throws IOException {
        if (log != null) log.error(message);
    }
    private static void complete(RunLog log, String status) throws IOException {
        if (log != null) log.complete(status);
    }

    private static String help() {
        return """
            Usage:
              java -jar jlinalg-<version>.jar --pheno FILE --id COLUMN
                --formula "y ~ covariates + <omics>" [--omics FILE] --out FILE
              java -jar jlinalg-<version>.jar ld-db list
              java -jar jlinalg-<version>.jar ld-db download
                --database NAME [--location DIRECTORY]
              java -jar jlinalg-<version>.jar mr-instruments search --trait TEXT
              java -jar jlinalg-<version>.jar mr-instruments download
                --study GCST... --out FILE [--p-threshold 5e-8]
              java -jar jlinalg-<version>.jar mr-instruments format
                --input FILE --out FILE [--map TARGET=SOURCE,...]
              java -jar jlinalg-<version>.jar clump --database DIRECTORY
                --instrument FILE --ld-threshold 0.001 --output FILE
              java -jar jlinalg-<version>.jar mr-xwas --exposure FILE
                --outcome FILE --output FILE --p-threshold X

            Core options:
              --omics FILE                 CSV/TSV, VCF, BCF, or BGEN matrix
              --pheno FILE                 Observation-by-variable CSV/TSV
              --id COLUMN                  Phenotype sample-ID column
              --formula FORMULA            R-style fixed/random formula
              --model auto|ols|lmm|glm|glmm|cox
              --family gaussian|binomial|poisson|gamma
              --grm FILE|PREFIX           Labeled dense matrix or GCTA prefix
              --individual-id COLUMN      Phenotype-to-GRM ID (defaults to --id)
              --ties efron|breslow         Cox tied-event method
              --df auto|satterth|kr
              --out FILE                   Log defaults to FILE.log

            Streaming and reproducibility:
              --block-size auto|N           Default uses JVM heap headroom
              --threads N
              --backend POLICY              preferred (default), cholmod, gpu,
                                            cuda, opencl, vulkan, onemkl,
                                            openblas, auto, or cpu
              --resume
              --checkpoint-every N
              --log FILE
              --no-log
              --dry-run
              --explain
              --overwrite

            Filtering and processing:
              --min-maf X --min-mac X --max-marker-missing X --min-info X
              --transform "<omics>=winsor(p=.01)|log1p()|zscore()"
              --transform-plugin JAR       Trusted Java transform provider
              --omics-type auto|gwas|ewas|expression|proteomics|generic
              --annot FILE --annot-id COLUMN --annot-cols c1,c2,...
              --case-value VALUE --control-value VALUE

            LD reference databases:
              ld-db list                  Show freely available choices
              ld-db download              Download and install a choice
              --database NAME             Required database identifier
              --location DIRECTORY        Install directory; default is .

            MR instrument preparation:
              mr-instruments search       Find downloadable GWAS by trait
              mr-instruments download     Stream significant MR candidates
              mr-instruments format       Normalize user GWAS/QTL columns

            LD clumping:
              clump --database DIRECTORY --instrument FILE --output FILE
              --population PANEL          AFR, AMR, EAS, EUR (default), or SAS
              --clump-kb N                Window in kilobases; default 10000
              --ld-threshold X            LD r-squared cutoff; default 0.001
              --p-threshold X             Index p-value cutoff; default 1
              --overwrite                 Replace an existing output file

            Parallel xWAS MR:
              mr-xwas --exposure FILE --outcome FILE --output FILE
              --p-threshold X             Retain screening p <= X
              --threads N                 Bounded exposure-outcome workers
              --pair-block-size N         Maximum resident pair evaluations
              Run mr-xwas --help for log-scale thresholds and table columns.
            """;
    }
}
