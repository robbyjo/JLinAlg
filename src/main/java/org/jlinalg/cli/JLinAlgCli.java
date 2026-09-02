/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.IOException;
import java.nio.file.Files;

/** Executable, bounded-memory statistical association command line. */
public final class JLinAlgCli {
    private JLinAlgCli() { }

    public static void main(String[] arguments) {
        int status = run(arguments);
        if (status != 0) System.exit(status);
    }

    static int run(String[] arguments) {
        CliOptions options;
        try {
            options = CliOptions.parse(arguments);
            if (options.help) {
                System.out.println(help());
                return 0;
            }
            if (options.version) {
                System.out.println(version());
                return 0;
            }
            options.validateForRun();
        } catch (RuntimeException exception) {
            System.err.println("jlinalg: " + exception.getMessage());
            System.err.println("Use --help for usage.");
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
            System.err.println("jlinalg: " + exception.getMessage());
            if (log == null) exception.printStackTrace(System.err);
            return 1;
        } finally {
            if (log != null) {
                try {
                    log.close();
                } catch (IOException exception) {
                    System.err.println(
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
            """;
    }
}
