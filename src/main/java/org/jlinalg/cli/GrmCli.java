/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.genetics.GenomicRelationshipMatrix;
import org.jlinalg.genetics.GenomicRelationshipOptions;
import org.jlinalg.pipeline.DataFormat;
import org.jlinalg.pipeline.VariantSource;
import org.jlinalg.pipeline.VariantSources;

/** Block-streamed additive GRM construction with labeled, reusable output. */
final class GrmCli {
    private GrmCli() { }

    static int run(String[] args, PrintStream console, PrintStream error) {
        try {
            if (Arrays.asList(args).contains("--help")) { console.println(help()); return 0; }
            Map<String, String> options = new HashMap<>();
            for (int i = 0; i < args.length; i++) {
                String key = args[i];
                if (!Set.of("--genotypes", "--out", "--maf", "--call-rate",
                        "--block-size", "--backend", "--bgen-samples").contains(key)
                        || i + 1 == args.length || options.put(key, args[++i]) != null)
                    throw new IllegalArgumentException("unknown, duplicate, or incomplete option: " + key);
            }
            for (String required : List.of("--genotypes", "--out"))
                if (!options.containsKey(required))
                    throw new IllegalArgumentException("required option: " + required);
            Path input = Path.of(options.get("--genotypes")).toAbsolutePath().normalize();
            Path output = Path.of(options.get("--out")).toAbsolutePath().normalize();
            String filename = output.getFileName().toString().toLowerCase(Locale.ROOT);
            if (!filename.endsWith(".tsv") && !filename.endsWith(".csv"))
                throw new IllegalArgumentException("GRM output must end in .tsv or .csv");
            var filters = new GenomicRelationshipOptions(
                Double.parseDouble(options.getOrDefault("--maf", "0.01")),
                Double.parseDouble(options.getOrDefault("--call-rate", "0.95")));
            int blockSize = Integer.parseInt(options.getOrDefault("--block-size", "256"));
            if (blockSize < 1) throw new IllegalArgumentException("block size must be positive");
            BackendPolicy backend = BackendPolicy.valueOf(
                options.getOrDefault("--backend", "cpu").toUpperCase(Locale.ROOT));
            DataFormat format = DataFormat.infer(input);
            if (options.containsKey("--bgen-samples") && format != DataFormat.BGEN)
                throw new IllegalArgumentException("--bgen-samples requires BGEN input");
            Path log = Path.of(output + ".log");
            PipelinePaths.requireFreshOutputs(output, log);
            try (RunLog journal = RunLog.openPlain(log, false)) {
                journal.metadata("command=grm\ninput=" + input + "\nformat=" + format
                    + "\nmaf=" + filters.minimumMinorAlleleFrequency()
                    + "\ncall_rate=" + filters.minimumCallRate() + "\nblock_size=" + blockSize
                    + "\nbackend_requested=" + backend + "\nmissing=mean imputation"
                    + "\nnormalization=mean standardized dosage outer product; common variant denominator\n");
                List<String> external = options.containsKey("--bgen-samples")
                    ? SampleFiles.read(Path.of(options.get("--bgen-samples"))) : null;
                if (external != null) journal.info("bgen_samples=" + options.get("--bgen-samples"));
                VariantSource source = VariantSources.open(input, format, external);
                int n = source.metadata().sampleIds().size();
                // Labeled output must round-trip through the existing --grm reader.
                for (String id : source.metadata().sampleIds())
                    if (!id.equals(id.trim()) || id.contains("\n") || id.contains("\r"))
                        throw new IllegalArgumentException("sample IDs must be single-line without surrounding spaces");
                journal.info("samples=" + n + "\ndense_matrix_bytes=" + Math.multiplyExact(8L, (long) n * n));
                GenomicRelationshipMatrix grm = GenomicRelationshipMatrix.fromSource(
                    source, filters, backend, blockSize);
                journal.metadata("variants_considered=" + grm.variantsConsidered()
                    + "\nvariants_used=" + grm.variantsUsed() + "\nvariants_excluded=" + grm.variantsExcluded()
                    + "\nbackend=" + grm.computationBackend().orElseThrow() + "\n");
                Path temporary = Files.createTempFile(output.getParent(), ".grm-", ".tmp");
                try {
                    write(grm, temporary, filename.endsWith(".csv") ? ',' : '\t');
                    Files.move(temporary, output);
                } finally { Files.deleteIfExists(temporary); }
                journal.complete("complete");
                console.println("GRM: " + output + "; samples=" + n + "; variants=" + grm.variantsUsed());
            }
            return 0;
        } catch (IOException | RuntimeException exception) {
            error.println("jlinalg: " + exception.getMessage());
            return 2;
        }
    }

    private static void write(GenomicRelationshipMatrix grm, Path path, char delimiter) throws IOException {
        List<String> ids = grm.sampleIds();
        String rowLabel = "sample_id";
        while (ids.contains(rowLabel)) rowLabel = "_" + rowLabel;
        try (BufferedWriter writer = Files.newBufferedWriter(path)) {
            writer.write(quote(rowLabel));
            for (String id : ids) { writer.write(delimiter); writer.write(quote(id)); }
            writer.newLine();
            for (String row : ids) {
                writer.write(quote(row));
                for (String column : ids) {
                    writer.write(delimiter);
                    writer.write(Double.toString(grm.relationship(row, column)));
                }
                writer.newLine();
            }
        }
    }

    private static String quote(String value) { return "\"" + value.replace("\"", "\"\"") + "\""; }

    static String help() {
        return """
            Usage: jlinalg grm --genotypes FILE --out grm.tsv
              [--maf 0.01] [--call-rate 0.95] [--block-size 256] [--backend cpu]
              [--bgen-samples cohort.sample]
            Input: VCF, VCF.gz, BCF, BGEN, or variant-by-sample dosage CSV/TSV.
            Output: labeled dense TSV/CSV usable directly with --grm, plus OUT.log.
            Uses all source samples, per-variant sqrt(2p(1-p)) standardization,
            mean-imputed missing dosages, and one common retained-variant denominator.
            Monomorphic variants are excluded. Supply LD-pruned diploid variants;
            this command does not prune LD or construct LOCO matrices.
            Genotypes are block streamed; the dense GRM needs O(samples^2) memory.
            Existing outputs are never overwritten. See docs/grm-cli.md.
            """;
    }
}
