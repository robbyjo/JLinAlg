/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

/** Presentation only. Never use this formatter for machine-readable result exports. */
final class StudyReportNumbers {
    private StudyReportNumbers() { }
    static String field(String key, Object value) {
        return value instanceof Number || java.util.Set.of("beta", "se", "p_value", "BH_FDR", "q_value", "fdr", "FDR",
            "i_squared", "mean", "sd", "percent", "male_percent", "odds_ratio", "fold_enrichment_unweighted",
            "adjusted_p_value", "precision", "recall", "cohens_d", "heterogeneity_p_value").contains(key)
                ? format(value) : String.valueOf(value == null ? "Not reported" : value);
    }
    static String format(Object value) {
        if (value == null) return "Not reported";
        String raw = value.toString();
        if (!raw.matches("[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?")) return raw;
        // Integer counts, coordinates and identifiers retain all digits.
        if (!raw.contains(".") && !raw.toLowerCase(Locale.ROOT).contains("e")) return raw;
        BigDecimal number = new BigDecimal(raw), magnitude = number.abs();
        if (number.signum() != 0 && (magnitude.compareTo(new BigDecimal("0.01")) < 0
                || magnitude.compareTo(new BigDecimal("1000000")) >= 0)) {
            DecimalFormat scientific = new DecimalFormat("0.00E0", DecimalFormatSymbols.getInstance(Locale.ROOT));
            scientific.setRoundingMode(RoundingMode.HALF_UP);
            return scientific.format(number).replace('E', 'e');
        }
        return number.setScale(2, RoundingMode.HALF_UP).toPlainString();
    }
}
