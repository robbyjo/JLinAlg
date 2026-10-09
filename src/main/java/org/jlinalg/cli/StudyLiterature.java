/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import static org.jlinalg.cli.StudyPublicEvidence.*;

/** Source-grounded excerpts, kept distinct from interpretation or endorsement. */
final class StudyLiterature {
    private StudyLiterature() { }

    static Map<String, Object> abstractEvidence(Map<String, Object> article, String gene, String phenotype, String audit) {
        Map<String, Object> result = new LinkedHashMap<>();
        for (String key : List.of("id", "source", "pmid", "pmcid", "doi", "title", "authorString",
                "pubYear", "journalInfo", "authorList", "pageInfo", "pubTypeList", "isOpenAccess", "commentCorrectionList"))
            if (article.containsKey(key)) result.put(key, article.get(key));
        result.put("source_audit", audit);
        result.put("assessment", "unreviewed_candidate_not_hypothesis_support");
        String source = text(article.get("abstractText"));
        // Keep raw provider markup in the response cache. Offsets below refer to
        // this explicitly normalized text, never to an unavailable full article.
        String plain = source.replaceAll("<[^>]*>", " ").replaceAll("\\s+", " ").trim();
        if (plain.isEmpty()) {
            result.put("excerpt_status", "abstract_unavailable");
            return result;
        }
        String lower = plain.toLowerCase(Locale.ROOT);
        int start = lower.indexOf(gene.toLowerCase(Locale.ROOT));
        if (start < 0) start = lower.indexOf(phenotype.toLowerCase(Locale.ROOT));
        if (start < 0) start = 0;
        // Never begin inside a word and never output more than twelve words.
        while (start > 0 && !Character.isWhitespace(plain.charAt(start - 1))) start--;
        Matcher words = Pattern.compile("\\S+").matcher(plain);
        words.region(start, plain.length());
        int end = start, count = 0;
        while (count < 12 && words.find()) { end = words.end(); count++; }
        String quote = plain.substring(start, end);
        result.put("excerpt", quote);
        result.put("excerpt_status", "verified_literal_substring_of_normalized_abstract");
        result.put("excerpt_word_count", count);
        result.put("excerpt_start", start); result.put("excerpt_end", end);
        result.put("excerpt_normalization", "Replace markup tags with spaces; collapse whitespace; trim");
        result.put("excerpt_scope", "Abstract excerpt only; may truncate a sentence and omit qualifications. Read the cached source before interpretation.");
        return result;
    }
}
