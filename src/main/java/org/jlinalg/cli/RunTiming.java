/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.time.Instant;
import java.util.Locale;

/** UTC wall timestamps with a monotonic elapsed-time clock. */
final class RunTiming {
    private final Instant started = Instant.now();
    private final long startedNanos = System.nanoTime();

    String started() { return started.toString(); }
    End finish() {
        long milliseconds = Math.max(0, (System.nanoTime() - startedNanos) / 1_000_000);
        return new End(Instant.now().toString(), milliseconds, readable(milliseconds));
    }
    record End(String finished, long elapsedMillis, String elapsed) { }

    static String readable(long milliseconds) {
        if (milliseconds < 0) throw new IllegalArgumentException("elapsed time must be nonnegative");
        long days = milliseconds / 86_400_000;
        long hours = milliseconds / 3_600_000 % 24;
        long minutes = milliseconds / 60_000 % 60;
        long seconds = milliseconds / 1000 % 60;
        StringBuilder result = new StringBuilder();
        if (days > 0) result.append(days).append(" d ");
        if (days > 0 || hours > 0) result.append(hours).append(" h ");
        if (days > 0 || hours > 0 || minutes > 0) result.append(minutes).append(" min ");
        return result.append(String.format(Locale.ROOT, "%d.%03d s", seconds, milliseconds % 1000)).toString();
    }
}
