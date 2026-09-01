/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

/** Validated right-censored or counting-process survival response. */
public final class CoxSurvivalData {
    private final double[] start;
    private final double[] stop;
    private final boolean[] event;
    private final int[] strata;

    public CoxSurvivalData(
            double[] start, double[] stop, boolean[] event, int[] strata) {
        if (start == null || stop == null || event == null
                || start.length == 0 || stop.length != start.length
                || event.length != start.length
                || (strata != null && strata.length != start.length))
            throw new IllegalArgumentException(
                "survival arrays must be nonempty and have equal lengths");
        this.start = start.clone();
        this.stop = stop.clone();
        this.event = event.clone();
        this.strata = strata == null ? new int[start.length] : strata.clone();
        int events = 0;
        for (int row = 0; row < start.length; row++) {
            if (!Double.isFinite(this.start[row])
                    || !Double.isFinite(this.stop[row])
                    || this.start[row] < 0
                    || !(this.stop[row] > this.start[row])
                    || this.strata[row] < 0)
                throw new IllegalArgumentException(
                    "survival rows require 0 <= start < stop and nonnegative strata");
            if (this.event[row]) events++;
        }
        if (events == 0)
            throw new IllegalArgumentException(
                "at least one observed event is required");
    }

    public static CoxSurvivalData rightCensored(
            double[] time, boolean[] event) {
        return rightCensored(time, event, null);
    }

    public static CoxSurvivalData rightCensored(
            double[] time, boolean[] event, int[] strata) {
        if (time == null)
            throw new IllegalArgumentException("event/censoring times are required");
        return new CoxSurvivalData(new double[time.length], time, event, strata);
    }

    public int observations() { return stop.length; }
    public double[] start() { return start.clone(); }
    public double[] stop() { return stop.clone(); }
    public boolean[] event() { return event.clone(); }
    public int[] strata() { return strata.clone(); }
    double[] startView() { return start; }
    double[] stopView() { return stop; }
    boolean[] eventView() { return event; }
    int[] strataView() { return strata; }
}
