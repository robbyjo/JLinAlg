/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import java.util.Arrays;
import java.util.TreeSet;

/** Dynamic scaled risk moments with rebuilds after cancellation-prone removals. */
final class CoxRiskMoments {
    private final double[] design, eta;
    private final int columns;
    private final TreeSet<Integer> active = new TreeSet<>();
    final double[] first, second;
    double maximum = Double.NEGATIVE_INFINITY, sum;
    private double removed;

    CoxRiskMoments(double[] design, int columns, double[] eta) {
        this.design=design; this.columns=columns; this.eta=eta;
        first=new double[columns]; second=new double[columns*columns];
    }

    void add(int row) {
        active.add(row);
        if (eta[row]>maximum) {
            double factor=Double.isFinite(maximum)?Math.exp(maximum-eta[row]):0;
            sum*=factor; removed*=factor;
            for(int j=0;j<first.length;j++)first[j]*=factor;
            for(int j=0;j<second.length;j++)second[j]*=factor;
            maximum=eta[row];
        }
        accumulate(row,Math.exp(eta[row]-maximum));
    }

    void remove(int row) {
        active.remove(row);
        double weight=Math.exp(eta[row]-maximum);
        accumulate(row,-weight); removed+=weight;
    }

    void refresh() {
        // Bound subtractive error relative to the remaining mass, and recover
        // risks that underflowed under a maximum which is no longer active.
        if (!(sum>0) || removed>.25*sum) {
            maximum=Double.NEGATIVE_INFINITY;
            for(int row:active)maximum=Math.max(maximum,eta[row]);
            sum=0; removed=0; Arrays.fill(first,0); Arrays.fill(second,0);
            for(int row:active)accumulate(row,Math.exp(eta[row]-maximum));
        }
        if (!(sum>0) || !Double.isFinite(sum))
            throw new IllegalArgumentException("Cox risk-set denominator is nonpositive or nonfinite");
    }

    private void accumulate(int row,double weight) {
        sum+=weight;
        for(int j=0;j<columns;j++) {
            double value=design[row*columns+j]; first[j]+=weight*value;
            for(int k=0;k<=j;k++)second[j*columns+k]+=weight*value*design[row*columns+k];
        }
    }
}
