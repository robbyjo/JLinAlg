/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;
final class HardyWeinberg {
    private HardyWeinberg() { }
    static Result calculate(double[] dosage,int[] groups,int selectedGroup) {
        var r=org.jlinalg.genetics.HardyWeinberg.calculate(dosage,groups,selectedGroup);
        return new Result(r.pValue(),r.samples(),r.method());
    }
    record Result(double pValue,int samples,String method) { }
}
