/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

/** Exact biallelic autosomal HWE test for hard additive calls. */
final class HardyWeinberg {
    private HardyWeinberg() { }

    static Result calculate(double[] dosage, int[] groups, int selectedGroup) {
        int homRef = 0;
        int heterozygous = 0;
        int homAlt = 0;
        for (int index = 0; index < dosage.length; index++) {
            if (groups != null && selectedGroup >= 0
                    && groups[index] != selectedGroup) continue;
            double value = dosage[index];
            if (!Double.isFinite(value)) continue;
            long call = Math.round(value);
            if (Math.abs(value - call) > 1e-8 || call < 0 || call > 2)
                return new Result(Double.NaN, 0, "non-hard-call");
            if (call == 0) homRef++;
            else if (call == 1) heterozygous++;
            else homAlt++;
        }
        int samples = homRef + heterozygous + homAlt;
        if (samples == 0) return new Result(Double.NaN, 0, "no-calls");
        return new Result(exact(homRef, heterozygous, homAlt),
            samples, "exact-hard-calls");
    }

    private static double exact(int firstHomozygotes, int heterozygotes,
            int secondHomozygotes) {
        int rareHomozygotes = Math.min(firstHomozygotes, secondHomozygotes);
        int rareCopies = 2 * rareHomozygotes + heterozygotes;
        int genotypes = firstHomozygotes + heterozygotes + secondHomozygotes;
        double[] probabilities = new double[rareCopies + 1];
        int midpoint = rareCopies * (2 * genotypes - rareCopies)
            / (2 * genotypes);
        if ((midpoint & 1) != (rareCopies & 1)) midpoint++;
        int heterozygote = midpoint;
        int rare = (rareCopies - midpoint) / 2;
        int common = genotypes - heterozygote - rare;
        probabilities[midpoint] = 1.0;
        double sum = 1.0;
        while (heterozygote > 1) {
            double value = probabilities[heterozygote]
                * heterozygote * (heterozygote - 1.0)
                / (4.0 * (rare + 1.0) * (common + 1.0));
            probabilities[heterozygote - 2] = value;
            sum += value;
            heterozygote -= 2;
            rare++;
            common++;
        }
        heterozygote = midpoint;
        rare = (rareCopies - midpoint) / 2;
        common = genotypes - heterozygote - rare;
        while (heterozygote <= rareCopies - 2) {
            double value = probabilities[heterozygote]
                * 4.0 * rare * common
                / ((heterozygote + 2.0) * (heterozygote + 1.0));
            probabilities[heterozygote + 2] = value;
            sum += value;
            heterozygote += 2;
            rare--;
            common--;
        }
        for (int index = 0; index < probabilities.length; index++)
            probabilities[index] /= sum;
        double observed = probabilities[heterozygotes];
        double p = 0.0;
        for (double probability : probabilities)
            if (probability <= observed + 1e-12) p += probability;
        return Math.min(1.0, p);
    }

    record Result(double pValue, int samples, String method) { }
}
