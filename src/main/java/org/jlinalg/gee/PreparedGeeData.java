/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

import java.util.Arrays;

/**
 * Immutable, cluster-sorted GEE inputs reusable across working structures,
 * covariance estimators, and warm-started fits.
 */
public final class PreparedGeeData {
    private final double[] response;
    private final double[] design;
    private final double[] weights;
    private final double[] offset;
    private final double[] scaleDesign;
    private final int scaleColumns;
    private final int rows;
    private final int columns;
    private final int[] cluster;
    private final int[] waves;
    private final int[] starts;
    private final int maximumWave;
    private final int[] retainedRows;
    private final int[] outputPosition;
    private final int originalRows;

    PreparedGeeData(
            double[] response,
            double[] design,
            double[] weights,
            double[] offset,
            double[] scaleDesign,
            int scaleColumns,
            int rows,
            int columns,
            int[] cluster,
            int[] waves,
            int[] starts,
            int maximumWave,
            int[] retainedRows,
            int[] outputPosition,
            int originalRows) {
        this.response = response;
        this.design = design;
        this.weights = weights;
        this.offset = offset;
        this.scaleDesign = scaleDesign;
        this.scaleColumns = scaleColumns;
        this.rows = rows;
        this.columns = columns;
        this.cluster = cluster;
        this.waves = waves;
        this.starts = starts;
        this.maximumWave = maximumWave;
        this.retainedRows = retainedRows;
        this.outputPosition = outputPosition;
        this.originalRows = originalRows;
    }

    public int observations() { return rows; }
    public int originalObservations() { return originalRows; }
    public int omittedObservations() { return originalRows - rows; }
    public int columns() { return columns; }
    public int clusters() { return starts.length - 1; }
    public int maximumWaveCount() { return maximumWave; }
    public int[] retainedRows() { return retainedRows.clone(); }

    public int minimumClusterSize() {
        int result = Integer.MAX_VALUE;
        for (int index = 0; index < clusters(); index++) {
            result = Math.min(result, starts[index + 1] - starts[index]);
        }
        return result;
    }

    public int maximumClusterSize() {
        int result = 0;
        for (int index = 0; index < clusters(); index++) {
            result = Math.max(result, starts[index + 1] - starts[index]);
        }
        return result;
    }

    double[] response() { return response; }
    double[] design() { return design; }
    double[] weights() { return weights; }
    double[] offset() { return offset; }
    double[] scaleDesign() { return scaleDesign; }
    int scaleColumns() { return scaleColumns; }
    int rows() { return rows; }
    int[] cluster() { return cluster; }
    int[] waves() { return waves; }
    int[] starts() { return starts; }
    int maximumWave() { return maximumWave; }
    int[] outputPosition() { return outputPosition; }
    int originalRows() { return originalRows; }

    double[] output(double[] sorted) {
        double[] result = new double[rows];
        for (int index = 0; index < rows; index++) {
            result[outputPosition[index]] = sorted[index];
        }
        return result;
    }

    PreparedGeeData withoutCluster(int omittedCluster) {
        if (omittedCluster < 0 || omittedCluster >= clusters()) {
            throw new IllegalArgumentException("cluster index is out of range");
        }
        int omittedStart = starts[omittedCluster];
        int omittedEnd = starts[omittedCluster + 1];
        int omitted = omittedEnd - omittedStart;
        int kept = rows - omitted;
        double[] keptResponse = removeRows(response, 1, omittedStart, omittedEnd);
        double[] keptDesign = removeRows(design, columns, omittedStart, omittedEnd);
        double[] keptWeights = removeRows(weights, 1, omittedStart, omittedEnd);
        double[] keptOffset = removeRows(offset, 1, omittedStart, omittedEnd);
        double[] keptScale = scaleDesign == null ? null
            : removeRows(scaleDesign, scaleColumns, omittedStart, omittedEnd);
        int[] keptCluster = removeRows(cluster, omittedStart, omittedEnd);
        int[] keptWaves = removeRows(waves, omittedStart, omittedEnd);
        int[] keptStarts = new int[starts.length - 1];
        int destination = 0;
        for (int index = 0; index < starts.length; index++) {
            if (index == omittedCluster + 1) continue;
            int value = starts[index];
            keptStarts[destination++] = value <= omittedStart ? value : value - omitted;
        }
        int[] identity = new int[kept];
        for (int row = 0; row < kept; row++) identity[row] = row;
        return new PreparedGeeData(keptResponse, keptDesign, keptWeights,
            keptOffset, keptScale, scaleColumns, kept, columns, keptCluster,
            keptWaves, keptStarts, maximumWave, identity, identity.clone(), kept);
    }

    private static double[] removeRows(
            double[] source, int width, int start, int end) {
        double[] result = new double[source.length - (end - start) * width];
        int first = start * width;
        int second = end * width;
        System.arraycopy(source, 0, result, 0, first);
        System.arraycopy(source, second, result, first, source.length - second);
        return result;
    }

    private static int[] removeRows(int[] source, int start, int end) {
        int[] result = Arrays.copyOf(source, source.length - (end - start));
        System.arraycopy(source, end, result, start, source.length - end);
        return result;
    }
}
