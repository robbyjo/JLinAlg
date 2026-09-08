/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

/** Scale-invariant rank checks on the analytic distribution Jacobian. */
final class SemInformation {
    private SemInformation() { }
    /** Pivoted, twice-reorthogonalized MGS of normalized parameter columns.
     * Near dependencies below 1e-8 are conservatively unavailable inference. */
    static boolean fullColumnRank(double[][] columns) {
        int k=columns.length;if(k==0)return true;
        int rows=columns[0].length;if(k>rows)return false;
        double[][] work=new double[k][];
        for(int i=0;i<k;i++) {
            work[i]=columns[i].clone();double norm=norm(work[i]);
            if(!(norm>0) || !Double.isFinite(norm))return false;
            for(int r=0;r<rows;r++)work[i][r]/=norm;
        }
        for(int j=0;j<k;j++) {
            int pivot=j;double maximum=norm(work[j]);
            for(int i=j+1;i<k;i++){double value=norm(work[i]);if(value>maximum){maximum=value;pivot=i;}}
            if(!(maximum>1e-8))return false;
            double[] tmp=work[j];work[j]=work[pivot];work[pivot]=tmp;
            for(int r=0;r<rows;r++)work[j][r]/=maximum;
            for(int i=j+1;i<k;i++)for(int pass=0;pass<2;pass++) {
                double projection=SemOptimizer.dot(work[i],work[j]);
                for(int r=0;r<rows;r++)work[i][r]-=projection*work[j][r];
            }
        }
        return true;
    }
    static double norm(double[] vector){double norm=0;for(double x:vector)norm=Math.hypot(norm,x);return norm;}
    static double[] unavailable(int k){double[] a=new double[k*k];java.util.Arrays.fill(a,Double.NaN);return a;}
    static boolean finiteCovariance(double[] covariance,int k) {
        for(double x:covariance)if(!Double.isFinite(x))return false;
        for(int i=0;i<k;i++)if(covariance[i*k+i]<0)return false;
        return true;
    }
}
