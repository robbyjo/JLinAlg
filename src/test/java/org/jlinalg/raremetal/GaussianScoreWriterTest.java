/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.raremetal;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.Path;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gwas.RemlAssociationScanner;
import org.jlinalg.ols.OlsOptions;
import org.jlinalg.pipeline.*;
import org.jlinalg.reml.*;
import org.jlinalg.settest.LinearSetTestNullModel;

class GaussianScoreWriterTest {
    @TempDir Path dir;
    @Test void phenotypeUnitsPreserveScoresAndCovarianceForBothProjections() throws Exception {
        for(boolean related:new boolean[]{false,true}) {
            double[][] base=export(related,1);
            for(double scale:new double[]{1e-12,1e12}) {
                double[][] scaled=export(related,scale);
                for(int i=0;i<2;i++) {
                    assertEquals(base[i][0],scaled[i][0]*scale,2e-8);
                    assertEquals(base[i][1],scaled[i][1]*scale,2e-8);
                    assertEquals(base[i][2],scaled[i][2]/scale,2e-8);
                    assertEquals(base[i][3],scaled[i][3],2e-10);
                }
                for(int i=0;i<4;i++)assertEquals(base[2][i],scaled[2][i]*scale*scale,2e-8);
            }
        }
    }
    private double[][] export(boolean related,double scale) throws Exception {
        int n=40;double[] y=new double[n];double[][] x=new double[n][2],g=new double[4][n];
        for(int i=0;i<n;i++) {
            x[i]=new double[]{1,i%2};y[i]=scale*(i%3+Math.sin(1.3*i));
            g[0][i]=i%3;g[1][i]=(i/3)%3;g[2][i]=1;g[3][i]=i%2;
        }
        // Include missing calls in an informative variant and in the constant one.
        g[1][0]=Double.NaN;g[2][0]=Double.NaN;
        VariantSource source=new VariantSource() {
            public VariantSourceMetadata metadata(){throw new UnsupportedOperationException();}
            public VariantBlockReader open(int[] order){return new VariantBlockReader() {
                int i;
                public VariantBlock read(int maximum){if(i==g.length)return null;int j=i++;
                    return new VariantBlock(j,List.of(new VariantRecord("v"+j,"1",10+j,"A","G",g[j],Double.NaN)));}
                public void close(){}
            };}
        };
        Path scores=dir.resolve(related+"-"+scale+".score.gz"),cov=dir.resolve(related+"-"+scale+".cov.gz");
        var ols=LinearSetTestNullModel.prepare(y,x,OlsOptions.defaults(),BackendPolicy.CPU);
        if(related) {
            // Identity covariance gives an independent analytic comparator for P.
            var options=RemlOptions.builder().initialVariances(ols.residualVariance())
                .varianceBounds(1e-10*scale*scale,1e10*scale*scale).build();
            var reml=RemlAssociationScanner.prepare(y,x,List.of(VarianceComponent.identity("residual",n)),options,BackendPolicy.CPU);
            GaussianScoreWriter.write(source,null,reml,100,10,"test",scores,cov);
        } else GaussianScoreWriter.write(source,null,ols,100,10,"test",scores,cov);
        List<String[]> rows=new ArrayList<>();
        try(var reader=RareMetalStudy.open(scores)) {
            for(String line;(line=reader.readLine())!=null;)if(!line.startsWith("#"))rows.add(line.split("\t"));
        }
        double[][] result=new double[3][4];
        for(int i=0;i<2;i++)for(int j=0;j<4;j++)result[i][j]=Double.parseDouble(rows.get(i)[13+j]);
        assertTrue(result[0][1]>0);assertTrue(result[0][3]>0&&result[0][3]<1);
        for(int i=2;i<4;i++) {
            assertEquals(0,Double.parseDouble(rows.get(i)[13]));assertEquals(0,Double.parseDouble(rows.get(i)[14]));
            assertEquals("NA",rows.get(i)[16]);
        }
        // Read raw disk covariance as the reader's absolute reconciliation tolerance
        // is a separate contract; no tolerance should conceal an exported zero.
        try(var reader=RareMetalStudy.open(cov)) {
            int row=0;
            for(String line;(line=reader.readLine())!=null;)if(!line.startsWith("#")) {
                String[] values=line.split("\t")[3].split(",");
                if(row==0){result[2][0]=Double.parseDouble(values[0]);result[2][1]=Double.parseDouble(values[1]);}
                if(row==1){result[2][2]=Double.parseDouble(values[0]);result[2][3]=Double.parseDouble(values[1]);}
                row++;
            }
        }
        return result;
    }
}
