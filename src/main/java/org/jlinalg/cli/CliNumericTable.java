/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Numeric view over the shared CLI delimited-table reader. */
final class CliNumericTable {
    private final DelimitedData data;
    private CliNumericTable(DelimitedData data) { this.data = data; }
    static CliNumericTable read(Path path) throws IOException { return new CliNumericTable(DelimitedData.read(path)); }
    int rows() { return data.rows().size(); }
    List<String> header() { return data.header(); }
    double[] column(String name) { int index=data.column(name);double[] result=new double[rows()];for(int i=0;i<result.length;i++)result[i]=number(data.rows().get(i)[index],i,name);return result; }
    double[][] matrix(List<String> names) { int[] index=names.stream().mapToInt(data::column).toArray();double[][] result=new double[rows()][index.length];for(int i=0;i<result.length;i++)for(int j=0;j<index.length;j++)result[i][j]=number(data.rows().get(i)[index[j]],i,names.get(j));return result; }
    double[][] design(List<String> names,boolean intercept){double[][] source=matrix(names),result=new double[rows()][names.size()+(intercept?1:0)];for(int i=0;i<result.length;i++){int offset=0;if(intercept)result[i][offset++]=1.0;System.arraycopy(source[i],0,result[i],offset,source[i].length);}return result;}
    int[] integerGroups(String name){int index=data.column(name);Map<String,Integer> values=new LinkedHashMap<>();int[] result=new int[rows()];for(int i=0;i<result.length;i++){String value=data.rows().get(i)[index].trim();if(value.isEmpty())throw new IllegalArgumentException("blank cluster at row "+(i+2));result[i]=values.computeIfAbsent(value,ignored->values.size());}return result;}
    double[][] withColumnSet(double[][] design,int designColumn,double value){double[][] result=copy(design);for(double[] row:result)row[designColumn]=value;return result;}
    int designColumn(List<String> names,String name,boolean intercept){int index=names.indexOf(name);if(index<0)throw new IllegalArgumentException("column is not in --predictors: "+name);return index+(intercept?1:0);}
    private static double number(String value,int row,String name){try{double result=Double.parseDouble(value.trim());if(!Double.isFinite(result))throw new NumberFormatException("not finite");return result;}catch(NumberFormatException failure){throw new IllegalArgumentException("non-numeric value at row "+(row+2)+", column "+name,failure);}}
    static List<String> columns(String value){List<String> result=Arrays.stream(value.split(",",-1)).map(String::trim).filter(v->!v.isBlank()).toList();if(result.isEmpty())throw new IllegalArgumentException("column list must not be empty");return result;}
    static double[][] copy(double[][] matrix){double[][] result=new double[matrix.length][];for(int i=0;i<matrix.length;i++)result[i]=matrix[i].clone();return result;}
}
