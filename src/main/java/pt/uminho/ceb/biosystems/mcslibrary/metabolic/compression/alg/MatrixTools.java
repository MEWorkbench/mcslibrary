package pt.uminho.ceb.biosystems.mcslibrary.metabolic.compression.alg;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringTokenizer;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.SingularValueDecomposition;
import org.ejml.ops.CommonOps;
import org.ejml.ops.NormOps;
import org.ejml.ops.SingularOps;

import pt.uminho.ceb.biosystems.mcslibrary.utilities.Pair;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;
import pt.uminho.ceb.biosystems.mew.utilities.java.StringUtils;

public class MatrixTools {
	public static double EPSILON = Math.pow(2, -52);

	public static double[][] computeKernel(double[][] matrix){
		DenseMatrix64F A = new DenseMatrix64F(matrix);
		SingularValueDecomposition<DenseMatrix64F> svd = DecompositionFactory.svd(A.numRows, A.numCols, true, true, false);
		svd.decompose(A);
		DenseMatrix64F nullspace = SingularOps.nullSpace(svd, null, EPSILON*A.numCols);
		System.out.println("Kernel dimension:"+nullspace.numRows+"by"+nullspace.numCols);
		return convertFromD64(nullspace);
	}

	private static double[][] convertFromD64(DenseMatrix64F matrix){
		double[][] reskernel = new double[matrix.numRows][matrix.numCols];
		for (int i = 0; i < reskernel.length; i++) {
			for (int j = 0; j < reskernel[0].length; j++) {
				reskernel[i][j] = matrix.get(i, j);
			}
		}
		return reskernel;
	}

	public static double[][] subMatrixCols(double[][] matrix, List<Integer> toKeepRows, List<Integer> toKeepCols){
		if (toKeepRows == null) {
			toKeepRows = new ArrayList<Integer>();
		}
		if (toKeepCols == null) {
			toKeepCols = new ArrayList<Integer>();
		}
		
		double[][] res = new double[toKeepRows.size()][toKeepCols.size()];
		for (int i = 0; i < res.length; i++) {
			for (int j = 0; j < res[0].length; j++) {
				res[i][j] = matrix[toKeepRows.get(i)][toKeepCols.get(j)];
			}
		}
		return res;
	}
	
	public static ArrayList<Integer> findNonZeroIdx(double[][] matrix, int col){
		ArrayList<Integer> idxs = new ArrayList<Integer>();
		for (int i = 0; i < matrix.length; i++) {
			if (Math.abs(matrix[i][col]) > Utilities.EPSILON) {
				idxs.add(i);
			}
		}
		return idxs;
	}
	
	public static int[] findNonZeroIdxArray(double[][] matrix, int col){
		ArrayList<Integer> idxs = new ArrayList<Integer>();
		for (int i = 0; i < matrix.length; i++) {
			if (Math.abs(matrix[i][col]) > Utilities.EPSILON) {
				idxs.add(i);
			}
		}
		
		int[] res = new int[idxs.size()];
		for (int i = 0; i < res.length; i++)
			res[i] = idxs.get(i);
		return res;
	}
	
	public static boolean hasAnyNonZero(double[][] rd, int col){
		return findNonZeroIdx(rd, col).size() > 0;
	}
	
	public static double[] subtract(double[] A, double[] B) {
		double[] res = new double[A.length];
		Arrays.fill(res, Double.NaN);
		if(A.length == B.length){
			for (int i = 0; i < res.length; i++) {
				res[i] = A[i] - B[i];
			}
		}
		return res;
	}
	
	public static double[] subtract(double[] A, double B) {
		double[] res = new double[A.length];
		for (int i = 0; i < res.length; i++) {
			res[i] = A[i] - B;
		}
		return res;
	}
	
	public static double[] abs(double[] vector){
		double[] res = new double[vector.length];
		for (int i = 0; i < vector.length; i++) {
			res[i] = Math.abs(vector[i]);
		}
		return res;
	}
	
	public static double[][] transpose(double[][] A){
		DenseMatrix64F mat = new DenseMatrix64F(A);
		CommonOps.transpose(mat);
		return convertFromD64(mat);
	}
	public static double[] divide(double[] num, double denom) {
		double[] res = new double[num.length];
		for (int i = 0; i < res.length; i++) {
			if (denom == 0) {
				res[i] = Double.NaN;
			} else { 
				res[i] = num[i] / denom;
			}
		}
		return res;
		
	}
	public static int[] findNonZeroIdx(boolean[] array){
		int dim = sum(array);
		int idx = 0;
		int[] indexes = new int[dim];
		for (int i = 0; i < array.length; i++) {
			if (array[i]) {
				indexes[idx] = i;
				idx++;
			}
		}
		return indexes;
	}
	
	public static int[] findNonZeroIdx(double[] array){
		int dim = countNZ(array);
		int idx = 0;
		int[] indexes = new int[dim];
		for (int i = 0; i < array.length; i++) {
			if (Math.abs(array[i]) > Utilities.EPSILON) {
				indexes[idx] = i;
				idx++;
			}
		}
		return indexes;
	}
	
	public static void writeCSV(double[][] matrix, String filename) throws IOException {
		BufferedWriter b = new BufferedWriter(new FileWriter(filename+".csv"));
		for (int i = 0; i < matrix.length; i++) {
			String line = "\n";
			if (i == 0) {
				line = "";
			}
			for (int j = 0; j < matrix[0].length; j++) {
				line = line + matrix[i][j] + ";";
			}
			b.write(line.substring(0, line.length()-1));
		}
		b.flush();
		b.close();
	}	
	
	public static void writeCSV(double[][] matrix, String filename, String header) throws IOException {
		BufferedWriter b = new BufferedWriter(new FileWriter(filename+".csv"));
		b.write(header);
		for (int i = 0; i < matrix.length; i++) {
			String line = "\n";
			for (int j = 0; j < matrix[0].length; j++) {
				if (j == 0) {
					line = line + matrix[i][j];
				} else {
					line = line + "," + matrix[i][j];
				}
			}
			b.write(line);
		}
		b.flush();
		b.close();
	}
	
	public static void writeCSV(double[][] matrix, String filename, List<String> header, List<String> rows) throws IOException {
		BufferedWriter b = new BufferedWriter(new FileWriter(filename+".csv"));
		b.write(",");
		b.write(StringUtils.concat(",", header));
		for (int i = 0; i < matrix.length; i++) {
			String line = "\n"+rows.get(i)+",";
			for (int j = 0; j < matrix[0].length; j++) {
				if (j == 0) {
					line = line + matrix[i][j];
				} else {
					line = line + "," + matrix[i][j];
				}
			}
			b.write(line);
		}
		b.flush();
		b.close();
	}
	public static double[][] outerProduct(double[][] matrix) {
		DenseMatrix64F mat = new DenseMatrix64F(matrix);
		DenseMatrix64F tmat = mat.copy();
		CommonOps.transpose(tmat);
		DenseMatrix64F resmat = new DenseMatrix64F(mat.numRows,mat.numRows);
		CommonOps.mult(mat, tmat, resmat);
		double[][] res = new double[resmat.numRows][resmat.numCols];
		for (int i = 0; i < res.length; i++) {
			for (int j = 0; j < res[0].length; j++) {
				res[i][j] = resmat.get(i, j);
			}
		}
		return res;
	}

	public static double getVectorMean(double[] vector) {
		double sum = 0;
		for (int i = 0; i < vector.length; i++) {
			sum = sum + vector[i];
		}
		return sum/vector.length;
	}

	public static int getMinimumIdx(double[] vector) {
		double min = Double.POSITIVE_INFINITY;
		int idx = 0;
		for (int i = 0; i < vector.length; i++) {
			if (vector[i] < min) {
				min = vector[i];
				idx = i;
			}
		}
		return idx;
	}
	
	public static int getMaximumIdx(double[] vector) {
		double min = Double.NEGATIVE_INFINITY;
		int idx = 0;
		for (int i = 0; i < vector.length; i++) {
			if (vector[i] > min) {
				min = vector[i];
				idx = i;
			}
		}
		return idx;
	}
	
	public static double getVectorNorm(double[] vector) {
		DenseMatrix64F v = new DenseMatrix64F(1,vector.length);
		for (int i = 0; i < vector.length; i++) {
			v.set(0, i, vector[i]);
		}
		double norm = NormOps.normP(v, 2);
		return norm;
	}
	
	public static double getMatrixNorm(double[][] matrix) {
		DenseMatrix64F mat = new DenseMatrix64F(matrix);
		double norm = NormOps.normP2(mat);
		return norm;
	}

	public static double[][] multiplyMatrix(double[][] a, double[][] b) {
		System.out.println("A: "+a.length+"x"+a[0].length);
		System.out.println("B: "+b.length+"x"+b[0].length);
		DenseMatrix64F da = new DenseMatrix64F(a);
		DenseMatrix64F db = new DenseMatrix64F(b);
		DenseMatrix64F dc = new DenseMatrix64F(a.length, b[0].length);
		CommonOps.mult(da, db, dc);
		double[][] c = new double[dc.numRows][dc.numCols];
		for (int i = 0; i < c.length; i++) {
			for (int j = 0; j < c[0].length; j++) {
				c[i][j] = dc.get(i, j);
			}
		}
		return c;

	}
	public static void printMatrix(Object[][] matrix){
		for (int i = 0; i < matrix.length; i++) {
			String str = "";
			for (int j = 0; j < matrix[0].length; j++) {
				str = str + matrix[i][j] + "\t";
			}
			System.out.println(str);
		}
	}
	public static double[][] readCSV(String filename, boolean header) throws IOException{
		ArrayList<double[]> res = new ArrayList<double[]>();
		BufferedReader br = new BufferedReader(new FileReader(filename));
		int mlength = 0;

		if (header) {
			br.readLine();
		}

		while (br.ready()) {
			String line = br.readLine();
			StringTokenizer tk = new StringTokenizer(line,",");
			int length = tk.countTokens();
			mlength = length > mlength ? length : mlength;
			double[] arr = new double[length];
			int i = 0;
			while (tk.hasMoreTokens()) {
				arr[i] = Double.parseDouble(tk.nextToken());
			}
			res.add(arr);
		}
		br.close();
		double[][] d = new double[res.size()][mlength];
		for (int i = 0; i < res.size(); i++) {
			d[i] = res.get(i);
		}
		return d;
	}
	public static Pair<double[][],int[]> removeZeroRows(double[][] matrix){
		ArrayList<Integer> nzrs = new ArrayList<Integer>();
 		for (int i = 0; i < matrix.length; i++) {
 			int zrCount = 0;
 			for (int j = 0; j < matrix[0].length; j++) {
				if (Math.abs(matrix[i][j]) < EPSILON) {
					zrCount ++;
				} else {
					break;
				}
			}
 			if (zrCount != matrix[0].length) {
				nzrs.add(i);
			}
		}
 		double[][] finalmatrix = new double[nzrs.size()][matrix[0].length];
 		for (int i = 0; i < nzrs.size(); i++) {
 			finalmatrix[i] = matrix[nzrs.get(i)];
		}
 		int[] narr = new int[nzrs.size()];
 		for (int i = 0; i < narr.length; i++) {
			narr[i] = nzrs.get(i);
		}
 		return new Pair<double[][], int[]>(finalmatrix, narr);
	}

	public static int[] getRrefPivots(double[][] matrix) {
		DenseMatrix64F mat = new DenseMatrix64F(matrix);
		DenseMatrix64F rrefmat = CommonOps.rref(mat, -1, null);
		double[][] matx = new double[rrefmat.numRows][rrefmat.numCols];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				matx[i][j] = rrefmat.get(i, j);
			}
		}
		Pair<double[][], int[]> res = removeZeroRows(matx);
		return res.getB();
	}
	
	public static double[][] cutMatrixRows(double[][] matrix, boolean[] toKeep){
		double[][] res = new double[sum(toKeep)][matrix[0].length];
		int idx = 0;
		for (int i = 0; i < matrix.length; i++) {
			if (toKeep[i]) {
				for (int j = 0; j < matrix[0].length; j++) {
					res[idx][j] = matrix[i][j];
				}
				idx++;
			} else {
				continue;
			}
		}
		return res;
	}
	
	public static double[][] cutMatrixColumns(double[][] matrix, boolean[] toKeep){
		double[][] res = new double[matrix.length][sum(toKeep)];
		int idx = 0;
		for (int i = 0; i < matrix[0].length; i++) {
			if (toKeep[i]) {
				for (int j = 0; j < matrix.length; j++) {
					res[j][idx] = matrix[j][i];
				}
				idx++;
			} else {
				continue;
			}
		}
		return res;
	}
	
	
	public static int sum(boolean[] array){
		int res = 0;
		for (int i = 0; i < array.length; i++) {
			if (array[i]) {
				res++;
			}
		}
		
		return res;
	}
	
	public static int countNZ(double[] array){
		int res = 0;
		for (int i = 0; i < array.length; i++) {
			if (Math.abs(array[i]) > Utilities.EPSILON) {
				res++;
			}
		}
		
		return res;
	}

	public static void describe(double[][] sub) {
		String str = "Matrix size is "+sub.length+" by "+sub[0].length;
		System.out.println(str);
	}
	
	


}
