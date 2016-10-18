/*******************************************************************************
 * Copyright 2016
 * CEB Centre of Biological Engineering
 * University of Minho
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This code is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this code. If not, see http://www.gnu.org/licenses/
 *
 * Created inside the BIOSYSTEMS Research Group
 * (http://www.ceb.uminho.pt/biosystems)
 *******************************************************************************/
package pt.uminho.ceb.biosystems.mcslibrary.utilities;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.StringTokenizer;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.solution.SolutionUtilities;
import pt.uminho.ceb.biosystems.mew.utilities.java.StringUtils;
/**
 * Class with static methods that are used in several other classes of this package.
 * @author V�tor
 *
 */
public final class Utilities {
	public static final double EPSILON = Math.pow(2, -52);
	public static final double PRECISION = 1e-10;
	public static final double INF = Double.MAX_VALUE;

	public static ArrayList<int[]> intArrayDecompressor(ArrayList<int[]> options) {
		ArrayList<int[]> container = new ArrayList<int[]>();
		_intArrayDecompressor(new int[]{}, options, container);
		return container;
	}


	public static List<String> readLines(String filename) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(filename));
		ArrayList<String> res = new ArrayList<String>();
		while (br.ready()) {
			res.add(br.readLine());
		}
		br.close();
		return res;
	}

	public static void writeLines(String filename, String delimiter, List<List<String>> data) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
		for (List<String> list : data) {
			bw.write(StringUtils.concat(delimiter, list));
		}
		bw.flush();
		bw.close();
	}
	
//	public static void writeLines(String filename, String delimiter, Set<Set<String>> data) throws IOException{
//		BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
//		for (Set<String> list : data) {
//			bw.write(StringUtils.concat(delimiter, list)+"\n");
//		}
//		bw.flush();
//		bw.close();
//	}
	
	public static void writeLines(String filename, String delimiter, Set<String> data) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
		for (String list : data) {
			bw.write(list+"\n");
		}
		bw.flush();
		bw.close();
	}



	private static void _intArrayDecompressor(int[] array, ArrayList<int[]> remaining, ArrayList<int[]> container){
		if (remaining.size() == 0) {
			container.add(array);
		} else {

			int[] next = remaining.get(remaining.size()-1);
			@SuppressWarnings("unchecked")
			ArrayList<int[]> newremaining = (ArrayList<int[]>) remaining.clone();
			newremaining.remove(remaining.size()-1);
			for (int i = 0; i < next.length; i++) {
				int[] newarray = new int[array.length+1];
				for (int j = 0; j < array.length; j++) {
					newarray[j] = array[j];
				}
				newarray[array.length] = next[i];
				_intArrayDecompressor(newarray, newremaining, container);
			}
		}
	}

	public static void convertSolutionToCSV(String fileIn, String fileOut) throws IOException {
//		String fileIn = "files/vitor/R_EX_succ_e__1428957900659mcs.txt";
//		String fileOut = "files/vitor/out.csv";

		BufferedReader br = new BufferedReader(new FileReader(fileIn));
		BufferedWriter bw = new BufferedWriter(new FileWriter(fileOut));

		bw.write("ID,SOLUTION");
		int i=0;
		while(br.ready()){
			String line = br.readLine();
			String [] tokens = line.split(";");

			bw.write("\nsol_"+i+","+StringUtils.concat(" ", tokens));
			i++;
		}

		bw.flush();
		bw.close();
		br.close();
	}
	
	public static void convertSolutionToCSV(List<List<String>> solutions, String fileOut) throws IOException {
//		String fileIn = "files/vitor/R_EX_succ_e__1428957900659mcs.txt";
//		String fileOut = "files/vitor/out.csv";

		BufferedWriter bw = new BufferedWriter(new FileWriter(fileOut));

		bw.write("ID,SOLUTION");
		for (int j = 0; j < solutions.size(); j++) {
			List<String> tokens = solutions.get(j);
			bw.write("\nsol_"+j+","+StringUtils.concat(" ", tokens));
			j++;
		}

		bw.flush();
		bw.close();
	}


	public static void printArray(int[] array){
		for (int i = 0; i < array.length; i++) {
			System.out.print(array[i]+", ");
		}
		System.out.println();
	}

	public static Reaction[] toReacArrayFromInt(DefaultMetabolicNetwork metaNet, int[] strs){
		int size = strs.length;
		int i = 0;
		Reaction[] res = new Reaction[size];
		for (int str : strs){
			Reaction r = metaNet.getReaction(str);
			res[i] = r;
			i++;
		}
		return res;
	}
	public static Reaction[] toReacArrayFromString(DefaultMetabolicNetwork metaNet, Set<String> strs){
		int size = strs.size();
		int i = 0;
		Reaction[] res = new Reaction[size];
		for (String str : strs){
			Reaction r = metaNet.getReaction(metaNet.getReactionIndex(str));
			res[i] = r;
			i++;
		}
		return res;
	}
	
	public static Reaction[] toReacArrayFromString(DefaultMetabolicNetwork metaNet, List<String> strs){
		int size = strs.size();
		int i = 0;
		Reaction[] res = new Reaction[size];

		for (int j = 0; j < size; j++) {
			String str = strs.get(j);
			Reaction r = metaNet.getReaction(metaNet.getReactionIndex(str));
			res[i] = r;
			i++;
		}
		return res;
	}


	public static void printNonZeroArray(double[] ds) {
		String str = "";
		for (int i = 0; i < ds.length; i++) {
			if (Math.abs(ds[i]) > PRECISION) {
				str = str + i + " = " + ds[i] + "; ";
			}
		}
		System.out.println(str);
	}

	public static List<String> getAllStringTokens(String s, String delimiter) {
		StringTokenizer t = new StringTokenizer(s, delimiter);
		List<String> res = new ArrayList<String>();
		while (t.hasMoreTokens()) {
			res.add(t.nextToken());
		}
		return res;

	}

	public static <K, V> void printMap(Map<K,V> map) {
		for (Entry<K, V> entry : map.entrySet()) {
			System.out.println(entry.getKey().toString()+" -> "+entry.getValue().toString());
		}
	}

	public static <K> Map<K,Integer> filterFromMap(Map<K,Integer> map, int threshold) {
		HashMap<K, Integer> filtered = new HashMap<K,Integer>();
		for (Entry<K, Integer> entry : map.entrySet()) {
			if (entry.getValue() > threshold) {
				filtered.put(entry.getKey(), entry.getValue());
			}
		}
		return SolutionUtilities.sortByValues(filtered);
	}
}
