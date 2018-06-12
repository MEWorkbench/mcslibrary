package pt.uminho.ceb.biosystems.mcslibrary.utilities;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.StringTokenizer;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.AbstractMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
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
	public static final double RELEVANCE_THRES = 1e-7;
	public static final double INF = Double.MAX_VALUE;

	
	
	public static ArrayList<int[]> intArrayDecompressor(ArrayList<int[]> options) {
		ArrayList<int[]> container = new ArrayList<int[]>();
		_intArrayDecompressor(new int[]{}, options, container);
		return container;
	}
	
	public static List<double[]> doubleArrayDecompressor(ArrayList<double[]> options) {
		ArrayList<double[]> container = new ArrayList<double[]>();
		_doubleArrayDecompressor(new double[]{}, options, container);
		return (List<double[]>) container;
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
			bw.write(StringUtils.concat(delimiter, list)+"\n");
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
	
	public static void writeLines(String filename, String delimiter, Collection<String> data) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
		for (String list : data) {
			bw.write(list+"\n");
		}
		bw.flush();
		bw.close();
	}
	
	public static Object[] concat(Object[] e1, Object[]e2){
		Object[] res = new Object[e1.length + e2.length];
		for (int i = 0; i < e1.length; i++) {
			res[i] = e1[i];
		}
		for (int i = 0; i < e2.length; i++) {
			res[i+e1.length] = e2[i];
		}
		return res;
	}
	
	
	public static <T extends Collection<S>, S extends Collection<String>> void writeDataset(T dataset, String filename, String delimiter) throws IOException{
		Collection<String> toWrite = new ArrayList<String>();
		for (Collection<String> collection : dataset) {
			String str = StringUtils.concat(delimiter, collection).replace("[", "").replace("]", "");
			toWrite.add(str);
		}
		writeLines(filename, delimiter, toWrite);
	}
	
	public static void write(String string, String filename) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
		bw.write(string);
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
	
	private static void _doubleArrayDecompressor(double[] array, ArrayList<double[]> remaining, ArrayList<double[]> container){
		if (remaining.size() == 0) {
			container.add(array);
		} else {

			double[] next = remaining.get(remaining.size()-1);
			@SuppressWarnings("unchecked")
			ArrayList<double[]> newremaining = (ArrayList<double[]>) remaining.clone();
			newremaining.remove(remaining.size()-1);
			for (int i = 0; i < next.length; i++) {
				double[] newarray = new double[array.length+1];
				for (int j = 0; j < array.length; j++) {
					newarray[j] = array[j];
				}
				newarray[array.length] = next[i];
				_doubleArrayDecompressor(newarray, newremaining, container);
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
	public static Reaction[] toReacArrayFromString(AbstractMetabolicNetwork metaNet, Set<String> strs){
		int size = strs.size();
		int i = 0;
		Reaction[] res = new Reaction[size];
		for (String str : strs){
			Reaction r = metaNet.getReaction(str);
			res[i] = r;
			i++;
		}
		return res;
	}
	
	public static Reaction[] toReacArrayFromString(AbstractMetabolicNetwork metaNet, List<String> strs){
		int size = strs.size();
		int i = 0;
		Reaction[] res = new Reaction[size];

		for (int j = 0; j < size; j++) {
			String str = strs.get(j);
			Reaction r = metaNet.getReaction(str);
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
			res.add(t.nextToken().trim());
		}
		return res;

	}
	
	public static List<String> getAllUntrimmedStringTokens(String s, String delimiter) {
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
	
	public static <T> T[] removeItemFromArray(T[] array, int index){
		T[] res = (T[]) new Object[array.length-1];
		for (int i = 0; i < index; i++) {
			res[i] = array[i];
		}
		for (int i = index; i < array.length; i++) {
			res[i] = array[i+1];
		}
		
		return res;
	}
	
	public static double[] removeItemFromArray(double[] array, int index){
		double[] res = new double[array.length-1];
		for (int i = 0; i < index; i++) {
			res[i] = array[i];
		}
		for (int i = index; i < array.length; i++) {
			res[i] = array[i+1];
		}
		
		return res;
	}
	
	public static <T extends Object> T[] appendItemToArray(T[] array, T item){
		@SuppressWarnings("unchecked")
		T[] res = (T[]) new Object[array.length+1];
		for (int i = 0; i < array.length; i++) {
			res[i] = array[i];
		}
		res[array.length] = item;
		return res;
	}
	
	public static double[] appendItemToArray(double[] array, double item){
		double[] res = new double[array.length+1];
		for (int i = 0; i < array.length; i++) {
			res[i] = array[i];
		}
		res[array.length] = item;
		return res;
	}
	
	public static Map<String,String> parseMap(String filename) throws IOException{
		Map<String,String> parameters = new HashMap<>();
		List<List<String>> tokens = SolutionUtilities.tokenizeFile(filename, "=");
		for (int i = 0; i < tokens.size(); i++) {
			parameters.put(tokens.get(i).get(0), tokens.get(i).get(1));
		}
		return parameters;
		
		
	}

	public static Pair<String,FluxBound[]> readBounds(String boundsPath, DefaultMetabolicNetwork dmn) throws IOException {
		List<String> lines = Utilities.readLines(boundsPath);
		FluxBound[] bounds = new FluxBound[lines.size()-1];
		for (int i = 1; i < lines.size(); i++) {
			String line = lines.get(i);
			List<String> tokens = Utilities.getAllStringTokens(line, ",");
			bounds[i-1] = new FluxBound(dmn.getReaction(tokens.get(0)), Double.parseDouble(tokens.get(1)), Double.parseDouble(tokens.get(2)));
		}
		return new Pair<String, FluxBound[]>(lines.get(0), bounds);
	}
	
	public static Map<String, String> readPropertyMap(String mapPath) throws IOException{
		HashMap<String, String> map = new HashMap<>();
		for (String line : readLines(mapPath)) {
			List<String> tokens = Utilities.getAllStringTokens(line, "=");
//			System.out.println(tokens);
			map.put(tokens.get(0), tokens.get(1));
		}
		return map;
		
	}
	
	public static Map<String, String> readMap(String mapPath, String delimiter) throws IOException{
		HashMap<String, String> map = new HashMap<>();
		for (String line : readLines(mapPath)) {
			List<String> tokens = Utilities.getAllStringTokens(line, delimiter);
//			System.out.println(tokens);
			if (tokens.size() == 1) {
				map.put(tokens.get(0), tokens.get(0));
			} else {
				map.put(tokens.get(0), tokens.get(1));
			}
		}
		return map;
		
	}
	
	public static List<Set<String>> toHashSets(List<List<String>> solutions){
		List<Set<String>> ls = new ArrayList<Set<String>>();
		for (int i = 0; i < solutions.size(); i++) {
			ls.add(new HashSet<>(solutions.get(i)));
		}
		return ls;
		
	}
	
	public static <K,V extends Object> void writeMap(Map<K,V> map, String path, String delimiter) throws IOException{
		List<String> lines = new ArrayList<>();
		for (K string : map.keySet()) {
			lines.add(string+delimiter+map.get(string));
		}
		Utilities.writeLines(path, delimiter, lines);
		
	}
}
