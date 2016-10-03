package pt.uminho.ceb.biosystems.mcslibrary.solution;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.TreeMap;
import org.paukov.combinatorics.CombinatoricsVector;
import org.paukov.combinatorics.Factory;
import org.paukov.combinatorics.Generator;
import org.paukov.combinatorics.ICombinatoricsVector;

import pt.uminho.ceb.biosystems.mcslibrary.enumeration.implementation.DefaultEnumerationResult;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Pair;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;

public class SolutionUtilities {
	public static List<List<String>> tokenizeFile(String filename, String delimiter) throws IOException{
		List<List<String>> res = new ArrayList<List<String>>();
		List<String> strs = Utilities.readLines(filename);
		for (String string : strs) {
			res.add(tokenize(string, delimiter));
		}
		return res;
	}

	public static <E extends Object> boolean isEqual(List<E> from, List<E> to){
		return (getCommonElements(from, to).size() == from.size()) && (from.size() == to.size());
	}
	public static <E extends Object> double[][] getDistanceMatrix(List<List<E>> from, List<List<E>> to) {
        double[][] distances = new double[from.size()][to.size()];
        for (int i = 0; i < from.size(); i++) {
			for (int j = 0; j < to.size(); j++) {
				distances[i][j] = overlapDistance(from.get(i), to.get(j));
			}
		}
		return distances;

	}


	public static <E extends Object> boolean hasElement(List<E> array, E element) {
		boolean res = false;
		for (int i = 0; i < array.size(); i++) {
			if (array.get(i).equals(element)) {
				res = true;
				break;
			}
		}
		return res;
	}

	public static <E extends Object> boolean containsSubset(List<E> superset, List<E> subset){
		return getCommonElements(superset, subset).size() == subset.size();
	}

	public static <E extends Object> boolean containsElement(List<List<E>> list, List<E> element){
		boolean found = false;
		for (List<E> e : list) {
			if (isEqual(element, e)) {
				found = true;
				break;
			}
		}
		return found;
	}
	public static <E extends Object> int overlapDistance(List<E> from, List<E> to) {
		int dist = 0;
		for (int i = 0; i < from.size(); i++) {
			if (!hasElement(to, from.get(i))) {
				dist++;
			}
		}

		return dist;
	}

	public static <E extends Object> List<E> getCommonElements(List<E> from, List<E> to) {
		List<E> res = new ArrayList<E>();
		for (int i = 0; i < from.size(); i++) {
			if (hasElement(to, from.get(i))) {
				res.add(from.get(i));
			}
		}
		return res;
	}

	public static <E extends Object> List<E> getExceedingElements(List<E> from, List<E> to) {
		List<E> res = new ArrayList<E>();
		for (int i = 0; i < from.size(); i++) {
			if (!hasElement(to, from.get(i))) {
				res.add(from.get(i));
			}
		}
		return res;
	}


	public static List<String> tokenize(String string, String delimiter) {
		List<String> res = new ArrayList<String>();
		StringTokenizer tok = new StringTokenizer(string,delimiter);
		while (tok.hasMoreTokens()) {
			res.add(tok.nextToken().trim());
		}
		return res;

	}

	public static List<List<String>> fromResultToList(DefaultEnumerationResult res) {
        List<List<String>> kns = new ArrayList<List<String>>();
        for (String[] list : res.toStringArrays()) {
			ArrayList<String> kn = new ArrayList<String>();
			for (String string : list) {
				kn.add(string);
			}
			kns.add(kn);
		}
		return kns;
	}

	public static <K, V extends Comparable<V>> Map<K, V> sortByValues(final Map<K, V> map) {
	    Comparator<K> valueComparator =  new Comparator<K>() {
	        public int compare(K k1, K k2) {
	            int compare = map.get(k2).compareTo(map.get(k1));
	            if (compare == 0) return 1;
	            else return compare;
	        }
	    };
	    Map<K, V> sortedByValues = new TreeMap<K, V>(valueComparator);
	    sortedByValues.putAll(map);
	    return sortedByValues;
	}

	public static Map<String, Integer> getReactionFrequencies(List<List<String>> sls) {
		Map<String, Integer> map = new TreeMap<String, Integer>();
		for (List<String> list : sls) {
			for (String string : list) {
				if (map.containsKey(string)) {
					map.put(string, map.get(string) + 1);
				} else {
					map.put(string, 1);
				}
			}
		}
		Map<String, Integer> finalmap = sortByValues(map);
		return finalmap;
	}

	public static Map<List<String>,Integer> getKmerFrequencies(List<List<String>> sls, int size) {
		Map<List<String>, Integer> map = new HashMap<List<String>,Integer>();
		for (List<String> solution : sls) {
			CombinatoricsVector<String> vector = new CombinatoricsVector<String>(solution);
			Generator<String> gen = Factory.createSimpleCombinationGenerator(vector, size);
			for (ICombinatoricsVector<String> combination : gen) {
				List<String> list = combination.getVector();
				if (map.containsKey(list)) {
					map.put(list, map.get(list)+1);
				} else {
					map.put(list, 1);
				}
			}
		}
		return sortByValues(map);
	}

	public static List<Pair<Boolean,List<String>>> getReactionAdditionCombinations(List<String> fromSol, List<String> toSol) {
		List<String> core = getCommonElements(fromSol, toSol);
		List<String> exc = getExceedingElements(toSol, fromSol);
		exc.addAll(getExceedingElements(fromSol, toSol));

		CombinatoricsVector<String> possibilityVector = new CombinatoricsVector<String>(exc);
		List<Pair<Boolean,List<String>>> res = new ArrayList<Pair<Boolean,List<String>>>();

		for (int i = 0; i < possibilityVector.getSize(); i++) {
			Generator<String> g = Factory.createSimpleCombinationGenerator(possibilityVector, i+1);
			for (ICombinatoricsVector<String> v : g) {
				List<String> sol = new ArrayList<String>();
				sol.addAll(core);
				sol.addAll(v.getVector());
				res.add(new Pair<Boolean, List<String>>(v.getSize() == toSol.size(),sol));
			}
		}
		return res;
	}
	
	public static int[][] getReactionKnockoutMatrix(List<List<String>> sols){
		int maxsize = 0;
		int nSols = sols.size();
		ArrayList<String> knockouts = new ArrayList<String>(getReactionFrequencies(sols).keySet());
		int[][] res = new int[nSols][maxsize];
		for (int i = 0; i < res.length; i++) {
			for (int j = 0; j < res[0].length; j++) {
				if (sols.get(i).contains(knockouts.get(j))) {
					res[i][j] = 1;
				}
			}
		}
		return res;
	}



}
