package pt.uminho.ceb.biosystems.mcslibrary.solution;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TreeMap;

import org.paukov.combinatorics.CombinatoricsVector;
import org.paukov.combinatorics.Factory;
import org.paukov.combinatorics.Generator;
import org.paukov.combinatorics.ICombinatoricsVector;

import pt.uminho.ceb.biosystems.mcslibrary.enumeration.implementation.DefaultEnumerationResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.AbstractMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.ReactionGroup;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.compression.alg.MatrixTools;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.CompressedMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.MapUtils;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Pair;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;
import pt.uminho.ceb.biosystems.mew.utilities.java.StringUtils;

public class SolutionUtilities {

	public static final String maxBiomass = "PFBA_BIOMASS";
	public static final String pfbaBPCY = "PFBA_BPCY";
	public static final String pfbaCYIELD = "PFBA_CYIELD";
	public static final String pfbaPYIELD = "PFBA_PYIELD";
	public static final String pfbaSUPT = "PFBA_SUBSTRATE";
	public static final String pfbaPFLUX = "PFBA_PRODUCT";
	public static final String pointOfCoupling = "FVA_POC";

	public static List<List<String>> tokenizeFile(String filename, String delimiter) throws IOException {
		List<List<String>> res = new ArrayList<List<String>>();
		List<String> strs = Utilities.readLines(filename);
		for (String string : strs) {
			res.add(tokenize(string, delimiter));
		}
		return res;
	}

	public static <E extends Object> boolean isEqual(List<E> from, List<E> to) {
		return (getCommonElements(from, to).size() == from.size()) && (from.size() == to.size());
	}

	public static <E extends Object> double[][] getDistanceMatrix(List<List<E>> from, List<List<E>> to) {
		double[][] distances = new double[from.size()][to.size()];
		for (int i = 0; i < from.size(); i++) {
			for (int j = 0; j < to.size(); j++) {
				distances[i][j] = (double) overlapDistance(from.get(i), to.get(j));
			}
		}
		return distances;

	}

	public static <E extends Object> int containsSameSet(List<E> what, List<List<E>> where) {
		for (int i = 0; i < where.size(); i++) {
			if (getCommonElements(what, where.get(i)).size() == what.size()) {
				return i;
			}
		}
		return -1;
	}

	public static <E extends Object> int isSameSet(List<E> what, List<List<E>> where) {
		for (int i = 0; i < where.size(); i++) {
			if (getCommonElements(what, where.get(i)).size() == what.size() && what.size() == where.get(i).size()) {
				return i;
			}
		}
		return -1;
	}

	public static FluxBound[] asFluxBound(List<String> s, AbstractMetabolicNetwork mn) {
		FluxBound[] fx = null;
		if (mn.getClass() == DefaultMetabolicNetwork.class) {
			return asFluxBound(s, (DefaultMetabolicNetwork) mn);
		} else if (mn.getClass() == CompressedMetabolicNetwork.class) {
			return asFluxBound(s, (CompressedMetabolicNetwork) mn);
		}
		return fx;
	}

	private static FluxBound[] asFluxBound(List<String> s, DefaultMetabolicNetwork mn) {
		FluxBound[] fb = new FluxBound[s.size()];
		for (int i = 0; i < s.size(); i++) {
			fb[i] = new FluxBound(mn.getReaction(s.get(i)), 0, 0);
		}
		return fb;
	}

	public static FluxBound[] asFluxBound(String[] s, DefaultMetabolicNetwork mn) {
		FluxBound[] fb = new FluxBound[s.length];
		for (int i = 0; i < s.length; i++) {
			fb[i] = new FluxBound(mn.getReaction(s[i]), 0, 0);
		}
		return fb;
	}

	private static FluxBound[] asFluxBound(List<String> s, CompressedMetabolicNetwork mn) {
		FluxBound[] fb = new FluxBound[s.size()];
		for (int i = 0; i < s.size(); i++) {
			ReactionGroup rg = mn.getReactionGroupFromReaction(s.get(i));
			Reaction r = rg.getReaction(mn.containsReaction(s.get(i)));
			fb[i] = new FluxBound(r, 0, 0);
		}
		return fb;
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

	public static <E extends Object> boolean containsSubset(List<E> superset, List<E> subset) {
		return getCommonElements(superset, subset).size() == subset.size();
	}

	public static <E extends Object> boolean containsElement(List<List<E>> list, List<E> element) {
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
		// return getExceedingElements(from, to).size() +
		// getExceedingElements(to, from).size();
		return Math.max(getExceedingElements(from, to).size(), getExceedingElements(to, from).size());
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

	public static <E extends Object> List<E> getCommonElements(List<List<E>> sets) {
		List<E> set = getCommonElements(sets.get(0), sets.get(1));
		for (int i = 2; i < sets.size(); i++) {
			set = getCommonElements(set, sets.get(i));
		}
		return set;
	}

	public static <E extends Object> List<E> getCommonElements(Collection<E> from, Collection<E> to) {
		List<E> res = new ArrayList<E>();
		for (E f : from) {
			for (E t : to) {
				if (f.equals(t)) {
					res.add(f);
				}
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
		StringTokenizer tok = new StringTokenizer(string, delimiter);
		while (tok.hasMoreTokens()) {
			res.add(tok.nextToken());
		}
		return res;

	}

	public static Map<String, Map<String, List<List<String>>>> folderImport(String rootFolder, String endPattern,
			String delimiter) throws IOException {
		File f = new File(rootFolder);
		File[] children = f.listFiles();
		Map<String, Map<String, List<List<String>>>> fmap;
		if (f.isDirectory()) {
			fmap = new HashMap<>();
			for (int i = 0; i < children.length; i++) {
				File child = children[i];
				
				if (child.isDirectory()) {
					String product = child.getName();
					HashMap<String, List<List<String>>> map = new HashMap<>();
					File[] solutions = child.listFiles();

					for (int j = 0; j < solutions.length; j++) {
						File sfile = solutions[j];
						String substr = sfile.getName();
						if (substr.length() < endPattern.length()) {
							continue;
						}
						String id = substr.substring(substr.length()-endPattern.length(), substr.length());
						String name = substr.substring(0, substr.length()-endPattern.length());
//						System.out.println(id);
						if (sfile.isFile() && id.contains(endPattern)) {
							List<List<String>> solution = SolutionUtilities.tokenizeFile(sfile.getAbsolutePath(),
									delimiter);
							map.put(name, solution);
						}
					}
					
					fmap.put(product, map);
				}
			}
		} else {
			System.out.println("Not a directory!");
			return null;
		}
		return fmap;

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
		Comparator<K> valueComparator = new Comparator<K>() {
			public int compare(K k1, K k2) {
				int compare = map.get(k2).compareTo(map.get(k1));
				if (compare == 0)
					return 1;
				else
					return compare;
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
//		Map<String, Integer> finalmap = sortByValues(map);
		return map;
	}

	public static Map<List<String>, Integer> getKmerFrequencies(List<List<String>> sls, int size) {
		Map<List<String>, Integer> map = new HashMap<List<String>, Integer>();
		for (List<String> solution : sls) {
			CombinatoricsVector<String> vector = new CombinatoricsVector<String>(solution);
			Generator<String> gen = Factory.createSimpleCombinationGenerator(vector, size);
			for (ICombinatoricsVector<String> combination : gen) {
				List<String> list = combination.getVector();
				if (map.containsKey(list)) {
					map.put(list, map.get(list) + 1);
				} else {
					map.put(list, 1);
				}
			}
		}
		return sortByValues(map);
	}

	public static List<Pair<Boolean, List<String>>> getReactionAdditionCombinations(List<String> fromSol,
			List<String> toSol) {
		List<String> core = getCommonElements(fromSol, toSol);
		List<String> exc = getExceedingElements(toSol, fromSol);
		exc.addAll(getExceedingElements(fromSol, toSol));

		CombinatoricsVector<String> possibilityVector = new CombinatoricsVector<String>(exc);
		List<Pair<Boolean, List<String>>> res = new ArrayList<Pair<Boolean, List<String>>>();

		for (int i = 0; i < possibilityVector.getSize(); i++) {
			Generator<String> g = Factory.createSimpleCombinationGenerator(possibilityVector, i + 1);
			for (ICombinatoricsVector<String> v : g) {
				List<String> sol = new ArrayList<String>();
				sol.addAll(core);
				sol.addAll(v.getVector());
				res.add(new Pair<Boolean, List<String>>(v.getSize() == toSol.size(), sol));
			}
		}
		return res;
	}

	public static int[][] getReactionKnockoutMatrix(List<List<String>> sols) {
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

	public static List<List<String>> getKSubsets(List<String> sol, int k) {
		CombinatoricsVector<String> possibilityVector = new CombinatoricsVector<String>(sol);
		Generator<String> g = Factory.createSimpleCombinationGenerator(possibilityVector, k);
		List<ICombinatoricsVector<String>> all = g.generateAllObjects();
		List<List<String>> combs = new ArrayList<>();
		for (ICombinatoricsVector<String> list : all) {
			combs.add(list.getVector());
		}
		return combs;

	}

	public static Set<Set<String>> getKSubsets(Set<String> sol, int k) {
		CombinatoricsVector<String> possibilityVector = new CombinatoricsVector<String>(sol);
		Generator<String> g = Factory.createSimpleCombinationGenerator(possibilityVector, k);
		List<ICombinatoricsVector<String>> all = g.generateAllObjects();
		Set<Set<String>> combs = new HashSet<>();
		for (ICombinatoricsVector<String> list : all) {
			combs.add(new HashSet<String>(list.getVector()));
		}
		return combs;

	}

	public static List<List<String>> getKSubsets(List<String> sol, int from, int to) {
		CombinatoricsVector<String> possibilityVector = new CombinatoricsVector<String>(sol);
		List<List<String>> combs = new ArrayList<>();
		for (int i = from; i <= to; i++) {
			Generator<String> g = Factory.createSimpleCombinationGenerator(possibilityVector, i);
			List<ICombinatoricsVector<String>> all = g.generateAllObjects();
			for (ICombinatoricsVector<String> list : all) {
				combs.add(list.getVector());
			}
		}

		return combs;

	}
	
	public static void writeFrequencyDataset(String filename, Map<String, List<List<String>>> solMap, boolean writeTotal) throws IOException{
		Map<String,Integer> general = new TreeMap<>();
		
		Map<String,Map<String, Integer>> sep = new TreeMap<>();
		for (String column : solMap.keySet()) {
			Map<String, Integer> map = getReactionFrequencies(solMap.get(column));
//			System.out.println(map);
			for (String rx : map.keySet()) {
				if (!general.containsKey(rx)) {
					general.put(rx, map.get(rx));
				} else {
					general.put(rx, map.get(rx) + general.get(rx));
				}
			}
			if (writeTotal) {
//				map.put("TOTAL", value);
				if (!general.containsKey("TOTAL")) {
					general.put("TOTAL", solMap.get(column).size());
				} else {
					general.put("TOTAL", solMap.get(column).size() + general.get("TOTAL"));
				}
				map.put("TOTAL", solMap.get(column).size());
				System.out.println(map.get("TOTAL"));
			}
			sep.put(column, map);
		}
		sep.put("TOTAL", general);
		

		
		MapUtils.sortByValue(general);

		double[][] data = new double[general.keySet().size()][sep.size()];
		List<String> rows = new ArrayList<String>(general.keySet());
		List<String> cols = new ArrayList<String>(sep.keySet());
		
		
//		System.out.println(rows);
//		System.out.println(cols);
		for (int i = 0; i < rows.size(); i++) {
//			System.out.println(sep.get(rows.get(i)));
			for (int j = 0; j < cols.size(); j++) {
				Integer res = sep.get(cols.get(j)).get(rows.get(i));
				data[i][j] = (res == null) ? 0 : res;
			}
		}
		MatrixTools.writeCSV(data, filename, cols, rows);
	}
	
	public static void writeSubsetFrequencyDataset(Map<String, List<List<String>>> solMap, String filename, int k) throws IOException{
		Set<String> tot = new HashSet<String>();
		tot.add("TOTAL");
		int sum = 0;
		Map<Set<String>,Integer> subsets = new LinkedHashMap<>();
		Map<String, Map<Set<String>, Integer>> map = new LinkedHashMap<>();
		for (String solsetID : solMap.keySet()) {
			Map<Set<String>, Integer> smap = new LinkedHashMap<>();
			for (List<String> set : solMap.get(solsetID)) {
				List<List<String>> subsubsets = getKSubsets(set,k);
				for (List<String> list : subsubsets) {
					Set<String> sset = new HashSet<>(list);
					if (!subsets.containsKey(sset)) {
						subsets.put(sset, 1);
						smap.put(sset, 1);
					} else {
						subsets.put(sset, subsets.get(sset)+1);
						if (!smap.containsKey(sset)) {
							smap.put(sset, 1);
						} else {
							smap.put(sset, smap.get(sset)+1);

						}
					}
				}
			}
			smap.put(tot, solMap.get(solsetID).size());
			sum += solMap.get(solsetID).size();
			System.out.println(tot+" "+solsetID+" "+solMap.get(solsetID).size());
			map.put(solsetID, smap);
		}
		subsets.put(tot, sum);
		map.put("TOTAL", subsets);
//		System.out.println(map.keySet());

//		MapUtils.sortByValue(map);
		MapUtils.sortByValue(subsets);
		
		double[][] data = new double[subsets.keySet().size()][map.size()];
		List<String> rows = new ArrayList<String>();
		List<Set<String>> krows = new ArrayList<Set<String>>();
		for (Set<String> string : subsets.keySet()) {
			rows.add(StringUtils.concat(" ", string));
			krows.add(string);
		}
		List<String> cols = new ArrayList<String>(map.keySet());
		
		for (int i = 0; i < rows.size(); i++) {
			for (int j = 0; j < cols.size(); j++) {
				Integer res = map.get(cols.get(j)).get(krows.get(i));
				data[i][j] = (res == null) ? 0 : res;
			}
		}
		MatrixTools.writeCSV(data, filename, cols, rows);
	}

}
