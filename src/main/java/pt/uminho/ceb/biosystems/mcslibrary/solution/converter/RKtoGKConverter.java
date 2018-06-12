package pt.uminho.ceb.biosystems.mcslibrary.solution.converter;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import pt.uminho.ceb.biosystems.mew.biocomponents.container.Container;

public class RKtoGKConverter {
	private Container cont;
	private ConvertReactionToGenesKos converter;

	public RKtoGKConverter(Container cont) {
		this.cont = cont;
		this.converter = new ConvertReactionToGenesKos(cont);

	}
	
	public TreeSet<String> getMinimal(Collection<String> rkColl) throws Exception{
//		Map<String, Collection<Set<String>>> gkos = converter.convertToGeneKos(rkColl);
		TreeSet<String> minimal = converter.minimalSet(rkColl);
		return minimal;
	}
	
	public List<List<String>> getMinimalGeneKOs(List<List<String>> solutions) throws Exception{
		List<List<String>> res = new ArrayList<List<String>>();
		for (int i = 0; i < solutions.size(); i++) {
			try {
				ArrayList<String> gko = new ArrayList<>(getMinimal(solutions.get(i)));
//				System.out.println(solutions.get(i)+"="+gko);
				res.add(gko);
			} catch (Exception e) {
//				System.out.println("error");
				res.add(new ArrayList<String>());
			}

		}
		return res;
	}
	
	public Set<Set<String>> getAllGeneKOs(List<String> solution) throws Exception{
//		System.out.println("Processing solution..."+solution);
		try {
			Map<String, Collection<Set<String>>> converted = converter.convertToGeneKos(solution);
//			System.out.println("Mapping completed... Printing result");
//			System.out.println("Given: "+solution+"; Converted: "+converted);
			return multiplyCombinations(converted);
		} catch (Exception e) {
			// TODO: handle exception
		}
		return null;

	}
	
	public Set<Set<String>> getMultipleGeneKOs(List<List<String>> solution) throws Exception{
		Set<Set<String>> master = new HashSet<Set<String>>();
//		int sum = 0;
		for (int i = 0; i < solution.size(); i++) {
			Set<Set<String>> kos = getAllGeneKOs(solution.get(i));
			if (kos != null) {
				master.addAll(kos);
			}
//			sum += kos.size();
		}
//		System.out.println("Repeated: "+sum+"; Set result: "+master.size());
		return master;
	}
	
	public Set<Set<String>> multiplyCombinations(Map<String,Collection<Set<String>>> master){
		Set<Set<String>> s = new HashSet<>();
//		System.out.println("Printing combination map");
//		System.out.println(master);
		List<String> options = new ArrayList<String>(master.keySet());
		_multiplyCombinations(master, options, new HashSet<String>(), s);
		return s;
	}
	
	private void _multiplyCombinations(Map<String,Collection<Set<String>>> master, List<String> options, Set<String> toAppend, Set<Set<String>> toStore){
		if (options.isEmpty()) {
//			System.out.println("No more options");
			toStore.add(toAppend);
		} else {
			String currentOption = options.get(0);
//			System.out.println("Unfolding option "+currentOption);
			Collection<Set<String>> sets = master.get(currentOption);
			List<String> theseOptions = new ArrayList<>(options);
			theseOptions.remove(0);
//			System.out.println(theseOptions.size()+" remaining options. "+theseOptions);
//			System.out.println(sets.size()+" possible knockout sets for option.");
			for (Set<String> set : sets) {
				Set<String> newSet = new HashSet<>(toAppend);
				newSet.addAll(set);
//				System.out.println(set);
				_multiplyCombinations(master, theseOptions, newSet, toStore);
			}
		}
	}
	
	public Set<Set<String>> filterCriticals(Set<Set<String>> gkSolutions, Set<String> criticals, int maxSize){
		Set<Set<String>> s = new HashSet<Set<String>>();
		for (Set<String> solution : gkSolutions) {
			if (solution.size() > maxSize) {
				continue;
			}
			boolean valid = true;
			for (String critical : criticals) {
				valid &= solution.contains(critical) ? false : true;
			}
			if (valid) {
				s.add(solution);
			}
		}
		return s;
	}
	
}
