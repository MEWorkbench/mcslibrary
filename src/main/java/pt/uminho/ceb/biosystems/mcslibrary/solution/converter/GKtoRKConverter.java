package pt.uminho.ceb.biosystems.mcslibrary.solution.converter;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.Container;
import pt.uminho.ceb.biosystems.mew.utilities.grammar.syntaxtree.AbstractSyntaxTree;
import pt.uminho.ceb.biosystems.mew.utilities.grammar.syntaxtree.AbstractSyntaxTreeNode;
import pt.uminho.ceb.biosystems.mew.utilities.grammar.syntaxtree.Environment;
import pt.uminho.ceb.biosystems.mew.utilities.grammar.syntaxtree.IEnvironment;
import pt.uminho.ceb.biosystems.mew.utilities.math.language.mathboolean.BooleanValue;
import pt.uminho.ceb.biosystems.mew.utilities.math.language.mathboolean.DataTypeEnum;
import pt.uminho.ceb.biosystems.mew.utilities.math.language.mathboolean.IValue;

public class GKtoRKConverter {
	private DefaultMetabolicNetwork dmn;
	private Container cont;

	public GKtoRKConverter(DefaultMetabolicNetwork metaNet, Container cont) {
		this.dmn = metaNet;
		this.cont = cont;
	}
	
	public List<List<String>> convertGKtoRK(List<List<String>> geneSolutions){
		List<List<String>> result = new ArrayList<>();
		for (int i = 0; i < geneSolutions.size(); i++) {
			IEnvironment<IValue> env = new Environment<>();
			for (int j = 0; j < geneSolutions.get(i).size(); j++) {
				env.associate(geneSolutions.get(i).get(j), new BooleanValue(false));
			}
			List<String> koRules = new ArrayList<>();
			System.out.println(geneSolutions.get(i));
			for (String rx : cont.getReactions().keySet()) {
				AbstractSyntaxTree<DataTypeEnum, IValue> tree = cont.getReaction(rx).getGeneRule();
//				System.out.println(tree);
//				System.out.println(env);
				if (tree != null) {
					try {
						Boolean eval = (Boolean) tree.evaluate(env).getValue();
						if (!eval) {
//							System.out.println(eval);
							koRules.add(rx);
						}	
					} catch (Exception e) {
						// TODO: handle exception
					}
				}
			}
			result.add(koRules);
		}
		return result;
	}
	
	public Map<Set<String>, List<String>> mapGKtoRK(Set<Set<String>> geneSolutions){
		Map<Set<String>, List<String>> map = new HashMap<Set<String>, List<String>>();
		for (Set<String> set : geneSolutions) {
			map.put(set, convertSingleGKtoRK(set));
		}
		return map;
	}
	
	public List<String> convertSingleGKtoRK(Set<String> geneSolutions){
		IEnvironment<IValue> env = new Environment<>();
		for (String string : cont.getGenes().keySet()) {
			if (geneSolutions.contains(string)) {
				env.associate(string, new BooleanValue(false));
			} else {
				env.associate(string, new BooleanValue(true));
			}
		}

		List<String> koRules = new ArrayList<>();
		for (String rx : cont.getReactions().keySet()) {
			AbstractSyntaxTree<DataTypeEnum, IValue> tree = cont.getReaction(rx).getGeneRule();
			if (tree != null) {
				try {
					Boolean eval = (Boolean) tree.evaluate(env).getValue();

					if (!eval) {
						koRules.add(rx);
					}	
				} catch (Exception e) {
					if (tree.getRootNode() != null) {
						e.printStackTrace();
					}
				}
			}
		}
		return koRules;
	}

}
