package pt.uminho.ceb.biosystems.mcslibrary.solution.converter;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.stream.XMLStreamException;

import org.xml.sax.SAXException;

import pt.uminho.ceb.biosystems.mcslibrary.solution.SolutionUtilities;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.Container;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.components.ReactionCI;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.io.readers.ErrorsException;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.io.readers.JSBMLReader;
import pt.uminho.ceb.biosystems.mew.biocomponents.validation.io.JSBMLValidationException;
import pt.uminho.ceb.biosystems.mew.utilities.datastructures.collection.CollectionUtils;
import pt.uminho.ceb.biosystems.mew.utilities.datastructures.map.MapUtils;
import pt.uminho.ceb.biosystems.mew.utilities.grammar.syntaxtree.AbstractSyntaxTree;
import pt.uminho.ceb.biosystems.mew.utilities.grammar.syntaxtree.AbstractSyntaxTreeNode;
import pt.uminho.ceb.biosystems.mew.utilities.math.language.mathboolean.DataTypeEnum;
import pt.uminho.ceb.biosystems.mew.utilities.math.language.mathboolean.IValue;
import pt.uminho.ceb.biosystems.mew.utilities.math.language.mathboolean.node.And;
import pt.uminho.ceb.biosystems.mew.utilities.math.language.mathboolean.node.Or;

public class ConvertReactionToGenesKos {
	
	Map<String, AbstractSyntaxTree<DataTypeEnum, IValue>> rules;
	Map<String, Set<String>> reactionGenes;
	
	Set<String> reactionsWithoutGPRs;
	
	
	public ConvertReactionToGenesKos(Container cont){
		initFromContainer(cont);
	}
	
	
	public Map<String, Collection<Set<String>>> convertToGeneKos(Collection<String> reactions) throws Exception{
		
		Map<String, Collection<Set<String>>> fRes = new HashMap<String, Collection<Set<String>>>();
		
		for(String rid : reactions){
			
			if(!rules.containsKey(rid)) {
				throw new Exception("Reaction " + rid + " does not have gene rule!");
//				return null;				
			}
			Set<Set<String>> values = getVariablesToSencetree(rules.get(rid), false);
			fRes.put(rid, values);
		}
		
		return fRes;	
	}
	
	
	public TreeSet<String> minimalSet(Collection<String> reactionsToKo) throws Exception{
		
		Map<String, Collection<Set<String>>> fRes = convertToGeneKos(reactionsToKo);
		Map<String, Integer> c = MapUtils.countMapCollectionSize(fRes);
		
		Map<Integer, Set<String>> stats = MapUtils.revertMap(c);
		
//		MapUtils.prettyPrint(stats);
		TreeSet<String> ole = new TreeSet<String>();
		
		for(Integer key: stats.keySet()){
			Set<String> reactions = stats.get(key);
			
			for(String r : reactions){
				
				Iterator<Set<String>> it = fRes.get(r).iterator();
				Set<String> genes = it.next();
				int common = CollectionUtils.getIntersectionValues(ole, genes).size();
				while(it.hasNext()){
					Set<String> genesaux = it.next();
					int commonAux = CollectionUtils.getIntersectionValues(ole, genesaux).size();
					if(commonAux > common){
						common = commonAux;
						genes = genesaux;
					}
				}
				
				ole.addAll(genes);
			}
			
		}
		
		return ole;
	}
	
	
	private void initFromContainer(Container cont){
		rules = new HashMap<String, AbstractSyntaxTree<DataTypeEnum, IValue>>();
		reactionGenes = new HashMap<String, Set<String>>();
		reactionsWithoutGPRs = new HashSet<String>();
		
		for(ReactionCI r : cont.getReactions().values()){
			
			Set<String> gene = r.getGenesIDs();
			AbstractSyntaxTree<DataTypeEnum, IValue> rule = r.getGeneRule();
			String id = r.getId();
			
			if(gene.size()>0){
				
				rules.put(id, rule);
				reactionGenes.put(id, gene);				
			}else{
				reactionsWithoutGPRs.add(id);
			}
			
			
		}
	}
	
	
	private Set<Set<String>> getVariablesToSencetree(
			AbstractSyntaxTree<DataTypeEnum, IValue> geneRule, boolean sence) throws Exception {
		
		AbstractSyntaxTreeNode<DataTypeEnum, IValue> node = geneRule.getRootNode();
		return getVariablesToSenceNode(node, sence);
	}
	
	static public Set<Set<String>> getVariablesToSenceNode(AbstractSyntaxTreeNode<DataTypeEnum, IValue> node, boolean sence) throws Exception{
		
		Set<Set<String>> ret = new HashSet<Set<String>>();
		
		if(node instanceof And){
			if(node.getNumberOfChildren()>2)
				throw new Exception("ERRO");
			
			Set<Set<String>> res1 = getVariablesToSenceNode(node.getChildAt(0), sence);
			Set<Set<String>> res2 = getVariablesToSenceNode(node.getChildAt(1), sence);
			
			if(sence)
				ret = propagate(res1, res2);
			else
				ret = aglumerate(res1, res2);
						
		}else if(node instanceof Or){
			if(node.getNumberOfChildren()>2)
				throw new Exception("ERRO");
			
			Set<Set<String>> res1 = getVariablesToSenceNode(node.getChildAt(0), sence);
			Set<Set<String>> res2 = getVariablesToSenceNode(node.getChildAt(1), sence);
			
			if(sence)
				ret = aglumerate(res1, res2);
			else
				ret = propagate(res1, res2);
			
			
//		}else if(node instanceof Not){
//			ret = getVariablesToSenceNode(node.getChildAt(0), !sence);
//		
		}else if(node.isLeaf()){
		
			Set<String> set = new HashSet<String>();
			set.add(node.toString());
			ret.add(set);
		}
//		System.out.println(node.getClass().getSimpleName()+" " + ret);
		return ret;
	}
	
	
	private static Set<Set<String>> aglumerate(Set<Set<String>> res1,
			Set<Set<String>> res2) {
		
		Set<Set<String>> ret = new HashSet<Set<String>>();
		
		ret.addAll(res1);
		ret.addAll(res2);
		
		return ret;
	}


	private static Set<Set<String>> propagate(Set<Set<String>> res1,
			Set<Set<String>> res2) {
		Set<Set<String>> ret = new HashSet<Set<String>>();
		
		
		for(Set<String> r1 : res1)
			for(Set<String> r2: res2)
				ret.add(CollectionUtils.getReunionValues(r1, r2));
		
		return ret;
	}
	
	public static void main(String[] args) throws Exception {
		String sFile = "/home/skapur/myCloud/Projectos/DeYeast/Solutions/SolutionAnalysis/TEST/MCS_ATP.txt";
		String mFile = "/home/skapur/myCloud/Projectos/DeYeast/FeistReplication/Models/iMM904/iMM904_peroxisome.xml";
		List<List<String>> sols = SolutionUtilities.tokenizeFile(sFile, ",");
		Container cont = new Container(new JSBMLReader(mFile,"a",false));
		ConvertReactionToGenesKos crtgk = new ConvertReactionToGenesKos(cont);
		RKtoGKConverter c = new RKtoGKConverter(cont);
		Set<Set<String>> kos = c.getMultipleGeneKOs(sols);
		System.out.println(kos.size());
		
	}
	
	
}
