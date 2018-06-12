package pt.uminho.ceb.biosystems.mcslibrary.solution.analysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.ctc.wstx.util.StringUtil;

import ilog.concert.IloException;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.compression.alg.MatrixTools;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.solution.analysis.simulation.components.EnvelopeProperties;
import pt.uminho.ceb.biosystems.mcslibrary.solution.analysis.simulation.components.EnvelopeResult;
import pt.uminho.ceb.biosystems.mcslibrary.solution.analysis.simulation.components.VariabilityAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.solution.converter.GKtoRKConverter;
import pt.uminho.ceb.biosystems.mcslibrary.solution.converter.RKtoGKConverter;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.Container;
import pt.uminho.ceb.biosystems.mew.utilities.java.StringUtils;

public class SolutionAnalysisPipeline {
	private DefaultMetabolicNetwork dmn;
	private List<List<String>> solutions;
	private Map<List<String>, SolutionAnalysisResult> res;
	private VariabilityAnalysis var;
	private FluxBound[][] envConds;
	private String name;
	private String[] envNames;
	
	public final static String NAME_PFBA = "pFBA_";
	public final static String NAME_ENVELOPE = "envelope_";
	public final static String NAME_VARIABILITY = "fluxLimits_";
	public final static String NAME_INDEX = "solutions_";
	
	public SolutionAnalysisPipeline(DefaultMetabolicNetwork dmn, List<List<String>> solutions, EnvelopeProperties[] envProp, FluxBound[][] envConds, String name, String[] envNames) throws IloException {
		this.dmn = dmn;
		this.solutions = solutions;
		this.var = new VariabilityAnalysis(dmn, envProp);
		this.envConds = envConds;
		this.name = name;
		this.envNames = envNames;
	}
	
	public Map<List<String>, SolutionAnalysisResult> filterAnalysisPipeline(RobustnessCriteria[] rcrit, FluxBound[] criteria, String critName,String writeRoot) throws IloException, IOException{
		SolutionFilter sf = new SolutionFilter(dmn);
		List<List<String>> valid = new ArrayList<List<String>>(solutions);
		if (rcrit != null) {
			valid = sf.doubleFilter(criteria, rcrit, valid);
		} else {
			valid = sf.filterBounds(criteria, valid);
		}
		System.out.println("Filtering complete. "+valid.size()+" out of "+solutions.size()+" were kept.");
		Map<List<String>, SolutionAnalysisResult> res = new LinkedHashMap<List<String>, SolutionAnalysisResult>();
		for (int i = 0; i < valid.size(); i++) {
			List<String> solution = valid.get(i);
			res.put(solution, new SolutionAnalysisResult(solution, dmn, var, envConds));
		}
		if (writeRoot != null) {
			outputToFiles(writeRoot, res, critName, name);
		}
		return res;
	}
	
	public Map<List<String>,SolutionAnalysisResult> geneKnockoutAnalysis(List<List<String>> solutions, String name, String[] envNames, String criteriaName, FluxBound[] criteria, String writeRoot, Container cont) throws Exception{
//		Map<List<String>, SolutionAnalysisResult> map = new LinkedHashMap<>();
		System.out.println(solutions.size()+" in for conversion to gene knockouts");
//		GKtoRKConverter gkrk = new GKtoRKConverter(dmn, cont);
//		List<List<String>> converted = gkrk.convertGKtoRK(solutions);
		RKtoGKConverter converter = new RKtoGKConverter(cont);
		List<List<String>> converted = converter.getMinimalGeneKOs(solutions);
		System.out.println("Converted="+converted.size());
		Map<List<String>, List<String>> gkrkMap = new LinkedHashMap<>();
		
		GKtoRKConverter reconverter = new GKtoRKConverter(dmn, cont);
		List<List<String>> reconverted = reconverter.convertGKtoRK(converted);
		System.out.println("Reconverted="+reconverted.size());
		for (int i = 0; i < converted.size(); i++) {
			if (reconverted.get(i).size() != 0) {
				gkrkMap.put(converted.get(i), reconverted.get(i));
			}
		}
		
		List<List<String>> finalreconverted = new ArrayList<List<String>>();
		List<List<String>> finalconverted = new ArrayList<List<String>>();
		for (List<String> list : gkrkMap.keySet()) {
			finalreconverted.add(gkrkMap.get(list));
			finalconverted.add(list);
		}
		
		SolutionAnalysisPipeline pp = new SolutionAnalysisPipeline(dmn, finalreconverted, this.var.getEnvelopeProperties(), envConds, name, envNames);
		Map<List<String>, SolutionAnalysisResult> filteredAnalysis = pp.filterAnalysisPipeline(null, criteria, criteriaName, null);
		
		LinkedHashMap<List<String>, SolutionAnalysisResult> sarmap = new LinkedHashMap<>();
		for (int i = 0; i < finalconverted.size(); i++) {
			if (filteredAnalysis.containsKey(finalreconverted.get(i))) {
				sarmap.put(finalconverted.get(i), filteredAnalysis.get(finalreconverted.get(i)));
			}
		}
		System.out.println(sarmap);
		if (sarmap.keySet().size() > 0) {
			outputToFiles(writeRoot, sarmap, criteriaName, name);
		}
		return sarmap;
	}
	
	private int outputToFiles(String root, Map<List<String>, SolutionAnalysisResult> analysis, String filterCriteria, String groupName) throws IOException{
//		Utilities.printMap(analysis);
		List<List<String>> indexedSolutions = new ArrayList<>(analysis.keySet());
		if (indexedSolutions.size() == 0) {
			return 0;
		}
		// WRITE INDEX DATASET
		String indexHeader = "id, solution, size\n";
		String indexDataset = "";
		List<String> index = new ArrayList<>();
		for (int i = 0; i < indexedSolutions.size(); i++) {
			indexDataset += "sol_"+i+","+StringUtils.concat(" ", indexedSolutions.get(i))+","+indexedSolutions.get(i).size()+"\n";
			index.add("sol_"+i);
		}
		Utilities.write(indexHeader+indexDataset.substring(0, indexDataset.length()-1), StringUtils.concat("#", new String[]{root,groupName,NAME_INDEX,filterCriteria,".csv"}));
		
		// WRITE FLUXDIST DATASET
		for (int i = 0; i < envConds.length; i++) {
			double[][] fluxDistributions = new double[indexedSolutions.size()][dmn.getNumOfReactions()];
			System.out.println("Flux distribution matrix "+fluxDistributions.length+" by "+fluxDistributions[0].length);
			for (int j = 0; j < fluxDistributions.length; j++) {
//				for (int k = 0; k < fluxDistributions[0].length; k++) {
				fluxDistributions[j] = analysis.get(indexedSolutions.get(j)).getSimulationResult(envConds[i]).getValues();
//				}
			}
			String[] filename = new String[]{root, groupName, filterCriteria, "_", SolutionAnalysisPipeline.NAME_PFBA, "env_"+envNames[i]};
			MatrixTools.writeCSV(fluxDistributions, StringUtils.concat("#", filename),dmn.getReactionNameList(),index);
		}
		
		// WRITE ENVELOPE DATASET
		EnvelopeProperties[] ep = this.var.getEnvelopeProperties();
		for (int i = 0; i < envConds.length; i++) {
			for (int j = 0; j < ep.length; j++) {
				int steps = ep[j].getSteps();
				List<String> header = new ArrayList<String>();
				for (int k = 0; k <= steps; k++) {
					header.add(String.valueOf(k)+"%_"+ep[j].getPivot());
				}
				double[][] envelopeMin = new double[indexedSolutions.size()][];
				double[][] envelopeMax = new double[indexedSolutions.size()][];
				for (int k = 0; k < indexedSolutions.size(); k++) {
					EnvelopeResult env = analysis.get(indexedSolutions.get(k)).getVariabilityAnalysis(envConds[i]).getEnvelope(j);
//					System.out.println("Envelope with "+env.length()+" points.");
					envelopeMin[k] = new double[env.length()];
					envelopeMax[k] = new double[env.length()];
					for (int m = 0; m < env.length(); m++) {
						envelopeMin[k][m] = env.getMinAtPoint(m);
						envelopeMax[k][m] = env.getMaxAtPoint(m);
					}
				}
				String[] filename1 = new String[]{root, groupName, filterCriteria, "_", SolutionAnalysisPipeline.NAME_ENVELOPE, "min", "env_"+envNames[i], "prop_"+ep[j].getConfName()};
				MatrixTools.writeCSV(envelopeMin, StringUtils.concat("#", filename1),header,index);
				
				String[] filename2 = new String[]{root, groupName, filterCriteria, "_", SolutionAnalysisPipeline.NAME_ENVELOPE, "max", "env_"+envNames[i], "prop_"+ep[j].getConfName()};
				MatrixTools.writeCSV(envelopeMax, StringUtils.concat("#", filename2),header,index);
			}
		}
		
		// WRITE FLUX LIMITS
		for (int i = 0; i < envConds.length; i++) {
			double[][] limitMin = new double[indexedSolutions.size()][dmn.getNumOfReactions()];
			double[][] limitMax = new double[indexedSolutions.size()][dmn.getNumOfReactions()];

			for (int j = 0; j < limitMin.length; j++) {
				for (int k = 0; k < limitMin[0].length; k++) {
					limitMin[j] = analysis.get(indexedSolutions.get(j)).getVariabilityAnalysis(envConds[i]).getMinLimits();
					limitMax[j] = analysis.get(indexedSolutions.get(j)).getVariabilityAnalysis(envConds[i]).getMaxLimits();
				}
			}
			String[] filename1 = new String[]{root, groupName, filterCriteria, "_", SolutionAnalysisPipeline.NAME_VARIABILITY, "env_"+envNames[i]};
			MatrixTools.writeCSV(limitMin, StringUtils.concat("#", filename1),dmn.getReactionNameList(),index);
			
			String[] filename2 = new String[]{root, groupName, filterCriteria, "_", SolutionAnalysisPipeline.NAME_VARIABILITY, "env_"+envNames[i]};
			MatrixTools.writeCSV(limitMax, StringUtils.concat("#", filename2),dmn.getReactionNameList(),index);
		}
		return 1;
	}
	
	
}
