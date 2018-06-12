package pt.uminho.ceb.biosystems.mcslibrary.solution.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ilog.concert.IloException;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba.CPLEXFluxBalanceAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba.CPLEXParsimoniousFluxBalanceAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.solution.SolutionUtilities;
import pt.uminho.ceb.biosystems.mcslibrary.solution.converter.GKtoRKConverter;
import pt.uminho.ceb.biosystems.mcslibrary.solution.converter.RKtoGKConverter;
import pt.uminho.ceb.biosystems.mcslibrary.solution.scoring.SolutionScorer;
import pt.uminho.ceb.biosystems.mcslibrary.solution.scoring.scoreitems.BPCYScoreItem;
import pt.uminho.ceb.biosystems.mcslibrary.solution.scoring.scoreitems.CarbonYieldScoreItem;
import pt.uminho.ceb.biosystems.mcslibrary.solution.scoring.scoreitems.IScoreItem;
import pt.uminho.ceb.biosystems.mcslibrary.solution.scoring.scoreitems.PFBAFluxValueItem;
import pt.uminho.ceb.biosystems.mcslibrary.solution.scoring.scoreitems.RobustnessScoreItem;
import pt.uminho.ceb.biosystems.mcslibrary.solution.scoring.scoreitems.SolutionSizeItem;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.MCSPipeline;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.SteadyStateModelReader;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.Container;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.io.readers.JSBMLReader;

public class SolutionFilterPipeline {
	private RKtoGKConverter rxToGene;
	private SolutionFilter sf;
	private GKtoRKConverter conv;

	public SolutionFilterPipeline(DefaultMetabolicNetwork dmn, Container cont, SolutionFilter sf) {
		this.rxToGene = new RKtoGKConverter(cont);
		this.conv = new GKtoRKConverter(dmn, cont);
		this.sf = sf;
	}
	
	
	public List<List<String>> filterSolutions(List<List<String>> solutions, RobustnessCriteria[] rcrit, FluxBound[] fcrit) throws IloException{
		
		System.out.println("Filtering solutions... (environment compliant)");
		List<List<String>> boundCompliant = sf.filterBounds(fcrit, solutions);
		System.out.println("Filtering "+boundCompliant.size()+" solutions... (robustness compliant)");
		List<List<String>> robustnessCompliant = sf.filterRobust(fcrit, rcrit, boundCompliant);
		System.out.println("Done");
		return robustnessCompliant;
	}
	
	public Set<Set<String>> generateAndFilterGK(List<List<String>> solutions, int maxGK, Set<String> criticals) throws Exception{
		Set<Set<String>> conv = rxToGene.getMultipleGeneKOs(solutions);
		Set<Set<String>> critFilt = rxToGene.filterCriticals(conv, criticals, maxGK);
		return critFilt;
	}
	
	public Map<Set<String>, List<String>> filterGKOs(Set<Set<String>> geneKOs, RobustnessCriteria[] rcrit, FluxBound[] fcrit) throws IloException{
		Map<Set<String>, List<String>> solutionMap = conv.mapGKtoRK(geneKOs);
		List<List<String>> rxSolutions = new ArrayList<List<String>>(solutionMap.values());
		List<List<String>> valid = filterSolutions(rxSolutions, rcrit, fcrit);
		System.out.println(solutionMap.values().size()+" solutions for "+solutionMap.keySet().size()+" gene knockouts.");
		System.out.println(valid.size()+" reaction knockouts.");
//		List<Set<String>> toRemove = new ArrayList<Set<String>>();
//		for (Set<String> list : solutionMap.keySet()) {
//			List<String> curRKO = solutionMap.get(list);
//			boolean invalid = SolutionUtilities.containsSameSet(curRKO, valid) < 0;
//			if (invalid) {
//				toRemove.add(list);
//			}
//		}
//		System.out.println("Ruling out "+toRemove.size()+" solutions");
		Map<Set<String>, List<String>> finalMap = new HashMap<>();

		for (Set<String> list : solutionMap.keySet()) {
			List<String> rko = solutionMap.get(list);
			boolean found = false;
			for (List<String> string : valid) {
				if (SolutionUtilities.isEqual(string, rko)) {
					found = true;
					break;
				}
			}
			if (found) {
				finalMap.put(list, rko);
			}
			
		}
		return finalMap;
	}
	
	public Map<Set<String>, List<String>> filterPipeline(List<List<String>> solutions, RobustnessCriteria[] rcrit, FluxBound[] fcrit, int maxGK, Set<String> criticals) throws Exception{
		System.out.println("\tDirect filter");
		List<List<String>> filtered = filterSolutions(solutions, rcrit, fcrit);
		System.out.println("\tGene KO step 1");
		Set<Set<String>> GKstep1 = generateAndFilterGK(filtered, maxGK, criticals);
		System.out.println("\tGene KO step 2");
		Map<Set<String>, List<String>> GKstep2 = filterGKOs(GKstep1, rcrit, fcrit);
		return GKstep2;
	}
	
	public static void main(String[] args) throws Exception {
		String rootFolder = "C:/Users/vitor/MEOCloud/Projectos/DeYeast/Solutions/ASAP/";
		String[] products = {"R_EX_succ_e_","R_EX_mal_L_e_","R_EX_fum_e_"};
//		String[] products = {"R_EX_fum_e_","R_EX_mal_L_e_"};
		String mFile = "C:/Users/vitor/MEOCloud/Projectos/DeYeast/FeistReplication/Models/iMM904/iMM904_peroxisome.xml";
		
		Container cont = new Container(new JSBMLReader(mFile, "a", false));
		MCSPipeline p = new MCSPipeline(cont, "");
		DefaultMetabolicNetwork meta = p.getMetabolicNetwork();
		SteadyStateModelReader.updateFormulae(cont, meta);
		Map<String,String> geneDict = Utilities.readMap("C:/Users/vitor/MEOCloud/Projectos/DeYeast/Models/SCEnsemblKey.txt", ",");

		Reaction substrate = meta.getReaction("R_EX_glc_e_");
		Reaction biomass = meta.getReaction("R_biomass_SC5_notrace");
		
		SolutionFilter sft = new SolutionFilter(meta);
		SolutionFilterPipeline sfp = new SolutionFilterPipeline(meta, cont, sft);
		SolutionScorer s = new SolutionScorer();
		
		FluxBound[] fcrit = new FluxBound[]{
				new FluxBound(meta.getReaction("R_EX_glc_e_"), -1.15, -1.15),
				new FluxBound(meta.getReaction("R_biomass_SC5_notrace"), 0.0001, Utilities.INF),
				new FluxBound(meta.getReaction("R_ATPM"), 1, 1)
//				null
		};
		
		FluxBound[] fcritATP = new FluxBound[]{
				new FluxBound(meta.getReaction("R_EX_glc_e_"), -1.15, -1.15),
				new FluxBound(meta.getReaction("R_biomass_SC5_notrace"), 0.0001, Utilities.INF),
				new FluxBound(meta.getReaction("R_ATPM"), 0, 1)
//				null
		};
		
		RobustnessCriteria[] rcrit = new RobustnessCriteria[]{
				null
		};
		
		int maxGK = 10;
		Set<String> criticals = new HashSet<String>(Utilities.readLines("C:/Users/vitor/MEOCloud/Projectos/DeYeast/Models/nontargets#GK#[aerobic#glucose].txt"));
		List<String> possible = new ArrayList<String>();
//		possible.add("MCS_FIXED");
//		possible.add("MCS_ATP");
		possible.add("MCS_WC");
		
		CPLEXFluxBalanceAnalysis fba = new CPLEXFluxBalanceAnalysis(meta);
		CPLEXParsimoniousFluxBalanceAnalysis pfba = new CPLEXParsimoniousFluxBalanceAnalysis(meta, 0.999999);
		
		System.out.println("Filtering...");
		for (int i = 0; i < products.length; i++) {
			String product = products[i];
			File searchFolder = new File(rootFolder+product+"/");
			File[] files = searchFolder.listFiles();
//			fcrit[fcrit.length-1] = new FluxBound(meta.getReaction(product), 0.0001, Utilities.INF);
			double rob = 0.9;
			if (!product.equals("R_EX_mal_L_e_")) {
				rcrit[rcrit.length-1] = new RobustnessCriteria("R_biomass_SC5_notrace", product, 0.01, "max", "min", "someCrit");
			} else {
				rcrit[rcrit.length-1] = new RobustnessCriteria("R_biomass_SC5_notrace", product, rob, "max", "min", "someCrit");
			}
			for (int j = 0; j < files.length; j++) {
				int f2 = files[j].getName().length();
				int f1 = f2-4;
				if (files[j].getName().substring(f1,f2).contains("txt") && possible.contains(files[j].getName().substring(0, f1))) {
					List<List<String>> solutions = SolutionUtilities.tokenizeFile(files[j].getAbsolutePath(), ",");
					List<List<String>> filtered = sfp.filterSolutions(solutions, rcrit, fcritATP);
					Map<Set<String>, List<String>> gkoFiltered = sfp.filterPipeline(filtered, rcrit, fcrit, maxGK, criticals);
					Utilities.writeMap(gkoFiltered, files[j].getAbsolutePath()+".filtered", "=");
				}
			}
		}
		
		IScoreItem[] scorers = new IScoreItem[]{
				new SolutionSizeItem(),
				new PFBAFluxValueItem(meta, pfba, fcrit, biomass, biomass, "max"),
				null,
				new PFBAFluxValueItem(meta, pfba, fcrit, substrate, biomass, "max"),
				null,
				null,
				null,
				null
		};
		
		
		Map<String, Set<Set<String>>> mothermap = new HashMap<>();
		Map<Set<String>, Map<IScoreItem, Double>> superScoreMap = new HashMap<>();
		Map<String,Map<String,IScoreItem>> prodScorers = new LinkedHashMap<>();
		
		System.out.println("Creating chassis...");
		for (int i = 0; i < products.length; i++) {
			Set<Set<String>> productSet = new HashSet<Set<String>>();
			String product = products[i];
			Reaction productRx = meta.getReaction(product);
			Map<String,IScoreItem> psMap = new LinkedHashMap<>();
			
			scorers[2] = new PFBAFluxValueItem(meta, pfba, fcrit, productRx, biomass, "max");
			scorers[4] = new BPCYScoreItem(meta, pfba, fcrit, biomass, productRx, substrate);
			scorers[5] = new RobustnessScoreItem(meta, fba, fcrit, biomass, productRx, 0.1, false);
			scorers[6] = new RobustnessScoreItem(meta, fba, fcrit, biomass, productRx, 0.9, false);
			scorers[7] = new CarbonYieldScoreItem(meta, pfba, fcrit, biomass, productRx, substrate);
			
			psMap.put("Biomass(PFBA)", scorers[1]);
			psMap.put("MinProductFlux@10%", scorers[5]);
			psMap.put("MinProductFlux@90%", scorers[6]);
			psMap.put("CYIELD", scorers[7]);
			
			prodScorers.put(product, psMap);
			
			File searchFolder = new File(rootFolder+product+"/");
			File[] files = searchFolder.listFiles();
			for (int j = 0; j < files.length; j++) {
//				System.out.println(files[j].getName());
				if (files[j].getName().contains(".filtered")) {
					Map<String, String> map = Utilities.readPropertyMap(files[j].getAbsolutePath());
					Map<Set<String>, List<String>> gtkmap = new HashMap<>();
					for (String sol : map.keySet()) {
						Set<String> gs = new HashSet<String>(Utilities.getAllStringTokens(sol.replace("[", "").replace("]", ""), ", "));
						List<String> ks = Utilities.getAllStringTokens(map.get(sol).replace("[", "").replace("]", ""), ", ");
						gtkmap.put(gs, ks);
						productSet.add(gs);
					}
					Map<Set<String>, Map<IScoreItem, Double>> scoreMap = s.getDataset(gtkmap, scorers, searchFolder+files[j].getName()+"_filteredDataset.csv");
					superScoreMap.putAll(scoreMap);
				}
			}
			mothermap.put(product, productSet);
		}
		
	}
}
