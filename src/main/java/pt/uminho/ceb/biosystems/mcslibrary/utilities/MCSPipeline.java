package pt.uminho.ceb.biosystems.mcslibrary.utilities;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;

import ilog.concert.IloException;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.AbstractEnumerationResult;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.EnumerationProblem;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.Solution;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.implementation.CPLEXEnumerationFilter;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.implementation.CPLEXEnumerationSolver;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.implementation.CPLEXIntegratedEnumerator;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.implementation.CompressedEnumerationResults;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.implementation.DefaultEnumerationResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.compression.NullspaceNetworkCompressor;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.ReactionConstraint;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.YieldConstraint;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fva.CPLEXFluxVariabilityAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fva.FluxVariabilityAnalysisResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.CompressedMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.Container;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.components.ReactionCI;
import pt.uminho.ceb.biosystems.mew.core.model.converters.ContainerConverter;
import pt.uminho.ceb.biosystems.mew.core.model.exceptions.InvalidSteadyStateModelException;
import pt.uminho.ceb.biosystems.mew.core.model.steadystatemodel.ISteadyStateModel;

public class MCSPipeline {

	private double M = 9999;
	private double INF = Utilities.INF;
	private ArrayList<String> exchange;
	private ArrayList<String> inflows;
	private ArrayList<String> outflows;
	private ArrayList<String> spontaneous;
	private DefaultMetabolicNetwork model;
	private ArrayList<FluxBound> undesired_targets;
	private ArrayList<FluxBound> desired_targets;
	private ArrayList<YieldConstraint> undesired_yields;
	private ArrayList<YieldConstraint> desired_yields;
	private ArrayList<String> excluded;
	private ArrayList<FluxBound> fvaconditions;
	private ArrayList<Solution> excludedSols;
	private ArrayList<String> optsingle;
	private CompressedMetabolicNetwork curcomp;
	private HashSet<Solution> seeds;
	private HashSet<String> forced;
	private AlgorithmType aty;
	private String boundaryString;
	private int threads;
	private boolean allowSolverLogging;
	private boolean initializesOnEnumeration;
	private HashMap<String, Object> samplingParams;
	private EnumerationProblem ep;
	
	public MCSPipeline(Container cont, String spontaneous_handle) {
		// TODO: allow any boundary string
		this.initializesOnEnumeration = true;
		this.allowSolverLogging = false;
		this.boundaryString = "_b";
		this.threads = 0;
		this.spontaneous = identifySpontaneousReactions(cont,
				spontaneous_handle);
		try {
			this.model = convertModel(cont, boundaryString);
		} catch (InvalidSteadyStateModelException e1) {
			e1.printStackTrace();
		}
		
		try {
			this.exchange = identifyExchangeReactions(cont, boundaryString).size() != 0 ? identifyExchangeReactions(cont, boundaryString)
					: identifyExchangeReactions(model);
			System.out.println("Exchange reactions :"+this.exchange.size());
		} catch (InvalidSteadyStateModelException e) {
			e.printStackTrace();
		}
		
		this.inflows = new ArrayList<String>();
		this.outflows = new ArrayList<String>();
		this.excluded = new ArrayList<String>();
		this.fvaconditions = new ArrayList<FluxBound>();
		this.excludedSols = new ArrayList<Solution>();
		this.optsingle = new ArrayList<String>();
		for (String exReac : this.exchange) {
			Reaction r = this.model.getReaction(exReac);
			double lb = r.getBounds().getLower();
			double ub = r.getBounds().getUpper();
			if (lb == 0) {
				outflows.add(exReac);
			} else if (ub == 0) {
				inflows.add(exReac);
			}
		}
		this.undesired_targets = new ArrayList<>();
		this.undesired_yields = new ArrayList<>();
		this.desired_targets = new ArrayList<>();
		this.desired_yields = new ArrayList<>();
		this.seeds = new HashSet<Solution>();
		this.forced = new HashSet<String>();
		this.aty = AlgorithmType.DEFAULT;
	}

	public void addFVABound(String reaction, double lb, double ub) {
		Reaction r = model.getReaction(reaction);
		FluxBound fb = new FluxBound(r, lb, ub);
		fvaconditions.add(fb);
	}

	public void setAlgorithmType(AlgorithmType t){
		this.aty = t;
	}
	public ArrayList<String> identifySingleReactions() {
		ArrayList<String> singleReactions = new ArrayList<String>();
		singleReactions.addAll(exchange);
		singleReactions.addAll(spontaneous);
		singleReactions.addAll(optsingle);
		// singleReactions.addAll(excluded);
		return singleReactions;
	}
	
	public CompressedMetabolicNetwork getCompressedMetabolicNetwork() {
		return curcomp;
	}
	
	public void addSingleReaction(String reaction){
		optsingle.add(reaction);
	}
	public void correctCapacities() {
		for (int i = 0; i < model.getNumOfReactions(); i++) {
			Reaction r = model.getReaction(i);
			double lb = r.getBounds().getLower();
			double ub = r.getBounds().getUpper();
			if (ub >= M) {
				ub = INF;
			}

			if (lb <= -M) {
				lb = -INF;
			}
			r.setBounds(new ReactionConstraint(lb, ub));
			model.setReaction(r, i);
		}
	}

	public void correctInOutflows() {
		for (String id : inflows) {
			Reaction r = model.getReaction(id);
			r.setBounds(new ReactionConstraint(r.getBounds().getLower(), 0));
		}
		for (String id : outflows) {
			Reaction r = model.getReaction(id);
			r.setBounds(new ReactionConstraint(0, r.getBounds().getUpper()));
		}
	}

	public void correctExchangeReactions() {
		for (String reaction : exchange) {
			Reaction r = model.getReaction(reaction);
			int rId = model.getReactionIndex(reaction);

			double lb = model.getLowerBound(rId);
			double ub = model.getUpperBound(rId);

			boolean xoutflow = (lb == 0);
			boolean xinflow = (ub == 0);

			if (xoutflow) {
				lb = 0;
			}

			if (xinflow) {
				ub = 0;
			}

			r.setBounds(new ReactionConstraint(lb, ub));
//			model.setReaction(r, rId);
		}
	}

	public void addFluxBound(String reaction, double lb, double ub,
			boolean block) {
		Reaction r = model.getReaction(reaction);
		FluxBound fb = new FluxBound(r, lb, ub);
		if (block) {
			undesired_targets.add(fb);
		} else {
			desired_targets.add(fb);
		}
	}
	
	public void addOverridingFluxBound(String reaction, double lb, double ub,
			boolean block) {
		Reaction r = model.getReaction(reaction);
		FluxBound fb = new FluxBound(r, lb, ub);
		fb.setOverride(true, true);
		if (block) {
			undesired_targets.add(fb);
		} else {
			desired_targets.add(fb);
		}
	}

	public void ignoreReaction(String reaction) {
		if (model.containsReaction(reaction) > -1)
			excluded.add(reaction);
	}

	public ArrayList<?> getConstraints(String ctype, boolean toBlock) {
		switch (ctype) {
		case "flux": {
			return toBlock ? undesired_targets : desired_targets;
		}
		case "yield": {
			return toBlock ? undesired_yields : desired_yields;
		}
		default:
			return null;
		}
	}

	public ArrayList<?> addConstraint(String ctype, boolean toBlock, Object constraint) {
		switch (ctype) {
		case "flux": {
			if (toBlock)
				undesired_targets.add((FluxBound)constraint);
			else
				desired_targets.add((FluxBound)constraint);
		}
		case "yield": {
			if (toBlock)
				undesired_yields.add((YieldConstraint)constraint);
			else
				desired_yields.add((YieldConstraint)constraint);
		}
		default:
			return null;
		}
	}


	public void resetConstraints() {
		desired_targets.clear();
		desired_yields.clear();
		undesired_targets.clear();
		undesired_yields.clear();
//		fvaconditions.clear();
//		excluded.clear();
		forced.clear();
		
	}

	public void addLowerYieldConstraint(String uptake, String product,
			double ratio, boolean block) {
		Reaction u = model.getReaction(uptake);
		Reaction p = model.getReaction(product);
		YieldConstraint yc = new YieldConstraint(u, p, ratio);
		yc.setAsLower(true);

		if (block) {
			undesired_yields.add(yc);
		} else {
			desired_yields.add(yc);
		}
	}

	public void addLowerYieldConstraint(String uptake, String product,
			double ratio, boolean block, double deviation) {
		Reaction u = model.getReaction(uptake);
		Reaction p = model.getReaction(product);
		YieldConstraint yc = new YieldConstraint(u, p, ratio);
		yc.setDeviation(deviation);
		yc.setAsLower(true);

		if (block) {
			undesired_yields.add(yc);
		} else {
			desired_yields.add(yc);
		}
	}
	
	
	public void addUpperYieldConstraint(String uptake, String product,
			double ratio, boolean block) {
		Reaction u = model.getReaction(uptake);
		Reaction p = model.getReaction(product);
		YieldConstraint yc = new YieldConstraint(u, p, ratio);
		yc.setAsLower(false);

		if (block) {
			undesired_yields.add(yc);
		} else {
			desired_yields.add(yc);
		}
	}

	public void addUpperYieldConstraint(String uptake, String product,
			double ratio, boolean block, double deviation) {
		Reaction u = model.getReaction(uptake);
		Reaction p = model.getReaction(product);
		YieldConstraint yc = new YieldConstraint(u, p, ratio);
		yc.setDeviation(deviation);
		yc.setAsLower(false);

		if (block) {
			undesired_yields.add(yc);
		} else {
			desired_yields.add(yc);
		}
	}
	public void addSolutionSeed(Solution sol){
		this.seeds.add(sol);
	}
	public static DefaultMetabolicNetwork convertModel(Container cont, String boundaryString)
			throws InvalidSteadyStateModelException {
		cont.removeMetabolites(identifyExternalMetabolites(cont, boundaryString));
		String bmid = cont.getBiomassId();
		ISteadyStateModel model = ContainerConverter.convert(cont);
		SteadyStateModelReader reader = new SteadyStateModelReader(model);
		DefaultMetabolicNetwork newmodel = reader.convertModel();
		newmodel.setBiomassReaction(bmid);
		return newmodel;
	}

	public DefaultMetabolicNetwork getMetabolicNetwork() {
		return this.model;
	}

	public void applyGeneRegulation(List<String> toBlock) {
		for (String string : toBlock) {
			int rId = model.getReactionIndex(string);
			Reaction r = model.getReaction(rId);
			r.setBounds(new ReactionConstraint(0, 0));
			model.setReaction(r, rId);
		}
		System.out.println(toBlock.size() + " reactions assumed as disabled.");
	}

	public static Set<String> identifyBlockedReactions(
			DefaultMetabolicNetwork dmodel, FluxBound[] undesiredFluxes,
			YieldConstraint[] undesiredYields) throws IloException {

		FluxVariabilityAnalysisResult fva = new CPLEXFluxVariabilityAnalysis(
				dmodel, undesiredFluxes, undesiredYields).solveFVA();
		return fva.getBlockedReactionNames();
	}
	
	public static Set<String> identifyExternalMetabolites(Container cont, String boundaryName) {
		return cont.identifyMetabolitesIdByPattern(Pattern.compile(".*"+boundaryName));
	}

	public ArrayList<String> identifyExchangeReactions(DefaultMetabolicNetwork model) {
		ArrayList<String> exchangeReactions = new ArrayList<String>();
		for (int i = 0; i < model.getNumOfReactions(); i++) {
			int metabs = 0;
			for (int j = 0; j < model.getNumOfMetabolites(); j++) {
				if (model.getStoichCoef(j, i) != 0) {
					metabs++;
				}
			}
			if (metabs == 1) {
				exchangeReactions.add(model.getReaction(i).getName());
			}
		}
		return exchangeReactions;
	}

	public static ArrayList<String> identifyExchangeReactions(Container cont, String boundaryString)
			throws InvalidSteadyStateModelException {
		// identifying external metabolites

		ISteadyStateModel fullmodel = ContainerConverter.convert(cont);
		ArrayList<String> exchangeReactionNames = new ArrayList<String>();
		for (String externalMetabolite : identifyExternalMetabolites(cont, boundaryString)) {
			Integer midx = fullmodel.getMetaboliteIndex(externalMetabolite);
			for (int i = 0; i < fullmodel.getNumberOfReactions(); i++) {
				if (Math.abs(fullmodel.getStoichiometricValue(midx, i)) > Utilities.PRECISION) {
					// if a reaction involves production or consumption of an
					// external metabolite, it is added to the exchange
					// reactions
					exchangeReactionNames.add(fullmodel.getReactionId(i));
				}
			}
		}
		return exchangeReactionNames;

	}

	public static ArrayList<String> identifySpontaneousReactions(Container cont,
			String spontaneous_handle) {
		ArrayList<String> spontaneousReactionNames = new ArrayList<String>(); // will
																		// store
																		// spontaneous
																		// reaction
																		// names
		if (spontaneous_handle != null) {
			for (ReactionCI rId : cont.getReactions().values()) {
				Set<String> genes = rId.getGenesIDs();
				if ((genes.contains(spontaneous_handle))) {
					spontaneousReactionNames.add(rId.getId());
				}
			}
		}
		return spontaneousReactionNames;
	}

	public void overrideBounds(String filename) throws IOException {
		model.overrideBounds(filename);
	}
	
	public void setSolverLogging(boolean allowSolverLogging){
		this.allowSolverLogging = allowSolverLogging;
	}
	
	public void initializesOnEnumeration(boolean initializes){
		this.initializesOnEnumeration = initializes;
	}
	public DefaultEnumerationResult enumerate(int maxParam, boolean constrained) {
		if (initializesOnEnumeration) {
			ep = initialize();
		}
		DefaultEnumerationResult mcs = null;
		
		
		
		try {
			AbstractEnumerationResult reslt = null;
			
			if (aty == AlgorithmType.INTEGRATED) {
				CPLEXIntegratedEnumerator enm = new CPLEXIntegratedEnumerator(ep);
				reslt = enm.solve(maxParam);
			} else if (aty == AlgorithmType.DEFAULT) {
				CPLEXEnumerationSolver enm = new CPLEXEnumerationSolver(ep);
				enm.setThreads(threads);
				enm.allowSolverLogging(allowSolverLogging);
				reslt = enm.solve(maxParam);
			} else if (aty == AlgorithmType.SAMPLING) {
				CPLEXIntegratedEnumerator enm = new CPLEXIntegratedEnumerator(ep);
				reslt = enm.iterativeSolve((3600*12), (3600*12), false, 1);
			}
			
			
			if (reslt.getClass() == CompressedEnumerationResults.class) {
				System.out.println("Results will be decompressed");
				mcs = ((CompressedEnumerationResults) reslt).decompressResult();
			} else {
				System.out.println("Results already uncompressed");
				mcs = (DefaultEnumerationResult) reslt;
			}
			
			if (aty == AlgorithmType.DEFAULT && constrained) {
				System.out.println("Filtering");
				System.out.println("Desired fluxes: ");
				for (FluxBound fluxBound : desired_targets) {
					System.out.println("\t"+fluxBound);
				}
				for (YieldConstraint yieldConstraint : desired_yields) {
					System.out.println("\t"+yieldConstraint);
				}
				mcs = new CPLEXEnumerationFilter(ep, mcs).calculateFilteredResults();
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		return mcs;

	}
	
	public EnumerationProblem initialize(){
		
		FluxBound[] uF = undesired_targets
				.toArray(new FluxBound[undesired_targets.size()]);
		FluxBound[] dF = desired_targets.toArray(new FluxBound[desired_targets
				.size()]);
		YieldConstraint[] uY = undesired_yields
				.toArray(new YieldConstraint[undesired_yields.size()]);
		YieldConstraint[] dY = desired_yields
				.toArray(new YieldConstraint[desired_yields.size()]);
		List<Reaction> single = new ArrayList<Reaction>();
		List<Reaction> blocked = new ArrayList<Reaction>();
		for (int i = 0; i < uF.length; i++) {
			Reaction r = uF[i].getReac();
			r.setBounds(uF[i].getBounds());
//			int idx = model.getReactionIndex(r.getName());
//			model.setReaction(r, idx);
		}

		ArrayList<FluxBound> fvaConditions = (ArrayList<FluxBound>) undesired_targets
				.clone();
		fvaConditions.addAll(fvaconditions);
		FluxBound[] fvauF = fvaConditions
				.toArray(new FluxBound[undesired_targets.size()
						+ fvaconditions.size()]);

		Reaction[] singleArray = Utilities.toReacArrayFromString(model,
				identifySingleReactions());
		Reaction[] blockedArray = null;
		try {
			blockedArray = Utilities.toReacArrayFromString(model,
					identifyBlockedReactions(model, fvauF, uY));
		} catch (IloException e) {
			e.printStackTrace();
			return null;
		}

		for (Reaction reaction : blockedArray) {
			blocked.add(reaction);
		}
		TreeSet<String> blockedOrd = new TreeSet<>();
		for (Reaction r : blocked)
			blockedOrd.add(r.getName());

		for (Reaction reaction : singleArray) {
			single.add(reaction);
		}

		System.out.println("Statistics:");
		System.out.println("\t" + blocked.size() + " blocked reactions.");
		System.out.println("\t" + single.size() + " single reactions.");
		System.out.println("\t" + exchange.size() + " exchange reactions.");
		System.out.println("\t" + spontaneous.size()
				+ " spontaneous reactions.");
		System.out.println("\t" + model.getNumberOfReversibleReactions()
				+ " reversible reactions");
		System.out.println("\t" + inflows.size() + " inflows.");
		System.out.println("\t" + outflows.size() + " outflows.");
		System.out.println("\t" + uF.length + " undesired flux vectors.");
		for (int i = 0; i < uF.length; i++) {
			System.out.println("\t\t" + uF[i]);
		}
		System.out.println("\t" + uY.length + " undesired yield constraints.");
		for (int i = 0; i < uY.length; i++) {
			System.out.println("\t\t" + uY[i]);
		}

		NullspaceNetworkCompressor red = new NullspaceNetworkCompressor(this.model);
		CompressedMetabolicNetwork comp = null;
		try {
			comp = red.reduceModel(blocked, single);
		} catch (IOException e) {
			e.printStackTrace();
			return null;
		}
		
		curcomp = comp;
		
		for (String reaction : excluded) {
			single.add(this.model.getReaction(reaction));
		}

		EnumerationProblem epLocal = new EnumerationProblem(comp, uF, dF, uY, dY,
				single.toArray(new Reaction[single.size()]));
		if (excludedSols != null) {
			epLocal.setExcludedSubsets(excludedSols
					.toArray(new Solution[excludedSols.size()]));
		}
		
		epLocal.setForcedReactions(Utilities.toReacArrayFromString(model, forced));
		this.ep = epLocal;
		return epLocal;
		
	}

	public DefaultEnumerationResult filter(String filename) {
		FluxBound[] uF = undesired_targets
				.toArray(new FluxBound[undesired_targets.size()]);
		FluxBound[] dF = desired_targets.toArray(new FluxBound[desired_targets
				.size()]);
		YieldConstraint[] uY = undesired_yields
				.toArray(new YieldConstraint[undesired_yields.size()]);
		YieldConstraint[] dY = desired_yields
				.toArray(new YieldConstraint[desired_yields.size()]);
		EnumerationProblem ep = new EnumerationProblem(model, uF, dF, uY, dY,
				Utilities.toReacArrayFromString(model,
						identifySingleReactions()));
		DefaultEnumerationResult unconstrained;
		try {
			unconstrained = DefaultEnumerationResult
					.fromFile(filename, ep);
			
			return new CPLEXEnumerationFilter(ep,
					unconstrained).calculateFilteredResults();
		} catch (IOException e) {
			e.printStackTrace();
			return null;
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}

	}
	
	public DefaultEnumerationResult filter(DefaultEnumerationResult res) {
		FluxBound[] uF = undesired_targets
				.toArray(new FluxBound[undesired_targets.size()]);
		FluxBound[] dF = desired_targets.toArray(new FluxBound[desired_targets
				.size()]);
		YieldConstraint[] uY = undesired_yields
				.toArray(new YieldConstraint[undesired_yields.size()]);
		YieldConstraint[] dY = desired_yields
				.toArray(new YieldConstraint[desired_yields.size()]);
		EnumerationProblem ep = new EnumerationProblem(model, uF, dF, uY, dY,
				Utilities.toReacArrayFromString(model,
						identifySingleReactions()));
		try {
			return new CPLEXEnumerationFilter(ep,
					res).calculateFilteredResults();
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}

	}

	public boolean isConstrainedProblem() {
		return (desired_targets.size() + desired_yields.size()) > 0;
	}

	public void addExcludedSolution(Solution solution) {
		this.excludedSols.add(solution);
	}
	
	public void addExcludedSolution(List<List<String>> solution) {
		for (List<String> list : solution) {
			int[] indexes = new int[list.size()];
			for (int i = 0; i < list.size(); i++) {
				indexes[i] = this.model.containsReaction(list.get(i));
			}
			Solution s = new Solution(indexes, model);
			this.excludedSols.add(s);
		}
	}

	public boolean isSeedingProblem() {
		return this.seeds.size() > 0;
	}
	
	public int getNumberOfSeeds() {
		return this.seeds.size();
	}
	
	public List<String[]> getSeedList() {
		List<String[]> res = new ArrayList<String[]>();
		for (Solution string : this.seeds) {
		}
		return res;
	}

	public void addForcedReaction(String string) {
		forced.add(string);
	}

	public int getNumberOfForcedReactions() {
		return forced.size();
	}

	public void addNonTargets(Collection<String> readLines) {
		for (String string : readLines) {
			if (this.model.containsReaction(string) > -1) {
				this.ignoreReaction(string);
			}
		}
	}

	public void setThreads(int threads) {
		this.threads = threads;
		
	}
	
}
