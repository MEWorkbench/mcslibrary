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
package pt.uminho.ceb.biosystems.mew.core.strainoptimization.strainoptimizationalgorithms.pathwayanalysis;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import pt.uminho.ceb.biosystems.mcslibrary.enumeration.EnumerationProblem;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.implementation.DefaultEnumerationResult;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.MCSPipeline;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.Container;
import pt.uminho.ceb.biosystems.mew.core.model.components.EnvironmentalConditions;
import pt.uminho.ceb.biosystems.mew.core.model.components.ReactionConstraint;
import pt.uminho.ceb.biosystems.mew.core.model.steadystatemodel.ISteadyStateModel;
import pt.uminho.ceb.biosystems.mew.core.simulation.components.GeneticConditions;
import pt.uminho.ceb.biosystems.mew.core.simulation.components.ReactionChangesList;
import pt.uminho.ceb.biosystems.mew.core.simulation.components.SimulationProperties;
import pt.uminho.ceb.biosystems.mew.core.simulation.components.SimulationSteadyStateControlCenter;
import pt.uminho.ceb.biosystems.mew.core.simulation.components.SteadyStateSimulationResult;
import pt.uminho.ceb.biosystems.mew.core.strainoptimization.algorithm.AbstractStrainOptimizationAlgorithm;
import pt.uminho.ceb.biosystems.mew.core.strainoptimization.optimizationresult.solution.RKSolution;
import pt.uminho.ceb.biosystems.mew.core.strainoptimization.optimizationresult.solution.SolutionFactory;
import pt.uminho.ceb.biosystems.mew.core.strainoptimization.optimizationresult.solutionset.RKSolutionSet;
import pt.uminho.ceb.biosystems.mew.solvers.builders.CPLEX3SolverBuilder;
import pt.uminho.ceb.biosystems.mew.utilities.datastructures.pair.Pair;

public class MinimalCutSetEnumerationAlgorithm extends AbstractStrainOptimizationAlgorithm<McslibraryGenericConfiguration> {

	/**
	 * 
	 */
	private static final long serialVersionUID = -8811595896473127528L;

//	@Override
//	public void putAllProperties(McslibraryGenericConfiguration configuration) {
//		// TODO Auto-generated method stub
//		
//	}

	public MinimalCutSetEnumerationAlgorithm() {
		setAlgorithmConfiguration(new McslibraryGenericConfiguration());
	}
	@Override
	public RKSolutionSet<McslibraryGenericConfiguration> execute(
			McslibraryGenericConfiguration configuration) throws Exception {
		Container cont = configuration.getContainer();
		String spontaneous_handle = configuration.getSpontaneousHandle();
		MCSPipeline p = new MCSPipeline(cont, spontaneous_handle);
		p.initializesOnEnumeration(false);
		
		EnvironmentalConditions env = configuration.getEnvironmentalConditions();
		EnvironmentalConditions uFlux = configuration.getUndesiredFluxes();
		EnvironmentalConditions dFlux = configuration.getDesiredFluxes();
		Map<Pair<String, String>, Pair<Boolean, Double>> uYield = configuration.getUndesiredYields();
		Map<Pair<String, String>, Pair<Boolean, Double>> dYield = configuration.getDesiredYields();
		if (uFlux != null) {
			for (String reaction : uFlux.keySet()) {
				ReactionConstraint rcons = uFlux.getReactionConstraint(reaction);
				p.addFluxBound(reaction, rcons.getLowerLimit(), rcons.getUpperLimit(), true);
			}
		}
		if (dFlux != null) {
			for (String reaction : dFlux.keySet()) {
				ReactionConstraint rcons = dFlux.getReactionConstraint(reaction);
				p.addFluxBound(reaction, rcons.getLowerLimit(), rcons.getUpperLimit(), false);
			}
		}
		if (env != null) {
			for (String envCond : env.keySet()) {
				ReactionConstraint rcons = env.getReactionConstraint(envCond);
				p.addFluxBound(envCond, rcons.getLowerLimit(), rcons.getUpperLimit(), false);
				p.addFluxBound(envCond, rcons.getLowerLimit(), rcons.getUpperLimit(), true);
			}
		}
		
		if (uYield != null) {
			for (Pair<String, String> yieldPair : uYield.keySet()) {
				String rU = yieldPair.getA();
				String rP = yieldPair.getB();
				double ratio = uYield.get(yieldPair).getB();
				boolean isUpper = uYield.get(yieldPair).getA();
				if (isUpper) {
					p.addUpperYieldConstraint(rU, rP, ratio, true);
				} else {
					p.addLowerYieldConstraint(rU, rP, ratio, true);
				}
			}
		}
		
		if (dYield != null) {
			for (Pair<String, String> yieldPair : dYield.keySet()) {
				String rU = yieldPair.getA();
				String rP = yieldPair.getB();
				double ratio = uYield.get(yieldPair).getB();
				boolean isUpper = uYield.get(yieldPair).getA();
				if (isUpper) {
					p.addUpperYieldConstraint(rU, rP, ratio, false);
				} else {
					p.addLowerYieldConstraint(rU, rP, ratio, false);
				}
			}
			
		}

		
		List<String> excTargets = configuration.getExcludedReactions();
		List<List<String>> excSolutions = configuration.getExcludedSolutions();
		
		if (excSolutions != null) {
			p.addExcludedSolution(excSolutions);
		}
		
		
		if (excTargets != null) {
			for (String reaction : excTargets) {
				p.ignoreReaction(reaction);
			}
		}
		
		
		if (env != null) {
			for (String reaction : env.keySet()) {
				ReactionConstraint rcons = env.getReactionConstraint(reaction);
				p.addFVABound(reaction, rcons.getLowerLimit(), rcons.getUpperLimit());
			}
		}
		
		p.correctCapacities();
		p.correctInOutflows();
//		try {
//			p.addSingleReaction("R_ATPM");
//		} catch (Exception e) {
//			// TODO: handle exception
//		}
		int maxSolSize = configuration.getMaximumSolutionSize();
		boolean cd = ( (dFlux == null) ? false : ( (dFlux.size() > 0) ? true : false));
		boolean cy = ( (dYield == null) ? false : ( (dYield.size() > 0) ? true : false));

		EnumerationProblem prob = p.initialize();
		
		System.out.println("Yield constraints");
		for (int i = 0; i < prob.getUndesiredYieldConstraints().length; i++) {
			System.out.println(prob.getUndesiredYieldConstraints()[i]);
		}
		
		System.out.println("Flux constraints");
		for (int i = 0; i < prob.getUndesiredFluxes().length; i++) {
			System.out.println(prob.getUndesiredFluxes()[i]);
		}
		DefaultEnumerationResult result = p.enumerate(maxSolSize, cd || cy);
		
		return buildResults(result, configuration);
	}

	private RKSolutionSet<McslibraryGenericConfiguration> buildResults(DefaultEnumerationResult res, McslibraryGenericConfiguration conf) throws InstantiationException, IllegalAccessException, IllegalArgumentException, InvocationTargetException, NoSuchMethodException, SecurityException {
		RKSolutionSet<McslibraryGenericConfiguration> sol = new RKSolutionSet<McslibraryGenericConfiguration>(conf);
		SolutionFactory<RKSolution> sf = new SolutionFactory<>();
		
		System.out.println("BASE CFG");
		System.out.println(sol.getBaseConfiguration());
		ArrayList<String[]> strRes = res.toStringArrays();
		
		for (String[] strings : strRes) {
			ReactionChangesList reactionList = new ReactionChangesList(strings);
			GeneticConditions solutionGeneticConditions = new GeneticConditions(reactionList);
			SteadyStateSimulationResult minimumFVA = generateFVAMin(solutionGeneticConditions, conf.getEnvironmentalConditions(), conf.getSteadyStateModel(), conf.getTargetId(), SimulationProperties.FBA);
			SteadyStateSimulationResult pFBA = generateBiomass(solutionGeneticConditions, conf.getEnvironmentalConditions(), conf.getSteadyStateModel(), conf.getBiomassId(), SimulationProperties.PFBA);

			List<Double> fit = new ArrayList<Double>();
			Map<String, SteadyStateSimulationResult> sssrMap = new HashMap<>();
			fit.add(minimumFVA.getOFvalue());
			sssrMap.put(SimulationProperties.FBA, minimumFVA);
			fit.add(pFBA.getFluxValues().getValue(conf.getBiomassId()));
			sssrMap.put(SimulationProperties.PFBA, pFBA);
			fit.add(getBiomassProductCoupledYield(pFBA.getFluxValues().getValue(conf.getBiomassId()), pFBA.getFluxValues().getValue(conf.getTargetId()), pFBA.getFluxValues().getValue(conf.getSubstrateId())));
			sssrMap.put(SimulationProperties.PFBA, pFBA);
			
			sol.addSolution(sf.getInstance(RKSolution.class, solutionGeneticConditions, sssrMap, fit));
//			sol.addSolution(new RKSolution(solutionGeneticConditions, sssrMap, fit));
		}
		return sol;
	}
	
	private SteadyStateSimulationResult generateFVAMin(GeneticConditions genes, EnvironmentalConditions envCond, ISteadyStateModel model, String product, String methodType){
		SimulationSteadyStateControlCenter ssscc = new SimulationSteadyStateControlCenter(envCond, genes, model, methodType);
		ssscc.setMaximization(false);
		ssscc.setSolver(CPLEX3SolverBuilder.ID);
		ssscc.setFBAObjSingleFlux(product, 1.0);
		SteadyStateSimulationResult simresult = null;
		try {
			simresult = ssscc.simulate();
		} catch (Exception e) {
			
		}
		ssscc.forceSolverCleanup();
		return simresult;
		
	}
	
	private SteadyStateSimulationResult generateBiomass(GeneticConditions genes, EnvironmentalConditions envCond, ISteadyStateModel model, String biomass, String methodType){
		SimulationSteadyStateControlCenter ssscc = new SimulationSteadyStateControlCenter(envCond, genes, model, methodType);
		ssscc.setMaximization(true);
		ssscc.setSolver(CPLEX3SolverBuilder.ID);
		ssscc.setFBAObjSingleFlux(biomass, 1.0);
		SteadyStateSimulationResult simresult = null;
		try {
			simresult = ssscc.simulate();
		} catch (Exception e) {
			
		}
		ssscc.forceSolverCleanup();
		return simresult;
		
	}
	
	private static double getBiomassProductCoupledYield(double b, double p, double s){
		return Math.abs((double)(b*p)/s);
	}
}
