package pt.uminho.ceb.biosystems.mcslibrary.solution.analysis;

import java.util.List;
import java.util.Map;

import ilog.concert.IloException;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.SimulationResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.solution.analysis.simulation.components.VariabilityAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.solution.analysis.simulation.components.VariabilityAnalysisResult;

public class SolutionAnalysisResult {
	private List<String> solution;
	private DefaultMetabolicNetwork dmn;
	private Map<FluxBound[], VariabilityAnalysisResult> variability;
	private SimulationMap simulations;
	private FluxBound[][] envConds;

	public SolutionAnalysisResult(List<String> solution, DefaultMetabolicNetwork dmn, VariabilityAnalysis var, FluxBound[][] envConds) throws IloException {
		this.solution = solution;
		this.dmn = dmn;
		if (var != null) {
			this.variability = var.processSolutionMultipleConditions(solution, envConds);
		}
		this.simulations = new SimulationMap(dmn, envConds, solution);
		this.envConds = envConds;
		
	}
	
	public VariabilityAnalysisResult getVariabilityAnalysis(FluxBound[] fluxBound){
		return variability.get(fluxBound);
	}
	
	public SimulationResult getSimulationResult(FluxBound[] fluxBound){
		return simulations.get(fluxBound);
	}
	
	public DefaultMetabolicNetwork getMetabolicNetwork(){
		return dmn;
	}
	
	public List<String> getSolution(){
		return solution;
	}
	
	public FluxBound[][] getEnvironmentalConditions(){
		return this.envConds;
	}
	
	
}
