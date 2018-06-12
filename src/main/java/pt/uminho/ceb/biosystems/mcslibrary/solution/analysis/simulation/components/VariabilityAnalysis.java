package pt.uminho.ceb.biosystems.mcslibrary.solution.analysis.simulation.components;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ilog.concert.IloException;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.implementation.CPLEXFlexibleFVA;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.ReactionConstraint;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba.CPLEXFluxBalanceAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fva.FluxVariabilityAnalysisResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Pair;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;

public class VariabilityAnalysis {
	private DefaultMetabolicNetwork dmn;
	private EnvelopeProperties[] envProp;
	private CPLEXFlexibleFVA fva;
	private CPLEXFluxBalanceAnalysis fba;

	public VariabilityAnalysis(DefaultMetabolicNetwork dmn, EnvelopeProperties[] envProp) throws IloException {
		this.dmn = dmn;
		this.envProp = envProp;
		this.fva = new CPLEXFlexibleFVA(dmn);
		this.fba = new CPLEXFluxBalanceAnalysis(dmn);
	}
	
	public Map<FluxBound[], VariabilityAnalysisResult> processSolutionMultipleConditions(List<String> solution, FluxBound[][] env) throws IloException{
		Map<FluxBound[], VariabilityAnalysisResult> res = new HashMap<FluxBound[], VariabilityAnalysisResult>();
		for (int i = 0; i < env.length; i++)
			res.put(env[i],processSolution(solution, env[i]));

		return res;
	}
	
	// class for FVA related analysis
	// - envelopes
	// - blocked / essential (limits)
	public VariabilityAnalysisResult processSolution(List<String> solution, FluxBound[] env) throws IloException{
		Reaction[] knockouts = Utilities.toReacArrayFromString(dmn, solution);
		// limits
		FluxVariabilityAnalysisResult limits = fva.solveFVA(env, knockouts);
		EnvelopeResult[] envelopes = new EnvelopeResult[envProp.length];
		for (int i = 0; i < envProp.length; i++) {
			EnvelopeProperties prop = envProp[i];
			int pivotId = dmn.getReactionIndex(prop.getPivot());
			Pair<double[], double[]> min = fba.solvePivotFBA(
					env, 
					knockouts, 
					dmn.getReaction(prop.getObjective()), 
					"min", 
					dmn.getReaction(prop.getPivot()), 
					new ReactionConstraint(limits.minFlux(pivotId), limits.maxFlux(pivotId)), 
					prop.getSteps());
			Pair<double[], double[]> max = fba.solvePivotFBA(
					env, 
					knockouts, 
					dmn.getReaction(prop.getObjective()), 
					"max", 
					dmn.getReaction(prop.getPivot()), 
					new ReactionConstraint(limits.minFlux(pivotId), limits.maxFlux(pivotId)), 
					prop.getSteps());
			double[] scale = min.getA();
			double[] pivotMin = min.getB();
			double[] pivotMax = max.getB();
			envelopes[i] = new EnvelopeResult(pivotMin, pivotMax, scale, prop);
		}

		return new VariabilityAnalysisResult(dmn, solution, limits, envelopes, env);
	}
	
	public EnvelopeProperties[] getEnvelopeProperties(){
		return this.envProp;
	}
	
}
