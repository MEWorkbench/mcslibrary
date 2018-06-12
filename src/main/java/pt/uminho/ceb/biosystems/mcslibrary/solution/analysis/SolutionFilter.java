package pt.uminho.ceb.biosystems.mcslibrary.solution.analysis;

import java.util.ArrayList;
import java.util.List;

import ilog.concert.IloException;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.SimulationResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba.CPLEXFluxBalanceAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;

public class SolutionFilter {
	private DefaultMetabolicNetwork dmn;
	private CPLEXFluxBalanceAnalysis fba;

	public SolutionFilter(DefaultMetabolicNetwork dmn) throws IloException {
		this.dmn = dmn;
		this.fba = new CPLEXFluxBalanceAnalysis(dmn);
	}
	
	public List<List<String>> filterBounds(FluxBound[] criteria, List<List<String>> solutions) throws IloException{
		List<List<String>> res = new ArrayList<>();
		for (int i = 0; i < solutions.size(); i++) {
			SimulationResult solution = fba.solveReactionKnockoutFBA(criteria, Utilities.toReacArrayFromString(dmn, solutions.get(i)), dmn.getReaction(dmn.getBiomassReactionId()), "max");
//			System.out.println(fba.robustnessPoint(criteria, Utilities.toReacArrayFromString(dmn, solutions.get(i)), dmn.getReaction("R_EX_succ_e_"), dmn.getReaction(dmn.getBiomassReactionId()), 1e-9, 0.01));
			if (!solution.getStatus().equals("Infeasible")) {
				res.add(solutions.get(i));
			}
		}
		return res;
	}
	
	public List<List<String>> filterRobust(FluxBound[] original, RobustnessCriteria[] criteria, List<List<String>> solutions) throws IloException{
		List<List<String>> res = new ArrayList<>();
		for (int i = 0; i < solutions.size(); i++) {
			boolean valid = true;
			for (int j = 0; j < criteria.length; j++) {
				valid &= criteria[j].evaluate(dmn, fba, solutions.get(i), original);
			}
//			System.out.println(valid);
			if (valid) {
				res.add(solutions.get(i));
			}
		}
		return res;
	}
	
	public List<List<String>> doubleFilter(FluxBound[] original, RobustnessCriteria[] criteria, List<List<String>> solutions) throws IloException{
		return filterRobust(original, criteria, filterBounds(original, solutions));
	}
}
