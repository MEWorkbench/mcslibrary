package pt.uminho.ceb.biosystems.mcslibrary.solution.analysis;

import java.util.LinkedHashMap;
import java.util.List;

import ilog.concert.IloException;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.SimulationResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba.CPLEXParsimoniousFluxBalanceAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;

public class SimulationMap extends LinkedHashMap<FluxBound[], SimulationResult>{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private DefaultMetabolicNetwork mn;
	private FluxBound[][] flx;
	private CPLEXParsimoniousFluxBalanceAnalysis fba;

	public SimulationMap(DefaultMetabolicNetwork mn, FluxBound[][] flx, List<String> solution) throws IloException {
		this.mn = mn;
		this.flx = flx;
		this.fba = new CPLEXParsimoniousFluxBalanceAnalysis(mn, 0.99999); 
		for (int i = 0; i < flx.length; i++) {
			put(flx[i], fba.solveKnockoutFBA(flx[i], Utilities.toReacArrayFromString(mn, solution) ,mn.getReaction(mn.getBiomassReactionId()), "max"));
		}
	}

	public FluxBound[][] getFluxBounds() {
		return flx;
	}

	public DefaultMetabolicNetwork getMetabolicNetwork() {
		return mn;
	}
}
