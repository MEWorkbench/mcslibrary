package pt.uminho.ceb.biosystems.mcslibrary.solution.scoring.scoreitems;

import java.util.List;

import ilog.concert.IloException;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.SimulationResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba.CPLEXParsimoniousFluxBalanceAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;

public class PFBAFluxValueItem implements IScoreItem{
	
	private DefaultMetabolicNetwork dmn;
	private CPLEXParsimoniousFluxBalanceAnalysis pfba;
	private Reaction flux;
	private Reaction objective;
	private String sense;
	private FluxBound[] bounds;

	public PFBAFluxValueItem(DefaultMetabolicNetwork dmn, CPLEXParsimoniousFluxBalanceAnalysis pfba, FluxBound[] bounds ,Reaction flux, Reaction objective, String sense) {
		this.dmn = dmn;
		this.pfba = pfba;
		this.flux = flux;
		this.objective = objective;
		this.sense = sense;
		this.bounds = bounds;
	}
	
	@Override
	public double evaluateReactionKnockout(List<String> rk) {
		double r = Double.NaN;
		try {
			Reaction[] ko = Utilities.toReacArrayFromString(dmn, rk);
			SimulationResult kopfba = this.pfba.solveKnockoutFBA(bounds, ko, objective, sense);
			r = kopfba.getFluxValue(flux);
		} catch (IloException e) {
			e.printStackTrace();
		}
		return r;
	}

	@Override
	public String getItemName() {
		return "PFBA("+flux.getName()+")";
	}

}
