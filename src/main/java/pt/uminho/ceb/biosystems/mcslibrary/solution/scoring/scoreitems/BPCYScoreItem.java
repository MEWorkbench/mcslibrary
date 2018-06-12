package pt.uminho.ceb.biosystems.mcslibrary.solution.scoring.scoreitems;

import java.util.List;

import ilog.concert.IloException;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.SimulationResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba.CPLEXParsimoniousFluxBalanceAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;

public class BPCYScoreItem implements IScoreItem{

	private DefaultMetabolicNetwork dmn;
	private CPLEXParsimoniousFluxBalanceAnalysis pfba;
	private Reaction biomass;
	private Reaction product;
	private Reaction substrate;
	private FluxBound[] bounds;

	public BPCYScoreItem(DefaultMetabolicNetwork dmn, CPLEXParsimoniousFluxBalanceAnalysis pfba, FluxBound[] bounds, Reaction biomass, Reaction product, Reaction substrate) {
		this.dmn = dmn;
		this.pfba = pfba;
		this.biomass = biomass;
		this.product = product;
		this.substrate = substrate;
		this.bounds = bounds;
	}
	
	@Override
	public double evaluateReactionKnockout(List<String> rk) {
		double r = Double.NaN;
		try {
			SimulationResult optimal = pfba.solveKnockoutFBA(bounds, Utilities.toReacArrayFromString(dmn, rk), biomass, "max");
			double b = optimal.getFluxValue(biomass);
			double p = optimal.getFluxValue(product);
			double s = optimal.getFluxValue(substrate);
			r = Math.abs((b*p)/s);
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return r;
	}

	@Override
	public String getItemName() {
		// TODO Auto-generated method stub
		return "("+biomass.getName()+"*"+product.getName()+")/"+substrate.getName();
	}

}
