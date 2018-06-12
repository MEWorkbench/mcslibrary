package pt.uminho.ceb.biosystems.mcslibrary.solution.scoring.scoreitems;

import java.util.List;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Metabolite;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.SimulationResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba.CPLEXParsimoniousFluxBalanceAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;

public class CarbonYieldScoreItem implements IScoreItem{
	private DefaultMetabolicNetwork dmn;
	private CPLEXParsimoniousFluxBalanceAnalysis pfba;
	private FluxBound[] bounds;
	private Reaction product;
	private Reaction substrate;
	private int pCarbons;
	private int sCarbons;
	private Reaction biomass;

	public CarbonYieldScoreItem(DefaultMetabolicNetwork dmn, CPLEXParsimoniousFluxBalanceAnalysis pfba, FluxBound[] bounds, Reaction biomass, Reaction product, Reaction substrate) {
		this.dmn = dmn;
		this.pfba = pfba;
		this.bounds = bounds;
		this.biomass = biomass;
		
		this.product = product;
		List<Metabolite> productMetabs = dmn.getMetabolitesInvolvedInReaction(dmn.getReactionIndex(product.getName()));
//		System.out.println(productMetabs);
		pCarbons = productMetabs.get(0).getCarbonAtoms();
		
		this.substrate = substrate;
		List<Metabolite> substrateMetabs = dmn.getMetabolitesInvolvedInReaction(dmn.getReactionIndex(substrate.getName()));
		sCarbons = substrateMetabs.get(0).getCarbonAtoms();
	}

	@Override
	public double evaluateReactionKnockout(List<String> rk) {
		double r = Double.NaN;
		try {
			SimulationResult optimal = pfba.solveKnockoutFBA(bounds, Utilities.toReacArrayFromString(dmn, rk), biomass, "max");
			double p = optimal.getFluxValue(product)*pCarbons;
			double s = optimal.getFluxValue(substrate)*sCarbons;
			r = p/s;
		} catch (Exception e) {
			// TODO: handle exception
		}
		return Math.abs(r);
	}

	@Override
	public String getItemName() {
		// TODO Auto-generated method stub
		return "CYield("+product.getName()+"/"+substrate.getName()+")";
	}
	
}
