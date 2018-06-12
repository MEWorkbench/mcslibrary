package pt.uminho.ceb.biosystems.mcslibrary.solution.scoring.scoreitems;

import java.util.List;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba.CPLEXFluxBalanceAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;

public class RobustnessScoreItem implements IScoreItem{

	private DefaultMetabolicNetwork dmn;
	private CPLEXFluxBalanceAnalysis fba;
	private Reaction biomass;
	private Reaction product;
	private FluxBound[] bounds;
	private double fraction;
	private boolean fixed;

	public RobustnessScoreItem(DefaultMetabolicNetwork dmn, CPLEXFluxBalanceAnalysis fba, FluxBound[] bounds, Reaction biomass, Reaction product, double fraction, boolean fixed) {
		this.dmn = dmn;
		this.fba = fba;
		this.biomass = biomass;
		this.product = product;
		this.bounds = bounds;
		this.fraction = fraction;
		this.fixed = fixed;
	}
	@Override
	public double evaluateReactionKnockout(List<String> rk) {
		double r = Double.NaN;
		try {
			r = fba.solveReactionKnockoutFVA(bounds, Utilities.toReacArrayFromString(dmn, rk), product, "min", fraction, biomass, fixed).getFluxValue(product);
		} catch (Exception e) {
			// TODO: handle exception
		}
		return r;
	}
	@Override
	public String getItemName() {
		return "FVAmin("+product.getName()+")@ "+((int) (fraction*100))+"% biomass";
	}

}
