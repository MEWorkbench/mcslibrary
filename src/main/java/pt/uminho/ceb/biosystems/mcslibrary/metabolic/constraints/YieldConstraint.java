package pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints;

import java.io.Serializable;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;




public class YieldConstraint implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 3863649109892138908L;
	private Reaction uptakeReaction;
	private Reaction productReaction;
	private double ratio;
	private double deviation;
	private boolean isLower = true;
/**
 *
 * @param uptake - {@link Reaction} that represents uptake.
 * @param product - {@link Reaction} that represents product synthesis
 * @param ratio - Maximum value allowed for Uptake/Product. If you wish to impose a lower bound instead, simply invert Product and Uptake and use 1/ratio.
 */
	public YieldConstraint(Reaction uptake, Reaction product, double ratio) {
		this.uptakeReaction = uptake;
		this.productReaction = product;
		this.ratio = ratio;
		this.deviation = 0;
	}

	public Reaction getUptakeReaction() {
		return uptakeReaction;
	}

	public Reaction getProductReaction() {
		return productReaction;
	}

	public double getRatio() {
		return ratio;
	}

	public String toString() {
		String str = uptakeReaction.getName()+" / "+productReaction.getName()+" <= "+ratio;
		return str;
	}
	
	public void setAsLower(boolean isLower) {
		this.isLower = true;
	}
	
	
	public boolean isLower() {
		return this.isLower;
	}

	public double getDeviation() {
		return deviation;
	}

	public void setDeviation(double deviation) {
		this.deviation = deviation;
	}

}
