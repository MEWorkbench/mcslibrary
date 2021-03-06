package pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints;

import java.io.Serializable;

public class ReactionConstraint implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 667095566372788119L;
	private double lower;

	public double getLower() {
		return lower;
	}

	public void setLower(double lower) {
		this.lower = lower;
	}

	public double getUpper() {
		return upper;
	}

	public void setUpper(double upper) {
		this.upper = upper;
	}

	private double upper;

	public ReactionConstraint(double lowerbound, double upperbound) {
		this.lower = lowerbound;
		this.upper = upperbound;
	}
	
	public ReactionConstraint getSplitForward(){
		double lb = (lower > 0) ? lower : 0;
		double ub = (upper > 0) ? upper : 0;
		return new ReactionConstraint(lb, ub);
	}
	
	public ReactionConstraint getSplitReverse(){
		double lb = (lower < 0) ? lower : 0;
		double ub = (upper < 0) ? upper : 0;
		return new ReactionConstraint(Math.abs(ub), Math.abs(lb));
	}

}
