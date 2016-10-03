package pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints;

import java.io.Serializable;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;

/**
 * This class represents a generic flux bound, desired or undesired.
 * @author V�tor
 *
 */
public class FluxBound implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = -4136563562346082910L;
	private Reaction reac;
	private ReactionConstraint bounds;
	private boolean overrideLB;
	private boolean overrideUB;
/**
 *
 * @param reac - Instance of {@link Reaction} to which the constraint will be applied.
 * @param lowerbound - Minimum value for the flux.
 * @param upperbound - Maximum value for the flux.
 */
	public FluxBound(Reaction reac, double lowerbound, double upperbound) {
		this.reac = reac;
		this.bounds = new ReactionConstraint(lowerbound, upperbound);
		this.overrideLB = false;
		this.overrideUB = false;
	}
	
	public void setOverride(boolean lower, boolean upper){
		this.overrideLB = lower;
		this.overrideUB = upper;
	}
	
	public boolean overridesLower(){
		return this.overrideLB;
	}
	
	public boolean overridesUpper(){
		return this.overrideUB;
	}
	
	
	public ReactionConstraint getBounds() {
		return bounds;
	}

	public Reaction getReac() {
		return reac;
	}

	public String toString() {
		String str = bounds.getLower()+" <= "+reac.getName()+" <= "+bounds.getUpper();
		return str;
	}

}
