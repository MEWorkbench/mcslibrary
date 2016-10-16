/*******************************************************************************
 * Copyright 2016
 * CEB Centre of Biological Engineering
 * University of Minho
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This code is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this code. If not, see http://www.gnu.org/licenses/
 *
 * Created inside the BIOSYSTEMS Research Group
 * (http://www.ceb.uminho.pt/biosystems)
 *******************************************************************************/
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
