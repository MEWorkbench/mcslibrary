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
package pt.uminho.ceb.biosystems.mcslibrary.metabolic;

import java.io.Serializable;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.ReactionConstraint;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;

/**
 * Class representing a single reaction.
 * @author V�tor
 *
 */
public class Reaction implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = -8934053258671788377L;
	ReactionConstraint bounds;
	String id;
	private boolean isExchange;
	
	/**
	 *
	 * @param id - Reaction id/name.
	 * @param limits - {@link ReactionConstraint} indicating the limits imposed by the model
	 * @param isExchange - {@link Boolean} indicating whether this reaction is a drain reaction. Null unless accessed from a {@link DefaultMetabolicNetwork}.
	 */
	public Reaction(String id, ReactionConstraint limits) {
		this.id = id;
		this.bounds = limits;
	}

	public boolean isReversible(){
		boolean rev = true;
//		if ((Math.abs(this.bounds.getLower()) < Utilities.PRECISION) || (this.bounds.getLower() > 0)) {
//			rev = false;
//		}
		if (this.bounds.getLower() >= 0) {
			rev = false;
		}
		return rev;
	}

	public ReactionConstraint getBounds() {
		return bounds;
	}

	public void setBounds(ReactionConstraint bounds) {
		this.bounds = bounds;
	}

	public String getName() {
		return id;
	}

	public boolean isExchange() {
		return isExchange;
	}
	@Override
	public String toString() {
		String str = getName()+"\t \t LB: "+getBounds().getLower()+"\t UB:"+getBounds().getUpper();
		return str;

	}
	public void setExchange(boolean isExchange) {
		this.isExchange = isExchange;
	}
}
