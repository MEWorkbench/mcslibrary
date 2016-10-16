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
