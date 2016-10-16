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
import java.util.ArrayList;
import java.util.List;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.ReactionConstraint;
/**
 * Class that represents a reaction cluster in a compressed network.
 * @author V�tor
 *
 */ 
public class ReactionGroup implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 2390682943768219370L;
	private List<Reaction> reactions;
	private ReactionConstraint bounds;
	/**
	 *
	 * @param lowerbound - Minimum value for this flux in the compressed network.
	 * @param upperbound - Maximum value for this flux in the compressed network.
	 */
	public ReactionGroup(double lowerbound, double upperbound) {
		reactions = new ArrayList<Reaction>();
		this.bounds = new ReactionConstraint(lowerbound, upperbound);
	}

	public ReactionConstraint getBounds() {
		return bounds;
	}

	public void setBounds(ReactionConstraint bounds) {
		this.bounds = bounds;
	}

	public void addReaction(Reaction reac){
		reactions.add(reac);
	}

	public int size(){
		return reactions.size();
	}

	public int containsReaction(String name){
		int res = -1;
		for (int i = 0; i < reactions.size(); i++) {
			if (reactions.get(i).getName() == name) {
				res = i;
				break;
			}
		}
		return res;
	}

	@Override
	public String toString() {
		String str = "Reaction group -- LB:"+this.getBounds().getLower()+" UB:"+this.getBounds().getUpper();
		for (int i = 0; i < this.size(); i++) {
			str = str+"\n \t"+"-"+this.reactions.get(i).toString();
		}
		return str;
	}

	public int containsReaction(Reaction reac){
		int res = -1;
		for (int i = 0; i < size(); i++) {
			if (reactions.get(i).equals(reac)){
				res = i;
				break;
			}
		}
		return res;

	}

	public Reaction getReaction(int index){
		return this.reactions.get(index);
	}

	public boolean isReversible(){
		boolean res = true;
		if (this.bounds.getLower() >= 0) {
			res = false;
		}
		return res;
	}

	public boolean hasAnyExchange() {
		boolean res = false;
		for (int i = 0; i < this.size(); i++) {
			res |= getReaction(i).isExchange();
		}
		return res;
	}

}
