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
package pt.uminho.ceb.biosystems.mcslibrary.enumeration;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;

public class Solution {
	private DefaultMetabolicNetwork model;
	protected int[] index;
	
	public Solution(int[] indexes, DefaultMetabolicNetwork model) {
		this.model = model;
		this.index = indexes;
	}
	
	public int getIndex(int pos) {
		return this.index[pos];
	}
	public int getSize() {
		return this.index.length;
	}
	
	public Reaction getReactionFromIndex(int index) {
		return model.getReaction(this.index[index]);
	}
	
	public String getReactionIdFromIndex(int index) {
		return model.getReaction(this.index[index]).getName();
	}

	public int[] getIndexes() {
		return index;
	}
	
	public DefaultMetabolicNetwork getMetabolicNetwork() {
		return this.model;
	}
}
