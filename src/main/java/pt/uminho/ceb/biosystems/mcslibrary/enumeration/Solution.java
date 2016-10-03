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
