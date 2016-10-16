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

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;

public abstract class AbstractMetabolicNetwork implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = -4140234916380002460L;
	private Metabolite[] metabolites;
	private double[][] matrix;
	protected int biomass;


	public int getMetaboliteIndex(String metaName) {
		int index = -1;
		for (int i = 0; i < this.metabolites.length; i++) {
			if (this.metabolites[i].getName().equals(metaName)) {
				index = i;
				break;
			}
		}
		return index;
	}

	public abstract FluxBound[] getExtraConstraints();
	
		public Metabolite getMetabolite(int index) {
		return metabolites[index];
	}

	public int getNumOfMetabolites() {
		return metabolites.length;
	}

	public void setMetabolites (Metabolite[] metab){
		this.metabolites = metab;
	}

	public void setMatrix (double[][] matrix){
		this.matrix = matrix;
	}
	public int getNumberOfReversibleReactions() {
		int count = 0;
		for (int i = 0; i < this.getNumOfReactions(); i++) {
			if (isReversible(i)) {
				count ++;
			}
		}
		return count;
	}
	public abstract double getLowerBound(int reactionIndex);

	public abstract double getUpperBound(int reactionIndex);

	public double[][] getStoichMatrix(){
		return this.matrix;
	}
	public abstract boolean isReversible(int reactionIndex);

	public abstract int containsReaction(String reacName);

	public abstract int getNumOfReactions();

	public abstract int getReactionIndex(String reac);

	public  double getStoichCoef(Metabolite meta, Reaction react){
		return this.matrix[getMetaboliteIndex(meta.getName())][getReactionIndex(react.getName())];
	}

	public double getStoichCoef(int metaboliteIndex, int reactionIndex){
		return this.matrix[metaboliteIndex][reactionIndex];
	}
	public void saveMetabolites(String filename) throws IOException{
		BufferedWriter metabfile = new BufferedWriter(new FileWriter(filename+".metabs"));
		for (int i = 0; i < this.getNumOfMetabolites(); i++) {
			Metabolite m = this.getMetabolite(i);
			String toWrite = m.getName()+";"+m.isExternal();
			metabfile.write(toWrite+"\n");
		}
		metabfile.flush();
		metabfile.close();
	}
	public abstract void saveNetwork(String filename) throws IOException;

	public void setBiomassReaction(String rname) {
		this.biomass = getReactionIndex(rname);
		
	}
	public int getBiomassReactionId() {
		return this.biomass;
	}
	
	public abstract String toString();
}
