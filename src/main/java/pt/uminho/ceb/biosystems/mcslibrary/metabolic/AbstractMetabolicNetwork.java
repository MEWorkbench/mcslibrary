package pt.uminho.ceb.biosystems.mcslibrary.metabolic;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.ReactionConstraint;

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
	
	public Metabolite[] getMetabolites(){
		return metabolites;
	}
	
	public abstract ReactionConstraint getBound(int index);
	
	public abstract boolean isReversible(int reactionIndex);

	public abstract int containsReaction(String reacName);

	public abstract int getNumOfReactions();

	public abstract int getReactionIndex(String reac);

	public abstract Reaction getReaction(String reac);
	
	public  double getStoichCoef(Metabolite meta, Reaction react){
		return this.matrix[getMetaboliteIndex(meta.getName())][getReactionIndex(react.getName())];
	}

	public double getStoichCoef(int metaboliteIndex, int reactionIndex){
		return this.matrix[metaboliteIndex][reactionIndex];
	}
	
	public double getStoichCoef(String metaboliteId, String reactionId){
		return this.matrix[getMetaboliteIndex(metaboliteId)][getReactionIndex(reactionId)];
	}
	
	public void setStoichCoef(String metaboliteId, String reactionId, double coef){
		this.matrix[getMetaboliteIndex(metaboliteId)][getReactionIndex(reactionId)] = coef;
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
