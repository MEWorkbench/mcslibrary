package pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints;

import java.io.IOException;
import java.io.Serializable;
import java.util.List;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Pair;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;

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
	private String name;
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
	
	public FluxBound(Reaction reac, double lowerbound, double upperbound, String name) {
		this.reac = reac;
		this.bounds = new ReactionConstraint(lowerbound, upperbound);
		this.overrideLB = false;
		this.overrideUB = false;
		this.setName(name);
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
	
	public static Pair<String,FluxBound[]> fromFile(String file, DefaultMetabolicNetwork mn) throws IOException{
		List<String> lines = Utilities.readLines(file);
		String name = lines.get(0);
		FluxBound[] bounds = new FluxBound[lines.size() - 1];
		for (int i = 0; i < bounds.length; i++) {
			bounds[i] = fromLine(lines.get(i+1), mn);
		}
		return new Pair<String, FluxBound[]>(name, bounds);
	}
	
	public static FluxBound fromLine(String line, DefaultMetabolicNetwork mn){
		List<String> tokens = Utilities.getAllStringTokens(line, ",");
		return new FluxBound(mn.getReaction(tokens.get(0)), Double.parseDouble(tokens.get(1)), Double.parseDouble(tokens.get(2)));
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

}
