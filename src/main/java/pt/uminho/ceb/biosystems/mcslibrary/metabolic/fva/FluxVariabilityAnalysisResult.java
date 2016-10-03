package pt.uminho.ceb.biosystems.mcslibrary.metabolic.fva;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.AbstractMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.ReactionGroup;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.CompressedMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Pair;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;


public class FluxVariabilityAnalysisResult {
	private static double epsilon = 1e-9;
	private AbstractMetabolicNetwork metaNet;
	private HashMap<Integer, Pair<Double,Double>> results;
	private boolean[] blockedReactions;


	public FluxVariabilityAnalysisResult(AbstractMetabolicNetwork metaNet, HashMap<Integer, Pair<Double, Double>> results) {
		this.metaNet = metaNet;
		this.results = results;
		this.blockedReactions = getBlockedReactions();
	}

	private Pair<Double, Double> getRange(int reactionindex){
		return results.get(reactionindex);
	}

	public double minFlux(int reactionindex){
		return getRange(reactionindex).getA();
	}
	public double maxFlux(int reactionindex){
		return getRange(reactionindex).getB();
	}

	public AbstractMetabolicNetwork getMetabolicNetwork(){
		return this.metaNet;
	}
	
	public double[] minToArray(){
		double[] res = new double[metaNet.getNumOfReactions()];
		for (int i = 0; i < metaNet.getNumOfReactions(); i++) {
			res[i] = minFlux(i);
		}
		return res;
	}
	
	public double[] maxToArray(){
		double[] res = new double[metaNet.getNumOfReactions()];
		for (int i = 0; i < metaNet.getNumOfReactions(); i++) {
			res[i] = maxFlux(i);
		}
		return res;
	}
	
	public void printResults(){
		for (int i = 0; i < results.size(); i++) {
			if (blockedReactions[i]) {
				System.out.println("BLOCKED REACTION - "+i+": "+minFlux(i)+" < V"+i+" > "+maxFlux(i));
			} else {
				System.out.println(i+": "+minFlux(i)+" < V"+i+" > "+maxFlux(i));
			}
		}
	}
	
	public void printResults(DefaultMetabolicNetwork metaNet){
		for (int i = 0; i < results.size(); i++) {
			if (blockedReactions[i]) {
				System.out.println("BLOCKED REACTION - "+metaNet.getReaction(i).getName()+": "+minFlux(i)+" < V"+i+" > "+maxFlux(i));
			} else {
				System.out.println(metaNet.getReaction(i).getName()+": "+minFlux(i)+" < V"+i+" > "+maxFlux(i));
			}
		}
	}
	public void printBlockedReactions() {
		for (String r : getBlockedReactionNames()) {
			System.out.println(r);
		}
	}
	private boolean[] getBlockedReactions(){
		boolean[] res = new boolean[this.metaNet.getNumOfReactions()];
		for (int i = 0; i < metaNet.getNumOfReactions(); i++) {
			if (Math.abs(minFlux(i)) < epsilon && Math.abs(maxFlux(i)) < epsilon) {
				res[i] = true;
			} else {
				res[i] = false;
			}
		}
		return res;
	}
	public int countBlockedReactions(){
		int res = 0;
		for (int i = 0; i < blockedReactions.length; i++) {
			if (blockedReactions[i]){
				res++;
			}
		}
		return res;
	}
	public Set<String> getBlockedReactionNames(){
		Set<String> res = new HashSet<String>();
		if (metaNet.getClass() == DefaultMetabolicNetwork.class) {
			for (int i = 0; i < blockedReactions.length; i++) {
				if (blockedReactions[i]) {
					res.add(((DefaultMetabolicNetwork) metaNet).getReaction(i).getName());
				}
			}
		} else if (metaNet.getClass() == CompressedMetabolicNetwork.class){
			for (int i = 0; i < blockedReactions.length; i++) {
				if (blockedReactions[i]) {
					ReactionGroup r = ((CompressedMetabolicNetwork) metaNet).getReactionGroup(i);
					for (int j = 0; j < r.size(); j++) {
						res.add(r.getReaction(j).getName());
					}
				}
			}
		}

		return res;
	}
	
	public boolean isUndefined(int index) {
		if (minFlux(index) == Double.NaN || maxFlux(index) == Double.NaN) {
			return true;
		} else {
			return false;
		}
		
	}
	
	public boolean isUnbounded(int index){
		if (minFlux(index) == Double.NEGATIVE_INFINITY || maxFlux(index) == Double.POSITIVE_INFINITY) {
			return true;
		} else {
			return false;
		}			
	}
	
	public boolean isEssential(int index){
		if (maxFlux(index) < -Utilities.EPSILON || minFlux(index) > Utilities.EPSILON) {
			return true;
		} else {
			return false;
		}
	}
	public boolean isBlocked(int reactionindex){
		return blockedReactions[reactionindex];
	}
	
	public void writeToFile(String path) throws IOException{
		BufferedWriter bf = new BufferedWriter(new FileWriter(path));
		for (int i = 0; i < metaNet.getNumOfReactions(); i++) {
			String string = i+","+minFlux(i)+","+maxFlux(i);
			bf.write(string+"\n");
		}
		bf.flush();
		bf.close();
		
	}

	public Set<String> getEssentialReactionNames() {
		Set<String> res = new HashSet<String>();
		if (metaNet.getClass() == DefaultMetabolicNetwork.class) {
			for (int i = 0; i < metaNet.getNumOfReactions(); i++) {
				if (blockedReactions[i]) {
					res.add(((DefaultMetabolicNetwork) metaNet).getReaction(i).getName());
				}
			}
		} else if (metaNet.getClass() == CompressedMetabolicNetwork.class){
			for (int i = 0; i < metaNet.getNumOfReactions(); i++) {
				if (blockedReactions[i]) {
					ReactionGroup r = ((CompressedMetabolicNetwork) metaNet).getReactionGroup(i);
					for (int j = 0; j < r.size(); j++) {
						res.add(r.getReaction(j).getName());
					}
				}
			}
		}
		return res;
	}

}
