package pt.uminho.ceb.biosystems.mcslibrary.solution.analysis;

import java.util.ArrayList;
import java.util.List;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Metabolite;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;

public class AnalysisUtilities {
	public static List<Reaction> getAllReactionsProducing(Metabolite metabolite, DefaultMetabolicNetwork mn){
		List<Reaction> r = new ArrayList<Reaction>();
		int idx = mn.getMetaboliteIndex(metabolite.getName());
		for (int i = 0; i < mn.getNumOfReactions(); i++) {
			if (mn.getStoichCoef(idx, i) > 0 || (mn.getStoichCoef(idx, i) < 0 && mn.isReversible(i))) {
				r.add(mn.getReaction(i));
			}
		}
		return r;
	}
	public static List<Reaction> getAllReactionsConsuming(Metabolite metabolite, DefaultMetabolicNetwork mn){
		List<Reaction> r = new ArrayList<Reaction>();
		int idx = mn.getMetaboliteIndex(metabolite.getName());
		for (int i = 0; i < mn.getNumOfReactions(); i++) {
			if (mn.getStoichCoef(idx, i) < 0 || (mn.getStoichCoef(idx, i) > 0 && mn.isReversible(i))) {
				r.add(mn.getReaction(i));
			}
		}
		return r;
	}
	
	public static List<String> getAllReactionsProducing(String metaboliteName, DefaultMetabolicNetwork mn){
		List<String> r = new ArrayList<String>();
		int idx = mn.getMetaboliteIndex(metaboliteName);
		for (int i = 0; i < mn.getNumOfReactions(); i++) {
			if (mn.getStoichCoef(idx, i) > 0 || (mn.getStoichCoef(idx, i) < 0 && mn.isReversible(i))) {
				r.add(mn.getReaction(i).getName());
			}
		}
		return r;
	}
	
	public static List<String> getAllReactionsConsuming(String metaboliteName, DefaultMetabolicNetwork mn){
		List<String> r = new ArrayList<String>();
		int idx = mn.getMetaboliteIndex(metaboliteName);
		for (int i = 0; i < mn.getNumOfReactions(); i++) {
			if (mn.getStoichCoef(idx, i) < 0 || (mn.getStoichCoef(idx, i) > 0 && mn.isReversible(i))) {
				r.add(mn.getReaction(i).getName());
			}
		}
		return r;
	}
}
