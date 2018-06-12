package pt.uminho.ceb.biosystems.mcslibrary.solution.analysis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import cern.colt.Arrays;
import ilog.concert.IloException;
import pt.uminho.ceb.biosystems.mcslibrary.datastructures.AbstractTree;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba.CPLEXFluxBalanceAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Pair;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;

public class TurnoverDependencyChecker {
	
	public static final String PRODUCING = "PRODUCING";
	public static final String CONSUMING = "CONSUMING";
	private DefaultMetabolicNetwork dmn;
	private FluxBound[] fb;
	private CPLEXFluxBalanceAnalysis fba;

	public TurnoverDependencyChecker(DefaultMetabolicNetwork dmn, FluxBound[] fb) throws IloException {
		this.fba = new CPLEXFluxBalanceAnalysis(dmn);
		this.dmn = dmn;
		this.fb = fb;
	}
	
	public Map<String,AbstractTree<TurnoverScore>> checkMetabolite(String metaboliteName) throws IloException{
		int mId = dmn.getMetaboliteIndex(metaboliteName);
		List<Pair<Double,Reaction>> fullObjProd = new ArrayList<Pair<Double,Reaction>>();
		List<Pair<Double,Reaction>> fullObjCons = new ArrayList<Pair<Double,Reaction>>();
		for (int i = 0; i < dmn.getNumOfReactions(); i++) {
			double coef = dmn.getStoichCoef(mId, i);
			if (Math.abs(coef) > Utilities.PRECISION) {
				if (dmn.isReversible(i) && coef > 0) {
					fullObjProd.add(new Pair<Double, Reaction>(1.0, dmn.getReaction(i)));
					fullObjCons.add(new Pair<Double, Reaction>(-1.0, dmn.getReaction(i)));
				} else if (dmn.isReversible(i) && coef < 0){
					fullObjProd.add(new Pair<Double, Reaction>(-1.0, dmn.getReaction(i)));
					fullObjCons.add(new Pair<Double, Reaction>(1.0, dmn.getReaction(i)));
				} else if (coef > 0){
					fullObjProd.add(new Pair<Double, Reaction>(1.0, dmn.getReaction(i)));
				} else if (coef < 0){
					fullObjCons.add(new Pair<Double, Reaction>(1.0, dmn.getReaction(i)));
				}
			}
		}


		Pair<Double, Reaction>[] ocons = fullObjCons.toArray(new Pair[fullObjCons.size()]);
		Pair<Double, Reaction>[] oprod = fullObjProd.toArray(new Pair[fullObjProd.size()]);
		
		System.out.println("Full objective - consumption: "+Arrays.toString(ocons));
		System.out.println("Full objective - production: "+Arrays.toString(oprod));
		
		AbstractTree<TurnoverScore> dCons = new AbstractTree<TurnoverScore>(new TurnoverScore(ocons, fba.solve(fb, ocons, "min").getObjectiveValue()));
		AbstractTree<TurnoverScore> dProd = new AbstractTree<TurnoverScore>(new TurnoverScore(oprod, fba.solve(fb, oprod, "min").getObjectiveValue()));
		
		branchCurrentNode(dCons);
		branchCurrentNode(dProd);

		Map<String, AbstractTree<TurnoverScore>> res = new HashMap<>();
		res.put(TurnoverDependencyChecker.PRODUCING, dProd);
		res.put(TurnoverDependencyChecker.CONSUMING, dCons);

		return res;
	}
	
	public void branchCurrentNode(AbstractTree<TurnoverScore> tree) throws IloException{
		
		if (tree.getValue().objective.length > 0) {
			for (int i = 1; i < tree.getValue().objective.length; i++) {
				Pair<Double,Reaction>[] obj = new Pair[tree.getValue().objective.length - i];
				for (int j = i; j < tree.getValue().objective.length; j++) {
					obj[j-i] = tree.getValue().objective[j];
				}
				double value = fba.solve(fb, obj, "min").getObjectiveValue();
				if (value > Utilities.PRECISION) {
					tree.addChild(new AbstractTree<TurnoverScore>(new TurnoverScore(obj, value)));
				}
			}
			for (int i = 0; i < tree.getNumberOfChildren(); i++) {
				branchCurrentNode(tree.getChild(i));
			}
		}
	}

}
