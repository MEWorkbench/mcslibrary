package pt.uminho.ceb.biosystems.mcslibrary.metabolic.fva;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.concert.IloObjectiveSense;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.Status;
import java.util.HashMap;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.AbstractMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.YieldConstraint;
import pt.uminho.ceb.biosystems.mcslibrary.solution.SolutionUtilities;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Pair;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;

public class CPLEXFluxVariabilityAnalysis {
	private static final IloObjectiveSense[] directions = {IloObjectiveSense.Minimize,IloObjectiveSense.Maximize};
	private AbstractMetabolicNetwork metaNet;
	private FluxBound[] fluxrestrictions;
	private YieldConstraint[] yieldrestrictions;

	public CPLEXFluxVariabilityAnalysis(AbstractMetabolicNetwork metaNet, FluxBound[] fluxrestrictions, YieldConstraint[] yieldrestrictions) {
		this.metaNet = metaNet;
		this.fluxrestrictions = fluxrestrictions;
		this.yieldrestrictions = yieldrestrictions;
	}
	
	public CPLEXFluxVariabilityAnalysis(AbstractMetabolicNetwork metaNet, FluxBound[] fluxrestrictions, Reaction[] knockouts, YieldConstraint[] yieldrestrictions) {
		this.metaNet = metaNet;
		this.fluxrestrictions = new FluxBound[fluxrestrictions.length + knockouts.length];
		for (int i = 0; i < fluxrestrictions.length; i++) {
			this.fluxrestrictions[i] = fluxrestrictions[i];
		}
		for (int i = 0; i < knockouts.length; i++) {
			this.fluxrestrictions[i + fluxrestrictions.length] = new FluxBound(knockouts[i], 0, 0);
		}
		this.yieldrestrictions = yieldrestrictions;
	}

	public FluxVariabilityAnalysisResult solveFVA() throws IloException{
		HashMap<Integer, Pair<Double,Double>> res = new HashMap<Integer, Pair<Double,Double>>();
		IloNumVar[] variables = new IloNumVar[metaNet.getNumOfReactions()];
		IloRange[] constraints = new IloRange[metaNet.getNumOfMetabolites()];
		IloRange[] yields = new IloRange[this.yieldrestrictions.length];
		IloCplex cplex = new IloCplex();
		cplex.setParam(IloCplex.BooleanParam.NumericalEmphasis, true);

		cplex.setOut(null);
		int assumed = 0;
		for (int i = 0; i < metaNet.getNumOfReactions(); i++) {
			variables[i] = cplex.numVar(metaNet.getLowerBound(i), metaNet.getUpperBound(i),"R"+i);
			if (variables[i].getLB() == 0 && variables[i].getUB() == 0) {
				assumed++;
			}
		}
		
		for (int i = 0; i < metaNet.getNumOfMetabolites(); i++) {
			IloLinearNumExpr linExp = cplex.linearNumExpr();
			for (int j = 0; j < metaNet.getNumOfReactions(); j++) {
				linExp.addTerm(metaNet.getStoichCoef(i, j), variables[j]);
			}
			constraints[i] = cplex.eq(linExp, 0, metaNet.getMetabolite(i).getName());
		}
//		System.out.println("Environmental conditions for FVA:");
//		System.out.println("\t Assuming "+assumed+" reactions as being disabled.");
		if (this.fluxrestrictions != null) {
			for (int i = 0; i < this.fluxrestrictions.length; i++) {

				FluxBound rest = this.fluxrestrictions[i];
//				System.out.println(rest);
				int varIdx = metaNet.containsReaction(rest.getReac().getName());
				if (varIdx > -1) {
					double lb = rest.getBounds().getLower();
					double ub = rest.getBounds().getUpper();
					if (lb != -Utilities.INF) {
						variables[varIdx].setLB(lb);
					}

					if (ub != Utilities.INF) {
						variables[varIdx].setUB(ub);
					}
				} else {
					System.out.println("Did not apply flux constraint for "+varIdx);
				}

			}
		}

		if (this.yieldrestrictions != null) {
			for (int i = 0; i < this.yieldrestrictions.length; i++) {
				YieldConstraint rest = this.yieldrestrictions[i];
				int uIdx = metaNet.containsReaction(rest.getUptakeReaction().getName());
				int pIdx = metaNet.containsReaction(rest.getProductReaction().getName());
				if (uIdx > -1 && pIdx > -1) {
					double deviation = rest.getDeviation();
					double ratio = rest.getRatio();
//					System.out.println("\t"+rest.getProductReaction().getName()+" production higher than "+ratio*variables[uIdx].getLB());
					variables[pIdx].setLB(ratio*variables[uIdx].getLB());
					IloLinearNumExpr exp = cplex.linearNumExpr();
					if (rest.isLower()) {
						exp.addTerm(1, variables[pIdx]);
						exp.addTerm(-ratio, variables[uIdx]);
						constraints[i] = cplex.le(exp, deviation);
					} else {
						exp.addTerm(-1, variables[pIdx]);
						exp.addTerm(ratio, variables[uIdx]);
						constraints[i] = cplex.le(exp, -deviation);
					}
				} else {
					System.out.println("Did not apply yield constraint");
				}

			}
		}


		cplex.add(variables);
		cplex.add(constraints);
		cplex.add(yields);

		for (int i = 0; i < variables.length; i++) {
			IloLinearNumExpr linExp = cplex.linearNumExpr();
			linExp.addTerm(1, variables[i]);
			double[] mmaxpair = {variables[i].getLB(),variables[i].getUB()}; 
			for (int j = 0; j < directions.length; j++) {
				IloObjective obj = cplex.objective(directions[j],linExp);
				cplex.add(obj);
				try {
					cplex.solve();
					Status status = cplex.getStatus();
//					System.out.println(status);
					if (status != IloCplex.Status.Optimal) {
						mmaxpair[j] = (j == 0) ? -Utilities.INF : Utilities.INF;
					} else {
						mmaxpair[j] = cplex.getValue(variables[i]);
					}
				} catch (Exception e) {
				} finally {
					cplex.remove(obj);
				}
			}
			Pair<Double,Double> range = new Pair<Double, Double>(mmaxpair[0], mmaxpair[1]);
			res.put(i, range);
		}
		return new FluxVariabilityAnalysisResult(metaNet, res);

	}
	
}
