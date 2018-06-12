package pt.uminho.ceb.biosystems.mcslibrary.enumeration.implementation;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.concert.IloObjectiveSense;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.Status;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.AbstractMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fva.FluxVariabilityAnalysisResult;

public class CPLEXFlexibleFVA {
	private AbstractMetabolicNetwork metaNet;
	private IloCplex cplex;
	private IloNumVar[] reactions;
	private IloObjectiveSense[] senses = { IloObjectiveSense.Minimize, IloObjectiveSense.Maximize };
	private IloObjective[][] objectives;
	private IloLinearNumExpr[] expressions;

	public CPLEXFlexibleFVA(AbstractMetabolicNetwork metaNet) throws IloException {

		this.metaNet = metaNet;
		this.cplex = new IloCplex();
		this.cplex.setOut(null);
		buildProblem(cplex);
		
	}

	private void buildProblem(IloCplex cplex) throws IloException {
		expressions = new IloLinearNumExpr[metaNet.getNumOfReactions()];
		reactions = new IloNumVar[metaNet.getNumOfReactions()];
		objectives = new IloObjective[2][metaNet.getNumOfReactions()];
		for (int i = 0; i < reactions.length; i++) {
			reactions[i] = cplex.numVar(metaNet.getLowerBound(i), metaNet.getUpperBound(i));
			expressions[i] = cplex.linearNumExpr();
			expressions[i].addTerm(1, reactions[i]);
			objectives[0][i] = cplex.objective(senses[0], expressions[i]);
			objectives[1][i] = cplex.objective(senses[1], expressions[i]);

		}
		IloRange[] metabolites = new IloRange[metaNet.getNumOfMetabolites()];

		for (int i = 0; i < metabolites.length; i++) {
			IloLinearNumExpr metaboliteExpression = cplex.linearNumExpr();
			for (int j = 0; j < reactions.length; j++) {
				metaboliteExpression.addTerm(metaNet.getStoichCoef(i, j), reactions[j]);
			}
			metabolites[i] = cplex.eq(metaboliteExpression, 0);
		}
		cplex.add(metabolites);
//		this.cplex.tuneParam();
	}

	public FluxVariabilityAnalysisResult solveFVA(FluxBound[] fluxBound, Reaction[] knockouts) throws IloException {

		if (fluxBound != null) {
			for (int i = 0; i < fluxBound.length; i++) {
				int id = metaNet.getReactionIndex(fluxBound[i].getReac().getName());
				reactions[id].setLB(fluxBound[i].getBounds().getLower());
				reactions[id].setUB(fluxBound[i].getBounds().getUpper());
			}	
		}

		if (knockouts != null) {
			for (int i = 0; i < knockouts.length; i++) {
				int id = metaNet.getReactionIndex(knockouts[i].getName());
				reactions[id].setLB(0);
				reactions[id].setUB(0);
			}
		}


		cplex.add(reactions);

		double[][] result = new double[2][reactions.length];
		double start = System.currentTimeMillis();
		
		for (int i = 0; i < reactions.length; i++) {
			for (int j = 0; j < senses.length; j++) {
				cplex.add(objectives[j][i]);
				try {
					cplex.solve();
					if (cplex.getStatus() == Status.Optimal) {
						result[j][i] = cplex.getObjValue();
					} else if (cplex.getStatus() == Status.Unbounded) {
						result[j][i] = j == 1 ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
					}
				} catch (Exception e) {
					result[j][i] = Double.NaN;
				}
				cplex.remove(objectives[j][i]);
			}
		}
		
		double finish = System.currentTimeMillis();
		System.out.println((finish-start)/1000);
		if (fluxBound != null) {
			for (int i = 0; i < fluxBound.length; i++) {
				int id = metaNet.getReactionIndex(fluxBound[i].getReac().getName());
				reactions[id].setLB(metaNet.getLowerBound(id));
				reactions[id].setUB(metaNet.getUpperBound(id));
			}
		}

		if (knockouts != null) {
			for (int i = 0; i < knockouts.length; i++) {
				int id = metaNet.getReactionIndex(knockouts[i].getName());
				reactions[id].setLB(metaNet.getLowerBound(id));
				reactions[id].setUB(metaNet.getUpperBound(id));
			}
		}
		cplex.remove(reactions);

		
		return new FluxVariabilityAnalysisResult(metaNet, result);
	}
}
