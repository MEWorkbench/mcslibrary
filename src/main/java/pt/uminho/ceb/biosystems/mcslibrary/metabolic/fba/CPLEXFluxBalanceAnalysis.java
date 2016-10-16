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
package pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.concert.IloObjectiveSense;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import java.util.ArrayList;
import java.util.List;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.AbstractMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;

public class CPLEXFluxBalanceAnalysis {

	AbstractMetabolicNetwork metaNet;
	IloCplex cplex;
	IloNumVar[] vars;
	IloRange[] constraints;

	public CPLEXFluxBalanceAnalysis(AbstractMetabolicNetwork abs) throws IloException {
		this.metaNet = abs;
		cplex = new IloCplex();
		cplex.setOut(null);
		vars = new IloNumVar[abs.getNumOfReactions()];
		constraints = new IloRange[abs.getNumOfMetabolites()];
		for (int i = 0; i < vars.length; i++) {
			vars[i] = cplex.numVar(abs.getLowerBound(i), abs.getUpperBound(i), "R"+i);
		}

		double[][] smat = abs.getStoichMatrix();

		for (int i = 0; i < smat.length; i++) {
			IloLinearNumExpr row = cplex.linearNumExpr();
			for (int j = 0; j < smat[0].length; j++) {
				row.addTerm(smat[i][j], vars[j]);
			}
			IloRange cons = cplex.eq(row, 0);
			constraints[i] = cons;
		}
		cplex.add(vars);
		cplex.add(constraints);
	}

	public CPLEXFluxBalanceAnalysis(AbstractMetabolicNetwork abs, double M) throws IloException {
		this.metaNet = abs;
		cplex = new IloCplex();
		cplex.setOut(null);
		vars = new IloNumVar[abs.getNumOfReactions()];
		constraints = new IloRange[abs.getNumOfMetabolites()];
		for (int i = 0; i < vars.length; i++) {
			vars[i] = cplex.numVar(abs.getLowerBound(i) < -M ? -M : abs.getLowerBound(i), abs.getUpperBound(i) > M ? M : abs.getUpperBound(i), "R"+i);
		}

		double[][] smat = abs.getStoichMatrix();

		for (int i = 0; i < smat.length; i++) {
			IloLinearNumExpr row = cplex.linearNumExpr();
			for (int j = 0; j < smat[0].length; j++) {
				row.addTerm(smat[i][j], vars[j]);
			}
			IloRange cons = cplex.eq(row, 0);
			constraints[i] = cons;
		}
		cplex.add(vars);
		cplex.add(constraints);
	}

	public ArrayList<FluxBalanceAnalysisResult> applyKnockout(List<List<String>> knockout, Reaction objective, FluxBound[] conditions, String sense) throws IloException{
		ArrayList<FluxBalanceAnalysisResult> res = new ArrayList<FluxBalanceAnalysisResult>();
		for (List<String> list : knockout) {
			System.out.println("Simulating "+list);
			int size = list.size();
			Reaction[] ko = new Reaction[size];
			for (int i = 0; i < size; i++) {
				String rname = list.get(i);
				if (metaNet.containsReaction(rname) < 0) {
					System.out.println("Reaction "+rname+" not found!!");
				}
				Reaction react = ((DefaultMetabolicNetwork) this.metaNet).getReaction(rname);
				ko[i] = react;
			}
			res.add(solveReactionKnockoutFBA(conditions, ko, objective, sense));
		}
		return res;
	}
	
	public double getMinimumFlux(FluxBound[] bounds, Reaction[] knockouts, Reaction objective) throws IloException{
		return this.solveReactionKnockoutFBA(bounds, knockouts, objective, "min").getFluxValue(objective);
	}
	public FluxBalanceAnalysisResult applyKnockout(Reaction objective, FluxBound[] conditions, String sense, List<String> knockout) throws IloException{
		int size = knockout.size();
		Reaction[] ko = new Reaction[size];
			for (int i = 0; i < size; i++) {
				String rname = knockout.get(i);
				Reaction react = ((DefaultMetabolicNetwork) this.metaNet).getReaction(rname);
				ko[i] = react;
			}
		FluxBalanceAnalysisResult res = solveReactionKnockoutFBA(conditions, ko, objective, sense);
		return res;
	}
	
	public double[] solveFVAsteps(FluxBound[] bounds, Reaction[] knockouts, int breaks, String sense, Reaction objective, Reaction variable) throws IloException{
		double[] res = new double[breaks+1];
		FluxBalanceAnalysisResult wt = solveReactionKnockoutFBA(bounds, knockouts, objective, sense);
		double maxObj = wt.getFluxValue(objective);
//		System.out.println(maxObj);
		res[0] = maxObj;
		for (int i = 0; i < breaks; i++) {
			res[i+1] = solveReactionKnockoutFVA(bounds, knockouts, objective, sense, (double)(i+1) / breaks, variable).getFluxValue(objective);
		}
		return res;
	}
	public FluxBalanceAnalysisResult solveReactionKnockoutFBA(FluxBound[] bounds, Reaction[] knockouts, Reaction objective, String sense) throws IloException {
		FluxBound[] fbaBounds = new FluxBound[bounds.length+knockouts.length];
		for (int i = 0; i < bounds.length; i++) {
			fbaBounds[i] = bounds[i];
		}
		
		for (int i = 0; i < knockouts.length; i++) {
			fbaBounds[i+bounds.length] = new FluxBound(knockouts[i], 0, 0);
		}
		return solve(fbaBounds, objective, sense);
	}
	
	public FluxBalanceAnalysisResult solveReactionKnockoutFVA(FluxBound[] bounds, Reaction[] knockouts, Reaction objective, String sense, double biomassPerc, Reaction biomassReaction) throws IloException {
		FluxBalanceAnalysisResult preRes = solveReactionKnockoutFBA(bounds, knockouts, biomassReaction, "max");
		double maxBio = preRes.getFluxValue(biomassReaction);
		FluxBound[] newBounds = new FluxBound[bounds.length+knockouts.length+1];
		for (int i = 0; i < bounds.length; i++) {
			newBounds[i] = bounds[i];
		}
		for (int i = 0; i < knockouts.length; i++) {
			newBounds[i+bounds.length] = new FluxBound(knockouts[i],0,0);
		}
		newBounds[bounds.length+knockouts.length] = new FluxBound(biomassReaction, biomassPerc*maxBio, Utilities.INF);
		for (int i = 0; i < newBounds.length; i++) {
		}
		return solveReactionKnockoutFBA(newBounds, knockouts, objective, sense);
	}
	
	
	public FluxBalanceAnalysisResult determineSOC(FluxBound[] bounds, Reaction[] knockouts, Reaction objective, Reaction variant, String sense, int resolution) throws IloException{
		FluxBound[] totalBounds = null;
		if (knockouts == null && (bounds != null)) {
			totalBounds = bounds;
		} else {
			totalBounds = new FluxBound[knockouts.length + bounds.length];
			for (int i = 0; i < knockouts.length; i++) {
				totalBounds[i] = new FluxBound(knockouts[i], 0, 0);
			}
			for (int i = 0; i < bounds.length; i++) {
				totalBounds[knockouts.length + i] = bounds[i];
			}
		}
		// iterate
		double[] objectiveVector = new double[resolution+1];
		objectiveVector[0] = solve(totalBounds, objective, sense).getFluxValue(objective);
		for (int i = 1; i < resolution; i++) {
			FluxBound[] varTotalBounds = new FluxBound[totalBounds.length + 1];
			objectiveVector[i] = solve(varTotalBounds, objective, sense).getFluxValue(objective);
		}
		objectiveVector[resolution] = solve(totalBounds, variant, "max").getFluxValue(objective);
		
		
		//	(product yield ^ 2) / slope
		return null;
		
		
		
	}
	public FluxBalanceAnalysisResult solve(FluxBound[] bounds, Reaction objective, String sense) throws IloException {
		if (bounds != null) {
			for (int i = 0; i < bounds.length; i++) {
				FluxBound fbound = bounds[i];
				int index = metaNet.getReactionIndex(fbound.getReac().getName());
				vars[index].setLB(fbound.getBounds().getLower());
				vars[index].setUB(fbound.getBounds().getUpper());
			}
		}
		IloNumVar objectiveVar = vars[metaNet.getReactionIndex(objective.getName())];
		IloLinearNumExpr objectiveExpr = cplex.linearNumExpr();
		objectiveExpr.addTerm(1, objectiveVar);
		IloObjective objectiveFun;
		int vsense = -1;

		if (sense == "max")
			vsense = 1;
		else if (sense == "min")
			vsense = 0;

		switch (vsense) {
		case 1:
			objectiveFun = cplex.objective(IloObjectiveSense.Maximize);
			break;
		case 0:
			objectiveFun = cplex.objective(IloObjectiveSense.Minimize);
			break;
		default:
			return null;
		}
		double[] values = null;
		objectiveFun.setExpr(objectiveExpr);
		cplex.add(objectiveFun);
		try {
			cplex.solve();
			values = cplex.getValues(vars);
			if (cplex.getStatus() != IloCplex.Status.Optimal) {
//				System.out.println("Exception found: Status = "+cplex.getStatus());
				for (int i = 0; i < values.length; i++) {
					values[i] = Double.NaN;
				}
			}
		} catch (Exception e) {
			System.out.println("Error!");
			e.printStackTrace();
			return new FluxBalanceAnalysisResult(metaNet,new double[vars.length]);
		}
		cplex.remove(objectiveFun);

		if (bounds != null) {
			for (int i = 0; i < bounds.length; i++) {
				FluxBound fbound = bounds[i];
				int index = metaNet.getReactionIndex(fbound.getReac().getName());
				vars[index].setLB(metaNet.getLowerBound(index));
				vars[index].setUB(metaNet.getUpperBound(index));
			}
		}
		return new FluxBalanceAnalysisResult(metaNet,values);
	}

}
