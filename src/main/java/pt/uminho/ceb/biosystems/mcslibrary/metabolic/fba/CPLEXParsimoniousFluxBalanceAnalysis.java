package pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba;

import java.util.Arrays;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.concert.IloObjectiveSense;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.DoubleParam;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.AbstractMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.SimulationResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.ReactionConstraint;

public class CPLEXParsimoniousFluxBalanceAnalysis {
	private AbstractMetabolicNetwork metaNet;
	private CPLEXFluxBalanceAnalysis fba;
	private IloCplex cplex;
	private IloNumVar[] fvars;
	private IloNumVar[] rvars;
	private IloObjective objectiveFun;
	private IloObjective originalObjectiveFun;
	private double relaxCoefficient;

	public CPLEXParsimoniousFluxBalanceAnalysis(AbstractMetabolicNetwork abs, double relaxCoef) throws IloException {
		this.fba = new CPLEXFluxBalanceAnalysis(abs);
		this.metaNet = abs;
		this.relaxCoefficient = relaxCoef;
		cplex = new IloCplex();
		cplex.setOut(null);
		cplex.setWarning(null);
		cplex.setParam(DoubleParam.EpRHS, 1e-9);
		cplex.setParam(DoubleParam.EpLin, 1e-9);

		fvars = new IloNumVar[abs.getNumOfReactions()];
		rvars = new IloNumVar[abs.getNumOfReactions()];

		for (int i = 0; i < fvars.length; i++) {
			ReactionConstraint fwd = abs.getBound(i).getSplitForward();
			ReactionConstraint rev = abs.getBound(i).getSplitReverse();
			fvars[i] = cplex.numVar(fwd.getLower(), fwd.getUpper(), "Rf"+i);
			rvars[i] = cplex.numVar(rev.getLower(), rev.getUpper(), "Rr"+i);
		}
		
//		for (int i = 0; i < rvars.length; i++) {
//
//		}

//		for (int i = 0; i < fvars.length; i++) {
//			System.out.println("Reaction "+i+" is reversible? "+abs.isReversible(i)+". Bounds: FWD["+fvars[i].getLB()+","+fvars[i].getUB()+"] REV["+rvars[i].getLB()+","+rvars[i].getUB()+"]");
//		}
		double[][] smat = abs.getStoichMatrix();
		IloRange[] constraints = new IloRange[abs.getNumOfMetabolites()];
		for (int i = 0; i < smat.length; i++) {
			IloLinearNumExpr row = cplex.linearNumExpr();
			for (int j = 0; j < smat[0].length; j++) {
				row.addTerm(smat[i][j], fvars[j]);
				row.addTerm(-smat[i][j], rvars[j]);
			}
			IloRange cons = cplex.eq(row, 0);
			constraints[i] = cons;
		}
		
		cplex.add(fvars);
		cplex.add(rvars);
		cplex.add(constraints);
		double[] ones = new double[abs.getNumOfReactions()];
		Arrays.fill(ones, 1);
		
		IloLinearNumExpr lexp = cplex.linearNumExpr();
		lexp.addTerms(fvars, ones);
		lexp.addTerms(rvars, ones);

		objectiveFun = cplex.objective(IloObjectiveSense.Minimize, lexp);
		originalObjectiveFun = cplex.objective(IloObjectiveSense.Minimize, lexp);
		cplex.add(originalObjectiveFun);
	}
	
	public SimulationResult solve(FluxBound[] bounds, Reaction objective, String sense) throws IloException{
		SimulationResult nonParsimonious = fba.solve(bounds, objective, sense);
		
		int objId = metaNet.getReactionIndex(objective.getName());
		double objValue = nonParsimonious.getFluxValue(objective)*this.relaxCoefficient;
		IloLinearNumExpr objExp = (IloLinearNumExpr) objectiveFun.getExpr();
		
		cplex.remove(originalObjectiveFun);
		objExp.remove(fvars[objId]);
		objExp.remove(rvars[objId]);
		objectiveFun.setExpr(objExp);
		cplex.add(objectiveFun);
//		System.out.println("Objective ID: "+objId+"; Objective value: "+objValue);
//		System.out.println("Using objective function: "+objectiveFun);
		
		if (bounds != null) {
			for (int i = 0; i < bounds.length; i++) {
				FluxBound fbound = bounds[i];
				int index = metaNet.getReactionIndex(fbound.getReac().getName());
				ReactionConstraint fwd = getForwardBound(fbound.getBounds().getLower(), fbound.getBounds().getUpper());
				ReactionConstraint rev = getReverseBound(fbound.getBounds().getLower(), fbound.getBounds().getUpper());

				fvars[index].setLB(fwd.getLower());
				fvars[index].setUB(fwd.getUpper());
				rvars[index].setLB(rev.getLower());
				rvars[index].setUB(rev.getUpper());
//				System.out.println("Bound vars before solving:"+fvars[index]+";"+fvars[index].getLB()+"-"+fvars[index].getUB()+"|"+rvars[index]+";"+rvars[index].getLB()+"-"+rvars[index].getUB());
			}
		}
		
		if (objValue > 0) {
			fvars[objId].setLB(objValue);
			fvars[objId].setUB(objValue);
			rvars[objId].setLB(0);
			rvars[objId].setUB(0);
		} else if (objValue < 0) {
			fvars[objId].setLB(0);
			fvars[objId].setUB(0);
			rvars[objId].setLB(-objValue);
			rvars[objId].setUB(-objValue);
		}

		double[] values = new double[fvars.length];
		cplex.exportModel("pfba_before.lp");
		try {
			cplex.solve();
			double[] fvalues = cplex.getValues(fvars);
			double[] rvalues = cplex.getValues(rvars);
			for (int i = 0; i < values.length; i++) {
				values[i] = fvalues[i] - rvalues[i];
			}
//			System.out.println("Status = "+cplex.getStatus()+"; Primal objective value: "+objValue+"; Norm: "+cplex.getObjValue());
			if (cplex.getStatus() != IloCplex.Status.Optimal) {
				for (int i = 0; i < values.length; i++) {
					values[i] = Double.NaN;
				}
			}
		} catch (Exception e) {
			System.out.println("Error!");
			e.printStackTrace();
			return new SimulationResult(metaNet,new double[fvars.length]);
		}
		if (bounds != null) {
			for (int i = 0; i < bounds.length; i++) {
				FluxBound fbound = bounds[i];
				int index = metaNet.getReactionIndex(fbound.getReac().getName());
				ReactionConstraint fwd = getForwardBound(this.metaNet.getLowerBound(index), this.metaNet.getUpperBound(index));
				ReactionConstraint rev = getReverseBound(this.metaNet.getLowerBound(index), this.metaNet.getUpperBound(index));

				fvars[index].setLB(fwd.getLower());
				fvars[index].setUB(fwd.getUpper());
				rvars[index].setLB(rev.getLower());
				rvars[index].setUB(rev.getUpper());
//				System.out.println("Bound vars after solving:"+fvars[index]+";"+fvars[index].getLB()+"-"+fvars[index].getUB()+"|"+rvars[index]+";"+rvars[index].getLB()+"-"+rvars[index].getUB());
			}
		}
		
		ReactionConstraint obfrc = getForwardBound(this.metaNet.getLowerBound(objId), this.metaNet.getUpperBound(objId));
		ReactionConstraint obrrc = getReverseBound(this.metaNet.getLowerBound(objId), this.metaNet.getUpperBound(objId));
		objExp.remove(fvars[objId]);
		objExp.remove(rvars[objId]);

		fvars[objId].setLB(obfrc.getLower());
		fvars[objId].setUB(obfrc.getUpper());
		rvars[objId].setLB(obrrc.getLower());
		rvars[objId].setUB(obrrc.getUpper());
		
		cplex.remove(objectiveFun);
		cplex.add(originalObjectiveFun);
		cplex.exportModel("pfba_after.lp");

		return new SimulationResult(metaNet,values);
	}
	
	public SimulationResult solveKnockoutFBA(FluxBound[] bounds, Reaction[] knockouts, Reaction objective, String sense) throws IloException {
		FluxBound[] fbaBounds = new FluxBound[bounds.length+knockouts.length];
		for (int i = 0; i < bounds.length; i++) {
			fbaBounds[i] = bounds[i];
		}
		
		for (int i = 0; i < knockouts.length; i++) {
			fbaBounds[i+bounds.length] = new FluxBound(knockouts[i], 0, 0);
		}
		return solve(fbaBounds, objective, sense);
	}
	
	public ReactionConstraint getReverseBound(double LB, double UB){
		double nlb = 0;
		double nub = 0;
		if (UB < 0) {
			nlb = -UB;
		}
		if (LB < 0) {
			nub = -LB;
		}

		return new ReactionConstraint(nlb, nub > 10000 ? 999999 : nub);
	}
	
	public ReactionConstraint getForwardBound(double LB, double UB){
		double nlb = 0;
		double nub = 0;
		if (UB > 0) {
			nub = UB;
		}
		if (LB > 0) {
			nlb = LB;
		}

		return new ReactionConstraint(nlb, nub > 10000 ? 999999 : nub);
	}
}
