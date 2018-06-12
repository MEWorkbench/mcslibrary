package pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Pair;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;

public class CPLEXFluxBalanceAnalysis {

	AbstractMetabolicNetwork metaNet;
	IloCplex cplex;
	IloNumVar[] vars;
	IloRange[] constraints;

	public CPLEXFluxBalanceAnalysis(AbstractMetabolicNetwork abs) throws IloException {
		this.metaNet = abs;
		cplex = new IloCplex();
		cplex.setParam(DoubleParam.EpRHS, 1e-9);
		cplex.setParam(DoubleParam.BarEpComp, 1e-9);
		cplex.setParam(DoubleParam.WorkMem, 12000);
		cplex.setParam(DoubleParam.EpLin, 1e-9);
		cplex.setWarning(null);
		cplex.setOut(null);
		

//		try {
//			cplex.setOut(new FileOutputStream("cplex_fba.out"));
//		} catch (FileNotFoundException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
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

	public AbstractMetabolicNetwork getMetabolicNetwork(){
		return this.metaNet;
	}
	
	public ArrayList<SimulationResult> applyKnockout(List<List<String>> knockout, Reaction objective, FluxBound[] conditions, String sense) throws IloException{
		ArrayList<SimulationResult> res = new ArrayList<SimulationResult>();
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
	
	public double getMaximumFlux(FluxBound[] bounds, Reaction[] knockouts, Reaction objective) throws IloException{
		return this.solveReactionKnockoutFBA(bounds, knockouts, objective, "max").getFluxValue(objective);
	}
	
	public SimulationResult applyKnockout(Reaction objective, FluxBound[] conditions, String sense, List<String> knockout) throws IloException{
		int size = knockout.size();
		Reaction[] ko = new Reaction[size];
			for (int i = 0; i < size; i++) {
				String rname = knockout.get(i);
				Reaction react = ((DefaultMetabolicNetwork) this.metaNet).getReaction(rname);
				ko[i] = react;
			}
		SimulationResult res = solveReactionKnockoutFBA(conditions, ko, objective, sense);
		return res;
	}
	
	public double[] solveFVAsteps(FluxBound[] bounds, Reaction[] knockouts, int breaks, String sense, Reaction objective, Reaction variable) throws IloException{
		double[] res = new double[breaks+1];
		for (int i = 0; i <= breaks; i++) {
//			SimulationResult min = solve(bounds, variable, "max");
//			SimulationResult max = solve(bounds, variable, "max");
			res[i] = solveReactionKnockoutFVA(bounds, knockouts, objective, sense, ((double)i / breaks)*0.999999, variable).getObjectiveValue();
		}
		return res;
	}
	public SimulationResult solveReactionKnockoutFBA(FluxBound[] bounds, Reaction[] knockouts, Reaction objective, String sense) throws IloException {
		FluxBound[] fbaBounds = new FluxBound[bounds.length+knockouts.length];
		for (int i = 0; i < bounds.length; i++) {
			fbaBounds[i] = bounds[i];
		}
		
		for (int i = 0; i < knockouts.length; i++) {
			fbaBounds[i+bounds.length] = new FluxBound(knockouts[i], 0, 0);
		}
		return solve(fbaBounds, objective, sense);
	}
	
	public SimulationResult solveReactionKnockoutFBA(FluxBound[] bounds, Reaction[] knockouts, Pair<Double,Reaction>[] objective, String sense) throws IloException {
		FluxBound[] fbaBounds = new FluxBound[bounds.length+knockouts.length];
		for (int i = 0; i < bounds.length; i++) {
			fbaBounds[i] = bounds[i];
		}
		
		for (int i = 0; i < knockouts.length; i++) {
			fbaBounds[i+bounds.length] = new FluxBound(knockouts[i], 0, 0);
		}
		return solve(fbaBounds, objective, sense);
	}
	
	public SimulationResult solveReactionKnockoutFVA(FluxBound[] bounds, Reaction[] knockouts, Reaction objective, String sense, double biomassPerc, Reaction biomassReaction) throws IloException {
		SimulationResult preRes = solveReactionKnockoutFBA(bounds, knockouts, biomassReaction, "max");
		double maxBio = preRes.getFluxValue(biomassReaction);
		FluxBound[] newBounds = new FluxBound[bounds.length+knockouts.length+1];
		for (int i = 0; i < bounds.length; i++) {
			newBounds[i] = bounds[i];
		}
		for (int i = 0; i < knockouts.length; i++) {
			newBounds[i+bounds.length] = new FluxBound(knockouts[i],0,0);
		}
		newBounds[bounds.length+knockouts.length] = new FluxBound(biomassReaction, biomassPerc*maxBio, Utilities.INF);
//		System.out.println(newBounds[bounds.length+knockouts.length]);
		
		return solveReactionKnockoutFBA(newBounds, knockouts, objective, sense);
	}
	
	public SimulationResult solveReactionKnockoutFVA(FluxBound[] bounds, Reaction[] knockouts, Reaction objective, String sense, double biomassPerc, Reaction biomassReaction, boolean fixed) throws IloException {
		SimulationResult preRes = solveReactionKnockoutFBA(bounds, knockouts, biomassReaction, "max");
		double maxBio = preRes.getFluxValue(biomassReaction);
		FluxBound[] newBounds = new FluxBound[bounds.length+knockouts.length+1];
		for (int i = 0; i < bounds.length; i++) {
			newBounds[i] = bounds[i];
		}
		for (int i = 0; i < knockouts.length; i++) {
			newBounds[i+bounds.length] = new FluxBound(knockouts[i],0,0);
		}
		newBounds[bounds.length+knockouts.length] = new FluxBound(biomassReaction, biomassPerc*maxBio, fixed ? biomassPerc*maxBio : Utilities.INF);
		
		return solveReactionKnockoutFBA(newBounds, knockouts, objective, sense);
	}
	
	
	public SimulationResult determineSOC(FluxBound[] bounds, Reaction[] knockouts, Reaction objective, Reaction variant, String sense, int resolution) throws IloException{
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
	public SimulationResult solve(FluxBound[] bounds, Reaction objective, String sense) throws IloException {
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

		if (sense.equals("max"))
			vsense = 1;
		else if (sense.equals("min"))
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
		double objVal = Double.NaN;
		double[] values = null;
		objectiveFun.setExpr(objectiveExpr);
		cplex.add(objectiveFun);
		try {
			cplex.solve();
			
			values = cplex.getValues(vars);
//			System.out.println(cplex.getStatus());
			if (cplex.getStatus() != IloCplex.Status.Optimal) {
//				System.out.println("Exception found: Status = "+cplex.getStatus());
				for (int i = 0; i < values.length; i++) {
					values[i] = Double.NaN;
				}
			} else {
				objVal = cplex.getObjValue();
			}
		} catch (Exception e) {
			System.out.println("Error!");
			e.printStackTrace();
			values = new double[vars.length];
			for (int i = 0; i < values.length; i++) {
				values[i] = Double.NaN;
			}
			return new SimulationResult(metaNet,values);

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
		return new SimulationResult(metaNet,values, objVal, cplex.getStatus().toString());
	}
	
	public Pair<double[],double[]> solvePivotFBA(FluxBound[] bounds, Reaction[] knockouts, Reaction objective, String objectiveSense, Reaction pivot, ReactionConstraint pivotRange, int breaks) throws IloException{
		FluxBound[] finalBounds = new FluxBound[bounds.length + 1];
		for (int i = 0; i < bounds.length; i++) {
			finalBounds[i] = bounds[i];
		}
		Pair<double[], double[]> res = new Pair<double[], double[]>(new double[breaks+1], new double[breaks+1]);
		double diff = pivotRange.getUpper() - pivotRange.getLower();
		for (int i = 0; i <= breaks; i++) {
			double fixed = ((diff*i)/breaks) + pivotRange.getLower();
			FluxBound pivotBound = new FluxBound(pivot, fixed, fixed);
			finalBounds[finalBounds.length-1] = pivotBound;
			SimulationResult curSol = solveReactionKnockoutFBA(finalBounds, knockouts, objective, objectiveSense);
			res.getA()[i] = fixed;
			res.getB()[i] = curSol.getFluxValue(objective);
		}
		return res;	
	}
	
	public SimulationResult solve(FluxBound[] bounds, Pair<Double,Reaction>[] objective, String sense) throws IloException {
		if (bounds != null) {
			for (int i = 0; i < bounds.length; i++) {
				FluxBound fbound = bounds[i];
				int index = metaNet.getReactionIndex(fbound.getReac().getName());
				vars[index].setLB(fbound.getBounds().getLower());
				vars[index].setUB(fbound.getBounds().getUpper());
			}
		}
		
		IloLinearNumExpr objectiveExpr = cplex.linearNumExpr();

		for (int i = 0; i < objective.length; i++) {
			Double objVal = objective[i].getA();
			Reaction objRx = objective[i].getB();
			IloNumVar objectiveVar = vars[metaNet.getReactionIndex(objRx.getName())];
			objectiveExpr.addTerm(objVal, objectiveVar);

		}
		
		IloObjective objectiveFun;
		int vsense = -1;

		if (sense.equals("max"))
			vsense = 1;
		else if (sense.equals("min"))
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
		double objValue = Double.NaN;
		try {
			cplex.solve();
			objValue = cplex.getObjValue();
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
			return new SimulationResult(metaNet,new double[vars.length]);
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
		return new SimulationResult(metaNet,values, objValue);
	}
	
	public double robustnessPoint(FluxBound[] bounds, Reaction[] knockouts, Reaction objective, Reaction pivot, double eps, double tol) throws IloException{
		
		double pivotMax = solveReactionKnockoutFBA(bounds, knockouts, pivot, "max").getFluxValue(pivot);
		double pivotDistance = (pivotMax)/2;
		double pivotPoint = pivotDistance;
		boolean found = false;
		
		FluxBound[] flexBound = new FluxBound[bounds.length+1];
		
		for (int i = 0; i < bounds.length; i++)
			flexBound[i] = bounds[i];
		
		while (!found) {
			
			flexBound[flexBound.length-1] = new FluxBound(pivot, pivotPoint, Utilities.INF);
			double pivotValue = solveReactionKnockoutFBA(flexBound, knockouts, pivot, "min").getFluxValue(objective);
			pivotDistance = pivotDistance/2;
//			System.out.println(pivotDistance/pivotMax+" "+pivotPoint/pivotMax);
			if (pivotDistance > tol*pivotMax) {
				if (pivotValue > eps) {
					pivotPoint -= pivotDistance;
				} else {
					pivotPoint += pivotDistance;
				}
			} else {
				found = true;
			}

		}
//		System.out.println(pivotPoint+"/"+pivotMax+"="+(pivotPoint/pivotMax));
		long a = Math.round((pivotPoint/pivotMax)*1e9);
		long b = Math.round(tol*1e9);
		long c = a/b;
//		System.out.println(a+"/"+b+"="+(c*tol));
		return c*tol;
	}
	
	public void close(){
		cplex.end();
	}
	
//	public void printMemoryUsage(){
//		cplex.us
//	}
	

}
