package pt.uminho.ceb.biosystems.mcslibrary.enumeration.implementation;

import ilog.concert.IloConstraint;
import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearIntExpr;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.concert.IloObjectiveSense;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;

import pt.uminho.ceb.biosystems.mcslibrary.enumeration.AbstractEnumerationSolver;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.EnumerationProblem;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.Solution;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.YieldConstraint;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;
import pt.uminho.ceb.biosystems.mew.utilities.java.TimeUtils;

/**
 * Subclass of {@link AbstractEnumerationSolver} that uses CPLEX to calculate
 * MCS for a given {@link EnumerationProblem}.
 *
 * @author V�tor
 *
 */
public class CPLEXEnumerationSolver extends AbstractEnumerationSolver {
	/**
	 *
	 * @param eprob
	 *            - The {@link EnumerationProblem} you want to calculate MCSs
	 *            for. Note that any following methods won't consider desired
	 *            fluxes or yields. See {@link CPLEXEnumerationFilter}.
	 */

	public CPLEXEnumerationSolver(EnumerationProblem eprob) {
		super(eprob);
		this.threads = 0;
		populate = true;
		allowSolverLogging = true;
	}

	private static final int M = 9999;
	private IloLinearNumExpr sizeexpression;
	private boolean populate;
	private int threads;
	private boolean allowSolverLogging;
	private BufferedOutputStream fos;

	// cplex parameters are set here
	public void setCplexParams(IloCplex cplex) throws IloException {
		cplex.setParam(IloCplex.IntParam.ClockType, 1);
		cplex.setParam(IloCplex.DoubleParam.WorkMem, 16000);
		cplex.setParam(IloCplex.DoubleParam.EpInt, 1e-10);
		System.out.println("CPLEX work memory: "+cplex.getParam(IloCplex.DoubleParam.WorkMem));
//		cplex.setParam(IloCplex.BooleanParam.NumericalEmphasis, true);
//		cplex.setParam(IloCplex.DoubleParam.EpRHS, 1e-9);
//		cplex.setParam(IloCplex.IntParam., arg1);
		cplex.setParam(IloCplex.IntParam.MIPEmphasis, 1);
//		cplex.setParam(IloCplex.IntParam.AdvInd, 1);
		cplex.setParam(IloCplex.IntParam.FPHeur, 1);
		cplex.setParam(IloCplex.IntParam.PopulateLim, 10000000);
		cplex.setParam(IloCplex.IntParam.SolnPoolReplace, 1);
		cplex.setParam(IloCplex.DoubleParam.SolnPoolAGap, 0);
		cplex.setParam(IloCplex.IntParam.SolnPoolCapacity, M);
		cplex.setParam(IloCplex.IntParam.SolnPoolIntensity, 4);
		cplex.setParam(IloCplex.IntParam.ReInv, 0);
//		cplex.setParam(DoubleParam.EpOpt, 1e-9);
	}


	/**
	 * @throws IOException
	 * @see pt.uminho.ceb.biosystems.mcslibrary.enumeration.AbstractEnumerationSolver#solve(int)
	 */
	
	public void setPopulate(boolean populate){
		this.populate = populate;
	}
	
	public void setThreads(int threads){
		this.threads = threads;
	}
	public DefaultEnumerationResult solve(int maxsize) throws IloException, IOException {


		IloCplex cplex = new IloCplex();
		setCplexParams(cplex);
		
		if (allowSolverLogging){
			File logFolder = new File("cplexLogs/");
			logFolder.mkdirs();
			fos = new BufferedOutputStream(new FileOutputStream("cplexLogs/cplexOut"+System.currentTimeMillis()+".log"));
			cplex.setOut(fos);
		} else {
			cplex.setOut(null);
		}

		cplex.setParam(IloCplex.IntParam.Threads, this.threads);
		
		// finding out useful quantities we'll need later
		int metabolites = getProblem().getMetabolicNetwork()
				.getNumOfMetabolites();
		int reactions = getProblem().getMetabolicNetwork().getNumOfReactions();
		int targets = getProblem().getUndesiredFluxes().length;
		int yields = getProblem().getUndesiredYieldConstraints().length;
		FluxBound[] extraCons = getProblem().getMetabolicNetwork().getExtraConstraints();
		int extra = extraCons.length;

		ArrayList<IloNumVar> variables = new ArrayList<IloNumVar>();
		ArrayList<IloRange> constraints = new ArrayList<IloRange>();

		double max = Double.MAX_VALUE;

		// creating variables:

			// U (transposed stoich matrix)
		for (int i = 0; i < metabolites; i++) {
			IloNumVar newVar = cplex.numVar(-max, max, "U" + i);
			variables.add(newVar);
		}
			// VP (forward fluxes)
		for (int i = 0; i < reactions; i++) {
			IloNumVar newVar = cplex.numVar(0, max, "VP" + i);
			variables.add(newVar);
		}
			// VN (reverse fluxes)
		for (int i = 0; i < reactions; i++) {
			IloNumVar newVar = cplex.numVar(0, max, "VN" + i);
			variables.add(newVar);
		}
			// T (lower bounds defined in the undesired fluxes)
		for (int i = 0; i < targets; i++) {
			IloNumVar newVar = cplex.numVar(0, max, "WL" + i);
			variables.add(newVar);
		}
			// T (upper bounds defined in the undesired fluxes)
		for (int i = 0; i < targets; i++) {
			IloNumVar newVar = cplex.numVar(0, max, "WU" + i);
			variables.add(newVar);
		}
			// T (subs/prod yield constraints in the undesired fluxes)
		for (int i = 0; i < yields; i++) {
			IloNumVar newVar = cplex.numVar(0, max, "WY" + i);
			variables.add(newVar);
		}
			// additional T (additional constraints other than +/- infinity)
		for (int i = 0; i < extra*2; i++) {
			IloNumVar newVar = cplex.numVar(0, max, "VEX" + i);
			variables.add(newVar);
		}
		IloNumVar c = cplex.numVar(1, max,"C");
		variables.add(c);

			// T-transposed b > -c
		IloLinearNumExpr targetexpression = cplex.linearNumExpr();

		// creating the transposed stoich matrix in the dual system
		// rows represented in U, VP and VN
		for (int i = 0; i < reactions; i++) {

			IloLinearNumExpr expression = cplex.linearNumExpr();
			for (int j = 0; j < metabolites; j++) {
				if (Math.abs(getProblem().getMetabolicNetwork()
						.getStoichCoef(j, i)) > 0) {
					expression.addTerm(getProblem().getMetabolicNetwork()
							.getStoichCoef(j, i), variables.get(j));
				}
			}

			expression.addTerm(1, variables.get(metabolites + i));
			expression.addTerm(-1, variables.get(metabolites + reactions + i));

			if (getProblem().getMetabolicNetwork().isReversible(i)) {
				constraints.add(cplex.eq(expression, 0));
			} else {
				constraints.add(cplex.ge(expression, 0));
			}
		}

		// this might be useful... eventually. it will store all fluxes in the undesired system
		ArrayList<Integer> involvedfluxes = new ArrayList<Integer>();

		// adding the undesired fluxes to the system (T matrix)
		for (int j = 0; j < getProblem().getUndesiredFluxes().length; j++) {
			FluxBound fbound = getProblem().getUndesiredFluxes()[j];
			int rindex = getProblem().getMetabolicNetwork().getReactionIndex(
					fbound.getReac().getName());
			if (rindex < 0) {
				System.out.println("WARNING!!"+ "Constraint for "+fbound.getReac().getName()+" could not be applied. Reaction does not exist.");
				continue;
			}
			involvedfluxes.add(rindex);
			IloLinearNumExpr expression = (IloLinearNumExpr) constraints.get(
					rindex).getExpr();
			/**
			 *if the undesired flux is off the range already defined in the model,
			 *it won't be added. if there is any specific limit to this flux
			 *it will be added in the additional constraints
			*/

			if (fbound.getBounds().getLower() > -M) {
				expression.addTerm(-1,
						variables.get(metabolites + reactions * 2 + j));
				targetexpression.addTerm(-fbound.getBounds().getLower(),
						variables.get(metabolites + reactions * 2 + j));
			}

			if (fbound.getBounds().getUpper() < M) {
				expression.addTerm(1, variables.get(metabolites + reactions * 2 + targets + j));
				targetexpression
						.addTerm(
								fbound.getBounds().getUpper(),
								variables.get(metabolites + reactions * 2 + targets + j));
			}

			constraints.get(rindex).setExpr(expression);

		}

		/** yield constraints are always added as lower bounds by default
		 */
		for (int j = 0; j < getProblem().getUndesiredYieldConstraints().length; j++) {
			YieldConstraint ycons = getProblem().getUndesiredYieldConstraints()[j];

			int rsindex = getProblem().getMetabolicNetwork().getReactionIndex(
					ycons.getUptakeReaction().getName());
			int rpindex = getProblem().getMetabolicNetwork().getReactionIndex(
					ycons.getProductReaction().getName());
//			involvedfluxes.add(rsindex);
//			involvedfluxes.add(rpindex);
			int factor = ycons.isLower() ? 1 : -1;

			IloLinearNumExpr rsexpression = (IloLinearNumExpr) constraints.get(
					rsindex).getExpr();
			rsexpression.addTerm(-ycons.getRatio() * factor, variables.get(metabolites
					+ reactions * 2 + targets * 2 + j));
			constraints.get(rsindex).setExpr(rsexpression);

			IloLinearNumExpr rpexpression = (IloLinearNumExpr) constraints.get(
					rpindex).getExpr();
			rpexpression.addTerm(1 * factor, variables.get(metabolites
					+ reactions * 2 + targets * 2 + j));
			constraints.get(rpindex).setExpr(rpexpression);

//			targetexpression.addTerm(0, variables.get(metabolites + reactions
//					* 2 + targets * 2 + j));

		}
		/** additional constraints that belong to the model itself are added here
		 * e.g. if the model already sets well defined limits to a given reaction
		 * other than infinity
		 */
		for (int j = 0; j < extraCons.length; j++) {
			FluxBound fbound = extraCons[j];
			int rindex = getProblem().getMetabolicNetwork().getReactionIndex(
					fbound.getReac().getName());
			if (rindex == -1 || involvedfluxes.contains(rindex)) {
				continue;
			}
			System.out.println("Adding extra constraint: "+fbound);
			involvedfluxes.add(rindex);
			IloLinearNumExpr expression = (IloLinearNumExpr) constraints.get(
					rindex).getExpr();

			if (fbound.getBounds().getLower() >= fbound.getReac().getBounds()
					.getLower()
					&& fbound.getBounds().getLower() > -M) {
				expression.addTerm(-1,
						variables.get(metabolites + reactions * 2 + targets * 2 + yields + j*2));
				targetexpression.addTerm(-fbound.getBounds().getLower(),
						variables.get(metabolites + reactions * 2 + targets * 2 + yields + j*2));
			}

			if (fbound.getBounds().getUpper() <= fbound.getReac().getBounds()
					.getUpper()
					&& fbound.getBounds().getUpper() < M) {
				
				expression.addTerm(1, 
						variables.get(metabolites + reactions * 2 + targets * 2 + yields + (j*2 + 1)));
				
				targetexpression
						.addTerm(
								fbound.getBounds().getUpper(),
								variables.get(metabolites + reactions * 2 + targets * 2 + yields + (j*2 + 1)));
			}
			constraints.get(rindex).setExpr(expression);
		}

		ArrayList<IloIntVar> ivariables = new ArrayList<IloIntVar>();
		ArrayList<IloConstraint> iconstraints = new ArrayList<IloConstraint>();

		// TARGET EXPRESSION IS ADDED HERE. DON'T CHANGE THE TARGET PAST THIS POINT
		targetexpression.addTerm(1, c);
		IloConstraint target = cplex.eq(targetexpression, 0);
		cplex.add(target);

		/** Indicator variables are added here. For each indicator two constraints are
		 *	added, one for the i = 0 case and another one for i = 1. If a VP or VN variable
		 *	is greater than 1, i = 1, otherwise it is 0.
		 */
		for (int i = 0; i < reactions * 2; i++) {
			IloIntVar indicatorVar = cplex.boolVar("I"
					+ variables.get(metabolites + i));
			IloConstraint condition0 = cplex.eq(variables.get(metabolites + i),
					0);
			IloConstraint definition0 = cplex.eq(indicatorVar, 0);
			IloConstraint condition1 = cplex.ge(
					variables.get(metabolites + i), 1);
			IloConstraint definition1 = cplex.eq(indicatorVar, 1);
			IloConstraint indicatorConstraint0 = cplex.ifThen(definition0,condition0);
			IloConstraint indicatorConstraint1 = cplex.ifThen(definition1,condition1);


			ivariables.add(indicatorVar);
			iconstraints.add(indicatorConstraint0);
			iconstraints.add(indicatorConstraint1);

		}

		// Adds constraints so that either VP or VN can be active
		for (int i = 0; i < reactions; i++) {
			IloLinearIntExpr exp = cplex.linearIntExpr();
			exp.addTerm(1, ivariables.get(i));
			exp.addTerm(1, ivariables.get(i + reactions));
			IloRange ran = cplex.le(exp, 1);
			cplex.add(ran);
		}

		// Adds a constraint so that certain reactions are not accounted for in the enumeration
		IloLinearNumExpr expi = cplex.linearNumExpr();
		for (Reaction reaction : getProblem().getExcludedReactions()) {
			int idx = getProblem().getMetabolicNetwork().getReactionIndex(
					reaction.getName());
			if (idx > -1) {
				expi.addTerm(1, ivariables.get(idx));
				expi.addTerm(1, ivariables.get(idx + reactions));
			}
		}
		IloRange rani = cplex.eq(0, expi);
		cplex.add(rani);
		
		// Adds constraints that do not allow MCSs to be subsets of a given solution pool
		for (Solution sol : getProblem().getExcludedSolutions()) {
			IloLinearIntExpr exp = cplex.linearIntExpr();
			for (int i = 0; i < sol.getSize(); i++) {
				int index = getProblem().getMetabolicNetwork().getReactionIndex(sol.getReactionIdFromIndex(i));
				if (index > -1) {
					exp.addTerm(ivariables.get(index), 1);
				}
			}
			cplex.add(cplex.le(exp, sol.getSize()));
		}

		IloNumVar[] vars = new IloNumVar[variables.size()];
		IloIntVar[] inds = new IloIntVar[ivariables.size()];

		for (IloNumVar var : variables) {
			cplex.add(var);
		}
		for (IloRange ran : constraints)
			cplex.add(ran);

		for (IloIntVar ind : ivariables) {
			cplex.add(ind);
		}
		for (IloConstraint icns : iconstraints)
			cplex.add(icns);

		for (int i = 0; i < vars.length; i++) {
			vars[i] = variables.get(i);
		} 

		for (int i = 0; i < inds.length; i++) {
			inds[i] = ivariables.get(i);
		}
		
		sizeexpression = cplex.linearNumExpr();
		for (IloIntVar var : ivariables)
			sizeexpression.addTerm(1, var);


		ArrayList<Double> times = new ArrayList<Double>();

		ArrayList<int[]> res = new ArrayList<int[]>();
		

//		cplex.exportModel("modelold.lp");
//		boolean solved = cplex.solve();
//		System.out.println("Problem can be solved? "+solved);
//		double[] arrtest = cplex.getValues(inds);
//		
//		for (int j = 0; j < arrtest.length; j++) {
//			if (arrtest[j] != 0) {
//				System.out.println("\t Prime solution index "+j+" = "+arrtest[j]);
//			}
//		}
//		
		
//		cplex.add(cplex.le(sizeexpression, maxsize));
		if (populate) {
			for (int i = 1; i <= maxsize; i++) {
				System.out.println("Starting size " + i + ".");
				long starttime = System.currentTimeMillis();

				

				IloRange sizeconstraint = cplex.eq(sizeexpression, i);
				cplex.add(sizeconstraint);
				cplex.populate();
				int solutions = cplex.getSolnPoolNsolns();
				System.out.println("Size " + i + " - Solutions: " + solutions);
				
				
				ArrayList<IloRange> cuts = new ArrayList<IloRange>();
				for (int j = 0; j < solutions; j++) {
					
					double[] values = cplex.getValues(inds, j);
					int[] mcsarray = recoverArray(values);
					res.add(recoverSolution(mcsarray, i));
					IloLinearIntExpr mcs = cplex.linearIntExpr();
					for (int k = 0; k < mcsarray.length; k++) {
						mcs.addTerm(mcsarray[k], inds[k]);
						mcs.addTerm(mcsarray[k], inds[mcsarray.length + k]);
					}

					IloRange cut = cplex.le(mcs, arraySum(mcsarray) - 1);
					cuts.add(cut);
				}
				for (IloRange cut : cuts)
					cplex.add(cut);

				cplex.remove(sizeconstraint);

				long spent = System.currentTimeMillis() - starttime;
				times.add((double) spent);
				System.out.println("Size "+i+" finished in "+TimeUtils.formatMillis(spent));

			}	
		} else {
			boolean exceptionRaised = false;
			int nCuts = 0;
			IloObjective ob = cplex.objective(IloObjectiveSense.Minimize, sizeexpression);
			cplex.add(ob);
			IloRange sizeconstraint = cplex.ge(sizeexpression, 1);
			cplex.add(sizeconstraint);
			do {
				try {
					cplex.solve();
					double[] values = cplex.getValues(inds);
					int[] mcsarray = recoverArray(values);
					res.add(recoverSolution(mcsarray));
					IloLinearIntExpr mcs = cplex.linearIntExpr();
					for (int k = 0; k < mcsarray.length; k++) {
						mcs.addTerm(mcsarray[k], inds[k]);
						mcs.addTerm(mcsarray[k], inds[mcsarray.length + k]);
					}
					IloRange cut = cplex.le(mcs, arraySum(mcsarray) - 1);
					cplex.add(cut);
					nCuts++;

				} catch (Exception e) {
					exceptionRaised = true;
					System.out.println("Last solver status: "+cplex.getStatus());
					System.out.println("Exception: "+e.getMessage());
				}

			} while (!exceptionRaised && nCuts < maxsize);
			
		}
		
		

		if (allowSolverLogging) {
			fos.flush();
			fos.close();	
		}
		
		double timesum = 0;
		for (Double time : times)
			timesum = timesum + time;
		System.out.println(timesum);
		DefaultEnumerationResult enumresult =
				getProblem().isCompressedProblem() ?
						new CompressedEnumerationResults(getProblem(), res).decompressResult()
					  : new DefaultEnumerationResult(getProblem(), res);
		return enumresult;

	}

	private int arraySum(int[] array) {
		int res = 0;
		for (int i = 0; i < array.length; i++) {
			res = res + array[i];
		}
		return res;
	}

	private int[] recoverArray(double[] values) {
		int halfway = values.length / 2;
		int[] res = new int[halfway];
		for (int i = 0; i < halfway; i++) {
			res[i] = (int) (values[i] + values[halfway + i]);
		}
		return res;
	}

	private int[] recoverSolution(int[] values, int size) {
		int[] fsol = new int[size];
		int count = 0;
		for (int i = 0; i < values.length; i++) {
			if (values[i] != 0) {
				fsol[count] = i;
				count++;
			}
		}
		return fsol;
	}
	
	private int[] recoverSolution(int[] values){
		ArrayList<Integer> inter = new ArrayList<Integer>();
		for (int i = 0; i < values.length; i++) {
			if (values[i] > 0) {
				inter.add(i);
			}
		}
		int[] res = new int[inter.size()];
		for (int i = 0; i < res.length; i++) {
			res[i] = inter.get(i);
		}
		return res;
	}
	
	// debugging method
	private static String sol2String(IloIntVar[] var, double[] values, String solutionName) {
		String res = solutionName;
		for (int i = 0; i < values.length; i++) {
			double val = values[i];
			String varName = var[i].getName();
			if (val != 0) {
				res += "\n";
				res += varName+" = "+val+";";
			}
		}
		return res;
	}
	
	private static String sol2String(IloNumVar[] var, double[] values, String solutionName) {
		String res = solutionName;
		for (int i = 0; i < values.length; i++) {
			double val = values[i];
			String varName = var[i].getName();
			if (val != 0) {
				res += "\n";
				res += varName+" = "+val+";";
			}
		}
		return res;
	}


	public void allowSolverLogging(boolean allowSolverLogging) {
		this.allowSolverLogging = allowSolverLogging;
	}
}
