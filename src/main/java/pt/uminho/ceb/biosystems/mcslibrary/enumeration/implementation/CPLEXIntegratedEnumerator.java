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
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import pt.uminho.ceb.biosystems.mcslibrary.enumeration.AbstractEnumerationResult;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.AbstractEnumerationSolver;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.EnumerationProblem;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.Solution;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.cplex.callbacks.DefaultCallback;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.cplex.callbacks.IntermediateSolverCallback;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.logging.EnumeratorLog;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.AbstractMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.compression.alg.MatrixTools;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.YieldConstraint;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fva.CPLEXFluxVariabilityAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fva.FluxVariabilityAnalysisResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.CompressedMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Pair;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;

public class CPLEXIntegratedEnumerator extends AbstractEnumerationSolver {
	private IloCplex cplex;
	private static double inf = Utilities.INF;
	private static final int M = 1000000;
	private int currentSize;
	private AbstractEnumerationResult result;
	private IloLinearNumExpr objective;
	private EnumeratorLog log;
	public IloIntVar[] ZP;
	public IloIntVar[] ZN;
	public IloNumVar[] VP;
	public IloNumVar[] VN;
	public IloNumVar[] We;
	public IloNumVar[] Wf;
	public IloNumVar[] Wy;
	private BufferedOutputStream fos;
	private IloObjective objectiveFunction;
	
	
	public CPLEXIntegratedEnumerator(EnumerationProblem eprob) throws IloException, IOException {
		super(eprob);
		initialize();
		setCplexParams(cplex);
		setCurrentSize(0+eprob.getForcedReactions().length);
		if (eprob.isCompressedProblem()) {
			result = new CompressedEnumerationResults(eprob);
		} else { 
			result = new DefaultEnumerationResult(eprob);
		}
		log = new EnumeratorLog(this);
	}
	
	public void setCplexParams(IloCplex cplex) throws IloException {
		cplex.setParam(IloCplex.IntParam.ClockType, 1);
		cplex.setParam(IloCplex.DoubleParam.WorkMem, 3096);
		cplex.setParam(IloCplex.DoubleParam.EpInt, 1e-10);
//		cplex.setParam(IloCplex.BooleanParam.NumericalEmphasis, true);
//		cplex.setParam(IloCplex.DoubleParam.EpRHS, 1e-9);
		cplex.setParam(IloCplex.IntParam.MIPEmphasis, 2);
		cplex.setParam(IloCplex.IntParam.Probe, 3);
		cplex.setParam(IloCplex.IntParam.PopulateLim, M);
		cplex.setParam(IloCplex.IntParam.SolnPoolReplace, 1);
		cplex.setParam(IloCplex.DoubleParam.SolnPoolAGap, 0);
		cplex.setParam(IloCplex.IntParam.SolnPoolCapacity, M);
		cplex.setParam(IloCplex.IntParam.SolnPoolIntensity, 4);
//		cplex.setParam(IloCplex.IntParam.MIPDisplay, 5);
//		cplex.setParam(DoubleParam.EpOpt, 1e-9);
//		cplex.tuneParam();
//		int[] c = new int[ZP.length];
//		Arrays.fill(c, 3);
//		cplex.setPriorities(ZP,c);
//		cplex.setPriorities(ZN,c);
		try {
			fos = new BufferedOutputStream(new FileOutputStream("cplexOut.log"));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		cplex.setOut(fos);	
	} 

	public void calculateNextSize() throws IloException, IOException{
		fos = new BufferedOutputStream(new FileOutputStream("cplexOut.log"));
		cplex.setOut(fos);	
		currentSize += 1;
		IloRange sizeConstraint = cplex.eq(objective, currentSize);
		cplex.add(sizeConstraint);
		cplex.populate();
		int numSols = cplex.getSolnPoolNsolns();
		ArrayList<IloRange> cuts = new ArrayList<IloRange>();
		for (int i = 0; i < numSols; i++) {
			double[] arrayPos = cplex.getValues(ZP,i);
			double[] arrayNeg = cplex.getValues(ZN,i);
			int[] sol = recoverArray(arrayPos, arrayNeg);
			result.addSolution(sol);
			cuts.add(getIntegerCut(sol));
		}
		cplex.remove(sizeConstraint);
		
		for (IloRange iloRange : cuts) {
			cplex.add(iloRange);
		}
		fos.flush();
		fos.close();
		
	}
	public int[] findNextSolution() throws IloException, IOException {

		cplex.solve();
		double[] arrayPos = cplex.getValues(ZP);
		double[] arrayNeg = cplex.getValues(ZN);
		int[] sol = recoverArray(arrayPos, arrayNeg);
		System.out.println(Arrays.toString(sol));
//		result.addSolution(sol);
		cplex.add(getIntegerCut(sol));
		fos.flush();
		fos.close();
		return sol;
	}
	
	public void clearModel() throws IloException{
		cplex.clearModel();
	}
	public void calculateNextSolution() throws IloException, IOException {

		cplex.solve();
		double[] arrayPos = cplex.getValues(ZP);
		double[] arrayNeg = cplex.getValues(ZN);
		int[] sol = recoverArray(arrayPos, arrayNeg);
		System.out.println(Arrays.toString(sol));
		result.addSolution(sol);
		cplex.add(getIntegerCut(sol));
		fos.flush();
		fos.close();
	}
	
	private void initialize() throws IloException, IOException{
		cplex = new IloCplex();
		// useful constants
		AbstractMetabolicNetwork metaNet = getProblem().getMetabolicNetwork();
		int M = metaNet.getNumOfMetabolites();
		int N = metaNet.getNumOfReactions();
		
		// creating variables for enumerator
		IloNumVar[] U = createNumVarArray(-inf, inf, M, cplex, "U");
		VP = createNumVarArray(0, inf, N, cplex, "VP");
		VN = createNumVarArray(0, inf, N, cplex, "VN");
		ZP = createBinVarArray(N, cplex, "ZP");
		ZN = createBinVarArray(N, cplex, "ZN");
		IloNumVar c = cplex.numVar(1, inf,"C");
		
		// Dual reaction constraints
		IloRange[] dualc = new IloRange[N];
		System.out.println(metaNet);
		for (int i = 0; i < N; i++) {
			IloLinearNumExpr numexp = cplex.linearNumExpr();
			for (int j = 0; j < M; j++) {
				numexp.addTerm(metaNet.getStoichCoef(j, i), U[j]);
			}
			numexp.addTerm(1, VP[i]);
			numexp.addTerm(-1, VN[i]);
			boolean isRev = metaNet.isReversible(i);
			dualc[i] = isRev ? cplex.eq(numexp, 0) : cplex.ge(numexp, 0);
		}
		
		// ArrayList to store involved reaction ids
		Set<Integer> involved = new HashSet<Integer>();

		// Mapping target modes
		// MAP<PAIR<SENSE,REACTIONID>, BOUNDVALUE>
		HashMap<Pair<Integer,Integer>,Double> Ttemp = new HashMap<Pair<Integer,Integer>, Double>();
		
		FluxBound[] uflxs = getProblem().getUndesiredFluxes();
		for (int i = 0; i < uflxs.length; i++) {
			FluxBound flx = uflxs[i];
			int idx = metaNet.getReactionIndex(flx.getReac().getName());
			if (idx > -1) {
				if (flx.getBounds().getLower() > -inf) {
					Ttemp.put(new Pair<Integer, Integer>(-1, idx), flx.getBounds().getLower());
					involved.add(idx);
				}
				if (flx.getBounds().getUpper() < inf) {
					Ttemp.put(new Pair<Integer, Integer>(1, idx), flx.getBounds().getUpper());
					involved.add(idx);
				}
			} else {
				continue;
			}
		}
		
		System.out.println("[ENUMERATOR] "+Ttemp.size()+" valid target modes.");
		System.out.println(Ttemp);
		
		// Mapping target yields
		// MAP<PAIR<PRODUCTID,SUBSTRATEID>, PAIR<PAIR<ISLOWER, RATIO>,DEVIATION>>
		// A yield constraint is something like: Ratio*Product - Substrate =/</> Deviation
		HashMap<Pair<Integer,Integer>,Pair<Pair<Boolean,Double>, Double>> UYtemp = new HashMap<Pair<Integer,Integer>, Pair<Pair<Boolean,Double>, Double>>();
		YieldConstraint[] uylds = getProblem().getUndesiredYieldConstraints();
		for (int i = 0; i < uylds.length; i++) {
			YieldConstraint yld = uylds[i];
			int pidx = metaNet.getReactionIndex(yld.getProductReaction().getName());
			int sidx = metaNet.getReactionIndex(yld.getUptakeReaction().getName());
			if (pidx > -1 && sidx > -1) {
				UYtemp.put(new Pair<Integer, Integer>(pidx, sidx), new Pair<Pair<Boolean,Double>,Double>(new Pair<Boolean, Double>(yld.isLower(), yld.getRatio()), yld.getDeviation()));
			}
		}
		
		System.out.println("[ENUMERATOR] "+UYtemp.size()+" valid target yield constraints.");
		
		// Mapping extra constraints
		// MAP<PAIR<SENSE,REACTIONID>, BOUNDVALUE>
		System.out.println("Reaction IDs involved in yield and rate restrictions: "+involved);
		HashMap<Pair<Integer,Integer>,Double> EXtemp = new HashMap<Pair<Integer,Integer>, Double>();
		FluxBound[] extraflxs = metaNet.getExtraConstraints();
		for (int i = 0; i < extraflxs.length; i++) {
			FluxBound flx = extraflxs[i];
			int idx = metaNet.getReactionIndex(flx.getReac().getName());
			if (idx > -1 && !involved.contains(idx)) {
				if (flx.getBounds().getLower() > -inf) {
					EXtemp.put(new Pair<Integer, Integer>(-1, idx), flx.getBounds().getLower());
				}
				if (flx.getBounds().getUpper() < inf) {
					EXtemp.put(new Pair<Integer, Integer>(1, idx), flx.getBounds().getUpper());
				}
			} else {
				continue;
			}
		}
		System.out.println("[ENUMERATOR] "+EXtemp.size()+" valid extra target modes.");

		/** Indicator variables are added here. For each indicator two constraints are
		 *	added, one for the i = 0 case and another one for i = 1. If a VP or VN variable
		 *	is greater than 1, i = 1, otherwise it is 0.
		 */
		
		IloConstraint[] indPZero = new IloConstraint[N];
		IloConstraint[] indNZero = new IloConstraint[N];
		IloConstraint[] indPOne = new IloConstraint[N];
		IloConstraint[] indNOne = new IloConstraint[N];

		for (int i = 0; i < N; i++) {
			IloConstraint condition0 = cplex.eq(VP[i],0);
			IloConstraint definition0 = cplex.eq(ZP[i], 0);
			IloConstraint condition1 = cplex.ge(VP[i], 1);
			IloConstraint definition1 = cplex.eq(ZP[i], 1);
			
			IloConstraint indicatorConstraint0 = cplex.ifThen(definition0,condition0);
			IloConstraint indicatorConstraint1 = cplex.ifThen(definition1,condition1);

			indPZero[i] = indicatorConstraint0;
			indPOne[i] = indicatorConstraint1;
		}
		
		for (int i = 0; i < N; i++) {
			IloConstraint condition0 = cplex.eq(VN[i],0);
			IloConstraint definition0 = cplex.eq(ZN[i], 0);
			IloConstraint condition1 = cplex.ge(VN[i], 1);
			IloConstraint definition1 = cplex.eq(ZN[i], 1);
			
			IloConstraint indicatorConstraint0 = cplex.ifThen(definition0,condition0);
			IloConstraint indicatorConstraint1 = cplex.ifThen(definition1,condition1);

			indNZero[i] = indicatorConstraint0;
			indNOne[i] = indicatorConstraint1;
		}
		
		// Reaction sense restrictions
		IloRange[] senseCons = new IloRange[N];
		for (int i = 0; i < N; i++) {
			IloLinearIntExpr exp = cplex.linearIntExpr();
			exp.addTerm(1, ZP[i]);
			exp.addTerm(1, ZN[i]);
			IloRange ran = cplex.le(exp, 1);
			senseCons[i] = ran;
		}
		
		// Reaction exclusions
		Reaction[] excludedReactions = getProblem().getExcludedReactions();
		ArrayList<IloRange> tempExclusionCons = new ArrayList<IloRange>();
		for (int i = 0; i < excludedReactions.length; i++) {
			int idx = metaNet.getReactionIndex(excludedReactions[i].getName());
			if (idx > -1) {
				IloLinearNumExpr exp = cplex.linearNumExpr();
				exp.addTerm(1, ZP[idx]);
				exp.addTerm(1, ZN[idx]);
				tempExclusionCons.add(cplex.le(exp, 0));
			}
		}
		
		// Reactions forced into an MCS (experimental!!)
		ArrayList<Integer> indsx = new ArrayList<Integer>();
		Reaction[] forcedReactions = getProblem().getForcedReactions();
//		ArrayList<IloRange> tempForcingCons = new ArrayList<IloRange>();
		IloLinearNumExpr tfexp = cplex.linearNumExpr();
		for (int i = 0; i < forcedReactions.length; i++) {
			int idx = metaNet.getReactionIndex(forcedReactions[i].getName());
			if (idx > -1 && !indsx.contains(idx)) {
				tfexp.addTerm(1, ZP[idx]);
				tfexp.addTerm(1, ZN[idx]);
				indsx.add(idx);
			}
		}
		
		IloRange reactForcingCons = cplex.eq(tfexp, indsx.size());
		
		// Solution subset exclusions
		Solution[] excludedSolutions = getProblem().getExcludedSolutions();
		ArrayList<IloRange> tempSolExclusionCons = new ArrayList<IloRange>();
		for (Solution sol : excludedSolutions) {
			IloLinearIntExpr exp = cplex.linearIntExpr();
			for (int i = 0; i < sol.getSize(); i++) {
				int index = getProblem().getMetabolicNetwork().getReactionIndex(sol.getReactionIdFromIndex(i));
				if (index > -1) {
					exp.addTerm(ZP[i], 1);
					exp.addTerm(ZN[i], 1);
				}
			}
			tempSolExclusionCons.add(cplex.le(exp, sol.getSize()-1));
		}
		IloRange[] solExclusionCons = tempSolExclusionCons.toArray(new IloRange[tempSolExclusionCons.size()]);
		
		
		// Adding target modes as cplex objects
		
		Wf = createNumVarArray(0, inf, Ttemp.size(), cplex, "Wf");
		Wy = createNumVarArray(0, inf, UYtemp.size(), cplex, "Wy");
		We = createNumVarArray(0, inf, EXtemp.size(), cplex, "We");
		
		IloLinearNumExpr Tw = cplex.linearNumExpr();
		
		int i = 0;
		for (Pair<Integer, Integer> entry : Ttemp.keySet()) {
			Integer sense = entry.getA();
			Integer rId = entry.getB();
			Double b = Ttemp.get(entry);
			
			IloLinearNumExpr exp = (IloLinearNumExpr) dualc[rId].getExpr();
			exp.addTerm(sense, Wf[i]);
			dualc[rId].setExpr(exp);
			Tw.addTerm(b*sense, Wf[i]);
			involved.add(rId);
			i++;
		}
		
		i = 0;
		for (Pair<Integer, Integer> entry : UYtemp.keySet()) {
			Integer pId = entry.getA();
			Integer sId = entry.getB();
			Boolean isLower = UYtemp.get(entry).getA().getA();
			int factor = isLower ? 1 : -1;
			Double ratio = UYtemp.get(entry).getA().getB();
			Double dev = UYtemp.get(entry).getB();
			
			IloLinearNumExpr pexp = (IloLinearNumExpr) dualc[pId].getExpr();
			IloLinearNumExpr sexp = (IloLinearNumExpr) dualc[sId].getExpr();
			pexp.addTerm(factor, Wy[i]);
			sexp.addTerm(factor*-ratio, Wy[i]);
			dualc[pId].setExpr(pexp);
			dualc[sId].setExpr(sexp);
			Tw.addTerm(dev, Wy[i]);
			i++;
		}
		
		i = 0;
		for (Pair<Integer, Integer> entry : EXtemp.keySet()) {
			Integer sense = entry.getA();
			Integer rId = entry.getB();
			Double b = EXtemp.get(entry);
			
			IloLinearNumExpr exp = (IloLinearNumExpr) dualc[rId].getExpr();
			exp.addTerm(sense, We[i]);
			dualc[rId].setExpr(exp);
			Tw.addTerm(b*sense, We[i]);
			i++;
		}
		
		objective = cplex.linearNumExpr();
		
		for (int j = 0; j < N; j++) {
			objective.addTerm(1,ZP[j]);
			objective.addTerm(1,ZN[j]);
		}
		

		objectiveFunction = cplex.objective(IloObjectiveSense.Minimize,objective);
		System.out.println("Objective function created: "+objectiveFunction);
		Tw.addTerm(1, c);
		IloRange targetConstraint = cplex.eq(Tw, 0);
		
		cplex.add(U);
		cplex.add(VP);
		cplex.add(VN);
		cplex.add(We);
		cplex.add(Wf);
		cplex.add(Wy);
		cplex.add(ZP);
		cplex.add(ZN);

		
		cplex.add(c);
		cplex.add(dualc);
		cplex.add(targetConstraint);
		cplex.add(senseCons);
		cplex.add(solExclusionCons);

		

		// experimental
//		cplex.add(cplex.eq(7, objective));
		cplex.add(indNOne);
		cplex.add(indNZero);
		cplex.add(indPOne);
		cplex.add(indPZero);

		cplex.add(objectiveFunction);
		// desired modes
		// Mahadevan 2015
		// setting up the N matrix 
		IloNumVar[] Nr = new IloNumVar[N];
		IloConstraint[] Nm = new IloConstraint[M];
		
		int amax = M;
		
		System.out.println("Calculating flux limits for compressed network...");
		System.out.println("Done");
		
		for (int id = 0; id < Nr.length; id++) {
			Nr[id] = cplex.numVar(metaNet.getLowerBound(id),metaNet.getUpperBound(id),"R"+id);
		}
		
		for (int id = 0; id < Nm.length; id++) {
			IloLinearNumExpr ex = cplex.linearNumExpr();
			for (int j = 0; j < Nr.length; j++) {
				ex.addTerm(metaNet.getStoichCoef(id, j), Nr[j]);
			}
			Nm[id] = cplex.eq(ex, 0);
		}
		
		// debug
		
		// adding the D matrix
		List<IloConstraint> f = new ArrayList<IloConstraint>();
		FluxBound[] dfs = getProblem().getDesiredFluxes();
		for (int it = 0; it < dfs.length; it++) {
			System.out.println("Applying desired constraint "+it);
			FluxBound fb = dfs[it];
			System.out.println("\t"+fb.toString());
			int id = metaNet.getReactionIndex(fb.getReac().getName());
			if (id > -1) {
				double lb = fb.getBounds().getLower();
				double ub = fb.getBounds().getUpper();
				IloLinearNumExpr Lex = cplex.linearNumExpr();
				IloLinearNumExpr Uex = cplex.linearNumExpr();
				Uex.addTerm(1, Nr[id]);
				Lex.addTerm(-1, Nr[id]);
				boolean overLB = lb <= -amax;
				boolean overUB = ub >= amax;

				if (!overUB) {
					f.add(cplex.le(Uex, ub));
					System.out.println(f.get(f.size()-1));
				}
				if (!overLB) {
					f.add(cplex.le(Lex, -lb));
					System.out.println(f.get(f.size()-1));
				}
			}
		}
		YieldConstraint[] dys = getProblem().getDesiredYieldConstraints();
		for (int it = 0; it < dys.length; it++) {
			YieldConstraint dy = dys[it];
			int pId = getProblem().getMetabolicNetwork().getReactionIndex(dy.getProductReaction().getName());
			int uId = getProblem().getMetabolicNetwork().getReactionIndex(dy.getUptakeReaction().getName());
			double yield = dy.getRatio();
			IloLinearNumExpr linexp = cplex.linearNumExpr();
			int k = dy.isLower() ? 1 : -1;
			linexp.addTerm(-k, Nr[pId]);
			linexp.addTerm(k*yield, Nr[uId]);
			f.add(cplex.le(linexp, -k*dy.getDeviation()));
		}
		
		IloConstraint[] D = f.toArray(new IloConstraint[f.size()]);
		
		// adding the required constraints for the integrated formulation
		
		IloConstraint[] ubc = new IloConstraint[Nr.length];
		IloConstraint[] lbc = new IloConstraint[Nr.length];
		
		CPLEXFluxVariabilityAnalysis fvaN = new CPLEXFluxVariabilityAnalysis(((CompressedMetabolicNetwork)metaNet).getParentNetwork(), dfs, dys);
		CPLEXFluxVariabilityAnalysis fvaC = new CPLEXFluxVariabilityAnalysis(((CompressedMetabolicNetwork)metaNet), dfs, dys);
//		
		
		FluxVariabilityAnalysisResult fvaResC = fvaC.solveFVA();
		fvaResC.printResults((DefaultMetabolicNetwork)((CompressedMetabolicNetwork)metaNet).getParentNetwork());
		double[] fvalb = fvaResC.minToArray();
		double[] fvaub = fvaResC.maxToArray();
		
		double[][] bounds = validateFVABounds(((CompressedMetabolicNetwork)metaNet).getSubMatrix(), fvaN.solveFVA(), fvaResC);
		double[] minValidate = bounds[0];
		double[] maxValidate = bounds[1];
		
		for (int j = 0; j < maxValidate.length; j++) {
			double devlb = minValidate[j];
			double devub = maxValidate[j];
			if (devlb > 1e-5 || devub > 1e-5) {
				System.out.println("Reaction "+j+" bounds are unstable.");
			}
		}
		System.out.println("Max FVAMIN deviation: "+minValidate[MatrixTools.getMaximumIdx(minValidate)]);
		System.out.println("Max FVAMAX deviation: "+maxValidate[MatrixTools.getMaximumIdx(maxValidate)]);
		
		for (int i1 = 0; i1 < Nr.length; i1++) {
			if (fvaResC.isUnbounded(i1)) {
				System.out.println("Reaction "+i1+" is unbounded.");
			}
			if (fvaResC.isEssential(i1)) {
				IloLinearNumExpr exp = cplex.linearNumExpr();
				exp.addTerm(1, ZP[i1]);
				exp.addTerm(1, ZN[i1]);
				tempExclusionCons.add(cplex.eq(exp, 0));
				System.out.println("Reaction "+i1+" is essential.");
			}
			IloLinearNumExpr lexp = cplex.linearNumExpr();
			IloLinearNumExpr uexp = cplex.linearNumExpr();
			
			double lb = fvalb[i1];
			double ub = fvaub[i1];
			
//			double lb = metaNet.getLowerBound(i1);
//			double ub = metaNet.getUpperBound(i1);

			lb = lb < -amax ? -amax : lb;
			ub = ub > amax ? amax : ub;

			if ((lb > 0 || ub < 0) && metaNet.isReversible(i1)) {
				System.out.println("Reversibility mismatch!"+lb+" < Reaction"+i1+" < "+ub);
			}
			lb = Math.abs(lb) > Utilities.PRECISION ? lb : 0;
			ub = Math.abs(ub) > Utilities.PRECISION ? ub : 0;

			lexp.addTerm(1, Nr[i1]);
			lexp.addTerm(lb, ZP[i1]);
			lexp.addTerm(lb, ZN[i1]);
			
			lbc[i1] = cplex.ge(lexp, lb);
			
			uexp.addTerm(1, Nr[i1]);
			uexp.addTerm(ub, ZP[i1]);
			uexp.addTerm(ub, ZN[i1]);
			
			ubc[i1] = cplex.le(uexp, ub);
		}
		
		
		if (D.length > 0) {
			cplex.add(Nr);
			cplex.add(Nm);
			cplex.add(D);
			cplex.add(lbc);
			cplex.add(ubc);
		}
		
		IloRange[] exclusionCons = tempExclusionCons.toArray(new IloRange[tempExclusionCons.size()]);
		cplex.add(exclusionCons);
		cplex.add(reactForcingCons);
		
		cplex.exportModel("model.lp");	
	}
	
//	public int[] solveIntermediateSolution(double stopTime, boolean persist){
//		IntermediateSolverCallback c = new IntermediateSolverCallback(stopTime, persist, ZP, ZN);
//		try {
//			cplex.clearCallbacks();
//			cplex.use(c);
//		} catch (IloException e) {
//			System.out.println("Exception during callback setup.");
//		}
//		try {
//			cplex.solve();
//		} catch (IloException e) {
//			System.out.println("Exception during solving.");
//		}
//		return c.getFinalSol();
//	}
	
	public DefaultEnumerationResult iterativeSolve(int stopTime, int totalTime, boolean persist, int iterMaxPerSol) throws IloException, IOException{
		AbstractEnumerationResult r = null;
		boolean isCompressed = false;
		if (this.getProblem().getMetabolicNetwork().getClass() == DefaultMetabolicNetwork.class) {
			r = new DefaultEnumerationResult(getProblem());
		} else if (this.getProblem().getMetabolicNetwork().getClass() == CompressedMetabolicNetwork.class){
			r = new CompressedEnumerationResults(getProblem());
			isCompressed = true;
		}
		int solCount = 0;
		IloRange[] tempConstraints;
		long startTime = System.currentTimeMillis();
		while ((((System.currentTimeMillis() - startTime) / 1000) < totalTime)) {
			solCount++;
			System.out.println("---------------------------Solution "+solCount+"---------------------------");
			// explore step
			System.out.println("---------------------------Phase 1: Probing");
			IntermediateSolverCallback c = new IntermediateSolverCallback(stopTime, persist, ZP, ZN);
			cplex.clearCallbacks();
			cplex.use(c);
			int[] curSol = null;
			System.out.println("\tSolver started.");
			cplex.solve();
			System.out.println("Probing complete");
			curSol = c.getFinalSol();
			// simplify step
			System.out.println("---------------------------Phase 2: Simplifying");
			int curIter = 0;
			if (curSol != null) {
				do {
					System.out.println("Iteration "+(1+curIter));
					System.out.println();
					int[] tempSol = c.getFinalSol();
//					if (tempSol == curSol) {
//						System.out.println("No solutions found in the specified time");
//						break;
//					}
					tempConstraints = new IloRange[this.getProblem().getMetabolicNetwork().getNumOfReactions() - tempSol.length];
					int consCount = 0;
					for (int i = 0; i < this.getProblem().getMetabolicNetwork().getNumOfReactions(); i++) {
						boolean hasItem = false;
						int j = 0;
						do {
							hasItem |= (i == tempSol[j]);
							j++;
						} while (!hasItem && j < tempSol.length);
						if (!hasItem) {
							IloLinearIntExpr nexp = cplex.linearIntExpr();
							nexp.addTerm(ZP[i],1);
							nexp.addTerm(ZN[i],1);
							tempConstraints[consCount] = cplex.le(nexp,0);
//							cplex.remove(ZP[i]);
//							cplex.remove(ZN[i]);
							consCount++;
						}
					}
					IloLinearIntExpr tempObj = cplex.linearIntExpr();
					for (int i = 0; i < tempSol.length; i++) {
						tempObj.addTerm(1,ZP[tempSol[i]]);
						tempObj.addTerm(1,ZN[tempSol[i]]);
					}
//					IloObjective ofx = cplex.objective(IloObjectiveSense.Minimize,tempObj);
//					cplex.remove(objectiveFunction);
//					cplex.add(ofx);
					cplex.add(tempConstraints);
					c.setStartTime(c.getCplexTime());
					
					
					
					cplex.solve();
					curSol = c.getFinalSol();
					System.out.println(Arrays.toString(curSol));
//					cplex.remove(tempConstraints);
//					cplex.remove(ofx);
//					cplex.add(objectiveFunction);
//					cplex.add(ZP);
//					cplex.add(ZN);
					curIter++;
					System.out.println("Iteration completed. Solution is optimal? "+cplex.getStatus());
					
				} while (cplex.getStatus()!=IloCplex.Status.Optimal && (curIter < iterMaxPerSol));
				if (cplex.getStatus()==IloCplex.Status.Optimal) {
					cplex.add(getIntegerCut(curSol));
				}
				cplex.remove(tempConstraints);
				for (Reaction ra : getProblem().getExcludedReactions()) {
					int i = getProblem().getMetabolicNetwork().containsReaction(ra.getName());
					if (i > -1) {
						IloLinearIntExpr tempObj = cplex.linearIntExpr();
						tempObj.addTerm(1,ZP[i]);
						tempObj.addTerm(1,ZN[i]);
						cplex.add(cplex.eq(tempObj, 0));
					}
				}
				r.addSolution(curSol);	
			} else {
				break;
			}
			
		}
		return isCompressed ? ((CompressedEnumerationResults) r).decompressResult() : (DefaultEnumerationResult) r;
	}
	
	public void endCplex(){
		this.cplex.end();
	}
	@Override
	public AbstractEnumerationResult solve(int maxsize) throws Exception {
		int sols = this.result.countResults();
		log.addEntry("CHECKPOINT", "Enumeration starting");
//		while (currentSize < maxsize) {
//			calculateNextSize();

//		}
		DefaultCallback callback = new DefaultCallback();
		cplex.use(callback);
		for (int i = 0; i < maxsize; i++) {
			calculateNextSolution();
			sols = this.result.countResults();
		}
		log.addEntry("CHECKPOINT", "Enumeration finished");
		return getResult();
	}
	
	
	private IloRange getIntegerCut(int[] mcsarray) throws IloException{
		IloLinearIntExpr mcs = cplex.linearIntExpr();
		for (int k = 0; k < mcsarray.length; k++) {
			mcs.addTerm(1, ZP[mcsarray[k]]);
			mcs.addTerm(1, ZN[mcsarray[k]]);
		}
		IloRange cut = cplex.le(mcs, mcsarray.length - 1);
		return cut;
	}
	
	private IloIntVar[] createIntVarArray(int lb, int ub, int vars, IloCplex cplex, String prefix) throws IloException{
		IloIntVar[] res = new IloIntVar[vars];
		for (int i = 0; i < res.length; i++) {
			res[i] = cplex.intVar(lb, ub, prefix+i);
		}
		return res;
	}
	
	private IloIntVar[] createBinVarArray(int vars, IloCplex cplex, String prefix) throws IloException{
		IloIntVar[] res = new IloIntVar[vars];
		for (int i = 0; i < res.length; i++) {
			res[i] = cplex.boolVar(prefix+i);
		}
		return res;
	}
	
	private IloNumVar[] createNumVarArray(double lb, double ub, int vars, IloCplex cplex, String prefix) throws IloException{
		IloNumVar[] res = new IloNumVar[vars];
		for (int i = 0; i < res.length; i++) {
			res[i] = cplex.numVar(lb, ub, prefix+i);
		}
		return res;
	}


	public AbstractEnumerationResult getResult() {
		return result;
	}

	public int getCurrentSize() {
		return currentSize;
	}

	private void setCurrentSize(int currentSize) {
		this.currentSize = currentSize;
	}
	
	private int arraySum(int[] array) {
		int res = 0;
		for (int i = 0; i < array.length; i++) {
			res = res + array[i];
		}
		return res;
	}

	public static int[] recoverArray(double[] positive, double[] negative) {
		int size = 0;
		ArrayList<Integer> indexes = new ArrayList<Integer>();
		double[] pres = new double[negative.length];
		for (int i = 0; i < negative.length; i++) {
			pres[i] = (int) (negative[i] + positive[i]);
			if (pres[i] > 0) {
				size += 1;
				indexes.add(i);
			}
		}
		
		int[] res = new int[size];
		for (int i = 0; i < indexes.size(); i++) {
			res[i] = indexes.get(i);
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
	
	private static double[][] validateFVABounds(double[][] sub, FluxVariabilityAnalysisResult fvanormal, FluxVariabilityAnalysisResult fvacomp){
		int N = sub.length;
		double[] rlb = new double[N];
		double[] rub = new double[N];
		double[] fvaNub = fvanormal.maxToArray();
		double[] fvaNlb = fvanormal.minToArray();
		double[] fvaCub = fvacomp.maxToArray();
		double[] fvaClb = fvacomp.minToArray();
		
		for (int i = 0; i < N; i++) {
			int[] nz = MatrixTools.findNonZeroIdx(sub[i]);
			ArrayList<Integer> pos = new ArrayList<Integer>();
			ArrayList<Integer> neg = new ArrayList<Integer>();
			for (int j = 0; j < nz.length; j++) {
				if (sub[i][nz[j]] > 0) {
					pos.add(nz[j]);
				} else {
					neg.add(nz[j]);
				}
				
				double[] posVal = new double[pos.size()];
				for (int k = 0; k < posVal.length; k++) {
					Integer idx = pos.get(k);
					posVal[k] = Math.abs((fvaClb[i] * sub[i][idx]) - fvaNlb[idx]);
				}
				
				double[] negVal = new double[neg.size()];
				for (int k = 0; k < negVal.length; k++) {
					Integer idx = neg.get(k);
					negVal[k] = Math.abs((fvaClb[i] * sub[i][idx]) - fvaNlb[idx]);
				}
				double pv = posVal.length > 0 ? posVal[MatrixTools.getMaximumIdx(posVal)] : Double.NEGATIVE_INFINITY;
				double nv = negVal.length > 0 ? negVal[MatrixTools.getMaximumIdx(negVal)] : Double.NEGATIVE_INFINITY;
				rlb[i] = Math.max(pv,nv);
				
				posVal = new double[pos.size()];
				for (int k = 0; k < posVal.length; k++) {
					Integer idx = pos.get(k);
					posVal[k] = Math.abs((fvaCub[i] * sub[i][idx]) - fvaNub[idx]);
				}
				
				negVal = new double[neg.size()];
				for (int k = 0; k < negVal.length; k++) {
					Integer idx = neg.get(k);
					negVal[k] = Math.abs((fvaCub[i] * sub[i][idx]) - fvaNub[idx]);
				}
				
				pv = posVal.length > 0 ? posVal[MatrixTools.getMaximumIdx(posVal)] : Double.NEGATIVE_INFINITY;
				nv = negVal.length > 0 ? negVal[MatrixTools.getMaximumIdx(negVal)] : Double.NEGATIVE_INFINITY;
				rub[i] = Math.max(pv,nv);
			}
		}
		return new double[][]{rlb,rub};
	}
	
	private boolean isFeasibleFluxVector(AbstractMetabolicNetwork mn, FluxBound[] fb, YieldConstraint[] yc) throws IloException{
		
		CPLEXFluxVariabilityAnalysis fva = new CPLEXFluxVariabilityAnalysis(mn, fb, yc);
		FluxVariabilityAnalysisResult res = fva.solveFVA();
		boolean feasible = true;
		for (int i = 0; i < mn.getNumOfReactions(); i++) {
			feasible &= !res.isUnbounded(i);
		}
		return feasible;
	}
}
