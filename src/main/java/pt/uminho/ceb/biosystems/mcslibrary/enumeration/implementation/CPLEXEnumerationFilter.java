package pt.uminho.ceb.biosystems.mcslibrary.enumeration.implementation;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjectiveSense;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import java.util.ArrayList;
import java.util.HashSet;

import pt.uminho.ceb.biosystems.mcslibrary.enumeration.EnumerationProblem;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.AbstractMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.ReactionConstraint;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.YieldConstraint;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.CompressedMetabolicNetwork;

/**
 * Implementation of a minimal cut set filtering algorithm in CPLEX that allows
 * for enumeration of constrained minimal cut sets, provided that you have
 * already calculated the original set of MCS.
 *
 * @author V�tor
 *
 */
public class CPLEXEnumerationFilter {
	protected EnumerationProblem ep;
	private DefaultEnumerationResult result;

	/**
	 *
	 * @param ep
	 *            - The {@link EnumerationProblem} containing all the desired
	 *            fluxes and the metabolic network. Note: This method will only
	 *            consider desired fluxes, as it's assumed the undesired fluxes
	 *            were already established in the results you are to provide to
	 *            this method.
	 * @param result
	 *            - The {@link DefaultEnumerationResult} instance that contains
	 *            the pre-calculated results.
	 */
	public CPLEXEnumerationFilter(EnumerationProblem ep,
			DefaultEnumerationResult result) {
		this.ep = ep;
		this.result = result;
	}

	/**
	 *
	 * @return A boolean array where each index determines whether the MCS with
	 *         the same index is or isn't a cMCS
	 * @throws IloException
	 */
	public DefaultEnumerationResult calculateFilteredResults()
			throws IloException {
		boolean[] boolArray = solve();
		ArrayList<int[]> results = new ArrayList<int[]>();
		for (int i = 0; i < boolArray.length; i++) {
			if (boolArray[i]) {
				results.add(this.result.getResult(i));
			}
		}
		DefaultEnumerationResult res = new DefaultEnumerationResult(this.ep,
				results);
		return res;
	}

	public boolean[] solve() throws IloException {
		AbstractMetabolicNetwork metaNet;
		if (ep.getMetabolicNetwork().getClass() == CompressedMetabolicNetwork.class) {
			metaNet = ((CompressedMetabolicNetwork) ep.getMetabolicNetwork()).getParentNetwork();
		} else {
			metaNet = ep.getMetabolicNetwork();
		}
		int reactions = metaNet.getNumOfReactions();
		int metabolites = metaNet.getNumOfMetabolites();

		IloCplex cplex = new IloCplex();
		cplex.setOut(null);
		IloLinearNumExpr objective = cplex.linearNumExpr();

		IloNumVar[] fluxes = new IloNumVar[reactions];
		IloRange[] stoichiometry = new IloRange[metabolites];

		for (int i = 0; i < reactions; i++) {
			double lb = metaNet.getLowerBound(i);
			double ub = metaNet.getUpperBound(i);
			fluxes[i] = cplex.numVar(lb, ub, "V" + i);
		}

		for (int i = 0; i < metabolites; i++) {
			IloLinearNumExpr row = cplex.linearNumExpr();
			for (int j = 0; j < reactions; j++) {
				row.addTerm(metaNet.getStoichCoef(i, j), fluxes[j]);
			}
			stoichiometry[i] = cplex.eq(row, 0);
		}

		HashSet<Integer> involved = new HashSet<Integer>();

		for (FluxBound fb : ep.getDesiredFluxes()) {
			int idx = metaNet.containsReaction(fb.getReac().getName());
			involved.add(idx);

			double lb = fb.getBounds().getLower();
			double ub = fb.getBounds().getUpper();

			fluxes[idx].setLB(lb);

			fluxes[idx].setUB(ub);
		}

		for (YieldConstraint yc : ep.getDesiredYieldConstraints()) {
			int uidx = metaNet.containsReaction(yc.getUptakeReaction()
					.getName());
			int pidx = metaNet.containsReaction(yc.getProductReaction()
					.getName());
			involved.add(uidx);
			involved.add(pidx);
			double ratio = yc.getRatio();

			if (yc.isLower())
				fluxes[pidx].setLB(ratio * fluxes[uidx].getLB());
			else
				fluxes[pidx].setUB(ratio * fluxes[uidx].getUB());

		}
		int biomassid = this.ep.getMetabolicNetwork().getBiomassReactionId();
		objective.addTerm(1, fluxes[biomassid]);
		cplex.add(fluxes);
		cplex.add(stoichiometry);
		cplex.add(cplex.objective(IloObjectiveSense.Maximize, objective));

		boolean[] resultArray = new boolean[result.countResults()];

		for (int i = 0; i < result.countResults(); i++) {
			double perc = (double) i * 100 / (double) result.countResults();
			if (Math.IEEEremainder(i, 1000) == 0) {
				System.out.println(perc + "%");
			}
			int[] indmcs = result.getResult(i);
			ReactionConstraint[] originalconstraints = new ReactionConstraint[indmcs.length];

			for (int j = 0; j < indmcs.length; j++) {
				originalconstraints[j] = new ReactionConstraint(
						fluxes[indmcs[j]].getLB(), fluxes[indmcs[j]].getUB());
				fluxes[indmcs[j]].setLB(0);
				fluxes[indmcs[j]].setUB(0);
			}

			boolean isFeasible = cplex.solve();

			resultArray[i] = isFeasible;

			for (int j = 0; j < indmcs.length; j++) {
				fluxes[indmcs[j]].setLB(originalconstraints[j].getLower());
				fluxes[indmcs[j]].setUB(originalconstraints[j].getUpper());
			}
		}

		return resultArray;
	}
}
