package pt.uminho.ceb.biosystems.mcslibrary.utilities;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.ReactionConstraint;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.Container;
import pt.uminho.ceb.biosystems.mew.core.model.components.Metabolite;
import pt.uminho.ceb.biosystems.mew.core.model.components.Reaction;
import pt.uminho.ceb.biosystems.mew.core.model.steadystatemodel.ISteadyStateModel;



public class SteadyStateModelReader {
	private ISteadyStateModel model;
	private Integer nreacts;
	private Integer nmetabs;

	public SteadyStateModelReader(ISteadyStateModel model) {
		this.model = model;
		nreacts = model.getNumberOfReactions();
		nmetabs = model.getNumberOfMetabolites();
	}

	public pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction[] createReactions(){
		pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction[] res = new pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction[nreacts];
		for (int i = 0; i < nreacts; i++) {
			Reaction r = model.getReaction(i);
			ReactionConstraint limits = new ReactionConstraint(r.getConstraints().getLowerLimit(), r.getConstraints().getUpperLimit());
			res[i] = new pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction(r.getId(), limits);
		}
		return res;
	}

	public pt.uminho.ceb.biosystems.mcslibrary.metabolic.Metabolite[] createMetabolites(){
		pt.uminho.ceb.biosystems.mcslibrary.metabolic.Metabolite[] res = new pt.uminho.ceb.biosystems.mcslibrary.metabolic.Metabolite[nmetabs];
		for (int i = 0; i < nmetabs; i++) {
			Metabolite m = model.getMetabolite(i);
			res[i] = new pt.uminho.ceb.biosystems.mcslibrary.metabolic.Metabolite(m.getId(), m.isBoundaryCondition());
		}
		return res;
	}

	public double[][] createMatrix(){
		double[][] res = new double[nmetabs][nreacts];
		for (int i = 0; i < nmetabs; i++) {
			for (int j = 0; j < nreacts; j++) {
				res[i][j] = model.getStoichiometricValue(i, j);
			}
		}
		return res;
	}

	public DefaultMetabolicNetwork convertModel(){
		DefaultMetabolicNetwork mod = new DefaultMetabolicNetwork(createMetabolites(),createReactions(),createMatrix());
		return mod;
	}
	
	public static void updateFormulae(Container cont, DefaultMetabolicNetwork mn){
		for (int i = 0; i < mn.getNumOfMetabolites(); i++) {
			try {
				mn.getMetabolite(i).setFormula(cont.getMetabolite(mn.getMetabolite(i).getName()).getFormula());
			} catch (Exception e) {
				System.out.println("Error");
			}
		}
	}

	public static void convertRxIdToName(Container cont, DefaultMetabolicNetwork mn) {
		for (int i = 0; i < mn.getNumOfReactions(); i++) {
			try {
				mn.getReaction(i).setName(cont.getReaction(mn.getReaction(i).getName()).getName());
			} catch (Exception e) {
				System.out.println("Error");
			}
		}
	}
}

