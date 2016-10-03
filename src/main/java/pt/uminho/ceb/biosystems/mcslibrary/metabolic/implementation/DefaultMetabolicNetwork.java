package pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.StringTokenizer;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.AbstractMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Metabolite;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.compression.alg.MatrixTools;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.ReactionConstraint;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Pair;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;

/**
 * Subclass of {@link AbstractMetabolicNetwork}. Holds a typical metabolic network.
 * @author V�tor
 *
 */
public class DefaultMetabolicNetwork extends AbstractMetabolicNetwork {
	Reaction[] reactions;

	@Override
	public boolean isReversible(int reactionIndex) {
		return reactions[reactionIndex].isReversible();
	}

	/**
	 * Constructor for the metabolic network.
	 * @param metab - Array with {@link Metabolite}.
	 * @param reacts - Array with {@link Reaction}.
	 * @param matrix - The stoichiometric matrix for this metabolic network.
	 */
	public DefaultMetabolicNetwork(Metabolite[] metab, Reaction[] reacts, double[][] matrix) {
		try {
			assert metab.length == matrix.length;
			assert reacts.length == matrix[0].length;
		} catch (AssertionError e) {
			System.out.println("Metabolite and reaction data mismatch");
		}
		setMatrix(matrix);
		setMetabolites(metab);
		setReactions(reacts);

	}

	
	public void determineExchange(Reaction[] reacts, Metabolite[] metabs) {
		for (int i = 0; i < reacts.length; i++) {
			boolean external = false;
			for (int j = 0; j < metabs.length; j++) {
				if (Math.abs(j) < MatrixTools.EPSILON) {
					external = true;
					break;
				}
			}
			reacts[i].setExchange(external);
		}
	}
	@Override
	public int containsReaction(String reacName) {
		int res = -1;
		for (int i = 0; i < reactions.length; i++) {
			if (reactions[i].getName().equals(reacName)){
				res = i;
				break;
			}
		}
		return res;
	}

	@Override
	public int getNumOfReactions() {
		return reactions.length;
	}

	public Reaction getReaction(int index) {
		return reactions[index];
	}

	public Reaction getReaction(String name){
		return getReaction(getReactionIndex(name));
	}

	@Override
	public int getReactionIndex(String reac) {
		return containsReaction(reac);
	}

	public void setReactions(Reaction[] reacs) {
		this.reactions = reacs;
	}

	public void setReaction(Reaction reac, int index) {
		this.reactions[index] = reac;
	}

	public void printSize() {
		System.out.println(this.getNumOfMetabolites()+" metabolites and "+this.getNumOfReactions()+" reactions.");
	}
	

	@Override
	public double getLowerBound(int reactionIndex) {
		return this.getReaction(reactionIndex).getBounds().getLower();
	}

	@Override
	public double getUpperBound(int reactionIndex) {
		return this.getReaction(reactionIndex).getBounds().getUpper();
	}

	@Override
	public FluxBound[] getExtraConstraints() {

		HashMap<Integer,Pair<Boolean,Boolean>> constraints = new HashMap<Integer,Pair<Boolean,Boolean>>();

		for (int i = 0; i < reactions.length; i++) {
			boolean unchanged = true;
			double lb = this.getLowerBound(i);
			double ub = this.getUpperBound(i);
			Pair<Boolean,Boolean> cons = new Pair<Boolean, Boolean>(false, false);
			if (lb > -Utilities.INF && Math.abs(lb) > Utilities.PRECISION) {
				cons.setA(true);
				unchanged = false;
			}
			if (ub < Utilities.INF) {
				cons.setB(true);
				unchanged = false;
			}
			if (!unchanged) {
				constraints.put(i, cons);
			}
		}

		FluxBound[] fbarr = new FluxBound[constraints.size()];
		int i = 0;
		for (Integer ridx : constraints.keySet()) {
			double lbv = getReaction(ridx).getBounds().getLower();
			double ubv = getReaction(ridx).getBounds().getUpper();
			fbarr[i] = new FluxBound(getReaction(ridx), lbv, ubv);
			i++;
		}
		return fbarr;
	}

	@Override
	public void saveNetwork(String filename) throws IOException {
		MatrixTools.writeCSV(this.getStoichMatrix(), filename);

		//save reactions
		BufferedWriter reacfile = new BufferedWriter(new FileWriter(filename+".reacts"));
		for (int i = 0; i < this.getNumOfReactions(); i++) {
			String reactionId = this.getReaction(i).getName();
			String lb = Double.toString(this.getReaction(i).getBounds().getLower());
			String ub = Double.toString(this.getReaction(i).getBounds().getUpper());
			String toWrite = reactionId+";"+lb+";"+ub;
			reacfile.write(toWrite+"\n");
		}
		reacfile.flush();
		reacfile.close();

		//save metabolites
		this.saveMetabolites(filename);
	}

	public void overrideBounds(String boundsfile) throws IOException {
		List<String> lines = Utilities.readLines(boundsfile);
		for (String string : lines) {
			StringTokenizer tok = new StringTokenizer(string,";");
			String rName = tok.nextToken();
			String lbs = tok.nextToken();
			String ubs = tok.nextToken();

			double lb;
			try {
				lb = Double.parseDouble(lbs);
			} catch (Exception e) {
				System.out.println("Bad lower bound on "+rName+". Defaulting to negative infinity");
				lb = -Utilities.INF;
			}

			double ub;
			try {
				ub = Double.parseDouble(ubs);
			} catch (Exception e) {
				System.out.println("Bad upper bound on "+rName+". Defaulting to positive infinity");
				ub = Utilities.INF;
			}

			int rId = this.containsReaction(rName);

			reactions[rId].setBounds(new ReactionConstraint(lb,ub));


		}
	}
	public static DefaultMetabolicNetwork loadNetwork(String filename) throws IOException{
		double[][] stmat = MatrixTools.readCSV(filename,false);
		List<String> metabfile = Utilities.readLines(filename+".metabs");
		List<String> reacfile = Utilities.readLines(filename+".reacs");

		Reaction[] reactions = new Reaction[reacfile.size()];
		Metabolite[] metabolites = new Metabolite[metabfile.size()];

		//load reactions
		for (int i = 0; i < reactions.length; i++) {
			StringTokenizer stk = new StringTokenizer(reacfile.get(i), ";");
			String name = stk.nextToken();
			double lb = Double.parseDouble(stk.nextToken());
			double ub = Double.parseDouble(stk.nextToken());
			reactions[i] = new Reaction(name, new ReactionConstraint(lb, ub));

		}

		//load metabolites
		for (int i = 0; i < metabolites.length; i++) {
			StringTokenizer stk = new StringTokenizer(metabfile.get(i), ";");
			String name = stk.nextToken();
			boolean external = Boolean.parseBoolean(stk.nextToken());
			metabolites[i] = new Metabolite(name,external);
		}

		DefaultMetabolicNetwork loadedNetwork = new DefaultMetabolicNetwork(metabolites, reactions, stmat);
		return loadedNetwork;

	}
	
	public void correctCapacities(double M) {
		for (int i = 0; i < this.getNumOfReactions(); i++) {
			Reaction r = this.getReaction(i);
			double lb = r.getBounds().getLower();
			double ub = r.getBounds().getUpper();
			if (ub >= M) {
				ub = Utilities.INF;
			}

			if (lb <= -M) {
				lb = -Utilities.INF;
			}
			r.setBounds(new ReactionConstraint(lb, ub));
			this.setReaction(r, i);
		}
	}


	@Override
	public String toString() {
		
		String str = "Metabolic network with "+reactions.length+" reactions and "+this.getNumOfMetabolites()+" metabolites.";
		for (int i = 0; i < reactions.length; i++) {
			str += "Reaction "+i+": "+reactions[i].toString()+"\n";
		}
		
		for (int i = 0; i < getNumOfMetabolites(); i++) {
			str += "Metabolite "+i+": "+getMetabolite(i).getName()+"\n";
		}
		return str;
	}

}
