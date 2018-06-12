package pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation;

import java.io.IOException;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.AbstractMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Metabolite;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.ReactionGroup;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.ReactionConstraint;
/**
 * Subclass of {@link AbstractMetabolicNetwork}. Holds a metabolic network where reactions are clusters of other reactions belonging to another, uncompressed metabolic network.
 * @see AbstractMetabolicNetwork
 * @author V�tor
 *
 */
public class CompressedMetabolicNetwork extends AbstractMetabolicNetwork {
	private ReactionGroup[] reactions;
	private AbstractMetabolicNetwork parentNetwork;
	private double[][] subMatrix;
	/**
	 * Constructor for the compressed metabolic network.
	 * @param metab - Array with {@link Metabolite}.
	 * @param rgs - Array with {@link ReactionGroup}.
	 * @param matrix - The stoichiometric matrix for this metabolic network.
	 * @param parentNetwork - The {@link DefaultMetabolicNetwork} that was used to get the compressed matrix.
	 */
	public CompressedMetabolicNetwork(Metabolite[] metab, ReactionGroup[] rgs, double[][] matrix, AbstractMetabolicNetwork parentNetwork) {
		try {
			assert metab.length == matrix.length;
			assert rgs.length == matrix[0].length;
		} catch (AssertionError e) {
			System.out.println("Metabolite and reaction data mismatch");
		}
		this.setParentNetwork(parentNetwork);
		setMatrix(matrix);
		setMetabolites(metab);
		setReactionGroups(rgs);

	}
	@Override
	public boolean isReversible(int reactionIndex) {
		return reactions[reactionIndex].isReversible();
	}

	@Override
	public int containsReaction(String reacName) {
		int res = -1;
		for (int i = 0; i < reactions.length; i++) {
			for (int j = 0; j < reactions[i].size(); j++) {
				if (reactions[i].containsReaction(reacName) > -1) {
					res = j;
					break;
				}
			}
		}
		return res;
	}

	@Override
	public int getNumOfReactions(){
		return reactions.length;
	}
	
	public Reaction getReaction(String reacName){
		return this.reactions[getReactionIndex(reacName)].getReaction(containsReaction(reacName));
	}
	
	public ReactionGroup getReactionGroupFromReaction(String reacName){
		return this.reactions[getReactionIndex(reacName)];
	}

	public int getTotalNumOfReactions() {
		int sum = 0;
		for (int i = 0; i < reactions.length; i++) {
			sum = sum + reactions[i].size();
		}
		return sum;
	}

	public ReactionGroup getReactionGroup(int index) {
		return reactions[index];
	}

	@Override
	public int getReactionIndex(String reac) {
		int res = -1;
		for (int i = 0; i < reactions.length; i++) {
			if (reactions[i].containsReaction(reac) > -1) {
				res = i;
				break;
			}
		}
		return res;
	}


	public void setReactionGroups(ReactionGroup[] rgs) {
		this.reactions = rgs;
	}
	public AbstractMetabolicNetwork getParentNetwork() {
		return parentNetwork;
	}
	private void setParentNetwork(AbstractMetabolicNetwork parentNetwork) {
		this.parentNetwork = parentNetwork;
	}
	public int[] convertGroupToIntArray(int rgIndex) {
		ReactionGroup rg = this.getReactionGroup(rgIndex);
		int[] res = new int[rg.size()];
		for (int i = 0; i < rg.size(); i++) {
			res[i] = this.getParentNetwork().getReactionIndex(rg.getReaction(i).getName());
		}
		return res;
	}
	@Override
	public double getLowerBound(int reactionIndex) {
		return this.getReactionGroup(reactionIndex).getBounds().getLower();
	}
	@Override
	public double getUpperBound(int reactionIndex) {
		// TODO Auto-generated method stub
		return this.getReactionGroup(reactionIndex).getBounds().getUpper();
	}
	@Override
	public FluxBound[] getExtraConstraints() {
		return getParentNetwork().getExtraConstraints();
	}

	public String toString() {
		String str = "";
		for (int i = 0; i < reactions.length; i++) {
			str = str + i+": "+reactions[i]+"\n";
		}
		return str;
	}

	@Override
	public void saveNetwork(String filename) throws IOException {
		this.parentNetwork.saveNetwork(filename);
	}
	public double[][] getSubMatrix() {
		return subMatrix;
	}
	public void setSubMatrix(double[][] subMatrix) {
		this.subMatrix = subMatrix;
	}
	public ReactionGroup[] getReactionGroups() {
		return reactions;
	}
	@Override
	public ReactionConstraint getBound(int index) {
		return this.reactions[index].getBounds();
	}

}
