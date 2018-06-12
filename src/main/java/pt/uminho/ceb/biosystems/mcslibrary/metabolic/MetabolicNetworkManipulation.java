package pt.uminho.ceb.biosystems.mcslibrary.metabolic;

import java.util.Arrays;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.CompressedMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;

public class MetabolicNetworkManipulation {

	
	public static void addReaction(Reaction r, AbstractMetabolicNetwork abs, double[] stoichCoefs) {
		
		if (abs.getClass() == DefaultMetabolicNetwork.class) {
			DefaultMetabolicNetwork mn = (DefaultMetabolicNetwork) abs;
			Object[] newa = Utilities.appendItemToArray(mn.getReactions(), r);
			Reaction[] newRx = Arrays.copyOf(newa, abs.getNumOfReactions()+1, Reaction[].class);
			double[][] matrix = abs.getStoichMatrix();
			for (int i = 0; i < stoichCoefs.length; i++) {
				matrix[i] = Utilities.appendItemToArray(matrix[i], stoichCoefs[i]);
			}
			abs.setMatrix(matrix);
			mn.setReactions(newRx);
			
		} else if (abs.getClass() == CompressedMetabolicNetwork.class) {
			CompressedMetabolicNetwork cmn = (CompressedMetabolicNetwork) abs;
			Object[] newa = Utilities.appendItemToArray(
					cmn.getReactionGroups(), 
					new ReactionGroup(r.getBounds().getLower(), r.getBounds().getUpper()));
			ReactionGroup[] newRg = Arrays.copyOf(newa, abs.getNumOfReactions()+1, ReactionGroup[].class);
			newRg[cmn.getNumOfReactions()].addReaction(r);
			double[][] matrix = abs.getStoichMatrix();
			for (int i = 0; i < stoichCoefs.length; i++) {
				matrix[i] = Utilities.appendItemToArray(matrix[i], stoichCoefs[i]);
			}
			abs.setMatrix(matrix);
			cmn.setReactionGroups(newRg);
		}
	}
	
	public static void removeReaction(String rName, AbstractMetabolicNetwork abs) {
		
		if (abs.getClass() == DefaultMetabolicNetwork.class) {
			DefaultMetabolicNetwork mn = (DefaultMetabolicNetwork) abs;
			int id = mn.getReactionIndex(rName);
			Reaction[] newReactions = new Reaction[abs.getNumOfReactions() - 1];
			double[][] newMatrix = new double[abs.getNumOfMetabolites()][abs.getNumOfReactions()-1];
			for (int i = 0; i < id; i++) {
				newReactions[i] = mn.getReaction(i);
				for (int j = 0; j < newMatrix.length; j++) {
					newMatrix[j][i] = abs.getStoichCoef(j, i);
				}
			}
			for (int i = id; i < abs.getNumOfReactions()-1; i++) {
				newReactions[i] = mn.getReaction(i+1);
				for (int j = 0; j < newMatrix.length; j++) {
					newMatrix[j][i] = abs.getStoichCoef(j, i+1);
				}
			}
			mn.setReactions(newReactions);
			mn.setMatrix(newMatrix);
		} else if (abs.getClass() == CompressedMetabolicNetwork.class) {
			CompressedMetabolicNetwork cmn = (CompressedMetabolicNetwork) abs;
			ReactionGroup rg = cmn.getReactionGroupFromReaction(rName);
			if (rg.size() < 2) {
				int id = cmn.getReactionIndex(rName);
				ReactionGroup[] newReactions = new ReactionGroup[abs.getNumOfReactions() - 1];
				double[][] newMatrix = new double[abs.getNumOfMetabolites()][abs.getNumOfReactions()-1];
				for (int i = 0; i < id; i++) {
					newReactions[i] = cmn.getReactionGroup(i);
					for (int j = 0; j < newMatrix.length; j++) {
						newMatrix[j][i] = abs.getStoichCoef(j, i);
					}
				}
				for (int i = id; i < abs.getNumOfReactions()-1; i++) {
					newReactions[i] = cmn.getReactionGroup(i+1);
					for (int j = 0; j < newMatrix.length; j++) {
						newMatrix[j][i] = abs.getStoichCoef(j, i+1);
					}
				}
				cmn.setReactionGroups(newReactions);
				cmn.setMatrix(newMatrix);
			} else {
				rg.removeReaction(rg.containsReaction(rName));
			}
		}

	}
	
	public static void addMetabolite(Metabolite m, AbstractMetabolicNetwork abs, double[] stoichCoefs) {
		Object[] mArray = Utilities.appendItemToArray(abs.getMetabolites(), m);
//		for (int i = 0; i < mArray.length; i++) {
////			System.out.println(mArray[i].getClass());
//		}
		Metabolite[] mArrayM = Arrays.copyOf(mArray, mArray.length, Metabolite[].class);
		double[][] newMatrix = new double[abs.getNumOfMetabolites()+1][abs.getNumOfReactions()];
		for (int i = 0; i < abs.getNumOfMetabolites(); i++) {
			newMatrix[i] = abs.getStoichMatrix()[i];
		}
		newMatrix[abs.getNumOfMetabolites()] = stoichCoefs;
		abs.setMetabolites(mArrayM);
		abs.setMatrix(newMatrix);
	}
	
	public static void removeMetabolite(String mName, AbstractMetabolicNetwork abs) {
		int id = abs.getMetaboliteIndex(mName);
		double[][] newMatrix = new double[abs.getNumOfMetabolites()-1][abs.getNumOfReactions()];
		for (int i = 0; i < id; i++) {
			newMatrix[i] = abs.getStoichMatrix()[i];
		}
		for (int i = id; i < abs.getNumOfMetabolites()-1; i++) {
			newMatrix[i] = abs.getStoichMatrix()[i+1];
		}
		Object[] mArray = Utilities.removeItemFromArray(abs.getMetabolites(), id);
		Metabolite[] mArrayM = Arrays.copyOf(mArray, mArray.length, Metabolite[].class);
		abs.setMetabolites(mArrayM);
		abs.setMatrix(newMatrix);
	}
	
	//TODO: encontrar metabolitos de prod/cons 
	//TODO: (getReactionWithSpecifiedBranching)
}
