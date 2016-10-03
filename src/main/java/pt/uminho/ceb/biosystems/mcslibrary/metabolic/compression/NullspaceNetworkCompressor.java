package pt.uminho.ceb.biosystems.mcslibrary.metabolic.compression;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Metabolite;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.ReactionGroup;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.compression.alg.MatrixTools;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.ReactionConstraint;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.CompressedMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Pair;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;

public class NullspaceNetworkCompressor {

	private DefaultMetabolicNetwork originalMetaNet;

	public NullspaceNetworkCompressor(DefaultMetabolicNetwork metaNet) {
		this.originalMetaNet = metaNet;
	}
	public static DefaultMetabolicNetwork getSubnetworkExc(DefaultMetabolicNetwork metaNet, List<Reaction> toRemove){
		Reaction[] reacs = new Reaction[metaNet.getNumOfReactions()-toRemove.size()];
		double[][] submatrix = new double[metaNet.getNumOfMetabolites()][metaNet.getNumOfReactions()-toRemove.size()];
		int idx = 0;
		for (int i = 0; i < metaNet.getNumOfReactions(); i++) {
			if (!toRemove.contains(metaNet.getReaction(i))) {
				reacs[idx] = metaNet.getReaction(i);
				for (int j = 0; j < metaNet.getNumOfMetabolites(); j++) {
					submatrix[j][idx] = metaNet.getStoichCoef(j, i);
				}
				idx++;
			}
		}
		Metabolite[] metabs = new Metabolite[metaNet.getNumOfMetabolites()];
		for (int i = 0; i < metabs.length; i++) {
			metabs[i] = metaNet.getMetabolite(i);
		}
		DefaultMetabolicNetwork subnet = new DefaultMetabolicNetwork(metabs, reacs, submatrix);
		return subnet;

	}

	public DefaultMetabolicNetwork getSubnetworkInc(DefaultMetabolicNetwork metaNet, Set<Reaction> toKeep){
		Reaction[] reacs = new Reaction[toKeep.size()];
		double[][] submatrix = new double[metaNet.getNumOfMetabolites()][toKeep.size()];

		int idx = 0;
		for (Reaction reac : toKeep) {
			reacs[idx] = reac;
			int i = 0;
			for (int j = 0; j < metaNet.getNumOfMetabolites(); j++) {
				submatrix[i][idx] = metaNet.getStoichCoef(metaNet.getMetabolite(j), reac);
				i++;
			}
		}

		Metabolite[] metabs = new Metabolite[metaNet.getNumOfMetabolites()];
		for (int i = 0; i < metabs.length; i++) {
			metabs[i] = metaNet.getMetabolite(i);
		}
		DefaultMetabolicNetwork subnet = new DefaultMetabolicNetwork(metabs, reacs, submatrix);
		return subnet;

	}

	public List<Reaction> getNullspaceBlocked(double[][] kernel, double tol, DefaultMetabolicNetwork metaNet) {
		List<Reaction> blocked = new ArrayList<Reaction>();
		for (int i = 0; i < kernel.length; i++) {
			boolean broken = true;
			for (int j = 0; j < kernel[0].length; j++) {
				if (Math.abs(kernel[i][j]) > tol) {
					broken = false;
					break;
				}
			}
			if (broken) {
				blocked.add(metaNet.getReaction(i));
			} else {
				continue;
			}
		}
		return blocked;
	}
	public CompressedMetabolicNetwork reduceModel(List<Reaction> toRemove, List<Reaction> toKeepSingle) throws IOException{
		// get the subnetwork M by N excluding columns in toRemove
		System.out.println("Started:");
		long starttime = System.currentTimeMillis();
		System.out.println("Network loaded");
		DefaultMetabolicNetwork metaNet = this.originalMetaNet;
		metaNet.printSize();

		System.out.println("Stage 1: remove conservation relations before removing ");
//		metaNet = removeConservationRelations(metaNet);
		metaNet.printSize();
		System.out.println(System.currentTimeMillis()-starttime);
		if (toRemove != null && toKeepSingle != null) {
			toKeepSingle.removeAll(toRemove);
		}

		System.out.println("Stage 2: remove blocked reactions in FVA");
		DefaultMetabolicNetwork nonSuppressedReactSubMatrix = getSubnetworkExc(metaNet, toRemove);
		nonSuppressedReactSubMatrix.printSize();
		System.out.println(System.currentTimeMillis()-starttime);

		System.out.println("Stage 3: generate kernel");
		double[][] kernel = MatrixTools.computeKernel(nonSuppressedReactSubMatrix.getStoichMatrix());
		System.out.println(System.currentTimeMillis()-starttime);

		System.out.println("Stage 4: remove kernel blocked reactions");
		double nbtol = Utilities.EPSILON * nonSuppressedReactSubMatrix.getNumOfReactions(); // tolerance for nullspace blocked reactions
		List<Reaction> blocked = getNullspaceBlocked(kernel, nbtol, nonSuppressedReactSubMatrix);
		if (toRemove != null && toKeepSingle != null) {
			toKeepSingle.removeAll(blocked);
		}

		System.out.println(blocked.size()+" blocked reactions.");
		DefaultMetabolicNetwork NSupNBlockSubNet = blocked.size() > 0 ? getSubnetworkExc(nonSuppressedReactSubMatrix, blocked) : nonSuppressedReactSubMatrix;
		System.out.println(System.currentTimeMillis()-starttime);
		kernel = blocked.size() > 0 ? MatrixTools.computeKernel(NSupNBlockSubNet.getStoichMatrix()) : kernel;

		System.out.println("Stage 5: generate CR matrix");
		double crtol = kernel.length*Utilities.EPSILON;
		double[][] cr = generateSubsetMatrix(kernel, crtol);
		System.out.println(System.currentTimeMillis()-starttime);
		System.out.println("Stage 6: generate subsets");
		LinkedHashMap<ArrayList<Reaction>,double[]> subsetmap = new LinkedHashMap<ArrayList<Reaction>,double[]>(kernel.length,(float) 0.75,false);
		ArrayList<Reaction> reactionsInSubset = new ArrayList<Reaction>();

		for (Reaction rtks : toKeepSingle) {
			ArrayList<Reaction> subset = new ArrayList<Reaction>();
			subset.add(rtks);
			subsetmap.put(subset,new double[]{1});
			reactionsInSubset.add(rtks);
		}

		System.out.println("Ungrouped reactions: "+toKeepSingle.size());

		for (int i = 0; i < cr.length; i++) {
			int p = cr.length - 1 - i;
			ArrayList<Integer> indexes = MatrixTools.findNonZeroIdx(cr, p);
			ArrayList<Reaction> reactions = new ArrayList<Reaction>();
			for (Integer idx : indexes) {
				if (!reactionsInSubset.contains(NSupNBlockSubNet.getReaction(idx))) {
					reactions.add(NSupNBlockSubNet.getReaction(idx));
				}
			}
			if (reactions.size() == 0) {
				continue;
			} else {
				reactionsInSubset.addAll(reactions);
				if (reactions.size() == 1) {
					ArrayList<Reaction> ss = new ArrayList<Reaction>();
					ss.addAll(reactions);
					subsetmap.put(ss,new double[]{1});
				} else {
					double[] lengths = new double[reactions.size()];
					for (int j = 0; j < lengths.length; j++) {
						for (int k = 0; k < reactions.size(); k++) {
							Reaction reac = reactions.get(k);
							int idx = NSupNBlockSubNet.getReactionIndex(reac.getName());
							lengths[k] = MatrixTools.getVectorNorm(kernel[idx]);
						}
						double[] difflen = new double[reactions.size()];
						double lenmean = MatrixTools.getVectorMean(lengths);
						for (int g = 0; g < lengths.length; g++) {
							difflen[g] = Math.abs((lengths[g]) - lenmean);
						}
						int ind = MatrixTools.getMinimumIdx(difflen);
						double min = lengths[ind];
						for (int f = 0; f < lengths.length; f++) {
							lengths[f] = lengths[f]/min;
						}
						double[] finalvec = new double[reactions.size()];
						for (int k = 0; k < finalvec.length; k++) {
							finalvec[k] = lengths[k]*cr[NSupNBlockSubNet.getReactionIndex(reactions.get(k).getName())][p];
						}
						subsetmap.put(reactions, finalvec);

					}
				}
			}
		}
		System.out.println(reactionsInSubset.size());

		ArrayList<ArrayList<Reaction>> irrViolSubsets = new ArrayList<ArrayList<Reaction>>();
		for (ArrayList<Reaction> reacts : subsetmap.keySet()) {
			double[] nums = subsetmap.get(reacts);
			boolean invalid = false;
			for (int i = 0; i < nums.length; i++) {
				if (!reacts.get(i).isReversible() && nums[i] < 0){
					invalid = true;
					break;
				}
			}
			if (invalid) {
				irrViolSubsets.add(reacts);
			}
		}

		ArrayList<ArrayList<Reaction>> finalIrrViolSubsets = new ArrayList<ArrayList<Reaction>>();
		if (irrViolSubsets.size()!=0) {
		for (ArrayList<Reaction> subset : irrViolSubsets) {
			boolean right = true;
			double[] arr = new double[subsetmap.get(subset).length];
			for (int i = 0; i < subsetmap.get(subset).length; i++) {
				arr[i] = -subsetmap.get(subset)[i];
				if (!subset.get(i).isReversible()){
					if (arr[i] < 0) {
						right = false;
						break;
					}
				}
			}
			if (right) {
				subsetmap.remove(subset);
				subsetmap.put(subset, arr);
			} else {
				finalIrrViolSubsets.add(subset);
				subsetmap.remove(subset);
			}
		}
		}
		// actual reduction
		
		double[][] sub = new double[subsetmap.size()][metaNet.getNumOfReactions()];
		int idxs = 0;
		for (Entry<ArrayList<Reaction>, double[]> entry : subsetmap.entrySet()) {
			for (int i = 0; i < entry.getValue().length; i++) {
				int idxr = metaNet.getReactionIndex(entry.getKey().get(i).getName());
				sub[idxs][idxr] = entry.getValue()[i];
			}
			idxs++;
		}
		
		CompressedMetabolicNetwork a = reduceMatrix(subsetmap, NSupNBlockSubNet, finalIrrViolSubsets);
		a.setSubMatrix(sub);
		for (int i = 0; i < a.getNumOfReactions(); i++) {
			if (toKeepSingle.contains(a.getReactionGroup(i).getReaction(0))) {
				a.getReactionGroup(i).setBounds(new ReactionConstraint(a.getReactionGroup(i).getReaction(0).getBounds().getLower(), a.getReactionGroup(i).getReaction(0).getBounds().getUpper()));
			}
		}
		System.out.println("FINAL SIZE: "+a.getNumOfMetabolites()+" by "+a.getNumOfReactions());
		System.out.println(System.currentTimeMillis()-starttime);
		return a;
	}
	public static CompressedMetabolicNetwork removeConservationRelations(
			CompressedMetabolicNetwork metaNet) {
		double[][] st = metaNet.getStoichMatrix();
		int[] idxs = MatrixTools.getRrefPivots(st);
		double[][] stcr = new double[idxs.length][st[0].length];
		ReactionGroup[] reacts = new ReactionGroup[metaNet.getNumOfReactions()];
		for (int i = 0; i < idxs.length; i++) {
			for (int j = 0; j < reacts.length; j++) {
				stcr[i][j] = st[idxs[i]][j];
			}
		}
		for (int i = 0; i < reacts.length; i++) {
			reacts[i] = metaNet.getReactionGroup(i);
		}
		Metabolite[] metab = new Metabolite[idxs.length];
		for (int i = 0; i < metab.length; i++) {
			metab[i] = metaNet.getMetabolite(idxs[i]);
		}
		CompressedMetabolicNetwork compNet = new CompressedMetabolicNetwork(metab, reacts, stcr, metaNet.getParentNetwork());
		//compNet = removeConservationRelations(compNet);
		return compNet;
	}

	public static DefaultMetabolicNetwork removeConservationRelations(
			DefaultMetabolicNetwork metaNet) {
		double[][] st = metaNet.getStoichMatrix();
		int[] idxs = MatrixTools.getRrefPivots(st);
		double[][] stcr = new double[idxs.length][st[0].length];
		Reaction[] reacts = new Reaction[metaNet.getNumOfReactions()];
		for (int i = 0; i < idxs.length; i++) {
			for (int j = 0; j < reacts.length; j++) {
				stcr[i][j] = st[idxs[i]][j];
			}
		}
		for (int i = 0; i < reacts.length; i++) {
			reacts[i] = metaNet.getReaction(i);
		}
		Metabolite[] metab = new Metabolite[idxs.length];
		for (int i = 0; i < metab.length; i++) {
			metab[i] = metaNet.getMetabolite(idxs[i]);
		}
		return new DefaultMetabolicNetwork(metab, reacts, stcr);
	}

	public double[][] generateSubsetMatrix(double[][] kernel, double tol) {
		double[][] corrmat = MatrixTools.outerProduct(kernel);
		for (int i = 0; i < corrmat.length; i++) {
			for (int j = i+1; j < corrmat[0].length; j++) {
				corrmat[i][j] = corrmat[i][j]/(Math.sqrt(corrmat[i][i]*corrmat[j][j]));
			}
			corrmat[i][i] = 1;
		}
		// keep upper triangular
		for (int i = 0; i < corrmat.length; i++) {
			for (int j = 0; j < i; j++) {
				corrmat[i][j] = 0;
			}
		}

		for (int i = 0; i < corrmat.length; i++) {
			for (int j = 0; j < corrmat[0].length; j++) {
				if (Math.abs(Math.abs(corrmat[i][j])-1) >= tol) {
					corrmat[i][j] = 0;
				} else {
					corrmat[i][j] = Math.round(corrmat[i][j]/Math.abs(corrmat[i][j]));
				}
			}
		}

		return corrmat;
	}

	private CompressedMetabolicNetwork reduceMatrix(HashMap<ArrayList<Reaction>, double[]> subsetMap, DefaultMetabolicNetwork metaNet, ArrayList<ArrayList<Reaction>> irrViolSubsets) throws IOException{
		ReactionGroup[] rgs = new ReactionGroup[subsetMap.size()];
		double[][] st = metaNet.getStoichMatrix();
		double[][] subtransposed = new double[st[0].length][subsetMap.size()];
		int ssi = 0;
		int irrevtotal = 0;
		for (ArrayList<Reaction> ss : subsetMap.keySet()) {
			ReactionGroup rg = new ReactionGroup(-Utilities.INF, Utilities.INF);
			double ub = Utilities.INF;
			double lb = -Utilities.INF;
			double[] arr = subsetMap.get(ss);
			int irrev = 0;
			for (int i = 0; i < arr.length; i++) {
				rg.addReaction(ss.get(i));
				int idx = metaNet.getReactionIndex(ss.get(i).getName());
				ub = Math.min(ub, metaNet.getUpperBound(idx));
				lb = Math.max(lb, metaNet.getLowerBound(idx));
				subtransposed[idx][ssi] = arr[i];
				if (!ss.get(i).isReversible()) {
					irrev++;
				}
			}
			if (irrev > 0 ){ //&& !irrViolSubsets.contains(ss)) {
				irrevtotal++;
			}
			rg.setBounds(new ReactionConstraint(lb, ub));
			rgs[ssi] = rg;
			ssi++;
		}
		
		double[][] rd = MatrixTools.multiplyMatrix(st,subtransposed);
		for (int i = 0; i < rd.length; i++) {
			for (int j = 0; j < rd[0].length; j++) {
				if (Math.abs(rd[i][j]) < Utilities.PRECISION) {
					rd[i][j] = 0;
				}
			}
		}
		MatrixTools.writeCSV(rd, "DEBUG_RD");
		Pair<double[][], int[]> zeroRowCorrection = MatrixTools.removeZeroRows(rd);
		rd = zeroRowCorrection.getA();
		int[] metabsToKeep = zeroRowCorrection.getB();
		Metabolite[] metabs = new Metabolite[metabsToKeep.length];
		for (int i = 0; i < metabs.length; i++) {
			metabs[i] = metaNet.getMetabolite(metabsToKeep[i]);
		}
		System.out.println("Irreversible subsets:"+irrevtotal);
		CompressedMetabolicNetwork finalnetwork = new CompressedMetabolicNetwork(metabs, rgs, rd, this.originalMetaNet);

//		finalnetwork = removeConservationRelations(finalnetwork);
		return finalnetwork;

	}
}
