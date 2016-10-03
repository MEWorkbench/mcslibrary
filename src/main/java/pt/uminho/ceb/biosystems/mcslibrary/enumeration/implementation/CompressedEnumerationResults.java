package pt.uminho.ceb.biosystems.mcslibrary.enumeration.implementation;

import java.util.ArrayList;
import java.util.Arrays;

import pt.uminho.ceb.biosystems.mcslibrary.enumeration.AbstractEnumerationResult;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.EnumerationProblem;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.CompressedMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;
/**
 * Subclass of {@link AbstractEnumerationResult}. Contains the results with the indexes for the compressed network. 
 * @author V�tor
 *
 */
public class CompressedEnumerationResults extends AbstractEnumerationResult{
	
	public CompressedEnumerationResults(EnumerationProblem problem,
			ArrayList<int[]> results) {
		super(problem, results);
	}
	
	public CompressedEnumerationResults(EnumerationProblem problem) {
		super(problem, new ArrayList<int[]>());
	}
	/**
	 * 
	 * @return
	 * An instance of {@link DefaultEnumerationResult} with decompressed results to match the original metabolic network
	 */
	public DefaultEnumerationResult decompressResult() {
		CompressedMetabolicNetwork compNet = (CompressedMetabolicNetwork) problem.getMetabolicNetwork();
		ArrayList<int[]> decompressedResults = new ArrayList<int[]>();
		for (int i = 0; i < this.results.size(); i++) {
			ArrayList<int[]> optionArray = new ArrayList<int[]>();
			for (int j = 0; j < results.get(i).length; j++) {
				optionArray.add(compNet.convertGroupToIntArray(results.get(i)[j]));
			}
			decompressedResults.addAll(Utilities.intArrayDecompressor(optionArray));
		}
		problem = new EnumerationProblem(compNet.getParentNetwork(),this.problem.getUndesiredFluxes(),this.problem.getDesiredFluxes(),this.problem.getUndesiredYieldConstraints(),this.problem.getDesiredYieldConstraints(),this.problem.getExcludedReactions());
		DefaultEnumerationResult res = new DefaultEnumerationResult(problem, decompressedResults);
		return res;
	}
	@Override
	public ArrayList<String[]> toStringArrays() {
		DefaultEnumerationResult res = decompressResult();
		return res.toStringArrays();
	}
	

	@Override
	public int countResults() {
		int sum = 0;
		CompressedMetabolicNetwork comp = (CompressedMetabolicNetwork) this.problem.getMetabolicNetwork();
		for (int i = 0; i < results.size(); i++) {
			int presize = 1;
			for (int j = 0; j < results.get(i).length; j++) {
					presize = presize
					* (comp.getReactionGroup(results.get(i)[j]).size());
			}
			sum += presize;
		}
		return sum;
	}
	
	public void addSolution(int[] mcs){
		this.results.add(mcs);
	}
	
	public void printResults(){
		ArrayList<String[]> finalresults = decompressResult().toStringArrays();
		for (int i = 0; i < finalresults.size(); i++) {
			System.out.println(Arrays.toString(finalresults.get(i)));
		}
	}

}
