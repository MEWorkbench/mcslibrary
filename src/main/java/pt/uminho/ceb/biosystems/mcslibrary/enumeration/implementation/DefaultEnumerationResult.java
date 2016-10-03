package pt.uminho.ceb.biosystems.mcslibrary.enumeration.implementation;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringTokenizer;

import pt.uminho.ceb.biosystems.mcslibrary.enumeration.AbstractEnumerationResult;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.EnumerationProblem;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.CompressedMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;
/**
 * @see AbstractEnumerationResult
 * @author V�tor
 *
 */
public class DefaultEnumerationResult extends AbstractEnumerationResult {

	public DefaultEnumerationResult(EnumerationProblem problem,
			ArrayList<int[]> results) {
		super(problem, results);
	}
	
	public DefaultEnumerationResult(EnumerationProblem problem) {
		super(problem, new ArrayList<int[]>());
	}



	public ArrayList<int[]> getResults() {
		return this.results;
	}

	public EnumerationProblem getProblem(){
		return this.problem;
	}

	public boolean containsExchange(int[] mcs){
		DefaultMetabolicNetwork metaNet = (DefaultMetabolicNetwork) ((this.problem.getMetabolicNetwork().getClass().equals(DefaultMetabolicNetwork.class) ? this.problem.getMetabolicNetwork() : ((CompressedMetabolicNetwork) this.problem.getMetabolicNetwork()).getParentNetwork()));
		boolean res = false;
		for (int i = 0; i < mcs.length; i++) {
			if (metaNet.getReaction(mcs[i]).isExchange()) {
				res = true;
				break;
			}
		}
		return res;
	}
	
	public void addSolution(int[] mcs){
		this.results.add(mcs);
	}
	
	public int countResultsWithoutDrains(){
		int sum = 0;
		for (int[] mcs : this.results) {
			if (!containsExchange(mcs)) {
				sum ++;
			}
		}
		return sum;
	}
	public ArrayList<String[]> toStringArrays() {
		DefaultMetabolicNetwork metaNet = (DefaultMetabolicNetwork) (this.problem.getMetabolicNetwork().getClass() == DefaultMetabolicNetwork.class ? ((DefaultMetabolicNetwork)this.problem.getMetabolicNetwork()) : ((CompressedMetabolicNetwork)this.problem.getMetabolicNetwork()).getParentNetwork());
		ArrayList<String[]> res = new ArrayList<String[]>();
		for (int i = 0; i < results.size(); i++) {
			int[] result = results.get(i);
			String[] strarray = new String[result.length];
			for (int j = 0; j < strarray.length; j++) {
				strarray[j] = metaNet.getReaction(result[j]).getName();
			}
			res.add(strarray);
		}
		return res;
	}

	@Override
	public int countResults() {
		return this.results.size();
	}

	public void printResults(){
		ArrayList<String[]> finalresults = toStringArrays();
		for (int i = 0; i < finalresults.size(); i++) {
			System.out.println(Arrays.toString(finalresults.get(i))+" "+Arrays.toString(this.results.get(i)));
		}
	}

	public static DefaultEnumerationResult fromFile(String filename, EnumerationProblem ep) throws IOException {
		List<String> lines = Utilities.readLines(filename);
		DefaultMetabolicNetwork metaNet = (DefaultMetabolicNetwork) ep.getMetabolicNetwork();
		ArrayList<int[]> results = new ArrayList<int[]>();
		for (String string : lines) {
			String line = string.replace("[", "").replace("]", "");
			StringTokenizer tok = new StringTokenizer(line, ";");
			int length = tok.countTokens();
			int[] res = new int[length];
			for (int i = 0; i < res.length; i++) {
				res[i] = metaNet.containsReaction(tok.nextToken());
			}
			results.add(res);
		}
		return new DefaultEnumerationResult(new EnumerationProblem(metaNet, null, null, null, null, null), results);

	}
	public int[] getResult(int i) {
		return this.results.get(i);
	}

	
	public void writeToFile(String filename) throws IOException {
		ArrayList<String[]> finalresults = toStringArrays();
		BufferedWriter bfw = new BufferedWriter(new FileWriter(filename));
		for (int i = 0; i < finalresults.size(); i++) {
			String[] mcs = finalresults.get(i);
			String str = "";
			for (int j = 0; j < mcs.length; j++) {
				str += mcs[j] + ";";
			}
			if (i == 0) {
				bfw.write(str.substring(0, str.length()-1));
			} else {
				bfw.write("\n"+str.substring(0, str.length()-1));
			}
		}
		bfw.flush();
		bfw.close();
	}

}
