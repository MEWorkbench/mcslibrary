/*******************************************************************************
 * Copyright 2016
 * CEB Centre of Biological Engineering
 * University of Minho
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This code is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this code. If not, see http://www.gnu.org/licenses/
 *
 * Created inside the BIOSYSTEMS Research Group
 * (http://www.ceb.uminho.pt/biosystems)
 *******************************************************************************/
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
// TODO: Auto-generated Javadoc

/**
 * The Class DefaultEnumerationResult.
 *
 * @author V�tor
 * @see AbstractEnumerationResult
 */
public class DefaultEnumerationResult extends AbstractEnumerationResult {

	/**
	 * Instantiates a new DefaultEnumerationResult.
	 *
	 * @param problem the EnumerationProblem
	 * @param results ArrayList with integer arrays containing reaction indices
	 */
	public DefaultEnumerationResult(EnumerationProblem problem,
			ArrayList<int[]> results) {
		super(problem, results);
	}
	
	/**
	 * Instantiates a new DefaultEnumerationResult with no solutions.
	 *
	 * @param problem the EnumerationProblem
	 */
	public DefaultEnumerationResult(EnumerationProblem problem) {
		super(problem, new ArrayList<int[]>());
	}



	/**
	 * Returns the results in integer array format
	 *
	 * @return an ArrayList with integer arrays
	 */
	public ArrayList<int[]> getResults() {
		return this.results;
	}

	/**
	 * Gets the EnumerationProblem.
	 *
	 * @return the EnumerationProblem
	 */
	public EnumerationProblem getProblem(){
		return this.problem;
	}

	/**
	 * Returns a true or false if the minimal cut set contains any exchange reactions.
	 *
	 * @param mcs the mcs
	 * @return true, if the minimal cut set contains any exchange reactions
	 */
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
	
	/* (non-Javadoc)
	 * @see pt.uminho.ceb.biosystems.mcslibrary.enumeration.AbstractEnumerationResult#addSolution(int[])
	 */
	public void addSolution(int[] mcs){
		this.results.add(mcs);
	}
	
	/**
	 * Count results without drains.
	 *
	 * @return the number of results without drains.
	 */
	public int countResultsWithoutDrains(){
		int sum = 0;
		for (int[] mcs : this.results) {
			if (!containsExchange(mcs)) {
				sum ++;
			}
		}
		return sum;
	}
	
	/* (non-Javadoc)
	 * @see pt.uminho.ceb.biosystems.mcslibrary.enumeration.AbstractEnumerationResult#toStringArrays()
	 */
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

	/* (non-Javadoc)
	 * @see pt.uminho.ceb.biosystems.mcslibrary.enumeration.AbstractEnumerationResult#countResults()
	 */
	@Override
	public int countResults() {
		return this.results.size();
	}

	/* (non-Javadoc)
	 * @see pt.uminho.ceb.biosystems.mcslibrary.enumeration.AbstractEnumerationResult#printResults()
	 */
	public void printResults(){
		ArrayList<String[]> finalresults = toStringArrays();
		for (int i = 0; i < finalresults.size(); i++) {
			System.out.println(Arrays.toString(finalresults.get(i))+" "+Arrays.toString(this.results.get(i)));
		}
	}

	/**
	 * Instantiates a DefaultEnumerationResult from a file.
	 *
	 * @param filename Path pointing to the file containing the minimal cut sets.
	 * @param ep the EnumerationProblem
	 * @return the DefaultEnumerationResult
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
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
	
	/**
	 * Gets a single minimal cut set.
	 *
	 * @param i the index
	 * @return the integer array containing indices which code for the reactions of a minimal cut set
	 */
	public int[] getResult(int i) {
		return this.results.get(i);
	}

	
	/**
	 * Write a DefaultEnumerationResult to a file.
	 *
	 * @param filename Path pointing to the file containing the minimal cut sets.
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
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
