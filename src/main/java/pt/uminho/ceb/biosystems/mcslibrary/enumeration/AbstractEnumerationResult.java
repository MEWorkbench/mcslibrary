/*
 *
 */
package pt.uminho.ceb.biosystems.mcslibrary.enumeration;

import java.util.ArrayList;

/**
 * An abstract class that contains the results from the solver and the problem that originated them.
 *
 */
public abstract class AbstractEnumerationResult {

	protected EnumerationProblem problem;
	protected ArrayList<int[]> results;

	/**
	 *
	 * @param problem - The EnumerationProblem instance that was used to make these results
	 * @param results - An ArrayList with the Integer arrays that resulted from
	 *
	 * Constructor for the Enumeration Result and subclasses implementing its methods.
	 */
	public AbstractEnumerationResult(EnumerationProblem problem, ArrayList<int[]> results) {
		this.problem = problem;
		this.results = results;
	}
	public abstract int countResults();
	/**
	 *
	 * @return An ArrayList with the results as Strings with the reaction names instead of their index.
	 */
	public abstract ArrayList<String[]> toStringArrays();
	public abstract void printResults();
	public abstract void addSolution(int[] mcs);
}
