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
