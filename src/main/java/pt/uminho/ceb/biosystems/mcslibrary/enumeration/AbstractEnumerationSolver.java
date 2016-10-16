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
/**
 * Abstract class that receives an enumeration problem. Contains methods to calculate minimal cut sets for the provided problem
 * @author vvieira
 *
 *
 */
public abstract class AbstractEnumerationSolver {
	private EnumerationProblem eprob;

	public AbstractEnumerationSolver(EnumerationProblem eprob) {
		this.eprob = eprob;
	}
	
	public EnumerationProblem getProblem() {
		return this.eprob;
	}
	
	/**
	 *
	 * @param maxsize - An integer specifying the maximum size of the minimal cut sets to be calculated
	 * @return The enumeration result as an instance of a subclass of AbstractEnumerationResult
	 * @throws Exception
	 */
	public abstract AbstractEnumerationResult solve(int maxsize) throws Exception;
}
