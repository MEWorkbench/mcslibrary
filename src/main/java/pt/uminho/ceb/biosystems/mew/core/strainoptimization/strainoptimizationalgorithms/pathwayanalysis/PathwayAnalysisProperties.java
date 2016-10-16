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
package pt.uminho.ceb.biosystems.mew.core.strainoptimization.strainoptimizationalgorithms.pathwayanalysis;


public class PathwayAnalysisProperties {
	
	// general PA
	public static final String CONTAINER = "mcslibrary.container";
	public static final String STEADY_STATE_MODEL = "mcslibrary.ssmodel";
	
	// mcslibrary
	public static final String UNDESIRED_FLUXES = "mcslibrary.fluxbounds.undesired";
	public static final String DESIRED_FLUXES = "mcslibrary.fluxbounds.desired";
	public static final String UNDESIRED_YIELDS = "mcslibrary.yieldconstraints.undesired";
	public static final String DESIRED_YIELDS = "mcslibrary.yieldconstraints.desired";
	public static final String ENVIRONMENTAL_CONDITIONS = "mcslibrary.fluxbounds.fva";
	public static final String EXCLUDED_TARGETS = "mcslibrary.excludedreactions";
	public static final String EXCLUDED_SOLUTIONS = "mcslibrary.excludedsolutions";
	
	public static final String MAXIMUM_SOLUTION_SIZE = "mcslibrary.maxsolutionsize";
	public static final String MILP_POPULATE = "mcslibrary.milp_populate";
	public static final String MAXIMUM_POOL_SIZE = "mcslibrary.maxpoolsize";
	public static final String SPONTANEOUS_HANDLE = "mcslibrary.spontaneoushandle";
	public static final String OPTIMIZATION_STRATEGY = "jecoli.optimizationstrategy";
	public static final String OPTIMIZATION_ALGORITHM = "jecoli.optimizationalalgorithm";
	
	// algorithms
	
	public static final String MCSENUMERATOR = "MCS";
	public static final String FLUXES_TO_OPTIMIZE = "mcslibrary.fluxestooptimize";
}
