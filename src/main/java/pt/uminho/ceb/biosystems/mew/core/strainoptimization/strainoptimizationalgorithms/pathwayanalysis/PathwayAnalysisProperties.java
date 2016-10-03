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
