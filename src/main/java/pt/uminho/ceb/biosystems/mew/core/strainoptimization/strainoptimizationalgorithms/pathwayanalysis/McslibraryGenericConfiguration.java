package pt.uminho.ceb.biosystems.mew.core.strainoptimization.strainoptimizationalgorithms.pathwayanalysis;

import java.util.List;
import java.util.Map;

import pt.uminho.ceb.biosystems.mew.core.model.components.EnvironmentalConditions;
import pt.uminho.ceb.biosystems.mew.core.simulation.components.SimulationProperties;
import pt.uminho.ceb.biosystems.mew.core.strainoptimization.configuration.GenericOptimizationProperties;
import pt.uminho.ceb.biosystems.mew.utilities.datastructures.pair.Pair;

public class McslibraryGenericConfiguration extends AbstractPathwayAnalysisConfiguration{

	/**
	 * 
	 */
	private static final long serialVersionUID = 3859041650999239953L;
//	public static final String ALGORITHM_NAME = "EA";
	
	public McslibraryGenericConfiguration() {
		super();
		
		mandatoryPropertyMap.put(GenericOptimizationProperties.MAX_SET_SIZE, Integer.class);

		optionalPropertyMap.put(PathwayAnalysisProperties.FLUXES_TO_OPTIMIZE, List.class);
		optionalPropertyMap.put(PathwayAnalysisProperties.UNDESIRED_YIELDS, Map.class);
		optionalPropertyMap.put(PathwayAnalysisProperties.MILP_POPULATE, Boolean.class);
		optionalPropertyMap.put(PathwayAnalysisProperties.UNDESIRED_FLUXES, EnvironmentalConditions.class);
		optionalPropertyMap.put(SimulationProperties.ENVIRONMENTAL_CONDITIONS, EnvironmentalConditions.class);
		optionalPropertyMap.put(PathwayAnalysisProperties.EXCLUDED_SOLUTIONS, List.class);
		optionalPropertyMap.put(PathwayAnalysisProperties.EXCLUDED_TARGETS, List.class);
		optionalPropertyMap.put(PathwayAnalysisProperties.MAXIMUM_POOL_SIZE, List.class);
		optionalPropertyMap.put(PathwayAnalysisProperties.DESIRED_FLUXES, EnvironmentalConditions.class);
		optionalPropertyMap.put(PathwayAnalysisProperties.DESIRED_YIELDS, Map.class);
		optionalPropertyMap.put(PathwayAnalysisProperties.SPONTANEOUS_HANDLE, String.class);
	}

	public EnvironmentalConditions getUndesiredFluxes(){
		return (EnvironmentalConditions) propertyMap.get(PathwayAnalysisProperties.UNDESIRED_FLUXES);
	}
	
	public String getBiomassId(){
		return ((List<String>) propertyMap.get(PathwayAnalysisProperties.FLUXES_TO_OPTIMIZE)).get(0);	}
	
	public String getTargetId(){
		return ((List<String>) propertyMap.get(PathwayAnalysisProperties.FLUXES_TO_OPTIMIZE)).get(1);	}
	
	public String getSubstrateId(){
		return ((List<String>) propertyMap.get(PathwayAnalysisProperties.FLUXES_TO_OPTIMIZE)).get(2);
	}
	
	public EnvironmentalConditions getDesiredFluxes(){
		return (EnvironmentalConditions) propertyMap.get(PathwayAnalysisProperties.DESIRED_FLUXES);
	}
	
	@SuppressWarnings("unchecked")
	public Map<Pair<String,String>, Pair<Boolean,Double>> getUndesiredYields() {
		return (Map<Pair<String,String>, Pair<Boolean,Double>>) propertyMap.get(PathwayAnalysisProperties.UNDESIRED_YIELDS);
	}
	
	@SuppressWarnings("unchecked")
	public Map<Pair<String,String>, Pair<Boolean,Double>> getDesiredYields() {
		return (Map<Pair<String,String>, Pair<Boolean,Double>>) propertyMap.get(PathwayAnalysisProperties.DESIRED_YIELDS);
	}
	
	public EnvironmentalConditions getEnvironmentalConditions(){
		return (EnvironmentalConditions) propertyMap.get(SimulationProperties.ENVIRONMENTAL_CONDITIONS);
	}
	
	public boolean isPopulateProblem() {
		return (boolean) propertyMap.get(PathwayAnalysisProperties.MILP_POPULATE);
	}
	
	public int getMaximumSolutionSize() {
		return (int) propertyMap.get(GenericOptimizationProperties.MAX_SET_SIZE);
	}
	
	public int getMaximumPoolSize() {
		return (int) propertyMap.get(PathwayAnalysisProperties.MAXIMUM_POOL_SIZE);
	}
	
	@SuppressWarnings("unchecked")
	public List<String> getExcludedReactions() {
		return (List<String>) propertyMap.get(PathwayAnalysisProperties.EXCLUDED_TARGETS);
	}
	
	//TODO: Change the way solutions are stored. Change this method once that is done
	@SuppressWarnings("unchecked")
	public List<List<String>> getExcludedSolutions() {
		return (List<List<String>>) propertyMap.get(PathwayAnalysisProperties.EXCLUDED_SOLUTIONS);
	}
	
	public Map<String, Map<String, Object>> getSimulationConfiguration(){
		return (Map<String, Map<String, Object>>) propertyMap.get(GenericOptimizationProperties.SIMULATION_CONFIGURATION); 
	}
	
	public String getSpontaneousHandle(){
		return (String) propertyMap.get(PathwayAnalysisProperties.SPONTANEOUS_HANDLE);
	}
	
	public void setUndesiredFluxes(EnvironmentalConditions cond){
        propertyMap.put(PathwayAnalysisProperties.UNDESIRED_FLUXES,cond);
	}
	
	public void setDesiredFluxes(EnvironmentalConditions cond){
        propertyMap.put(PathwayAnalysisProperties.DESIRED_FLUXES,cond);
	}

	public void setUndesiredYields(Map<Pair<String,String>, Pair<Boolean,Double>> cond){
        propertyMap.put(PathwayAnalysisProperties.UNDESIRED_YIELDS,cond);
	}
	
	public void setDesiredYields(Map<Pair<String,String>, Pair<Boolean,Double>> cond) {
		propertyMap.put(PathwayAnalysisProperties.DESIRED_YIELDS, cond);
	}
	
	public void setIsPopulateProblem(boolean isPopulateProblem) {
		propertyMap.put(PathwayAnalysisProperties.MILP_POPULATE, isPopulateProblem);
	}
	
	public void setMaximumSolutionSize(int maxSolutionSize) {
		propertyMap.put(GenericOptimizationProperties.MAX_SET_SIZE, maxSolutionSize);
	}
	
	public void setMaximumPoolSize(int maxPoolSize) {
		propertyMap.put(PathwayAnalysisProperties.MAXIMUM_POOL_SIZE, maxPoolSize);
	}
	
	public void setExcludedSolutions(List<List<String>> excludedSolutions){
		propertyMap.put(PathwayAnalysisProperties.EXCLUDED_SOLUTIONS, excludedSolutions);
	}
	
	public void setExcludedReactions(List<String> excludedReactions) {
		propertyMap.put(PathwayAnalysisProperties.EXCLUDED_TARGETS, excludedReactions);
	}
	
	public void setSpontaneousHandle(String spontaneousHandle) {
		propertyMap.put(PathwayAnalysisProperties.SPONTANEOUS_HANDLE, spontaneousHandle);
	}

	public void setEnvironmentalConditions(
			EnvironmentalConditions environmentalConditions) {
		propertyMap.put(SimulationProperties.ENVIRONMENTAL_CONDITIONS, environmentalConditions);
	}
}

