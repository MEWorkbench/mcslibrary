package pt.uminho.ceb.biosystems.mew.core.strainoptimization.strainoptimizationalgorithms.pathwayanalysis;

import java.util.Map;

import pt.uminho.ceb.biosystems.mew.biocomponents.container.Container;
import pt.uminho.ceb.biosystems.mew.core.model.steadystatemodel.ISteadyStateModel;
import pt.uminho.ceb.biosystems.mew.core.strainoptimization.configuration.GenericConfiguration;
import pt.uminho.ceb.biosystems.mew.core.strainoptimization.configuration.GenericOptimizationProperties;
import pt.uminho.ceb.biosystems.mew.core.strainoptimization.configuration.IGenericConfiguration;
import pt.uminho.ceb.biosystems.mew.core.strainoptimization.configuration.ISteadyStateConfiguration;
import pt.uminho.ceb.biosystems.mew.core.strainoptimization.objectivefunctions.IObjectiveFunction;
import pt.uminho.ceb.biosystems.mew.core.strainoptimization.strainoptimizationalgorithms.jecoli.JecoliOptimizationProperties;
import pt.uminho.ceb.biosystems.mew.utilities.datastructures.map.indexedhashmap.IndexedHashMap;

public class AbstractPathwayAnalysisConfiguration extends GenericConfiguration implements IGenericConfiguration, ISteadyStateConfiguration{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -5609319229733633773L;

	public AbstractPathwayAnalysisConfiguration() {
		super();
		mandatoryPropertyMap.put(PathwayAnalysisProperties.CONTAINER, Container.class);
		mandatoryPropertyMap.put(GenericOptimizationProperties.SIMULATION_CONFIGURATION, Map.class);
		mandatoryPropertyMap.put(GenericOptimizationProperties.STEADY_STATE_MODEL, ISteadyStateModel.class);
		mandatoryPropertyMap.put(GenericOptimizationProperties.OPTIMIZATION_ALGORITHM, String.class);
		
	}
	
	public Container getContainer(){
		return (Container) propertyMap.get(PathwayAnalysisProperties.CONTAINER);
	}
	
	public ISteadyStateModel getSteadyStateModel(){
		return (ISteadyStateModel) propertyMap.get(GenericOptimizationProperties.STEADY_STATE_MODEL);
	}
	
	public void setContainer(Container cont){
		propertyMap.put(PathwayAnalysisProperties.CONTAINER, cont);
	}
	
	public void setOptimizationAlgorithm(String oa){
		propertyMap.put(GenericOptimizationProperties.OPTIMIZATION_ALGORITHM, oa);
	}
	
	public String getOptimizationAlgorithm(){
		return (String) propertyMap.get(GenericOptimizationProperties.OPTIMIZATION_ALGORITHM);
	}
	@Override
	public IndexedHashMap<IObjectiveFunction, String> getObjectiveFunctionsMap() {
		return null;
	}

	@Override
	public void setObjectiveFunctionsMap(
			IndexedHashMap<IObjectiveFunction, String> objectiveFunctionMap) {
		this.propertyMap.put(GenericOptimizationProperties.MAP_OF2_SIM, objectiveFunctionMap);
	}

	@Override
	public Map<String, Map<String, Object>> getSimulationConfiguration() {
		return (Map<String, Map<String, Object>>) this.getProperty(GenericOptimizationProperties.SIMULATION_CONFIGURATION);
	}

	@Override
	public void setSimulationConfiguration(
			Map<String, Map<String, Object>> simulationConfiguration) {
		this.propertyMap.put(GenericOptimizationProperties.SIMULATION_CONFIGURATION, simulationConfiguration);
	}

	@Override
	public void setModel(ISteadyStateModel model) {
		propertyMap.put(GenericOptimizationProperties.STEADY_STATE_MODEL, model);
		
	}

	@Override
	public String getOptimizationStrategy() {
		return "RK";
	}

	@Override
	public void setOptimizationStrategy(String optimizationStrategy) {
		propertyMap.put(GenericOptimizationProperties.OPTIMIZATION_STRATEGY, optimizationStrategy);
	}

	@Override
	public boolean getIsGeneOptimization() {
		return false;
	}

	@Override
	public boolean getIsOverUnderExpression() {
		return false;
	}
}
