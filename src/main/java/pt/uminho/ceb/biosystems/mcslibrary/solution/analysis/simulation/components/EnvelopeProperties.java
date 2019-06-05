package pt.uminho.ceb.biosystems.mcslibrary.solution.analysis.simulation.components;

import java.io.IOException;
import java.util.Map;

import pt.uminho.ceb.biosystems.mcslibrary.solution.SolutionUtilities;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;

public class EnvelopeProperties {
	private String pivot;
	private String objective;
	private String confName;
	private int steps;
	
	public final static String OBJECTIVE_NAME = "OBJECTIVE";
	public final static String PIVOT_NAME = "PIVOT";
	public final static String CONFNAME = "NAME";
	public final static String STEPS = "STEPS";


	public EnvelopeProperties(String objective, String pivot, String confName, int steps) {
		/*
		 * objective: Reaction ID representing the flux to maximize/minimize
		 * pivot: Reaction ID representing the flux to vary along its range
		 * confName: String with the name of this configuration
		 * steps: Number of steps to divide the admissible range for the pivot reaction
		 */
		this.objective = objective;
		this.pivot = pivot;
		this.confName = confName;
		this.steps = steps;
	}

	public String getPivot() {
		return pivot;
	}

	public String getObjective() {
		return objective;
	}

	public String getConfName() {
		return confName;
	}
	
	public int getSteps(){
		return steps;
	}
	
	public static EnvelopeProperties fromFile(String file) throws IOException{
		Map<String, String> map = Utilities.readPropertyMap(file);
		return new EnvelopeProperties(map.get(OBJECTIVE_NAME), map.get(PIVOT_NAME), map.get(CONFNAME), Integer.parseInt(map.get(STEPS)));
	}
	
	
}
