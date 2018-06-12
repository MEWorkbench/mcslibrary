package pt.uminho.ceb.biosystems.mcslibrary.solution.analysis;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import ilog.concert.IloException;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.SimulationResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba.CPLEXFluxBalanceAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;

public class RobustnessCriteria {
	private String objective;
	private double fractionOfMinimum;
	private String decisive;
	private String senseObj;
	private String senseDec;
	private String name;
	
	public static String objP = "OBJECTIVE";
	public static String fomP = "FRACTION";
	public static String decP = "DECISIVE";
	public static String snsoP = "OBJECTIVE_SENSE";
	public static String snsdP = "DECISIVE_SENSE";
	public static String nameP = "NAME";

	

	public RobustnessCriteria(String objective, String decisive, double fractionOfMinimum, String senseObj, String senseDec, String name) {
		this.objective = objective;
		this.fractionOfMinimum = fractionOfMinimum;
		this.decisive = decisive;
		this.senseObj = senseObj;
		this.senseDec = senseDec;
		this.name = name;
	}
	
	public boolean evaluate(DefaultMetabolicNetwork dmn, CPLEXFluxBalanceAnalysis fba, List<String> solution, FluxBound[] remainingCriteria) throws IloException{
		double pivotMin = fba.solveReactionKnockoutFVA(remainingCriteria, Utilities.toReacArrayFromString(dmn, solution), dmn.getReaction(decisive), senseDec, fractionOfMinimum ,dmn.getReaction(objective)).getFluxValue(dmn.getReaction(decisive));
		return pivotMin > 1e-9;
	}
	
	public String getName(){
		return this.name;
	}

	public static RobustnessCriteria fromFile(String absolutePath, DefaultMetabolicNetwork mn) throws IOException {
		Map<String, String> map = Utilities.readPropertyMap(absolutePath);
		return new RobustnessCriteria(map.get(objP), map.get(decP), Double.parseDouble(map.get(fomP)), map.get(snsoP), map.get(snsdP), map.get(nameP));
	}
}
