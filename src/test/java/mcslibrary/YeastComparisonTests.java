package mcslibrary;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.stream.XMLStreamException;

import org.xml.sax.SAXException;

import ilog.concert.IloException;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.implementation.DefaultEnumerationResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.ReactionConstraint;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba.CPLEXFluxBalanceAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba.CPLEXParsimoniousFluxBalanceAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.solution.SolutionUtilities;
import pt.uminho.ceb.biosystems.mcslibrary.solution.analysis.RobustnessCriteria;
import pt.uminho.ceb.biosystems.mcslibrary.solution.analysis.SolutionAnalysisPipeline;
import pt.uminho.ceb.biosystems.mcslibrary.solution.analysis.SolutionFilter;
import pt.uminho.ceb.biosystems.mcslibrary.solution.analysis.SolutionFilterPipeline;
import pt.uminho.ceb.biosystems.mcslibrary.solution.analysis.simulation.components.EnvelopeProperties;
import pt.uminho.ceb.biosystems.mcslibrary.solution.scoring.SolutionScorer;
import pt.uminho.ceb.biosystems.mcslibrary.solution.scoring.scoreitems.BPCYScoreItem;
import pt.uminho.ceb.biosystems.mcslibrary.solution.scoring.scoreitems.CarbonYieldScoreItem;
import pt.uminho.ceb.biosystems.mcslibrary.solution.scoring.scoreitems.IScoreItem;
import pt.uminho.ceb.biosystems.mcslibrary.solution.scoring.scoreitems.PFBAFluxValueItem;
import pt.uminho.ceb.biosystems.mcslibrary.solution.scoring.scoreitems.RobustnessScoreItem;
import pt.uminho.ceb.biosystems.mcslibrary.solution.scoring.scoreitems.SolutionSizeItem;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.MCSPipeline;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.SteadyStateModelReader;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.Container;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.io.readers.ErrorsException;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.io.readers.JSBMLReader;
import pt.uminho.ceb.biosystems.mew.biocomponents.validation.io.JSBMLValidationException;

public class YeastComparisonTests {
	static MCSPipeline mcs;
	static CPLEXFluxBalanceAnalysis fba;
	static DefaultMetabolicNetwork model;
	private static Reaction biomassReaction;
	private static Container cont;
	private static CPLEXParsimoniousFluxBalanceAnalysis pfba;
	
	// Defining important reaction names for use with the various methods below
	final static String atpm = "R_ATPM";
	final static String biomass = "R_biomass_SC5_notrace";
	final static String glucose = "R_EX_glc_e_";
	final static String o2 = "R_EX_o2_e_";
	final static String product = "R_EX_succ_e_";

	public static void setUpPipeline() throws FileNotFoundException, IOException, XMLStreamException, ErrorsException, ParserConfigurationException, SAXException, JSBMLValidationException, IloException{
		
		// Read the SBML model from the specified path
		String path = "src/test/resources/iMM904_RP/iMM904_corrected_201609001170165318.xml";
		cont = new Container(new JSBMLReader(path, "default", false));
		mcs = new MCSPipeline(cont, "");
		model = mcs.getMetabolicNetwork();
		
		biomassReaction = model.getReaction(biomass);
		
		// Set default bounds for maintenance ATP demand and oxygen/glucose uptakes.
		model.getReaction(atpm).setBounds(new ReactionConstraint(1, 1));
		model.getReaction(o2).setBounds(new ReactionConstraint(-Utilities.INF, Utilities.INF));
        model.getReaction(glucose).setBounds(new ReactionConstraint(-1.15, Utilities.INF));
		
		mcs.correctCapacities();
        // Correct drain reversibilities
        mcs.correctInOutflows();
        // Exclude critical reactions validated from literature from the enumeration
        // This speeds up the enumeration process 
        mcs.addNonTargets(Utilities.readLines("src/test/resources/iMM904_RP/SupportFiles/nontargets#[aerobic#glucose].txt"));
        
	}
	
	
	public static void setDesiredSpace(double glucoseUptake) {
		/*
		 * Set conditions that the MCSs are required to comply with
		 * to become constrained MCSs 
		 */
		mcs.addFluxBound(glucose, -glucoseUptake, Utilities.INF, false);
		mcs.addFluxBound(biomass, 0.001, Utilities.INF, false);
		mcs.addFluxBound(atpm, 1, Utilities.INF, true);
	}
	
	public static void setMCSFIXED(double glucoseUptake, double minimumYield) throws FileNotFoundException, IOException, XMLStreamException, ErrorsException, ParserConfigurationException, SAXException, JSBMLValidationException, IloException {
		/*
		 * Applies the MCS_FIXED formulation to the mcs object (MCSPipeline instance)
		 */
		setUpPipeline();
		mcs.addFluxBound(glucose, -glucoseUptake, -glucoseUptake, true);
		mcs.addLowerYieldConstraint(glucose, product, -minimumYield, true);
		setDesiredSpace(glucoseUptake);
	}
	
	public static void setMCSATP(double glucoseUptake, double minimumYield) throws FileNotFoundException, IOException, XMLStreamException, ErrorsException, ParserConfigurationException, SAXException, JSBMLValidationException, IloException {
		/*
		 * Applies the MCS_ATP formulation to the mcs object (MCSPipeline instance)
		 */
		setUpPipeline();
		mcs.addFluxBound(glucose, -glucoseUptake, Utilities.INF, true);
		mcs.addFluxBound(atpm, 1, Utilities.INF, true);
		mcs.addLowerYieldConstraint(glucose, product, -minimumYield, true);
		setDesiredSpace(glucoseUptake);

	}
	
	public static void main(String[] args) throws FileNotFoundException, IOException, XMLStreamException, ErrorsException, ParserConfigurationException, SAXException, JSBMLValidationException, IloException {
		// Define important constants
		double glcUp = 1.15; // glucose uptake in mmol/gDW/h
		double minYld = 0.001; // minimum product/substrate yield
		int maxSize = 6; // maximum solution size (modify with caution - avoid very large sizes)
		
		
		// Set the desired formulation - choose only one line to run! Do not uncomment both!
		setMCSATP(glcUp, minYld); String algName = "MCS_ATP";
		//setMCSFIXED(glcUp, minYld); String algName = "MCS_FIXED";

		// Enumerate the cMCSs and generate a List<List<String>> with knockouts as list of reaction IDs
		DefaultEnumerationResult results = mcs.enumerate(maxSize, true);
		List<List<String>> solutions = convertERtoList(results);
		
		// Create the filtering pipeline instances
		SolutionFilter sf = new SolutionFilter(model);
		SolutionFilterPipeline sfp = new SolutionFilterPipeline(model, cont, sf);
		
		// Define the robustness criteria
		// Check whether the product flux is above 0 when the biomass flux is at 90% of its maximum value
		RobustnessCriteria[] rcrit = new RobustnessCriteria[] {
				new RobustnessCriteria(biomass, product, 0.9, "max", "min", "90% Robustness")
		};
		
		// Environmental conditions definition
		FluxBound[] fcrit = new FluxBound[] {
				new FluxBound(model.getReaction(atpm), 1, 1),
				new FluxBound(model.getReaction(glucose), -glcUp, -glcUp),
				new FluxBound(model.getReaction(o2), -Utilities.INF, Utilities.INF)
		};
		
		// Filter solutions according to the environmental conditions and robustness criteria
		// Yields a new List of List<String>
		List<List<String>> filteredSolutions = sfp.filterSolutions(solutions, rcrit, fcrit);
		
		// Define the properties for FVA analysis
		// The product will be minimized and maximized from 0 to 100% of the maximum biomass
		// in 50 steps (incrementing each 2% of the maximum biomass)
		EnvelopeProperties[] envProp = new EnvelopeProperties[] {
				new EnvelopeProperties(product, biomass, "Product-by-biomass", 50)
		};
		
		FluxBound[][] envConds = new FluxBound[][] {fcrit};
		String[] envNames = new String[] {"aerobic"};
		
		// Create a solution analysis pipeline object with the parameters defined above
		SolutionAnalysisPipeline sap = new SolutionAnalysisPipeline(model, filteredSolutions, envProp, envConds, algName, envNames);

		// Create a SolutionScorer object to output the Metrics Overview datasets
		SolutionScorer ss = new SolutionScorer();
		
		// Check if results folder exists
		File f = new File(algName+"/");
		
		if (!f.exists()) {
			f.mkdir();
		}
		
		// Create FBA and PFBA solvers
		fba = new CPLEXFluxBalanceAnalysis(model);
		pfba = new CPLEXParsimoniousFluxBalanceAnalysis(model, 0.99999);
		
		// Update the formulas on the model object with information from the container
		// MCSPipeline does not do this by default to avoid errors!
		SteadyStateModelReader.updateFormulae(cont, model);
		
		// Create an array of IScoreItem to score the solutions
		IScoreItem[] scorers = new IScoreItem[] {
				new BPCYScoreItem(model, pfba, envConds[0], model.getReaction(biomass), model.getReaction(product), model.getReaction(glucose)), // calculate BPCY
				new RobustnessScoreItem(model, fba, fcrit,model.getReaction(biomass),model.getReaction(product), 0.1, false), // return minimum product flux value at 10% of the maximum biomass, etc...
				new RobustnessScoreItem(model, fba, fcrit,model.getReaction(biomass),model.getReaction(product), 0.5, false),
				new RobustnessScoreItem(model, fba, fcrit,model.getReaction(biomass),model.getReaction(product), 0.9, false),
				new PFBAFluxValueItem(model, pfba, fcrit,  model.getReaction(biomass), model.getReaction(biomass), "max"), // maximize biomass and return the biomass flux value
				new PFBAFluxValueItem(model, pfba, fcrit,  model.getReaction(product), model.getReaction(biomass), "max"), // maximize biomass and return the product flux value
				new CarbonYieldScoreItem(model, pfba, fcrit, model.getReaction(biomass), model.getReaction(product), model.getReaction(glucose)), // calculate the carbon yield according to the formulae
				new SolutionSizeItem() // outputs solution size
		};
		
		// Analyse and write the datasets to a specific folder
		ss.writeDataset(filteredSolutions, scorers, algName+"/MetricsOverview");
		sap.filterAnalysisPipeline(rcrit, fcrit, "Aerobic", f.getAbsolutePath()+"/");
	}
	
	public static List<List<String>> convertERtoList(DefaultEnumerationResult results) {
		List<List<String>> solutions = new ArrayList<List<String>>();
		for (String[] list : results.toStringArrays()) {
			List<String> sol = new ArrayList<String>();
			for (String string : list) {
				sol.add(string);
			}
			solutions.add(sol);
		}
		return solutions;
	}
}
