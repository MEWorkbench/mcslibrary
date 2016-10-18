package mcslibrary;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.stream.XMLStreamException;

import org.xml.sax.SAXException;

import cern.colt.Arrays;
import ilog.concert.IloException;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.implementation.DefaultEnumerationResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.ReactionConstraint;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba.CPLEXFluxBalanceAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba.FluxBalanceAnalysisResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.MCSPipeline;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.Container;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.io.readers.ErrorsException;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.io.readers.JSBMLReader;
import pt.uminho.ceb.biosystems.mew.biocomponents.validation.io.JSBMLValidationException;

public class YeastTests {
	static MCSPipeline mcs;
	static CPLEXFluxBalanceAnalysis fba;
	static DefaultMetabolicNetwork model;
	private static Reaction biomassReaction;
	
	final static String atpm = "R_ATPM";
	final static String biomass = "R_biomass_SC5_notrace";
	final static String glucose = "R_EX_glc_e_";
	final static String o2 = "R_EX_o2_e_";

	public static void setUpPipeline() throws FileNotFoundException, IOException, XMLStreamException, ErrorsException, ParserConfigurationException, SAXException, JSBMLValidationException, IloException{
		String path = "iMM904_RP/iMM904_corrected_201609001170165318.xml";
		Container cont = new Container(new JSBMLReader(path, "default", false));
		mcs = new MCSPipeline(cont, "");
		model = mcs.getMetabolicNetwork();
		biomassReaction = model.getReaction(biomass);
		
		model.getReaction(atpm).setBounds(new ReactionConstraint(1, 1));
		model.getReaction(o2).setBounds(new ReactionConstraint(-Utilities.INF, Utilities.INF));
        model.getReaction(glucose).setBounds(new ReactionConstraint(-1.15, Utilities.INF));
		
		mcs.correctCapacities();
        // correct drain reversibilities
        mcs.correctInOutflows();
	}
	
	public static void main(String[] args) throws FileNotFoundException, IOException, XMLStreamException, ErrorsException, ParserConfigurationException, SAXException, JSBMLValidationException, IloException {
		setUpPipeline();
		
		fba = new CPLEXFluxBalanceAnalysis(model);
		FluxBalanceAnalysisResult simulationWildType = fba.solve(null, biomassReaction, "max");
		System.out.println("Wild-type biomass: "+simulationWildType.getFluxValue(biomassReaction));
		
		String product = "R_EX_succ_e_";
		int maximumSize = 3;
		runYeastTest(product, maximumSize);
	}
	
	public static void runYeastTest(String product, int maximumSize) throws IOException, IloException{
        // set remove absolute bounds from unbounded reactions
		Reaction productReaction = model.getReaction(product);
        
        // reactions to exclude from compression
        mcs.addSingleReaction(biomass);
        mcs.addSingleReaction(atpm);
        mcs.addSingleReaction(glucose);
        mcs.addSingleReaction(o2);
        mcs.addSingleReaction(product);
        
        List<String> nontargets = Utilities.readLines("iMM904_RP/SupportFiles/nontargets#[aerobic#glucose].txt");
        mcs.addNonTargets(nontargets);
        
        // parameters
        double minSYield = 0.01; // minimum Product/Substrate yield
        double glucoseUptake = 1.15; // minimum Substrate yield
        double viableBiomass = 0.0001; // minimum Biomass for filtering
        
        // add a flux bound to the intervention problem
        mcs.addFluxBound(
        		glucose, // which flux
        		-glucoseUptake, // lower bound
        		Utilities.INF, // upper bound
        		true); // is it an undesired flux?
        
        
        // add a yield constraint to the intervention problem
        mcs.addLowerYieldConstraint(
        		glucose, // denominator flux (substrate)
        		product, // numerator flux (product)
        		-minSYield, // yield threshold (LowerYield -> Y >= x; UpperYield -> Y <= x)
        		true); // is it an undesired constraint?

        mcs.addFluxBound(atpm, 1, 1, true);
        
        mcs.addFluxBound(biomass, viableBiomass, Utilities.INF, false);
        mcs.addFluxBound(glucose, -glucoseUptake, -glucoseUptake, false);
        mcs.addUpperYieldConstraint(glucose, product, -minSYield, false);
        mcs.addFluxBound(atpm, 1, 1, false);

        
        // add non targets (list of reaction IDs to exclude from the enumeration - speeds up)
//        mcs.addNonTargets(Utilities.readLines("iMM904_files/nontargets.txt"));
        
        
        DefaultEnumerationResult r = mcs.enumerate(maximumSize, true);
        printSolutionDataset(r, productReaction);
        // any solution with minimum product of 0 is not strongly coupled
	}
	
	public static void printSolutionDataset(DefaultEnumerationResult r, Reaction productReaction) throws IloException{
		FluxBound[] environmentalConditions = new FluxBound[]{
				new FluxBound(model.getReaction(atpm), 1, 1), 
				new FluxBound(model.getReaction(glucose), -1.15, -1.15),
		};
		
		ArrayList<String[]> strings = r.toStringArrays();
		System.out.println("Solution,MaximumBiomass,MinimumProductFlux");
		for (int i = 0; i < r.countResults(); i++) {
			String solutionContent = Arrays.toString(strings.get(i)).replace(",", " ");
			Reaction[] knockouts = Utilities.toReacArrayFromInt(model, r.getResult(i));
			FluxBalanceAnalysisResult simulationMaxBiomass = fba.solveReactionKnockoutFBA(environmentalConditions, knockouts, biomassReaction, "max");
			FluxBalanceAnalysisResult simulationMinProduct = fba.solveReactionKnockoutFBA(environmentalConditions, knockouts, productReaction, "min");

			System.out.println(solutionContent+","+simulationMaxBiomass.getFluxValue(biomassReaction)+","+simulationMinProduct.getFluxValue(productReaction));
		}
	}
}
