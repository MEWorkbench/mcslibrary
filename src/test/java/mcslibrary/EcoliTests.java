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
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba.CPLEXFluxBalanceAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.MCSPipeline;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.Container;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.io.readers.ErrorsException;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.io.readers.JSBMLReader;
import pt.uminho.ceb.biosystems.mew.biocomponents.validation.io.JSBMLValidationException;

public class EcoliTests {
	
	static MCSPipeline mcs;
	static CPLEXFluxBalanceAnalysis fba;
	
	public static void setUpPipeline() throws FileNotFoundException, IOException, XMLStreamException, ErrorsException, ParserConfigurationException, SAXException, JSBMLValidationException, IloException{
		String path = "models/iAF1260_flux2.xml";
		Container cont = new Container(new JSBMLReader(path, "default", false));
		mcs = new MCSPipeline(cont, "s0001");
		fba = new CPLEXFluxBalanceAnalysis(mcs.getMetabolicNetwork());
	}
	
	public static void main(String[] args) throws FileNotFoundException, IOException, XMLStreamException, ErrorsException, ParserConfigurationException, SAXException, JSBMLValidationException, IloException {
		
		setUpPipeline();		// do not comment this line when running any of the examples
		
		/* 
		 * each of these functions take an integer to define the maximum solution size.
		 * feel free to edit these when running these examples or comment ones you don't need
		 * be advised: computation times scale exponentially with solution size.
		 * 
		 */
		
		runSyntheticLethalsProblem(1); // should not exceed 4
		runEthanolProductionProblem(3);
		runFumarateProductionProblem(6);
		runSerineProductionProblem(7);
	}
	

	public static void runSyntheticLethalsProblem(int maximumSize){
		
		mcs.resetConstraints(); // remove any previous constraints
		
        mcs.correctCapacities(); // change capacities over M (defined on MCS pipeline) to infinity
        mcs.correctInOutflows(); // correct reaction directionality on exchange reactions

        // exclude the biomass reaction from being compressed
        mcs.addSingleReaction("R_Ec_biomass_iAF1260_core_59p81M"); 
        
        // define the undesired space (3rd argument as true for solution space to block)
        mcs.addFluxBound("R_EX_glc_e_", -18.5, Utilities.INF, true); // substrate uptake
        mcs.addFluxBound("R_ATPM", 8.39, Utilities.INF, true); // maintenance ATP restrictions
        mcs.addFluxBound("R_Ec_biomass_iAF1260_core_59p81M", 0.001, Utilities.INF, true); // block vectors with biomass above 1e-3
        
        // define the desired space (3rd argument as false for solution space to allow)
        mcs.addFluxBound("R_EX_glc_e_", -18.5, Utilities.INF, false); // substrate uptake
        mcs.addFluxBound("R_ATPM", 8.39, Utilities.INF, false); // maintenance ATP restrictions

        // enumerate minimal cut sets 
        // 1st argument sets the maximum solution size
        // 2nd argument defines whether the solutions will be filtered according to the problem
        // the latter concerns the desired space
        DefaultEnumerationResult res = mcs.enumerate(maximumSize,false);
        
        printEnumerationResult(res);
	}
	
	public static void runEthanolProductionProblem(int maximumSize) throws IOException{
		
		mcs.resetConstraints(); // remove any previous constraints

        String substrate = "R_EX_glc_e_"; // glucose exchange reaction name
        String product = "R_EX_etoh_e_"; // ethanol exchange reaction name
        String oxygen = "R_EX_o2_e_"; // oxygen exchange reaction name
        String biomass = "R_Ec_biomass_iAF1260_core_59p81M"; // biomass reaction name
        String atpm = "R_ATPM"; // maintenance ATP reaction name
        
        double glcUpt = 20; // glucose (substrate) uptake
        double oxyUpt = 20; // oxygen uptake
        double mATP = 8.39; // minimum maintenance ATP flux
        double minYld = 0.5; // minimum product/substrate yield threshold
        double minBio = 0.001; // minimum biomass threshold
        
        // bounds specific for aerobic conditions (needed to replicate MCSEnumerator)
        mcs.overrideBounds("iAF1260_files/bounds.csv");
        
        // reactions that are not to be compressed
        String[] singleReactions = new String[]{substrate, product, oxygen, biomass, atpm};
        
        for (String reaction : singleReactions)
        	mcs.addSingleReaction(reaction);

        
        mcs.correctCapacities();
        
        // block reactions contained in the file
        // needed to replicate MCSEnumerator validation cases
		List<String> grm = Utilities.readLines("iAF1260_files/aerobic.grm");
        mcs.applyGeneRegulation(grm);
        
        mcs.correctInOutflows();
        
        // eliminates reactions that are blocked if this constraint is forced upon the model
        mcs.addFVABound("R_Ec_biomass_iAF1260_core_59p81M", 0.001, Utilities.INF);

        // define the undesired space
        mcs.addFluxBound(atpm, mATP, Utilities.INF, true);
		mcs.addFluxBound(substrate, -glcUpt, Utilities.INF, true); 
        mcs.addFluxBound(oxygen, -oxyUpt, Utilities.INF, true);
        mcs.addLowerYieldConstraint(substrate, product, -minYld, true); // block vectors with Substrate/Product yield below minYld
        // addLowerYieldConstraint for Y(P/S) > -minYld (glucose uptake carries negative flux so the sign is flipped)
        
        // define the desired space
        mcs.addFluxBound(biomass, minBio, Utilities.INF, false);
        mcs.addUpperYieldConstraint(substrate, product, -minYld, false);
        mcs.addFluxBound(oxygen, -oxyUpt, Utilities.INF, false);
        mcs.addFluxBound(substrate, -glcUpt, Utilities.INF, false);
        mcs.addFluxBound(atpm, mATP, Utilities.INF, false);
        // addUpperYieldConstraint for Y(P/S) > -minYld (glucose uptake carries negative flux so the sign is flipped)

        DefaultEnumerationResult res = mcs.enumerate(maximumSize,true);
        
        printEnumerationResult(res);
	}
	
	public static void runFumarateProductionProblem(int maximumSize) throws IOException{
		
		mcs.resetConstraints(); // remove any previous constraints

        mcs.overrideBounds("iAF1260_files/bounds.csv");
        mcs.correctCapacities();
        mcs.applyGeneRegulation(Utilities.readLines("iAF1260_files/aerobicoxidativestress.grm"));
        mcs.correctInOutflows();

        mcs.addFluxBound("R_EX_glc_e_", -20, Utilities.INF, true);
        mcs.addFluxBound("R_ATPM", 8.39, Utilities.INF, true);
        mcs.addFluxBound("R_EX_o2_e_", -20, Utilities.INF, true);
        mcs.addLowerYieldConstraint("R_EX_glc_e_","R_EX_fum_e_", -0.5, true);

        mcs.addFVABound("R_Ec_biomass_iAF1260_core_59p81M", 0.001, Utilities.INF);

        mcs.addFluxBound("R_Ec_biomass_iAF1260_core_59p81M", 0.001, Utilities.INF, false);
        mcs.addFluxBound("R_EX_glc_e_", -20, Utilities.INF, false);
        mcs.addFluxBound("R_ATPM", 8.39, Utilities.INF, false);
        mcs.addFluxBound("R_EX_o2_e_", -20, Utilities.INF, false);
        
        DefaultEnumerationResult res = mcs.enumerate(maximumSize,true);
        
        printEnumerationResult(res);
	}
	
	public static void runSerineProductionProblem(int maximumSize) throws IOException{
		
		mcs.resetConstraints(); // remove any previous constraints

        mcs.overrideBounds("iAF1260_files/bounds.csv");
        mcs.correctCapacities();
        mcs.applyGeneRegulation(Utilities.readLines("iAF1260_files/aerobicoxidativestress.grm"));
        mcs.correctInOutflows();

        mcs.addFluxBound("R_EX_glc_e_", -20, Utilities.INF, true);
        mcs.addFluxBound("R_ATPM", 8.39, Utilities.INF, true);
        mcs.addFluxBound("R_EX_o2_e_", -20, Utilities.INF, true);
        mcs.addLowerYieldConstraint("R_EX_glc_e_","R_EX_ser_L_e_", -0.5, true);

        mcs.addFVABound("R_Ec_biomass_iAF1260_core_59p81M", 0.001, Utilities.INF);

        mcs.addFluxBound("R_Ec_biomass_iAF1260_core_59p81M", 0.001, Utilities.INF, false);
        mcs.addFluxBound("R_EX_glc_e_", -20, Utilities.INF, false);
        mcs.addFluxBound("R_ATPM", 8.39, Utilities.INF, false);
        mcs.addFluxBound("R_EX_o2_e_", -20, Utilities.INF, false);
        
        DefaultEnumerationResult res = mcs.enumerate(maximumSize,true);
        
        printEnumerationResult(res);
	}
	
	public static void printEnumerationResult(DefaultEnumerationResult res){
        // get string arrays corresponding to the MCSs
        ArrayList<String[]> resultsAsStrings = res.toStringArrays();
        
        for (int i = 0; i < resultsAsStrings.size(); i++) {
			System.out.println("Solution "+i+" = "+Arrays.toString(resultsAsStrings.get(i)));
		}
	}
	
}
