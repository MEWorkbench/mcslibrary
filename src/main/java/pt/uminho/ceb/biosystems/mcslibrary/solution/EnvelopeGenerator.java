package pt.uminho.ceb.biosystems.mcslibrary.solution;

import ilog.concert.IloException;
import java.io.FileNotFoundException;
import java.io.IOException;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.stream.XMLStreamException;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.xml.sax.SAXException;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba.CPLEXFluxBalanceAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba.FluxBalanceAnalysisResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.ChartUtilities;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.MCSPipeline;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.ScatterPlot;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.Container;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.io.readers.ErrorsException;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.io.readers.JSBMLReader;
import pt.uminho.ceb.biosystems.mew.biocomponents.validation.io.JSBMLValidationException;
import pt.uminho.ceb.biosystems.mew.core.model.exceptions.InvalidSteadyStateModelException;

public class EnvelopeGenerator {
	private DefaultMetabolicNetwork model;
	private CPLEXFluxBalanceAnalysis fba;
	private FluxBound[] env;

	public EnvelopeGenerator(DefaultMetabolicNetwork model, FluxBound[] env) throws IloException {
		this.model = model;
		this.fba = new CPLEXFluxBalanceAnalysis(model);
		this.env = env;
	}
	
	public ScatterPlot getProductionEnvelope(String title, String product, String biomass, String[] ko, int breaks) throws IloException {
		FluxBound[] fbf = new FluxBound[ko.length+env.length+1];
		Reaction biomassR = model.getReaction(biomass);
		Reaction productR = model.getReaction(product);
		for (int i = 0; i < ko.length; i++) {
			fbf[i] = new FluxBound(model.getReaction(ko[i]),0,0);
		}
		
		for (int i = 0; i < env.length; i++) {
			fbf[i+ko.length] = env[i];
		}
		
		fbf[ko.length+env.length] = new FluxBound(biomassR, 0, Utilities.INF);

		XYSeries max = new XYSeries("Max");
		XYSeries min = new XYSeries("Min");
		XYSeriesCollection coll = new XYSeriesCollection();
		coll.addSeries(min);
		coll.addSeries(max);
		
		double maxBiomass = fba.solve(fbf, biomassR, "max").getFluxValue(biomassR);
		FluxBalanceAnalysisResult maxdn = fba.solve(fbf, productR, "max");
		FluxBalanceAnalysisResult mindn = fba.solve(fbf, productR, "min");
		double pmaxn = maxdn.getFluxValue(productR);
		double pminn = mindn.getFluxValue(productR);
		max.add(0, pmaxn);
		min.add(0, pminn);
		
		for (int i = 1; i < breaks; i++) {
			double thres = (i * maxBiomass)/breaks;
			fbf[ko.length+env.length] = new FluxBound(biomassR, thres, Utilities.INF);
			FluxBalanceAnalysisResult maxd = fba.solve(fbf, productR, "max");
			FluxBalanceAnalysisResult mind = fba.solve(fbf, productR, "min");
			double pmax = maxd.getFluxValue(productR);
			double pmin = mind.getFluxValue(productR);
			max.add(thres, pmax);
			min.add(thres, pmin);
		}
		
		FluxBalanceAnalysisResult maxdb = fba.solve(fbf, biomassR, "max");
		double pmaxb = maxdb.getFluxValue(productR);
		max.add(maxBiomass, pmaxb);
		min.add(maxBiomass, pmaxb);
		
		ScatterPlot sp = new ScatterPlot(title,"Production Envelope", biomass, product, coll);
		return sp;
	}
	
	public static void main(String[] args) throws FileNotFoundException, IOException, XMLStreamException, ErrorsException, ParserConfigurationException, SAXException, JSBMLValidationException, InvalidSteadyStateModelException, IloException {
		Container cont = new Container(new JSBMLReader("/home/skapur/Downloads/iMM904_sl.xml","a",false));
		MCSPipeline p = new MCSPipeline(cont, "");
		DefaultMetabolicNetwork metaNet = p.getMetabolicNetwork();
		EnvelopeGenerator g = new EnvelopeGenerator(metaNet, new FluxBound[]{});
		
		ScatterPlot penv = g.getProductionEnvelope("SAMPLE", "R_EX_succ_e_", "R_biomass_SC5_notrace", new String[]{"R_SUCD2_u6m","R_SUCD3_u6m","R_G6PDH2","R_PYRDC","R_ACOAHim"}, 100);
		ChartUtilities.displayPlot(penv);
	}
}

