package pt.uminho.ceb.biosystems.mcslibrary.solution.analysis;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.stream.XMLStreamException;

import org.xml.sax.SAXException;

import ilog.concert.IloException;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.solution.SolutionUtilities;
import pt.uminho.ceb.biosystems.mcslibrary.solution.analysis.simulation.components.EnvelopeProperties;
import pt.uminho.ceb.biosystems.mcslibrary.solution.converter.RKtoGKConverter;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.MCSPipeline;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Pair;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.Container;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.io.readers.ErrorsException;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.io.readers.JSBMLReader;
import pt.uminho.ceb.biosystems.mew.biocomponents.validation.io.JSBMLValidationException;
import pt.uminho.ceb.biosystems.mew.utilities.java.StringUtils;

public class AnalysisFileManager {
	
	public final String EXTENSION_FILTER = "fcrit";
	public final String EXTENSION_ENVCON = "econd";
	public final String EXTENSION_ENVPRO = "eprop";
	public final String EXTENSION_FILTRB = "frobu";


	public final String[] PARAM_EXTS = new String[]{EXTENSION_FILTER, EXTENSION_ENVCON, EXTENSION_ENVPRO, EXTENSION_FILTRB};
	
	private String root;
	private Map<String, File> filemap;
	private Map<String, File> solutionmap;
	private Map<String, File> parammap;
	
	private FluxBound[][] filterCriteria;
	private String[] filterNames;
	
	private FluxBound[][] environmentalConditions;
	private String[] envNames;

	private EnvelopeProperties[] envelopeProperties;
	private DefaultMetabolicNetwork mn;
	
	private boolean willConvert;
	private Container cont;
	private RobustnessCriteria[] robCriteria;
	
	public AnalysisFileManager(String rootFolder, String modelPath, boolean convertToGenes) throws FileNotFoundException, IOException, XMLStreamException, ErrorsException, ParserConfigurationException, SAXException, JSBMLValidationException {
		this.root = rootFolder;
		this.cont = new Container(new JSBMLReader(modelPath,"a",false));
		this.mn = new MCSPipeline(cont,"").getMetabolicNetwork();
		this.willConvert = convertToGenes;
		buildFileMap();
		addAnalysisFolder();
		retrieveParameters();
	}
	
	private void addAnalysisFolder(){
		File f = new File(root+"/Analysis");
		if (!f.exists()) {
			f.mkdirs();
		}
	}
	
	private void buildFileMap(){
		File f = new File(root);
		File[] files = f.listFiles();
		
		filemap = new HashMap<>();
		solutionmap = new HashMap<>();
		parammap = new HashMap<>();

		for (int i = 0; i < files.length; i++) {
			List<String> parts = Utilities.getAllStringTokens(files[i].getName(), ".");
			String extension = parts.size() > 1 ? parts.get(parts.size() - 1) : ".";
			System.out.println("File found with "+extension+" extension. ("+parts+").");
			if (!files[i].isDirectory()) {
				filemap.put(files[i].getName().replaceAll("[.].*", ""), files[i]);
			}
//			System.out.println(Arrays.binarySearch(PARAM_EXTS, extension)+" "+parts);
			boolean found = false;
			for (int j = 0; j < PARAM_EXTS.length; j++) {
				if (PARAM_EXTS[j].equals(extension)) {
					found = true;
					break;
				}
			}
			if (found) {
				parammap.put(parts.get(0), files[i]);
			} else {
				parts.remove(parts.get(parts.size()-1));
				if (parts.size() > 0) {
					solutionmap.put(StringUtils.concat("", parts), files[i]);
					System.out.println("Solution file added ("+parts+").");
				}

			}
		}
		
		System.out.println(filemap.size()+" total files. "+parammap.size()+" parameter files found. "+solutionmap.size()+" solution files found.");
	}
	
	private void retrieveParameters() throws IOException{
		List<EnvelopeProperties> eps = new ArrayList<>();
		List<Pair<String,FluxBound[]>> ecs = new ArrayList<>();
		List<Pair<String,FluxBound[]>> fcs = new ArrayList<>();
		List<RobustnessCriteria> rf = new ArrayList<>();

		for (String file : parammap.keySet()) {
			File f = parammap.get(file);
			List<String> parts = Utilities.getAllStringTokens(f.getName(), ".");
			String ext = parts.get(parts.size()-1);
			switch (ext) {
			case EXTENSION_ENVPRO:
			{
				System.out.println("Found envelope properties file at "+f.getAbsolutePath());
				eps.add(EnvelopeProperties.fromFile(f.getAbsolutePath()));
			}
				break;
			case EXTENSION_FILTER:
			{
				System.out.println("Found filter criteria file at "+f.getAbsolutePath());
				fcs.add(FluxBound.fromFile(f.getAbsolutePath(), mn));
			}
				break;
			case EXTENSION_ENVCON:
			{
				System.out.println("Found environmental conditions file at "+f.getAbsolutePath());
				ecs.add(FluxBound.fromFile(f.getAbsolutePath(), mn));
			}
				break;
			case EXTENSION_FILTRB:
			{
				rf.add(RobustnessCriteria.fromFile(f.getAbsolutePath(), mn));
			}
				break;
			default:
				break;
			}
		}
		filterCriteria = new FluxBound[fcs.size()][];
		environmentalConditions = new FluxBound[ecs.size()][];
		envelopeProperties = new EnvelopeProperties[eps.size()];
		robCriteria = new RobustnessCriteria[rf.size()];
		
		filterNames = new String[fcs.size()];
		envNames = new String[ecs.size()];
		
		for (int i = 0; i < fcs.size(); i++){
			filterCriteria[i] = fcs.get(i).getB();
			filterNames[i] = fcs.get(i).getA();
		}
		for (int i = 0; i < ecs.size(); i++){
			environmentalConditions[i] = ecs.get(i).getB();
			envNames[i] = ecs.get(i).getA();
		}
		for (int i = 0; i < eps.size(); i++)
			envelopeProperties[i] = eps.get(i);
		
		for (int i = 0; i < rf.size(); i++)
			robCriteria[i] = rf.get(i);
	}
	
	public void run() throws Exception{
		for (String fname : solutionmap.keySet()) {
			File file = solutionmap.get(fname);
//			RKtoGKConverter rkgk = new RKtoGKConverter(cont);
			SolutionAnalysisPipeline sap = new SolutionAnalysisPipeline(mn, SolutionUtilities.tokenizeFile(file.getAbsolutePath(), ","), envelopeProperties, environmentalConditions, fname, envNames);
			for (int i = 0; i < filterCriteria.length; i++) {
				System.out.println("Analysing... "+fname+" solution file with filtering criteria "+filterNames[i]+".");
				Map<List<String>, SolutionAnalysisResult> fap = sap.filterAnalysisPipeline(robCriteria, filterCriteria[i], filterNames[i], root+"/Analysis/");
				List<List<String>> currentSolutions = new ArrayList<>(fap.keySet());
				if (willConvert) {
					File f = new File(root+"/Analysis/GeneKO/");
					f.mkdirs();
//					List<List<String>> geneKOs = rkgk.getMinimalGeneKOs(currentSolutions);
//					System.out.println(geneKOs);
					sap.geneKnockoutAnalysis(currentSolutions, fname+"#GENEKO", envNames, filterNames[i], filterCriteria[i], root+"/Analysis/GeneKO/", cont);
				}
			}
		}
	}
	
	public void runFilter(boolean printDataset){
		for (String fname : solutionmap.keySet()) {
			File f = solutionmap.get(fname);
			
		}
	}
	
	
	
	
}
