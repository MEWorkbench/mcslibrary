package pt.uminho.ceb.biosystems.mcslibrary.utilities;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

public class BashQueueGenerator {
	
	public static String PROPERTY_FVABOUND = "FVAF";
	public static String PROPERTY_ALGORITHM_TYPE = "ATY";
	public static String PROPERTY_UNDESIRED_FLUX = "UF";
	public static String PROPERTY_DESIRED_FLUX = "DF";
	public static String PROPERTY_UNDESIRED_YIELD_LOWER = "UY";
	public static String PROPERTY_UNDESIRED_YIELD_UPPER = "UYU";
	public static String PROPERTY_DESIRED_YIELD_LOWER = "DY";
	public static String PROPERTY_DESIRED_YIELD_UPPER = "DYU";
	public static String PROPERTY_MAXSIZE = "S";
	public static String PROPERTY_NONTARGETS = "E";
	public static String PROPERTY_BOUNDOVERRIDE = "B";
	public static String PROPERTY_GENEREGULATION = "G";
	public static String PROPERTY_SINGLE = "SNG";
	public static String PROPERTY_SOLUTIONEXCLUSION = "SOL";
	public static String PROPERTY_SOLUTION_SEED = "SEED";
	
	public static String BASH_JAVAEXEC = "jexec";
	public static String BASH_VM = "jvmargs";
	public static String BASH_JARFILE = "jarfile";
	public static String BASH_THREADS = "threads";
	public static String BASH_QUEUECOMMAND = "queuecommand";

	public static String getBashString(Map<String,String> bashParams, String modelName, String problemName, String resultsName, String enumerationParamsName){
		String str = bashParams.get(BASH_JAVAEXEC) + " " + 
					 bashParams.get(BASH_VM) + " -jar " +
					 bashParams.get(BASH_JARFILE) + " ";
		str += modelName + " " + problemName + " " + resultsName + " " + bashParams.get(BASH_THREADS);
		
		return str;
	}
	
	public static String generateProblemString(Map<String,List<String>> problem){
		String str = "";
		for (String string : problem.keySet()) {
			List<String> params = problem.get(string);
			for (int i = 0; i < params.size(); i++) {
				String param = params.get(i);
				if (param.length() > 0) {
					str += string+"!"+param+"\n";
				}
			}
		}
		return str.substring(0, str.length()-1);
	}
	
	public static String generateModelFile(String modelPath, String spont){
		return modelPath + "\n" + spont;
	}
	
	public static void generateAllQueues(String rootdir, Map<String,String> bash, List<Map<String,List<String>>> problem, String modelPath, String spont, String resultsName, String enumerationParams) throws IOException{
		String model = generateModelFile(modelPath, spont);
		String modelName = modelPath+".mcsmodel";
		Utilities.write(model, rootdir+modelName);
		
		File proFile = new File(rootdir+"Problems");
		try {
			proFile.mkdir();
		} catch (Exception e) {
			// TODO: handle exception
		}
		
		String queueFile = "";

		for (int i = 0; i < problem.size(); i++) {
			
			String problemString = generateProblemString(problem.get(i));
			String problemName = "Problems/problem_"+i;
			Utilities.write(problemString, rootdir+problemName);
			
			String bashString = getBashString(bash, modelName, problemName, resultsName+"_problem_"+i+"_results", enumerationParams);
			String bashName = rootdir+"bash_"+i+".sh";
			queueFile += bash.get(BashQueueGenerator.BASH_QUEUECOMMAND)+" "+"bash_"+i+".sh\n";
			Utilities.write(bashString, bashName);

		}
		
		Utilities.write(queueFile, rootdir+resultsName+"_queue.sh");
		
	}
	
}
