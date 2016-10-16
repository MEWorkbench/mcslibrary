/*******************************************************************************
 * Copyright 2016
 * CEB Centre of Biological Engineering
 * University of Minho
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This code is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this code. If not, see http://www.gnu.org/licenses/
 *
 * Created inside the BIOSYSTEMS Research Group
 * (http://www.ceb.uminho.pt/biosystems)
 *******************************************************************************/
package pt.uminho.ceb.biosystems.mcslibrary.utilities;

import ilog.concert.IloException;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.StringTokenizer;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.stream.XMLStreamException;
import org.xml.sax.SAXException;

import pt.uminho.ceb.biosystems.mcslibrary.enumeration.Solution;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.implementation.DefaultEnumerationResult;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.Container;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.io.readers.ErrorsException;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.io.readers.JSBMLReader;
import pt.uminho.ceb.biosystems.mew.biocomponents.validation.io.JSBMLValidationException;
import pt.uminho.ceb.biosystems.mew.core.model.exceptions.InvalidSteadyStateModelException;

public class MainClass {
	public static void main(String[] args) throws FileNotFoundException, XMLStreamException, ErrorsException, IOException, InvalidSteadyStateModelException, IloException, ParserConfigurationException, SAXException, JSBMLValidationException {
		if (args.length == 0) {
			args = new String[4];
			args[0] = "/home/skapur/Dropbox/modeltest/tests/imm904.mcsmodel";
			args[1] = "/home/skapur/Dropbox/modeltest/tests/succinate_production.mcsproblem";
			args[2] = "iaf1260_Ethanol";
			args[3] = "8";
		}

		List<String> modelconfs = Utilities.readLines(args[0]);
		String modelpath = modelconfs.get(0);
		String spontaneousId = modelconfs.get(1);

		Container cont = new Container(new JSBMLReader(modelpath,"a",false));
		MCSPipeline pipeline = new MCSPipeline(cont, spontaneousId);

		String algorithmType = null;
		String boundOverride = null;
		String geneRegulation = null;
		String nonTargets = null;
		BufferedReader br = new BufferedReader(new FileReader(args[1]));
		int maxsize = 1;

		while (br.ready()) {
			String line = br.readLine();
			System.out.println(br.ready());
			StringTokenizer stk = new StringTokenizer(line,"!");
			String type = stk.nextToken();
			String values = stk.nextToken();
			System.out.println("[CONFIGURATION]: "+type+" += "+values);
			switch (type) {
			
			case "ATY":
					{
					algorithmType = values;
					break;
					}
			case "UF":
					{
					StringTokenizer stkv = new StringTokenizer(values,";");
					System.out.println(values);
					String name = stkv.nextToken();
					String lb = stkv.nextToken();
					System.out.println(lb);
					double newlb = 0;
					if (lb.equals("i")) {
						newlb = Utilities.INF;
					} else if (lb.equals("-i")){
						newlb = - Utilities.INF;
					} else {
						newlb = Double.parseDouble(lb);
					}
	
					String ub = stkv.nextToken();
					System.out.println(ub);
					double newub = 0;
					if (ub.equals("i")) {
						newub = Utilities.INF;
					} else if (ub.equals("-i")){
						newub = - Utilities.INF;
					} else {
						newub = Double.parseDouble(ub);
					}
	
					pipeline.addFluxBound(name, newlb, newub, true);
					break;
					}
					
			case "FVAF":
					{
						StringTokenizer stkv = new StringTokenizer(values,";");
						System.out.println(values);
						String name = stkv.nextToken();
		
						String lb = stkv.nextToken();
						System.out.println(lb);
						double newlb = 0;
						if (lb.equals("i")) {
							newlb = Utilities.INF;
						} else if (lb.equals("-i")){
							newlb = - Utilities.INF;
						} else {
							newlb = Double.parseDouble(lb);
						}
		
						String ub = stkv.nextToken();
						System.out.println(ub);
						double newub = 0;
						if (ub.equals("i")) {
							newub = Utilities.INF;
						} else if (ub.equals("-i")){
							newub = - Utilities.INF;
						} else {
							newub = Double.parseDouble(ub);
						}
		
						pipeline.addFVABound(name, newlb, newub);
						break;
					}
					
			case "DF":
					{
					StringTokenizer stkv = new StringTokenizer(values,";");
					System.out.println(values);
					String name = stkv.nextToken();
	
					String lb = stkv.nextToken();
					System.out.println(lb);
					double newlb = 0;
					if (lb.equals("i")) {
						newlb = Utilities.INF;
					} else if (lb.equals("-i")){
						newlb = - Utilities.INF;
					} else {
						newlb = Double.parseDouble(lb);
					}
	
					String ub = stkv.nextToken();
					System.out.println(ub);
					double newub = 0;
					if (ub.equals("i")) {
						newub = Utilities.INF;
					} else if (ub.equals("-i")){
						newub = - Utilities.INF;
					} else {
						newub = Double.parseDouble(ub);
					}
	
					pipeline.addFluxBound(name, newlb, newub, false);
					break;
					}

			case "UY":
				{
					StringTokenizer stkv = new StringTokenizer(values,";");
					String substrate = stkv.nextToken();
					String product = stkv.nextToken();
					double ratio = Double.parseDouble(stkv.nextToken());
					double deviation = 0;
					if (stkv.hasMoreTokens()) {
						deviation = Double.parseDouble(stkv.nextToken());
					}
					pipeline.addLowerYieldConstraint(substrate, product, ratio, true, deviation);
					break;
				}

			case "UYU":
				{
					StringTokenizer stkv = new StringTokenizer(values,";");
					String substrate = stkv.nextToken();
					String product = stkv.nextToken();
					double ratio = Double.parseDouble(stkv.nextToken());
					double deviation = 0;
					if (stkv.hasMoreTokens()) {
						deviation = Double.parseDouble(stkv.nextToken());
					}
					pipeline.addUpperYieldConstraint(substrate, product, ratio, true, deviation);
					break;
				}

			case "DY":
				{
					StringTokenizer stkv = new StringTokenizer(values,";");
					String substrate = stkv.nextToken();
					String product = stkv.nextToken();
					double ratio = Double.parseDouble(stkv.nextToken());
					double deviation = 0;
					if (stkv.hasMoreTokens()) {
						deviation = Double.parseDouble(stkv.nextToken());
					}
					pipeline.addLowerYieldConstraint(substrate, product, ratio, false, deviation);
					break;
				}

			case "DYU":
				{
					StringTokenizer stkv = new StringTokenizer(values,";");
					String substrate = stkv.nextToken();
					String product = stkv.nextToken();
					double ratio = Double.parseDouble(stkv.nextToken());
					double deviation = 0;
					if (stkv.hasMoreTokens()) {
						deviation = Double.parseDouble(stkv.nextToken());
					}
					pipeline.addUpperYieldConstraint(substrate, product, ratio, false, deviation);
					break;
				}
				
			case "S":
				{
					maxsize = Integer.parseInt(values);
					break;
				}
				
			case "E":
				{
					nonTargets = values;
//					pipeline.ignoreReaction(values);
					break;
				}
				
			case "B":
				{
					boundOverride = values;
					break;
				}
				
			case "G":
				{
					geneRegulation = values;
					break;
				}
				
			case "SNG":
				{
					StringTokenizer tok = new StringTokenizer(values,";");
					String r = tok.nextToken();
					pipeline.addSingleReaction(r);
					break;
				}
				
			case "SOL":
				{
					StringTokenizer tok = new StringTokenizer(values,";");
					int size = tok.countTokens();
					int[] sol = new int[size];
					for (int i = 0; i < size; i++) {
						sol[i] = pipeline.getMetabolicNetwork().containsReaction(tok.nextToken());
					}
					pipeline.addExcludedSolution(new Solution(sol, pipeline.getMetabolicNetwork()));
					break;
				}
				
			case "SEED":
				{
					StringTokenizer tok = new StringTokenizer(values,";");
					int size = tok.countTokens();
					int[] sol = new int[size];
					for (int i = 0; i < size; i++) {
						sol[i] = pipeline.getMetabolicNetwork().containsReaction(tok.nextToken());
					}
					pipeline.addSolutionSeed(new Solution(sol, pipeline.getMetabolicNetwork()));
					break;
				}
			}
		}
		br.close();

		if (boundOverride != null)
				pipeline.overrideBounds(boundOverride);

		pipeline.correctCapacities();
		pipeline.correctInOutflows();
		if (algorithmType != null) {
			if (algorithmType.contains("INTEGRATED")) {
				System.out.println("Algorithm set: Integrated");
				pipeline.setAlgorithmType(AlgorithmType.INTEGRATED);
			} else if (algorithmType.contains("DEFAULT")){
				System.out.println("Algorithm set: Default");
				pipeline.setAlgorithmType(AlgorithmType.DEFAULT);
			} else if (algorithmType.contains("SAMPLING")) {
				System.out.println("Algorithm set: Sampling");
				pipeline.setAlgorithmType(AlgorithmType.SAMPLING);
			}
		}
		
		if (nonTargets != null) {
			pipeline.addNonTargets(Utilities.readLines(nonTargets));
		}
		if (geneRegulation != null)
				pipeline.applyGeneRegulation(Utilities.readLines(geneRegulation));

		
		DefaultEnumerationResult re = null;
		int threads = Integer.parseInt(args[3]);
		System.out.println("Seeding problem? "+pipeline.isSeedingProblem());
		if (pipeline.isSeedingProblem()) {
//			List<String[]> seedList = pipeline.getSeedList();
//			SolutionSeedingPipeline ssp = new SolutionSeedingPipeline();
//			DefaultEnumerationResult[] rel = ssp.enumerate(seedList, pipeline, maxsize);
		} else {
			pipeline.setThreads(threads);
			re = pipeline.enumerate(maxsize, pipeline.isConstrainedProblem());
		}
		System.out.println(re.countResults());
		re.writeToFile(args[2]+"_"+System.currentTimeMillis()+"mcs.txt");
		System.out.println("FINISHED!");
		}
	}
