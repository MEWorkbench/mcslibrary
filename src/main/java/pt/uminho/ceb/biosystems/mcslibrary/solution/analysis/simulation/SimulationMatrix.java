package pt.uminho.ceb.biosystems.mcslibrary.solution.analysis.simulation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import ilog.concert.IloException;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.implementation.CPLEXFlexibleFVA;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.SimulationResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.compression.alg.MatrixTools;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba.CPLEXFluxBalanceAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba.CPLEXParsimoniousFluxBalanceAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fva.CPLEXFluxVariabilityAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fva.FluxVariabilityAnalysisResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;
import pt.uminho.ceb.biosystems.mew.utilities.java.StringUtils;

public class SimulationMatrix {
	
	private double[][] matrix;
	private List<List<String>> solutions;
	private DefaultMetabolicNetwork metaNet;
	private double[][] min;
	private double[][] max;
	private double[][] essentials;
	private double[][] blocked;
	private double[][] essentialFlux;

	
	public SimulationMatrix(List<List<String>> solutions, DefaultMetabolicNetwork metaNet) {
		this.solutions = solutions;
		this.matrix = new double[solutions.size()][metaNet.getNumOfReactions()];
		this.metaNet = metaNet;
	}
	
	public int numSolutions(){
		return solutions.size();
	}
	
	public int numRx(){
		return metaNet.getNumOfReactions();
	}
	
	public double getValue(int row, int col){
		return this.matrix[row][col];
	}
	public void fill(FluxBound[] fluxBound, String objective, String sense) throws IloException{
		if (objective == null) {
			fillFVA(fluxBound, sense);
		} else {
			fillFBA(fluxBound, objective, sense);
		}
	}
	
	
	
	public void fillParsimonious(FluxBound[] flx, Reaction objective, String sense) throws Exception{
		CPLEXParsimoniousFluxBalanceAnalysis pfba = new CPLEXParsimoniousFluxBalanceAnalysis(metaNet, 0.999999);
		for (int i = 0; i < solutions.size(); i++) {
			SimulationResult sol = pfba.solveKnockoutFBA(flx, Utilities.toReacArrayFromString(metaNet, solutions.get(i)), objective, sense);
			matrix[i] = sol.getValues();
		}
	}

	private void fillFVA(FluxBound[] fluxBound, String sense) throws IloException {
		essentials = new double[solutions.size()][metaNet.getNumOfReactions()];
		blocked = new double[solutions.size()][metaNet.getNumOfReactions()];
		
		for (int i = 0; i < solutions.size(); i++) {
			FluxBound[] newFlux = new FluxBound[fluxBound.length+solutions.get(i).size()];
			for (int j = 0; j < fluxBound.length; j++) {
				newFlux[j] = fluxBound[j];
			}
			for (int j = 0; j < solutions.get(i).size(); j++) {
				newFlux[fluxBound.length+j] = new FluxBound(metaNet.getReaction(solutions.get(i).get(j)), 0, 0);
			}
			CPLEXFluxVariabilityAnalysis fva = new CPLEXFluxVariabilityAnalysis(metaNet, newFlux, null);
			FluxVariabilityAnalysisResult r = fva.solveFVA();
			for (int j = 0; j < metaNet.getNumOfReactions(); j++) {
				matrix[i][j] = sense.equals("min") ? r.minFlux(j) : r.maxFlux(j);
				essentials[i][j] = r.isEssential(j) ? 1 : 0;
				blocked[i][j] = r.isBlocked(j) ? 1 : 0;
			}
		}
	}
	
	public void calcFVA(FluxBound[] fluxBound) throws IloException {
		essentials = new double[solutions.size()][metaNet.getNumOfReactions()];
		blocked = new double[solutions.size()][metaNet.getNumOfReactions()];
		max = new double[solutions.size()][metaNet.getNumOfReactions()];
		min = new double[solutions.size()][metaNet.getNumOfReactions()];
		essentialFlux = new double[solutions.size()][metaNet.getNumOfReactions()];
		
		CPLEXFlexibleFVA fva = new CPLEXFlexibleFVA(metaNet);
		for (int i = 0; i < solutions.size(); i++) {
			System.out.println("Solution "+i);
			FluxBound[] newFlux = new FluxBound[fluxBound.length+solutions.get(i).size()];
			for (int j = 0; j < fluxBound.length; j++) {
				newFlux[j] = fluxBound[j];
			}
			for (int j = 0; j < solutions.get(i).size(); j++) {
				newFlux[fluxBound.length+j] = new FluxBound(metaNet.getReaction(solutions.get(i).get(j)), 0, 0);
			}
//			CPLEXFluxVariabilityAnalysis fva = new CPLEXFluxVariabilityAnalysis(metaNet, newFlux, null);
			FluxVariabilityAnalysisResult r = fva.solveFVA(newFlux, null);
			for (int j = 0; j < metaNet.getNumOfReactions(); j++) {
				max[i][j] = r.maxFlux(j);
				min[i][j] = r.minFlux(j);
				essentials[i][j] = r.isEssential(j) ? 1 : 0;
				blocked[i][j] = r.isBlocked(j) ? 1 : 0;
				essentialFlux[i][j] = Math.min(Math.abs(r.minFlux(j)), Math.abs(r.maxFlux(j)));
			}
		}
	}

	private void fillFBA(FluxBound[] fluxBound, String objective, String sense) throws IloException {
		CPLEXFluxBalanceAnalysis fba = new CPLEXFluxBalanceAnalysis(metaNet);
		for (int i = 0; i < solutions.size(); i++) {
			matrix[i] = fba.solveReactionKnockoutFBA(fluxBound, Utilities.toReacArrayFromString(metaNet, solutions.get(i)), metaNet.getReaction(objective), sense).getValues();
		}
	}
	
	public void fillEssentialFluxValue(FluxBound[] fluxBound, double eps) throws IloException{
		this.matrix = new double[solutions.size()][metaNet.getNumOfReactions()];
		CPLEXFluxVariabilityAnalysis fvawt = new CPLEXFluxVariabilityAnalysis(metaNet, fluxBound, null);
		FluxVariabilityAnalysisResult wt = fvawt.solveFVA();
		List<Integer> nonBlocked = new ArrayList<Integer>();
		for (int i = 0; i < metaNet.getNumOfReactions(); i++) {
			if (!wt.isBlocked(i)) {
				nonBlocked.add(i);
			}
		}
		
		System.out.println(nonBlocked.size()+" candidates out of "+metaNet.getNumOfReactions());
		for (int i = 0; i < solutions.size(); i++) {
			Reaction[] ko = Utilities.toReacArrayFromString(metaNet, solutions.get(i));
			CPLEXFluxVariabilityAnalysis fva = new CPLEXFluxVariabilityAnalysis(metaNet, fluxBound, ko, null);
			FluxVariabilityAnalysisResult fba = fva.solveFVA();
			System.out.println("Solution "+i+"; "+solutions.get(i));
			for (int j = 0; j < nonBlocked.size(); j++) {
				if (fba.isEssential(nonBlocked.get(j), eps)) {
					matrix[i][nonBlocked.get(j)] = Math.min(Math.abs(fba.minFlux(nonBlocked.get(j))), Math.abs(fba.maxFlux(nonBlocked.get(j))));
				}
			}
		}
	}

	public void writecsv(String filename) throws IOException{
		List<String> names = new ArrayList<String>();
		for (int i = 0; i < metaNet.getNumOfReactions(); i++) {
			names.add(metaNet.getReaction(i).getName());
		}
		List<String> rows = new ArrayList<String>();
		for (int i = 0; i < solutions.size(); i++) {
			rows.add(StringUtils.concat(" ", solutions.get(i)));
		}
		MatrixTools.writeCSV(matrix, filename, names, rows);
	}
	
	
	public void writeExtra(String filenameEssential, String filenameBlocked, String filenameFVA) throws IOException{
		List<String> names = new ArrayList<String>();
		for (int i = 0; i < metaNet.getNumOfReactions(); i++) {
			names.add(metaNet.getReaction(i).getName());
		}
		MatrixTools.writeCSV(essentials, filenameEssential, StringUtils.concat(",", names));
		MatrixTools.writeCSV(blocked, filenameBlocked, StringUtils.concat(",", names));
		MatrixTools.writeCSV(max, filenameFVA+"_max", StringUtils.concat(",", names));
		MatrixTools.writeCSV(min, filenameFVA+"_min", StringUtils.concat(",", names));
		MatrixTools.writeCSV(essentialFlux, filenameEssential+"_flux", StringUtils.concat(",", names));
	}
}
