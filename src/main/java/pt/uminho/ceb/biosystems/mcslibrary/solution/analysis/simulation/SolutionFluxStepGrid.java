package pt.uminho.ceb.biosystems.mcslibrary.solution.analysis.simulation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import ilog.concert.IloException;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.SimulationResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.compression.alg.MatrixTools;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba.CPLEXFluxBalanceAnalysis;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;
import pt.uminho.ceb.biosystems.mew.utilities.java.StringUtils;

public class SolutionFluxStepGrid {
	private double[][] matrix;
	private List<List<String>> solutions;
	private DefaultMetabolicNetwork metaNet;
	private double minv;
	private double maxv;
	
	public SolutionFluxStepGrid(List<List<String>> solutions, DefaultMetabolicNetwork metaNet) {
		this.solutions = solutions;
		this.metaNet = metaNet;
	}
	
	public int numSolutions(){
		return solutions.size();
	}
	
	public int numSteps(){
		return metaNet.getNumOfReactions();
	}
	
	public double getValue(int row, int col){
		return this.matrix[row][col];
	}
	public void fill(FluxBound[] fluxBound, String objective, String variable, String sense, int steps) throws IloException{
		this.matrix = new double[solutions.size()][steps];
		CPLEXFluxBalanceAnalysis fba = new CPLEXFluxBalanceAnalysis(metaNet);
		for (int i = 0; i < solutions.size(); i++) {
			SimulationResult min = fba.solve(fluxBound, metaNet.getReaction(variable), "min");
			SimulationResult max = fba.solve(fluxBound, metaNet.getReaction(variable), "max");
			minv = min.getObjectiveValue();
			maxv = max.getObjectiveValue();
			double[] s = fba.solveFVAsteps(fluxBound, Utilities.toReacArrayFromString(metaNet, solutions.get(i)), steps-1, sense, metaNet.getReaction(objective), metaNet.getReaction(variable));
			matrix[i] = s;
		}
		for (int i = 0; i < fluxBound.length; i++) {
			
		}
	}
	
	
	public void writecsv(String filename) throws IOException{
		List<String> names = new ArrayList<String>();
		double diff = maxv - minv;
		for (int i = 0; i < matrix[0].length; i++) {
			names.add(Double.toString(minv + (diff * (double)i/matrix[0].length)));
		}
		MatrixTools.writeCSV(matrix, filename, StringUtils.concat(",", names));
	}
}
