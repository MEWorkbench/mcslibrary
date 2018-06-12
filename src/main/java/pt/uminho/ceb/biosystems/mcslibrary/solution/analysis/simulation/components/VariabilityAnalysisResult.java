package pt.uminho.ceb.biosystems.mcslibrary.solution.analysis.simulation.components;

import java.util.List;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.fva.FluxVariabilityAnalysisResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;

public class VariabilityAnalysisResult {
	private DefaultMetabolicNetwork dmn;
	private List<String> solution;
	private FluxVariabilityAnalysisResult limits;
	private EnvelopeResult[] envelope;
	private FluxBound[] bounds;

	public VariabilityAnalysisResult(DefaultMetabolicNetwork dmn, List<String> solution, FluxVariabilityAnalysisResult limits, EnvelopeResult[] envelopes, FluxBound[] bounds){
		this.dmn = dmn;
		this.solution = solution;
		this.limits = limits;
		this.envelope = envelopes;
		this.bounds = bounds;
	}
	
	public double[] getEssentials(){
		double[] ess = new double[dmn.getNumOfReactions()];
		for (int i = 0; i < ess.length; i++)
			ess[i] = limits.isEssential(i) ? 1 : 0;
		return ess;
	}
	
	public double[] getBlocked(){
		double[] blk = new double[dmn.getNumOfReactions()];
		for (int i = 0; i < blk.length; i++)
			blk[i] = limits.isBlocked(i) ? 1 : 0;
		return blk;
	}
	
	public double[] getMinLimits(){
		return this.limits.minToArray();
	}
	
	public double[] getMaxLimits(){
		return this.limits.maxToArray();
	}
	
	public int numberOfEnvelopes(){
		return envelope.length;
	}
	
	public EnvelopeResult getEnvelope(int index){
		return envelope[index];
	}
	
	public FluxBound[] getConditions(){
		return this.bounds;
	}
	
	public List<String> getSolution(){
		return this.solution;
	}
	
}
