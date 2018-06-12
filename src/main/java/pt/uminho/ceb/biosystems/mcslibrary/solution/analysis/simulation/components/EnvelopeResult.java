package pt.uminho.ceb.biosystems.mcslibrary.solution.analysis.simulation.components;

public class EnvelopeResult {
	private double[] min;
	private double[] max;
	private EnvelopeProperties prop;
	private double[] scale;
	private int points;

	public EnvelopeResult(double[] min, double[] max, double[] scale, EnvelopeProperties prop) {
		this.min = min;
		this.max = max;
		this.prop = prop;
		this.scale = scale;
		this.points = scale.length;
	}
	
	public EnvelopeProperties properties(){
		return this.prop;
	}
	
	public double getMinAtPoint(int i){
		return min[i];
	}
	
	public double getMaxAtPoint(int i){
		return max[i];
	}
	
	public int length(){
		return this.points;
	}
	
	public double scaleAt(int i){
		return this.scale[i];
	}
	
}
