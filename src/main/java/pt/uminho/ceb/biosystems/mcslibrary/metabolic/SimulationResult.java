package pt.uminho.ceb.biosystems.mcslibrary.metabolic;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.CompressedMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;

public class SimulationResult {

	private AbstractMetabolicNetwork metaNet;
	private double[] values;
	private double objectiveValue;
	private String status;

	public SimulationResult(AbstractMetabolicNetwork metaNet,
			double[] values) {
		this.metaNet = metaNet;
		this.values = values;
	}
	
	public SimulationResult(AbstractMetabolicNetwork metaNet,
			double[] values, double objectiveValue) {
		this.metaNet = metaNet;
		this.values = values;
		this.objectiveValue = objectiveValue;
	}
	
	public SimulationResult(AbstractMetabolicNetwork metaNet,
			double[] values, double objectiveValue, String status) {
		this.metaNet = metaNet;
		this.values = values;
		this.objectiveValue = objectiveValue;
		this.status = status;
	}

	public String getStatus(){
		return this.status;
	}
	public double getObjectiveValue(){
		return this.objectiveValue;
	}
	public double getFluxValue(Reaction r) {
		return values[metaNet.getReactionIndex(r.getName())];
	}
	
	public double[] getValues(){
		return this.values;
	}

	public void printValue(Reaction r) {
		System.out.println(r.getName()+" = "+getFluxValue(r));
	}

	public String toString(){
		Class<? extends AbstractMetabolicNetwork> klaz = metaNet.getClass();
		String res = "";
		if (klaz == CompressedMetabolicNetwork.class) {
			CompressedMetabolicNetwork mninstance = (CompressedMetabolicNetwork) metaNet;
			for (int i = 0; i < mninstance.getNumOfReactions(); i++) {
				ReactionGroup r = mninstance.getReactionGroup(i);
				for (int j = 0; j < r.size(); j++) {
					res += r.getReaction(j).getName()+";";
				}
				res += "\t"+values[i]+"\n";
				
			}
		} else if (klaz == DefaultMetabolicNetwork.class) {
			DefaultMetabolicNetwork mninstance = (DefaultMetabolicNetwork) metaNet;
			for (int i = 0; i < mninstance.getNumOfReactions(); i++) {
				res += mninstance.getReaction(i).getName()+"\t"+values[i]+"\n";
			}
		}
		return res;

	}

}
