package pt.uminho.ceb.biosystems.mcslibrary.solution.analysis;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Pair;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;

public class TurnoverScore {
	public Pair<Double, Reaction>[] objective;
	public Double score;

	public TurnoverScore(Pair<Double,Reaction>[] objective, Double score) {
		this.objective = objective;
		this.score = score;
	}
	
	public String toString(){
		String objStr = "";
		for (int i = 0; i < objective.length; i++) {
			objStr += objective[i].getA() + "*" + objective[i].getB().getName() + " + ";
		}
		String finalstr = objStr.substring(0, objStr.length() - 2) + " = " + score;
		return finalstr;
	}
	
	public boolean isEqual(TurnoverScore object){
		return ((this.score - object.score) < Utilities.PRECISION) && equalObjectives(object);
	}

	private boolean equalObjectives(TurnoverScore object) {
		boolean eq = true;
		for (int i = 0; i < objective.length; i++) {
			eq |= (objective[i].getA() == object.objective[i].getA()) && (objective[i].getB().equals(object.objective[i].getB()));
		}
		return eq;
	}
	
}
