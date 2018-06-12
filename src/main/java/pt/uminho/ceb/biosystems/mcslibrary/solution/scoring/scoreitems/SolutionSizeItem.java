package pt.uminho.ceb.biosystems.mcslibrary.solution.scoring.scoreitems;

import java.util.List;

public class SolutionSizeItem implements IScoreItem{

	@Override
	public double evaluateReactionKnockout(List<String> rk) {
		return rk.size();
	}

	@Override
	public String getItemName() {
		return "SolutionSize";
	}

}
