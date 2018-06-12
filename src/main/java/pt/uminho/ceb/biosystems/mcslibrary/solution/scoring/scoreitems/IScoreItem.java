package pt.uminho.ceb.biosystems.mcslibrary.solution.scoring.scoreitems;

import java.util.List;

public interface IScoreItem {
	public double evaluateReactionKnockout(List<String> rk);
	public String getItemName();
}

