package pt.uminho.ceb.biosystems.mcslibrary.solution.scoring;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.compression.alg.MatrixTools;
import pt.uminho.ceb.biosystems.mcslibrary.solution.scoring.scoreitems.IScoreItem;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;
import pt.uminho.ceb.biosystems.mew.utilities.java.StringUtils;

public class SolutionScorer {

	public SolutionScorer() {
//		this.cont = cont;
//		this.dmn = dmn;
	}
	
	public Map<Set<String>, Map<IScoreItem,Double>> getDataset(Map<Set<String>,List<String>> solutions, IScoreItem[] scorers, String writePath) throws IOException{
		List<String> header = new ArrayList<>();
		header.add("SolutionName");
		header.add("Solution");
		for (int i = 0; i < scorers.length; i++) {
			header.add(scorers[i].getItemName());
		}
		
		List<List<String>> data  = new ArrayList<>();
		data.add(header);
		
		Map<Set<String>, Map<IScoreItem,Double>> map = new HashMap<>();

		for (Set<String> list : solutions.keySet()) {
			List<String> rk = solutions.get(list);
			List<String> lst = new ArrayList<String>();
			lst.add(StringUtils.concat("+", list));
			lst.add(StringUtils.concat(" ", solutions.get(list)));
			Map<IScoreItem,Double> submap = new HashMap<>();
			for (int i = 0; i < scorers.length; i++) {
				double score = scorers[i].evaluateReactionKnockout(rk);
				lst.add(Double.toString(score));
				submap.put(scorers[i], score);
			}
			map.put(list, submap);
			data.add(lst);
		}

		Utilities.writeDataset(data, writePath, ",");
		return map;
		
	}
	
	public void writeDataset(List<List<String>> solutions, IScoreItem[] scorers, String writePath) throws IOException{
//		Map<List<String>, Map<IScoreItem,Double>> map = new HashMap<>();
		double[][] data = new double[solutions.size()][scorers.length];
		List<String> header = new ArrayList<String>();
		List<String> rows = new ArrayList<String>();
		for (int i = 0; i < scorers.length; i++) {
			IScoreItem scoreItem = scorers[i];
			header.add(scoreItem.getItemName());
			for (int j = 0; j < solutions.size(); j++) {
				data[j][i] = scoreItem.evaluateReactionKnockout(solutions.get(j));
				if (i == 0) {
					rows.add(StringUtils.concat(" ", solutions.get(j)));
				}
//				if (data[j][i] > 0) {
//					System.out.println(solutions.get(j));
//				}
			}
		}
		MatrixTools.writeCSV(data, writePath, header, rows);
	}
}
