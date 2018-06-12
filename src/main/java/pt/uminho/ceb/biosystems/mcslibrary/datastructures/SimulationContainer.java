package pt.uminho.ceb.biosystems.mcslibrary.datastructures;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.SimulationResult;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mew.utilities.java.StringUtils;

public class SimulationContainer {
	public Map<String,SimulationResult> map;
	private DefaultMetabolicNetwork dmn;
	
	public SimulationContainer(DefaultMetabolicNetwork dmn) {
		this.map = new HashMap<String,SimulationResult>();
		this.dmn = dmn;
	}
	
	public String toString(){
		String r = "Reaction\t";
		List<String> keys = new ArrayList<String>(this.map.keySet());
		r += StringUtils.concat("\t", keys); 
		
		for (int i = 0; i < dmn.getNumOfReactions(); i++) {
			r += "\n"+dmn.getReaction(i).getName();
			for (int j = 0; j < keys.size(); j++) {
				r += "\t"+map.get(keys.get(j)).getFluxValue(dmn.getReaction(i));
				
			}
		}
		return r;
	}
}
