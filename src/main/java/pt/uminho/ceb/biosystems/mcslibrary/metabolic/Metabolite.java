package pt.uminho.ceb.biosystems.mcslibrary.metabolic;

import java.io.Serializable;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import pt.uminho.ceb.biosystems.mcslibrary.utilities.RegexUtils;

/**
 * Class that represents a metabolite for use in metabolic networks
 * @author V�tor
 *
 */
public class Metabolite implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 3472441352248104475L;
	private String name;
	private boolean isExternal;
	private String formula;
/**
 * 
 * @param name - The metabolite's name.
 * @param isExternal - Boolean indicating whether this metabolite is in the external compartment
 */
	public Metabolite(String name, boolean isExternal) {
		this.name = name;
		this.isExternal = isExternal;
	}

	public String getName() {
		return this.name;
	}

	public boolean isExternal() {
		return isExternal;
	}
	
	public String getFormula(){
		return this.formula;
	}
	
	public void setFormula(String formula){
		this.formula = formula;
	}
	
	public int getCarbonAtoms(){
		List<List<String>> m = RegexUtils.findAll("C([0-9]+)", formula);
		if (m.size() == 0) {
			return 0;
		} else {
			List<String> c = m.get(0);
			String s = c.get(0).replace("C", "");
			return Integer.parseInt(s);
		}
	}
	
	public String toString(){
		return this.name + " = " + this.formula;
	}
}
