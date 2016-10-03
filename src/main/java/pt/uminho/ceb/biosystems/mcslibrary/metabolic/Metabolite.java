package pt.uminho.ceb.biosystems.mcslibrary.metabolic;

import java.io.Serializable;

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
	
}
