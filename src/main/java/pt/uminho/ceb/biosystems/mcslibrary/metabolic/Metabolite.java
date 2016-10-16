/*******************************************************************************
 * Copyright 2016
 * CEB Centre of Biological Engineering
 * University of Minho
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This code is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this code. If not, see http://www.gnu.org/licenses/
 *
 * Created inside the BIOSYSTEMS Research Group
 * (http://www.ceb.uminho.pt/biosystems)
 *******************************************************************************/
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
