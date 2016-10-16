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
package pt.uminho.ceb.biosystems.mcslibrary.metabolic.fba;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.AbstractMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.ReactionGroup;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.CompressedMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;

public class FluxBalanceAnalysisResult {

	private AbstractMetabolicNetwork metaNet;
	private double[] values;

	public FluxBalanceAnalysisResult(AbstractMetabolicNetwork metaNet,
			double[] values) {
		this.metaNet = metaNet;
		this.values = values;
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
