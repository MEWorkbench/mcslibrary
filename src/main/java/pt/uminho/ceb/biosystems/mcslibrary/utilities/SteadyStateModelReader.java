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
package pt.uminho.ceb.biosystems.mcslibrary.utilities;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.ReactionConstraint;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mew.core.model.components.Metabolite;
import pt.uminho.ceb.biosystems.mew.core.model.components.Reaction;
import pt.uminho.ceb.biosystems.mew.core.model.steadystatemodel.ISteadyStateModel;



public class SteadyStateModelReader {
	private ISteadyStateModel model;
	private Integer nreacts;
	private Integer nmetabs;

	public SteadyStateModelReader(ISteadyStateModel model) {
		this.model = model;
		nreacts = model.getNumberOfReactions();
		nmetabs = model.getNumberOfMetabolites();
	}

	public pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction[] createReactions(){
		pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction[] res = new pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction[nreacts];
		for (int i = 0; i < nreacts; i++) {
			Reaction r = model.getReaction(i);
			ReactionConstraint limits = new ReactionConstraint(r.getConstraints().getLowerLimit(), r.getConstraints().getUpperLimit());
			res[i] = new pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction(r.getId(), limits);
		}
		return res;
	}

	public pt.uminho.ceb.biosystems.mcslibrary.metabolic.Metabolite[] createMetabolites(){
		pt.uminho.ceb.biosystems.mcslibrary.metabolic.Metabolite[] res = new pt.uminho.ceb.biosystems.mcslibrary.metabolic.Metabolite[nmetabs];
		for (int i = 0; i < nmetabs; i++) {
			Metabolite m = model.getMetabolite(i);
			res[i] = new pt.uminho.ceb.biosystems.mcslibrary.metabolic.Metabolite(m.getId(), m.isBoundaryCondition());
		}
		return res;
	}

	public double[][] createMatrix(){
		double[][] res = new double[nmetabs][nreacts];
		for (int i = 0; i < nmetabs; i++) {
			for (int j = 0; j < nreacts; j++) {
				res[i][j] = model.getStoichiometricValue(i, j);
			}
		}
		return res;
	}

	public DefaultMetabolicNetwork convertModel(){
		DefaultMetabolicNetwork mod = new DefaultMetabolicNetwork(createMetabolites(),createReactions(),createMatrix());
		return mod;
	}
}
