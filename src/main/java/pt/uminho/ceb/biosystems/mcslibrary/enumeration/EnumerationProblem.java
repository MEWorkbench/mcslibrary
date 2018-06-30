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
package pt.uminho.ceb.biosystems.mcslibrary.enumeration;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.List;
import java.util.StringTokenizer;

import pt.uminho.ceb.biosystems.mcslibrary.metabolic.AbstractMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.Reaction;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.FluxBound;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.constraints.YieldConstraint;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.CompressedMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.metabolic.implementation.DefaultMetabolicNetwork;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;
/**
 * Class that contains all information required to proceed with the calculation of minimal cut sets.
 * @author V�tor
 *
 */
public class EnumerationProblem implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = -3464962484033153317L;
	private FluxBound[] undesiredfluxes;
	private AbstractMetabolicNetwork metanet;
	private YieldConstraint[] undesiredyieldbounds;
	private FluxBound[] desiredfluxes;
	private boolean isCompressedProblem;
	private YieldConstraint[] desiredyieldbounds;
	private Reaction[] excludedReactions;
	private Solution[] excludedSubsets;
	private Reaction[] forced;
/**
 * Constructor for the enumeration problem.
 * @param metanet - An instance of any of the subclasses of {@link AbstractMetabolicNetwork}.
 * @param undesiredFluxes - An array of {@link FluxBound} corresponding to the undesired fluxes
 * @param desiredFluxes - An array of {@link FluxBound} corresponding to the desired fluxes
 * @param undesiredYields - An array of {@link YieldConstraint} corresponding to undesired yield ratios
 * @param desiredYields - An array of {@link YieldConstraint} corresponding to desired yield ratios
 */
	public EnumerationProblem(AbstractMetabolicNetwork metanet, FluxBound[] undesiredFluxes, FluxBound[] desiredFluxes, YieldConstraint[] undesiredYields, YieldConstraint[] desiredYields, Reaction[] excReactions) {
		this.metanet = metanet;
		this.undesiredfluxes = undesiredFluxes;
		this.undesiredyieldbounds = undesiredYields;
		this.desiredfluxes = desiredFluxes;
		this.desiredyieldbounds = desiredYields;
		this.excludedSubsets = new Solution[]{};
		this.setCompressedProblem(false);
		if (metanet.getClass().equals(CompressedMetabolicNetwork.class)){
			this.setCompressedProblem(true);
		}
		this.excludedReactions = excReactions;
	}

	public AbstractMetabolicNetwork getMetabolicNetwork() {
		return metanet;
	}
	
	public void setExcludedSubsets(Solution[] sol) {
		this.excludedSubsets = sol;
	}

	public void setMetabolicNetwork(AbstractMetabolicNetwork metanet) {
		this.metanet = metanet;
	}

	public FluxBound[] getUndesiredFluxes() {
		return undesiredfluxes;
	}

	public void setFluxbounds(FluxBound[] fluxbounds) {
		this.undesiredfluxes = fluxbounds;
	}

	public YieldConstraint[] getUndesiredYieldConstraints(){
		return this.undesiredyieldbounds;
	}

	public YieldConstraint[] getDesiredYieldConstraints(){
		return this.desiredyieldbounds;
	}

	public boolean isCompressedProblem() {
		return isCompressedProblem;
	}

	private void setCompressedProblem(boolean isCompressedProblem) {
		this.isCompressedProblem = isCompressedProblem;
	}

	public FluxBound[] getDesiredFluxes() {
		return desiredfluxes;
	}

	public void setDesiredfluxes(FluxBound[] desiredfluxes) {
		this.desiredfluxes = desiredfluxes;
	}

	public Reaction[] getExcludedReactions() {
		return excludedReactions;
	}
	
	public Solution[] getExcludedSolutions() {
		return this.excludedSubsets;
	}

	public void setExcludedReactions(Reaction[] excludedReactions) {
		this.excludedReactions = excludedReactions;
	}
	
	public void setForcedReactions(Reaction[] forcedReactions){
		this.forced = forcedReactions;
	}
	public static void saveFluxBoundArray(FluxBound[] flx, String filename) throws IOException{
		BufferedWriter flxfile = new BufferedWriter(new FileWriter(filename));
		for (FluxBound fluxBound : flx) {
			String reactionId = fluxBound.getReac().getName();
			String lb = Double.toString(fluxBound.getBounds().getLower());
			String ub = Double.toString(fluxBound.getBounds().getUpper());
			String toWrite = reactionId+";"+lb+";"+ub;
			flxfile.write(toWrite+"\n");
		}
		flxfile.flush();
		flxfile.close();
	}
	
	public static void saveYieldConstraintArray(YieldConstraint[] flx, String filename) throws IOException{
		BufferedWriter flxfile = new BufferedWriter(new FileWriter(filename));
		for (YieldConstraint fluxBound : flx) {
			String substrate = fluxBound.getUptakeReaction().getName();
			String product = fluxBound.getProductReaction().getName();
			String ratio = Double.toString(fluxBound.getRatio());
			String toWrite = substrate+";"+product+";"+ratio;
			flxfile.write(toWrite+"\n");
		}
		flxfile.flush();
		flxfile.close();
	}
	
	public static void saveMetabolicNetwork(AbstractMetabolicNetwork abs, String filename) throws IOException{
		abs.saveNetwork(filename);
	}
	public void saveEnumerationProblem(String problemname) throws IOException{
		FluxBound[] undflux = this.getUndesiredFluxes();
		FluxBound[] desflux = this.getDesiredFluxes();
		YieldConstraint[] undyield = this.getUndesiredYieldConstraints();
		YieldConstraint[] desyield = this.getDesiredYieldConstraints();
		
		saveMetabolicNetwork(this.metanet, problemname);
		saveFluxBoundArray(undflux, problemname+".uf");
		saveFluxBoundArray(desflux, problemname+".df");
		saveYieldConstraintArray(undyield, problemname+".uy");
		saveYieldConstraintArray(desyield, problemname+".dy");
	}
	
	public static FluxBound[] loadFluxArray(String filename, DefaultMetabolicNetwork metaNet) throws IOException{
		
		List<String> fluxfile = Utilities.readLines(filename);
		FluxBound[] res = new FluxBound[fluxfile.size()];
		for (int i = 0; i < fluxfile.size(); i++) {
			StringTokenizer tok = new StringTokenizer(fluxfile.get(i),";");
			Reaction r = metaNet.getReaction(metaNet.getReactionIndex(tok.nextToken()));
			Double lb = Double.parseDouble(tok.nextToken());
			Double ub = Double.parseDouble(tok.nextToken());
			res[i] = new FluxBound(r, lb, ub);
		}
		return res;
	}
	public static YieldConstraint[] loadYieldArray(String filename, DefaultMetabolicNetwork metaNet) throws IOException{
		List<String> yldfile = Utilities.readLines(filename);
		YieldConstraint[] res = new YieldConstraint[yldfile.size()];
		for (int i = 0; i < yldfile.size(); i++) {
			StringTokenizer tok = new StringTokenizer(yldfile.get(i),";");
			Reaction s = metaNet.getReaction(metaNet.getReactionIndex(tok.nextToken()));
			Reaction p = metaNet.getReaction(metaNet.getReactionIndex(tok.nextToken()));
			Double ratio = Double.parseDouble(tok.nextToken());
			res[i] = new YieldConstraint(s, p, ratio);
		}
		return res;
	}
	public static EnumerationProblem loadEnumerationProblem(String problemname) throws IOException{
		DefaultMetabolicNetwork metaNet = DefaultMetabolicNetwork.loadNetwork(problemname);
		FluxBound[] uf = loadFluxArray(problemname, metaNet);
		FluxBound[] df = loadFluxArray(problemname, metaNet);
		YieldConstraint[] uy = loadYieldArray(problemname, metaNet);
		YieldConstraint[] dy = loadYieldArray(problemname, metaNet);
		return new EnumerationProblem(metaNet, uf, df, uy, dy, new Reaction[]{});
	}

	public Reaction[] getForcedReactions() {
		return forced;
	}
	
	public void serializeData(String filename) throws IOException{
	    FileOutputStream fos = new FileOutputStream(filename);
	    ObjectOutputStream oos = new ObjectOutputStream(fos);
	    oos.writeObject(this);
	    oos.close();
	}
	
	public static EnumerationProblem readSerializedData(String filename) throws ClassNotFoundException, IOException{
	   FileInputStream fin = new FileInputStream(filename);
	   ObjectInputStream ois = new ObjectInputStream(fin);
	   EnumerationProblem obj = (EnumerationProblem) ois.readObject();
	   ois.close();
	   return obj;
	}
}
