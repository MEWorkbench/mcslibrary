package pt.uminho.ceb.biosystems.mcslibrary.enumeration.cplex.callbacks;

import ilog.concert.IloException;
import ilog.cplex.IloCplex.MIPInfoCallback;
import pt.uminho.ceb.biosystems.mcslibrary.utilities.Utilities;

public class DefaultCallback extends MIPInfoCallback{
	int curObj;
	private int solSize;
	
	public DefaultCallback() {
		super();
		this.solSize = Integer.MAX_VALUE;

	}

	@Override
	protected void main() throws IloException {
		if (hasIncumbent()) {
			if (getIncumbentObjValue() < this.solSize) {
				this.solSize = (int) getIncumbentObjValue();
				System.out.println("CPLEX found a solution with size "+getIncumbentObjValue());
			}
		}
	}

}
