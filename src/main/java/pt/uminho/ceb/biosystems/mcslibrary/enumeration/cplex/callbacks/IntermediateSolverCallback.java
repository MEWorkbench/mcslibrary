package pt.uminho.ceb.biosystems.mcslibrary.enumeration.cplex.callbacks;

import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.cplex.IloCplex.MIPInfoCallback;
import pt.uminho.ceb.biosystems.mcslibrary.enumeration.implementation.CPLEXIntegratedEnumerator;

public class IntermediateSolverCallback extends MIPInfoCallback{
	int curObj;
	private int solSize;
	private double starttime;
	private double stopTime;
	private IloIntVar[] ZP;
	private IloIntVar[] ZN;
	private int[] mcs;
	private boolean continueIfNoSols;
	private boolean abortWasCalled;
	
	public IntermediateSolverCallback(double stopTime, boolean continueIfNoSols, IloIntVar[] ZP, IloIntVar[] ZN) throws IloException {
		super();
		this.starttime = System.currentTimeMillis();
		this.ZP = ZP;
		this.ZN = ZN;
		this.stopTime = stopTime;
		this.solSize = Integer.MAX_VALUE;
		this.continueIfNoSols = continueIfNoSols;
		this.abortWasCalled = false;
	}

	@Override
	protected void main() throws IloException {
		if (hasIncumbent()) {
			if (getIncumbentObjValue() < this.solSize) {
				this.solSize = (int) getIncumbentObjValue();
				System.out.println("CPLEX found a solution with size "+getIncumbentObjValue());
				double[] zpv = getIncumbentValues(ZP);
				double[] znv = getIncumbentValues(ZN);
				this.mcs = CPLEXIntegratedEnumerator.recoverArray(zpv, znv);
//				System.out.println(Arrays.toString(this.mcs));
//				System.out.println(Arrays.toString(zpv));
//				System.out.println(Arrays.toString(znv));

				
			}
		} 
		if (((System.currentTimeMillis()-this.starttime)/1000 > this.stopTime)){
//			System.out.println((getCplexTime()-this.starttime)+" "+this.stopTime);
			if(hasIncumbent()){
				System.out.println("Initial optimization ended");
				this.abort();
				abortWasCalled = true; 
			} else if (continueIfNoSols){
				
			} else {
				System.out.println("\tCPLEX was aborted");
				this.abort();
				abortWasCalled = true;
			}
		}


//		} else if (getCplexTime() > this.stopTime && !hasIncumbent() && continueIfNoSols){
//			System.out.println("Out of time");
//		} else if (getCplexTime() < this.stopTime && !hasIncumbent() && continueIfNoSols){
//		} else if (getCplexTime() < this.stopTime && !hasIncumbent()){
//		} else{
//			System.out.println("Optimizer aborted");
//			this.abort();
//			abortWasCalled = true;
//		}
	}
	
	public boolean hasSolution(){
		return this.mcs != null;
	}

	public int[] getFinalSol() {
		return this.mcs;
	}

	public boolean hasOptimalSolution(){
		return !abortWasCalled;
	}

	public void setStartTime(double cplexTime) {
		this.starttime = cplexTime;
		
	}
}
