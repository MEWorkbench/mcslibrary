package pt.uminho.ceb.biosystems.mcslibrary.enumeration;
/**
 * Abstract class that receives an enumeration problem. Contains methods to calculate minimal cut sets for the provided problem
 * @author vvieira
 *
 *
 */
public abstract class AbstractEnumerationSolver {
	private EnumerationProblem eprob;

	public AbstractEnumerationSolver(EnumerationProblem eprob) {
		this.eprob = eprob;
	}
	
	public EnumerationProblem getProblem() {
		return this.eprob;
	}
	
	/**
	 *
	 * @param maxsize - An integer specifying the maximum size of the minimal cut sets to be calculated
	 * @return The enumeration result as an instance of a subclass of AbstractEnumerationResult
	 * @throws Exception
	 */
	public abstract AbstractEnumerationResult solve(int maxsize) throws Exception;
}
