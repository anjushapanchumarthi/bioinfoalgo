/**
 * 
 */
package gui;


import java.util.Vector;


/**
 * This class is the superclass for the BioinfAlgo algorithm
 * GUI. It defines and provides access functions 
 * for a generic GUI creation and algorithm call.
 * 
 * @author Martin Mann - http://www.bioinf.uni-freiburg.de/~mmann/
 *
 */
public abstract class BioinfAlgorithm {
	

	  /**
	   * List of all needed parameters and their default values that can be
	   * used by subclasses.
	   */
	protected Vector<AlgorithmParameter> parameters = new Vector<AlgorithmParameter>();

	
	/**
	 * Access to the algorithm name.
	 * 
	 * @return	the name of the algorithm
	 */
	abstract public String getName();
	
	
	/**
	 * A short description of the algorithms functionality.
	 * 
	 * @return	the short description of the algorithm
	 */
	abstract public String getDescription();
	
	/**
	 * The returned parameter descriptions will be used to set up the
	 * parameter input interface in the GUI. The parameters are listed in
	 * the order given by the returned vector.
	 * 
	 * @return the necessary parameters for the algorithm
	 */
	abstract public Vector<AlgorithmParameter> getInputParameters();
	
	/**
	 * Starts the algorithm with the given parameters. The algorithm
	 * output will be returned as string. In error case the error messages
	 * are returned via String too.
	 * 
	 * @param params	the input parameters for the run
	 * @return			the output of the algorithm or an error message
	 */
	abstract public String run(Vector<AlgorithmParameter> params);
	
	/**
	 * Runs the given algorithm using its default parameter values as the input
	 * of the call. Output is written to System.out.
	 * 
	 * @param algorithm the algorithm to run
	 */
	static 
	public 
	void
	runAlgorithmDefaults( BioinfAlgorithm algorithm ) {
		
		  // announce algorithm
		System.out.println(
				"\n" +
				"\n============= ALGORITHM =============" +
				"\n '" + algorithm.getName() + "'" +
				"\n============ DESCRIPTION ============" +
				"\n '" + algorithm.getDescription() + "'" +
				"\n=============== INPUT ===============");
		
		  // copy default values as input values and report
		Vector<AlgorithmParameter> params = algorithm.getInputParameters();

		for (int i = 0; i < params.size(); i++) {
			params.elementAt(i).data = params.elementAt(i).defVal;
			System.out.println("\n "
					+ params.elementAt(i).name
					+ " = "
					+ params.elementAt(i).defVal.toString() );
		}
		
		System.out.println("\n============== RESULT ===============\n");

		System.out.println(algorithm.run(params));

		System.out.println("\n=============== DONE ================\n");
	}
}
