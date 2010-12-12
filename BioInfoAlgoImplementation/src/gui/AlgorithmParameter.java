/**
 * @author Martin Mann - http://www.bioinf.uni-freiburg.de/~mmann/
 */
package gui;



/**
 * This class provides a generic representation for the input parameters
 * of  BioinfAlgorithm. It is used to build up a GUI for user inputs and
 * to run the algorithm
 * 
 * @author Martin Mann - http://www.bioinf.uni-freiburg.de/~mmann/
 *
 */
public class AlgorithmParameter {
	
	/**
	 * A short name of the parameter.
	 */
	public String name = "";
	
	/**
	 * A detailed description of the parameter.
	 */
	public String description = "";

	/**
	 * The class of the Parameter. (Currently all Number derived classes 
	 * and String supported.)
	 */
	@SuppressWarnings("rawtypes")
	public Class type = Object.class;
	
	/**
	 * The value of the parameter submitted to run the algorithm.
	 */
	public Object data = null;
	
	/**
	 * The default value of the parameter.
	 */
	public Object defVal = null;

	
	/**
	 * Creates a parameter and initialises the internal datastructures.
	 * 
	 * @param name_		a short name of the parameter
	 * @param descr_	a detailed parameter description
	 * @param type_		the type of the parameter (String or Number derived classes)
	 */
	@SuppressWarnings("rawtypes")
	public AlgorithmParameter(String name_, String descr_, Class type_) {
		name = name_;
		description = descr_;
		type = type_;
	}

	/**
	 * Creates a parameter and initialises the internal datastructures.
	 * 
	 * @param name_		a short name of the parameter
	 * @param descr_	a detailed parameter description
	 * @param type_		the type of the parameter (String or Number derived classes)
	 * @param defaultData the default value of the parameter
	 */
	@SuppressWarnings("rawtypes")
	public AlgorithmParameter(String name_, String descr_, Class type_, Object defaultData) {
		this(name_, descr_, type_);
			// check if the default value can be casted
		if (defaultData != null && !defaultData.getClass().equals(type_))
			throw new ClassCastException();
		else
			defVal = defaultData;
	}
}
