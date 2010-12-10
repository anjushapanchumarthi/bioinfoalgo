package cthoelken;

import gui.AlgorithmParameter;
import gui.BioinfAlgorithm;
import gui.StringList;

import java.util.Vector;

 /**
  * Example implementation of the BioinfAlgorithm interface.
  * 
  * @author Martin Mann - http://www.bioinf.uni-freiburg.de/~mmann/
  *
  */
public class NeedlemanWunsch extends BioinfAlgorithm {



	 /**
	  * Constructs an example algorithm and initializes all allowed parameters.
	  */
	public NeedlemanWunsch() { 
		
		// create one example parameter for each type of parameter allowed
		
		super.parameters.add(new AlgorithmParameter(	
				"String-Parameter"
				, "This is a simple string that will be set within" +
					" the GUI via a single line TextField." 
				, String.class
				, "My Default String Value."));
		super.parameters.add(new AlgorithmParameter(	
				"StringList-Parameter"
				, "This is a string that allows for multiple lines" +
					" that will be set within" +
					" the GUI via a multi line TextArea." 
				, StringList.class
				, new StringList("My Default String\n List Value of\n multiple lines!")));
		super.parameters.add(new AlgorithmParameter(
				"Boolean-Parameter"
				, "This is a boolean that will be set within" +
					" the GUI via two RadioButtons." 
				, Boolean.class 
				, new Boolean(true)));
		super.parameters.add(new AlgorithmParameter(
				"Integer-Parameter"
				, "This is a integer that will be set within" +
					" the GUI via a dedicated TextField for Integers." 
				, Integer.class 
				, new Integer(12345)));
		super.parameters.add(new AlgorithmParameter(
				"Double-Parameter"
				, "This is a double that will be set within" +
					" the GUI via a dedicated TextField for Doubles." 
				, Double.class 
				, new Double(-1.2345)));

	}

	@Override
	public Vector<AlgorithmParameter> getInputParameters() {
		return super.parameters;
	}

	@Override
	public String getName() {
		return new String("Example of BioinfAlgorithm");
	}

	@Override
	public String getDescription() {
		return new String("This class exemplifies the available parameter" +
				" types supported by the GUI into which subclasses are" +
				" integrated.");
	}
	
	@Override
	public String run(Vector<AlgorithmParameter> params) {
		
		String retVal = new String("");
		
		  // ##########  PARSE INPUT PARAMETERS FOR ERRORS  ###########
		
		retVal += "\n mhh.. I assume my default values are fine ! ;) \n";
		
		  // ##########  RUN THE PROGRAM  ###########

		 // (JUST TO EXEMPLIFY SOME OUTPUT I REPORT THE INPUT PARAMETERS ...)
		for (int i = 0; i < params.size(); i++) {
			retVal	+= "\n Input : "
					+ params.elementAt(i).name
					+ " = "
					+ params.elementAt(i).defVal.toString();
		}
		
		 // (TO EXEMPLIFY SOME CALCULATION I DETERMINE THE NUMBER OF LINES
		 //  WITHIN THE StringList PARAMETER ...)
		int newLinePos = 0, newLineCount = 0;
		String multiLineParam = ((StringList)params.elementAt(1).data).toString();
		do  {
			newLinePos = (multiLineParam.indexOf('\n', newLinePos) + 1);
			newLineCount++;
		} while (newLinePos > 0);
		  // append to output
		retVal	+= "\n\n The StringList parameter contained " +
				+ newLineCount
				+ " lines!\n";
		
		
		  // returning the algorithm results ...
		return retVal;
	}

	/**
	 * Creates an instance of this class and calls the run method using the
	 * default parameters.
	 * 
	 * @param args program parameters (completely ignored)
	 */
	public static void main(String[] args) {
		  // create an instance of this class
		NeedlemanWunsch myInstance = new NeedlemanWunsch();
		  // run the example the instance with the default parameters
		BioinfAlgorithm.runAlgorithmDefaults( myInstance );
	}

}
