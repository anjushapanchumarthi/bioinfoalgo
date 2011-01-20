package cthoelken;

import gui.AlgorithmParameter;
import gui.BioinfAlgorithm;
import gui.StringList;

import java.util.LinkedList;
import java.util.Random;
import java.util.Vector;

 /**
  * Implementation of the Needleman Wunsch Algorithm for pairwise alignments
  * with linear gap costs.
  * 
  * @author Clemens Thoelken
  *
  */
public class SumOfPairs extends BioinfAlgorithm {

	protected Alignment sequences;
	protected boolean usePAM;
	protected double gapCosts;
	protected Alignment alignment;
	protected SubstitutionMatrix omega;

	 /**
	  * Constructor which generates an empty vector of parameters of the needed 
	  * types.
	  */
	public SumOfPairs() { 
		
		// create all needed parameters for the algorithm to work.
		
		super.parameters.add(new AlgorithmParameter(	
				"Sequence 1"
				, "A first sequence of the amino acid alphabet with the length " +
						"1-128 characters is required." 
				, StringList.class
				, "CLEMENS"));
		super.parameters.add(new AlgorithmParameter(
				"PAM / BLOSUM"
				, "Choose YES to use PAM or NO to use BLOSUM for scoring." 
				, Boolean.class 
				, new Boolean(true)));
		super.parameters.add(new AlgorithmParameter(
				"Gap costs"
				, "A double for the constant gap costs used for scoring."
				, Double.class 
				, new Double(-4.0)));

	}

	@Override
	public Vector<AlgorithmParameter> getInputParameters() {
		return super.parameters;
	}

	@Override
	public String getName() {
		return new String("Sum-Of-Pairs");
	}
	

	@Override
	public String getDescription() {
		return new String("Implemenation of the Sum-Of-Pairs approach " +
				"to multiple alignment with the help of dynamic programming.");
	}
	
	private Alignment calculateFirstColumn(Alignment algn) {
		double score = 0;
		for(int i = 0; i < algn.size(); i++) {
			for(int j = i+1; j < algn.size(); j++) {
				
				score += omega.getScore(algn.getSequence(i).charAt(0), algn.getSequence(j).charAt(0));
			}
		}
		return algn;
	}
	
	/**
	 * Main method of the algorithm.
	 * 
	 * @param params The filled out parameters are entered externally.
	 * 
	 * @return Output string containing the used parameters, the results and errors.
	 */
	@Override
	public String run(Vector<AlgorithmParameter> params) {
		
		String retVal = new String("");
		
		  // ##########  PARSE INPUT PARAMETERS FOR ERRORS  ###########
		
		//TODO Parse input for errors!!!
		String sequenceString = (String) params.elementAt(0).data;
		usePAM = (Boolean) params.elementAt(1).data;
		gapCosts = (Double) params.elementAt(2).data;
		
		omega = new SubstitutionMatrix(usePAM, gapCosts);
		
		retVal += "\n mhh.. I assume my default values are fine ! ;) \n";
		
		  // ##########  RUN THE PROGRAM  ###########

		 // (JUST TO EXEMPLIFY SOME OUTPUT I REPORT THE INPUT PARAMETERS ...)
		for (int i = 0; i < params.size(); i++) {
			retVal	+= "\n Input : "
					+ params.elementAt(i).name
					+ " = "
					+ params.elementAt(i).defVal.toString();
		}
				
		  // returning the algorithm results
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
		SumOfPairs myInstance = new SumOfPairs();
		  // run the example the instance with the default parameters
		BioinfAlgorithm.runAlgorithmDefaults( myInstance );
	}

}
