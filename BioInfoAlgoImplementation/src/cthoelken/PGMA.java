package cthoelken;

import gui.AlgorithmParameter;
import gui.BioinfAlgorithm;
import gui.StringList;

import java.util.LinkedList;
import java.util.Random;
import java.util.Vector;

 /**
  * Implementation of the Weighted and Unweighted Pair Group Method using the 
  * Arithmetic mean for multiple alignments.
  * 
  * @author Clemens Thoelken
  *
  */
public class PGMA extends BioinfAlgorithm {

	private String sequenceData;
	private Boolean usePAM;
	private Boolean weighted;
	private Integer gapCosts;


	/**
	 * Constructor which generates an empty vector of parameters of the needed 
	 * types.
	 */
	public PGMA() { 
		// create additional parameters for the algorithm to work.
		super.parameters.add(new AlgorithmParameter(
				"PAM / BLOSUM"
				, "Choose YES to use PAM or NO to use BLOSUM for scoring." 
				, StringList.class 
				, new StringList("\n\n;Kommentar 1 \nABADEREAGEERGAA\n>" +
						"Sequence 1 [Die Schwarmm�ücke] \"Plutonium Maximum\"" +
						"\n;Kommentar2\n\nABCDEFG\nABCDFGE\nABDGFEC\n>" +
						"Sequence 2\nABGEFGFEGEFCEG\nG*\n>Sequnce 3\n" +
						"ABGEGEGEGEGEGEGEG\n!Kommentar 3\nAGEGEGEGE**")));
		super.parameters.add(new AlgorithmParameter(
				"Weighted / Unweighted"
				, "Choose YES to use WEIGHTED or NO to use UNWEIGHTED pairing." 
				, Boolean.class 
				, new Boolean(true)));
		super.parameters.add(new AlgorithmParameter(
				"PAM / BLOSUM"
				, "Choose YES to use PAM or NO to use BLOSUM for scoring." 
				, Boolean.class 
				, new Boolean(true)));
		super.parameters.add(new AlgorithmParameter(
				"Gap costs"
				, "An integer for the constant gap costs used for scoring."
				, Integer.class 
				, new Integer(-1)));
	}

	@Override
	public String getName() {
		return new String("UPGMA / WPGMA");
	}

	@Override
	public String getDescription() {
		return new String("Implementation of the Weighted and Unweighted " +
				"Pair Group Method using the Arithmetic mean for multiple " +
				"alignments.");
	}
	
	@Override
	public Vector<AlgorithmParameter> getInputParameters() {
		return super.parameters;
	}
	
	public Alignment parseFasta(String fasta) {
		Alignment algn = new Alignment(2);
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
		Alignment sequences = parseFasta((String) params.elementAt(0).data);
		weighted = (Boolean) params.elementAt(1).data;
		usePAM = (Boolean) params.elementAt(2).data;
		gapCosts = (Integer) params.elementAt(3).data;
		
		  // ##########  RUN THE PROGRAM  ###########

		 // (JUST TO EXEMPLIFY SOME OUTPUT I REPORT THE INPUT PARAMETERS ...)
		for (int i = 0; i < params.size(); i++) {
			retVal	+= "\n Input : "
					+ params.elementAt(i).name
					+ " = "
					+ params.elementAt(i).defVal.toString();
		}
		
		CostMatrix D = new CostMatrix(sequences.size(), sequences.size());
		NeedlemanWunsch nw = new NeedlemanWunsch();
		int[] bestPair = new int[2];
		int[] ranks = new int[sequences.size()];
		int max = 0;
		
		for(int i = 0; i <= sequences.size(); i++) {
			ranks[i] = sequences.size();
			for(int j = 0; j <= sequences.size()-1; i++) {
				D.set(i, j, (int) nw.getScore(sequences.getSequence(i), 
						sequences.getSequence(j), usePAM, gapCosts));
				if(D.get(i, j) > max && i != j) {
					max = D.get(i, j);
					bestPair[0] = i; bestPair[1] = j;
				}
			}
		}
		ranks[bestPair[0]] = 0; ranks[bestPair[1]] = 0;
		for(int k = 0; k <= sequences.size(); k++) {
			//TODO Hab keinen Schimmer was ich hier noch machen will
		}
		
		boolean notFinished = true;
		while(notFinished) {
			
		}
		
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
		Gotoh myInstance = new Gotoh();
		  // run the example the instance with the default parameters
		BioinfAlgorithm.runAlgorithmDefaults( myInstance );
	}

}
