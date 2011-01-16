package cthoelken;

import gui.AlgorithmParameter;
import gui.BioinfAlgorithm;
import gui.StringList;
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

	private String seqenceData;
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
				, new StringList("\n\n!Kommentar 1 \nABADEREAGEERGAA\n>" +
						"Sequence 1 [Die Schwarmmücke] \"Plutonium Maximum\"" +
						"\n!Kommentar2\n\nABCDEFG\nABCDFGE\nABDGFEC\n>" +
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
		seqenceData = (String) params.elementAt(0).data;
		weighted = (Boolean) params.elementAt(1).data;
		usePAM = (Boolean) params.elementAt(2).data;
		gapCosts = (Integer) params.elementAt(3).data;
		
		retVal += "\n mhh.. I assume my default values are fine ! ;) \n";
		
		  // ##########  RUN THE PROGRAM  ###########

		 // (JUST TO EXEMPLIFY SOME OUTPUT I REPORT THE INPUT PARAMETERS ...)
		for (int i = 0; i < params.size(); i++) {
			retVal	+= "\n Input : "
					+ params.elementAt(i).name
					+ " = "
					+ params.elementAt(i).defVal.toString();
		}
		
		NeedlemanWunsch nw = new NeedlemanWunsch
		
		retVal += "\n Maximal Score: " + calculate();
		
		// print final cost matrix to console   VERBOSE!!!
		//TODO add to return string!!!
		System.out.println(M.toString());
		
		backtrack(seq1.length()-1, seq2.length()-1, new Alignment(2));
		
		// ### print alignments to return string ###
		for(int k=0; k<algnmts.size(); k++) 
			retVal += "\n########### Alignment "+k+":\n"+algnmts.get(k).toString();
		
		  // returning the algorithm results
		System.out.println("######################\n"+omega.getScore('b', 'b'));
		return retVal;
	}
	
	private void backtrack(int x, int y, Alignment algn) {
		
		System.out.println("backtrack: x="+x+" y="+y);
		
		if(x < 0 || y < 0) return;	// we are out of bounds, shouldn't happen
		
		char[] column = new char[2];
		char[] match = new char[1];
		
		if(x == 0 && y == 0) {	// we are done, great!
//			column[0] = seq1.charAt(x); column[1] = seq2.charAt(y);
//			match[0] = (seq1.charAt(x) == seq2.charAt(y)) ? '|' : '*';
//			algnmts.add(algn.addFirst(column, match));
			System.out.println("All is good!");
			algnmts.add(algn);
			return;
		}
		
		int[] psbl = {0, 0, 0};
		
		System.out.println("\n M(x,y) = "+M.get(x, y)+" - M(x,y-1) = "+M.get(x, y-1)+" - SEQ2: "+seq2.charAt(y)); //TODO
		
		// insertion
		if(M.get(x, y) == M.get(x-1, y) + omega.getScore(seq1.charAt(x), '_'))
			psbl[0] = 1;
		
		// deletion, start new thread from there
		if(M.get(x, y) == M.get(x, y-1) + omega.getScore('_', seq2.charAt(y)))
			psbl[1] = 1;

		// match/mismatch, lets carry on
		if(M.get(x, y) == M.get(x-1, y-1) + omega.getScore(seq1.charAt(x), seq2.charAt(y)))
			psbl[2] = 1;
		
		System.out.println("\n "+psbl[0]+""+psbl[1]+""+psbl[2]); //TODO
		
		int num_psbl = psbl[0]+psbl[1]+psbl[2];
		if(randomBackTrace && num_psbl > 1) {
			int random;
			if(num_psbl == 3) {
				psbl[0] = 0; psbl[1] = 0; psbl[2] = 0; psbl[new Random().nextInt(3)] = 1;
			}
			if(num_psbl == 2) {
				random = new Random().nextInt(2);
				if(psbl[0] == 1) {
					if(psbl[1] == 1) { 
						psbl[0] = random; psbl[1] = 1-random;
					} else { psbl[0] = random; psbl[2] = 1-random; }
				} else { psbl[1] = random; psbl[2] = 1-random; }
			}
			System.out.println("\n "+psbl[0]+""+psbl[1]+""+psbl[2]); //TODO
		}
		
		match[0] = ' ';
		
		if(psbl[0] == 1) {
			System.out.println("\n I go left");
			column[0] = seq1.charAt(x); column[1] = '_'; 
			backtrack(x-1, y, algn.addFirst(column, match));
		}
		if(psbl[1] == 1) {
			System.out.println("\n I go right");
			column[0] = '_'; column[1] = seq2.charAt(y);
			backtrack(x, y-1, algn.addFirst(column, match));
		}
		if(psbl[2] == 1) {
			System.out.println("\n I go diagonally");
			column[0] = seq1.charAt(x); column[1] = seq2.charAt(y);
			match[0] = (seq1.charAt(x) == seq2.charAt(y)) ? '|' : '*';
			backtrack(x-1, y-1, algn.addFirst(column, match));
		}
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
