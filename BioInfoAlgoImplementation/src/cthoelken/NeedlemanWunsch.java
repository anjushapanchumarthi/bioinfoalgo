package cthoelken;

import gui.AlgorithmParameter;
import gui.BioinfAlgorithm;
import java.util.LinkedList;
import java.util.Vector;

 /**
  * Implementation of the Needleman Wunsch Algorithm for pairwise alignments
  * with linear gap costs.
  * 
  * @author Clemens Thoelken
  *
  */
public class NeedlemanWunsch extends BioinfAlgorithm {

	protected String seq1;
	protected String seq2;
	protected boolean usePAM;
	protected double gapCosts;
	protected boolean randomBackTrace;
	protected SubstitutionMatrix omega;
	protected CostMatrix M;
	protected LinkedList<Alignment> algnmts = new LinkedList<Alignment>();

	 /**
	  * Constructor which generates an empty vector of parameters of the needed 
	  * types.
	  */
	public NeedlemanWunsch() { 
		
		// create all needed parameters for the algorithm to work.
		
		super.parameters.add(new AlgorithmParameter(	
				"Sequence 1"
				, "A first sequence of the amino acid alphabet with the length " +
						"1-128 characters is required." 
				, String.class
				, "MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAA" +
						"EGLVSVKVSDDFTI"));
		super.parameters.add(new AlgorithmParameter(	
				"Sequence 2"
				, "A second sequence of the amino acid alphabet with the length " +
						"1-128 characters is required." 
				, String.class
				, "SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLK" +
						"SQIATIGENLVVR"));
		super.parameters.add(new AlgorithmParameter(
				"use PAM (BLOSUM otherwise)"
				, "Choose YES to use PAM or NO to use BLOSUM for scoring." 
				, Boolean.class 
				, new Boolean(true)));
		super.parameters.add(new AlgorithmParameter(
				"single random alignment (all alignments otherwise)"
				, "Choose YES to output only one random optimal alignment or NO" +
						"  to output all optimal alignments at once." 
				, Boolean.class 
				, new Boolean(false)));
		super.parameters.add(new AlgorithmParameter(
				"Gap costs"
				, "A double for the constant gap costs used for scoring."
				, Double.class 
				, new Double(-1.0)));

	}

	@Override
	public Vector<AlgorithmParameter> getInputParameters() {
		return super.parameters;
	}

	@Override
	public String getName() {
		return new String("Needleman-Wunsch");
	}
	

	@Override
	public String getDescription() {
		return new String("Implemenation of the Needleman-Wunsch approach " +
				"to pairwise alignment with the help of dynamic programming.");
	}
	
	/** Get the score of two Strings from external algorithms
	 * @param s1 Sequence string 1
	 * @param s2 sequence string 2
	 * @param usePAM Use PAM for substitution
	 * @param gapCosts Gap costs for the alignment
	 * @return Score of the alignment
	 */
	public double getScore(String s1, String s2, boolean usePAM, double gapCosts) {
		seq1 = s1; seq2 = s2; this.usePAM = usePAM; this.gapCosts = gapCosts;
		omega = new SubstitutionMatrix(usePAM, gapCosts);
		M = new CostMatrix(seq1.length()+1, seq2.length()+1);
		return calculate();
	}
	
	/** Get the Alignment of two strings from external algorithms
	 * @param s1 Sequence string 1
	 * @param s2 sequence string 2
	 * @param usePAM Use PAM for substitution
	 * @param gapCosts Gap costs for the alignment
	 * @return The actual alignment
	 */
	public Alignment getAlignment(String s1, String s2, boolean usePAM, double gapCosts) {
		seq1 = s1; seq2 = s2; this.usePAM = usePAM; this.gapCosts = gapCosts;
		randomBackTrace = true;
		omega = new SubstitutionMatrix(usePAM, gapCosts);
		M = new CostMatrix(seq1.length()+1, seq2.length()+1);
		calculate(); backtrack(seq1.length()-1, seq2.length()-1, new Alignment(2));
		return algnmts.getFirst();
	}
	
	/** Calculates the costmatrix
	 * @return Returns the score of the alignment
	 */
	private double calculate() {
		
		seq1 = "#" + seq1; seq2 = "#" + seq2; //increase sequence length, disregarded afterwards
		
		// fill the cost matrix row-wise
		for(int i=0; i<seq1.length(); i++) {
			for(int j=0; j<seq2.length(); j++) {
				if(!(i==0 && j==0))					
					M.set(i, j, Util.maxValue(M.get(i, j-1) + gapCosts,
							M.get(i-1, j) + gapCosts, M.get(i-1, j-1)
							+ omega.getScore(seq1.charAt(i), seq2.charAt(j))));
			}
		}
		return M.score();
	}
	

	/** Backtracks the cost matrix for feasable paths
	 * @param x x-offset
	 * @param y y-offset
	 * @param algn The alginment so far, starts with empty alignment
	 */
	private void backtrack(int x, int y, Alignment algn) {
		
		if(x < 0 || y < 0) return;	// we are out of bounds, shouldn't happen
		
		char[] column = new char[2];
		char[] match = new char[1];
		
		if(x == 0 && y == 0) {	// we are done, great!
			algnmts.add(algn);
			return;
		}
		
		int[] psbl = {0, 0, 0};
		
		// insertion
		if(M.get(x, y) == M.get(x-1, y) + gapCosts) psbl[0] = 1;
		
		// deletion, start new thread from there
		if(M.get(x, y) == M.get(x, y-1) + gapCosts) psbl[1] = 1;
		
		// match/mismatch, lets carry on
		if(M.get(x, y) == M.get(x-1, y-1) + omega.getScore(seq1.charAt(x), seq2.charAt(y)))
			psbl[2] = 1;
		
		//int num_psbl = psbl[0]+psbl[1]+psbl[2];
		if(randomBackTrace == true) psbl = Util.chooseRandom(psbl);
//		if(randomBackTrace && num_psbl > 1) {
//			int random;
//			if(num_psbl == 3) {
//				random = new Random().nextInt(3);
//				psbl[0] = 0; psbl[1] = 0; psbl[2] = 0; psbl[random] = 1;
//			}
//			if(num_psbl == 2) {
//				random = new Random().nextInt(2);
//				if(psbl[0] == 1) {
//					if(psbl[1] == 1) { 
//						psbl[0] = random; psbl[1] = 1-random;
//					} else { psbl[0] = random; psbl[2] = 1-random; }
//				} else { psbl[1] = random; psbl[2] = 1-random; }
//			}
//		}
		
		match[0] = ' ';
		
		if(psbl[0] == 1) {
			column[0] = seq1.charAt(x); column[1] = '_'; 
			backtrack(x-1, y, algn.addFirst(column, match));
		}
		if(psbl[1] == 1) {
			column[0] = '_'; column[1] = seq2.charAt(y);
			backtrack(x, y-1, algn.addFirst(column, match));
		}
		if(psbl[2] == 1) {
			column[0] = seq1.charAt(x); column[1] = seq2.charAt(y);
			match[0] = (seq1.toUpperCase().charAt(x) == seq2.toUpperCase().charAt(y)) ? '|' : '*';
			backtrack(x-1, y-1, algn.addFirst(column, match));
		}
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

		seq1 = (String) params.elementAt(0).data;
		if(!Util.isValidSequence(seq1)) return "Sequence 1 is not valid!";
		
		seq2 = (String) params.elementAt(1).data;
		if(!Util.isValidSequence(seq2)) return "Sequence 2 is not valid!";
		
		usePAM = (Boolean) params.elementAt(2).data;
		randomBackTrace = (Boolean) params.elementAt(3).data;
		
		try{
			if(params.elementAt(4).data.getClass() == Double.class)
				gapCosts = (Double) params.elementAt(4).data;
			else return "Gap costs are not a valid decimal value!";
		} catch(Exception e) {return "Gap costs are not a valid decimal value!";}
		
		omega = new SubstitutionMatrix(usePAM, gapCosts);
		M = new CostMatrix(seq1.length()+1, seq2.length()+1);
		
		  // ##########  RUN THE PROGRAM  ###########
		
		retVal += "\n Maximal Score: " + calculate();
		
		backtrack(seq1.length()-1, seq2.length()-1, new Alignment(2));
		
		// ### print alignments to return string ###
		for(int k=0; k<algnmts.size(); k++) 
			retVal += "\n########### Alignment "+k+":\n"+algnmts.get(k).toString();
		
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
		NeedlemanWunsch myInstance = new NeedlemanWunsch();
		  // run the example the instance with the default parameters
		BioinfAlgorithm.runAlgorithmDefaults( myInstance );
	}

}
