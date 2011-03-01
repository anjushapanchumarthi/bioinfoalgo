package cthoelken;

import gui.AlgorithmParameter;
import gui.BioinfAlgorithm;
import gui.StringList;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.Vector;

 /**
  * Implementation of the Sum-Of-Pairs Algorithm for optimal multiple alignments
  * with linear gap costs.
  * 
  * @author Clemens Thoelken
  *
  */
public class SumOfPairs extends BioinfAlgorithm {

	protected String seq1;
	protected String seq2;
	protected boolean usePAM;
	protected double gapCosts;
	protected boolean randomBackTrace;
	protected SubstitutionMatrix omega;
	protected CostMatrix M;
	protected LinkedList<Alignment> algnmts = new LinkedList<Alignment>();
	protected String seq3;
	private Alignment sequences;

	 /**
	  * Constructor which generates an empty vector of parameters of the needed 
	  * types.
	  */
	public SumOfPairs() { 
		
		// create all needed parameters for the algorithm to work.
		
		super.parameters.add(new AlgorithmParameter(	
				"FASTA Input"
				, "Type or paste three sequences in FASTA format into the field." 
				, StringList.class
				, new StringList(">Sequenz 1;Kommentarzeile A" +
						"\nMTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKDT" +
						"\nLVSVKVSDDFTIAAMRPSYLSYEDLDMT" +
						"\n>Sequenz 2\n;Kommentarzeile B\n;Kommentarzeile C" +
						"\nSATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGV" +
						"\nATIGENLVVRRFATLKAGANGV" +
						"\n>Sequenz 3\n;Kommentarzeile B\n;Kommentarzeile C" +
						"\nGATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINDT" +
						"\nATIGENLVVRRFA")));
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
				, new Boolean(true)));
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
		return new String("Sum-Of-Pairs");
	}
	
	@Override
	public String getDescription() {
		return new String("Implemenation of the Sum-Of-Pairs approach " +
				"to multiple alignment with the help of dynamic programming.");
	}
	
	/** Scores a pre-aligned alignment based on SOP-Score
	 * @param algn Pre-aligned alignment of same lenght sequences
	 * @param usePAM Use PAM for substitution
	 * @param gapCosts Gap costs for the alignment
	 * @return Returns the total score
	 */
	public static double score(Alignment algn, boolean usePAM, double gapCosts) {
		SubstitutionMatrix omega = new SubstitutionMatrix(usePAM, gapCosts);
		double score = 0;
		for(int x = 0; x < algn.size(); x++) {
			if(algn.getSeq(x).length() != algn.getSeq(0).length())
				throw new IllegalArgumentException("Sequences in alignment are not of the same length!");
			else for(int y = x+1; y < algn.size(); y++)
				for(int z = 0; z < algn.getSeq(0).length(); z++)
					score += omega.getScore(algn.getSeq(x).charAt(z), algn.getSeq(y).charAt(z));
		}
		return score;
	}
	
	/** Fills the CostMatrix with values and returns the optimal score
	 * @return Optimal Score
	 */
	private double score() {
		seq1 = "#" + seq1; seq2 = "#" + seq2; seq3 = "#" + seq3; //increase sequence length, disregarded afterwards
		M = new CostMatrix(seq1.length(), seq2.length(), seq3.length());
		double[] score = new double[7];
		for(int x = 0; x < seq1.length(); x++) {
			for(int y = 0; y < seq2.length(); y++) {
				for(int z = 0; z < seq3.length(); z++) {
					if(x+y+z != 0) {
						score[0] = M.get(x-1, y-1, z-1) + omega.getScore(// 111
								seq1.charAt(x), 
								seq2.charAt(y), 
								seq3.charAt(z));
						score[1] = M.get(x-1, y-1, z) + omega.getScore(	// 11_
								seq1.charAt(x), seq2.charAt(y)) + 2*gapCosts;
						score[2] = M.get(x-1, y, z-1) + omega.getScore(	// 1_1
								seq1.charAt(x), seq3.charAt(z)) + 2*gapCosts;
						score[3] = M.get(x, y-1, z-1) + omega.getScore(	// _11
								seq2.charAt(y), seq3.charAt(z)) + 2*gapCosts;
						score[4] = M.get(x-1, y, z) + 2*gapCosts;		// 1__
						score[5] = M.get(x, y-1, z) + 2*gapCosts;		// _1_
						score[6] = M.get(x, y, z-1) + 2*gapCosts;		// __1

						M.set(x, y, z, new Double(Util.maxValue(score[0], score[1], score[2], score[3], score[4], score[5], score[6])));
					}
				}
			}
		}
		return M.score();
	}
	
	/** Backtracks the CostMatrix for feasable paths and adds them to the internal alignment list
	 * @param x x-offset
	 * @param y y-offset
	 * @param z z-offset
	 * @param algn Initial alignment, starts with an empty alignment
	 */
	private void backtrack(int x, int y, int z, Alignment algn) {
		if(x < 0 || y < 0 || z < 0) return;
		
		if(x == 0 && y == 0 && z == 0) {	// we are done, great!
			algnmts.add(algn);
			return;
		}
		
		int[] psbl = {0, 0, 0, 0, 0, 0, 0};
		if(M.get(x, y, z) == M.get(x-1, y-1, z-1) // perfect match
				+ omega.getScore(seq1.charAt(x), seq2.charAt(y), seq3.charAt(z)))
			psbl[0] = 1;
		if(M.get(x, y, z) == M.get(x-1, y-1, z)   // match in X and Y
				+ omega.getScore( seq1.charAt(x), seq2.charAt(y)) + 2*gapCosts)
			psbl[1] = 1;
		if(M.get(x, y, z) == M.get(x-1, y, z-1)   // match in X and Z
				+ omega.getScore( seq1.charAt(x), seq3.charAt(z)) + 2*gapCosts)
			psbl[2] = 1;
		if(M.get(x, y, z) == M.get(x, y-1, z-1)   // match in Y and Z
				+ omega.getScore( seq2.charAt(y), seq3.charAt(z)) + 2*gapCosts)
			psbl[3] = 1;
		if(M.get(x, y, z) == M.get(x-1, y, z) + 2*gapCosts) psbl[4] = 1; // insert in Y and Z
		if(M.get(x, y, z) == M.get(x, y-1, z) + 2*gapCosts) psbl[5] = 1; // insert in X and Z
		if(M.get(x, y, z) == M.get(x, y, z-1) + 2*gapCosts) psbl[6] = 1; // insert in Y and Z
		
		char[] column = new char[3];
		char[] match = new char[2];
		
		match[0] = ' '; match[1] = '_';
		
		if(randomBackTrace == true) psbl = Util.chooseRandom(psbl); //Random BT
		
		if(psbl[0] == 1) {
			match[0] = (seq1.toUpperCase().charAt(x) == seq2.toUpperCase().charAt(y)) ? '|' : '*';
			match[1] = (seq2.toUpperCase().charAt(y) == seq3.toUpperCase().charAt(z)) ? '|' : '*';
			column[0] = seq1.charAt(x);
			column[1] = seq2.charAt(y);
			column[2] = seq3.charAt(z);
			backtrack(x-1, y-1, z-1, algn.addFirst(column, match));
		}
		if(psbl[1] == 1) {
			match[0] = (seq1.toUpperCase().charAt(x) == seq2.toUpperCase().charAt(y)) ? '|' : '*';
			match[1] = ' ';
			column[0] = seq1.charAt(x);
			column[1] = seq2.charAt(y);
			column[2] = '_';
			backtrack(x-1, y-1, z, algn.addFirst(column, match));
		}
		if(psbl[2] == 1) {
			match[0] = ' ';
			match[1] = ' ';
			column[0] = seq1.charAt(x);
			column[1] = '_';
			column[2] = seq3.charAt(z);
			backtrack(x-1, y, z-1, algn.addFirst(column, match));
		}
		if(psbl[3] == 1) {
			match[0] = ' ';
			match[1] = (seq2.toUpperCase().charAt(y) == seq3.toUpperCase().charAt(z)) ? '|' : '*';
			column[0] = '_';
			column[1] = seq2.charAt(y);
			column[2] = seq3.charAt(z);
			backtrack(x, y-1, z-1, algn.addFirst(column, match));
		}
		if(psbl[4] == 1) {
			match[0] = ' ';
			match[1] = ' ';
			column[0] = seq1.charAt(x);
			column[1] = '_';
			column[2] = '_';
			backtrack(x-1, y, z, algn.addFirst(column, match));
		}
		if(psbl[5] == 1) {
			match[0] = ' ';
			match[1] = ' ';
			column[0] = '_';
			column[1] = seq2.charAt(y);
			column[2] = '_';
			backtrack(x, y-1, z, algn.addFirst(column, match));
		}
		if(psbl[6] == 1) {
			match[0] = ' ';
			match[1] = ' ';
			column[0] = '_';
			column[1] = '_';
			column[2] = seq3.charAt(z);
			backtrack(x, y, z-1, algn.addFirst(column, match));
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
		
		//TODO Parse input for errors!!!
		try {
			sequences = Util.parseFasta((StringList) params.elementAt(0).data);
		} catch(IllegalArgumentException e) { return ""+e; }
		if(sequences.size() != 3)
			return "Please enter 3 valid Sequences in the FASTA format above!";
		seq1 = sequences.getSeq(0);
		seq2 = sequences.getSeq(1);
		seq3 = sequences.getSeq(2);
		
		usePAM = (Boolean) params.elementAt(1).data;
		randomBackTrace = (Boolean) params.elementAt(2).data;
		
		try{
			if(params.elementAt(3).type == Double.class)
				gapCosts = (Double) params.elementAt(3).data;
			else return "Gap costs are not a valid decimal value!";
		} catch(Exception e) {return "Gap costs are not a valid decimal value!";}
		
		omega = new SubstitutionMatrix(usePAM, gapCosts);
		
		  // ##########  RUN THE PROGRAM  ###########
		
		System.out.println("\nScore: "+score());
		
		backtrack(seq1.length()-1, seq2.length()-1, seq3.length()-1, new Alignment(3));
		
		// ### print alignments to return string ###
		for(int k=0; k<algnmts.size(); k++) 
			retVal += "\n########### Alignment "+(k+1)+":\n"+algnmts.get(k).toString();
				
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
