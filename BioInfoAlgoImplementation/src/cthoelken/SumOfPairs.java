package cthoelken;

import gui.AlgorithmParameter;
import gui.BioinfAlgorithm;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.Vector;

 /**
  * Implementation of the Needleman Wunsch Algorithm for pairwise alignments
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
				, String.class
				, "ABCE"));
		super.parameters.add(new AlgorithmParameter(	
				"Sequence 1"
				, "A first sequence of the amino acid alphabet with the length " +
						"1-128 characters is required." 
				, String.class
				, "ABDE"));
		super.parameters.add(new AlgorithmParameter(	
				"Sequence 1"
				, "A first sequence of the amino acid alphabet with the length " +
						"1-128 characters is required." 
				, String.class
				, "BDEE"));
		super.parameters.add(new AlgorithmParameter(
				"PAM / BLOSUM"
				, "Choose YES to use PAM or NO to use BLOSUM for scoring." 
				, Boolean.class 
				, new Boolean(true)));
		super.parameters.add(new AlgorithmParameter(
				"Random / Exhaustive Backtracking"
				, "Choose YES to output only one random optimal alignment or NO" +
						"  to output all optimal alignments at once." 
				, Boolean.class 
				, new Boolean(false)));
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
	
	private double score() {
		M = new CostMatrix(seq1.length()-1, seq2.length()-1, seq3.length()-1);
		double[] score = new double[7];
		for(int x = 0; x < seq1.length()-1; x++) {
			for(int y = 0; y < seq2.length()-1; y++) {
				for(int z = 0; z < seq3.length()-1; z++) {
					if(x+y+z != 0) {
						score[0] = M.get(x-1, y-1, z-1) + omega.getScore(// 111
								seq1.charAt(x), 
								seq2.charAt(y), 
								seq3.charAt(z));
						score[1] = M.get(x-1, y-1, z) + omega.getScore(	// 11_
								seq1.charAt(x), 
								seq2.charAt(y)) + 2*gapCosts;
						score[2] = M.get(x-1, y, z-1) + omega.getScore(	// 1_1
								seq1.charAt(x), 
								seq3.charAt(z)) + 2*gapCosts;
						score[3] = M.get(x, y-1, z-1) + omega.getScore(	// _11
								seq2.charAt(y), 
								seq3.charAt(z)) + 2*gapCosts;
						score[4] = M.get(x-1, y, z) + 2*gapCosts;		// 1__
						score[5] = M.get(x, y-1, z) + 2*gapCosts;		// _1_
						score[6] = M.get(x, y, z-1) + 2*gapCosts;		// __1

						M.set(x, y, z, new Double(NeedlemanWunsch.maxValue(score[0], score[1], score[2], score[3], score[4], score[5], score[6])));
						System.out.println("M["+(x)+"]["+(y)+"]["+(z)+"]="+M.get(x, y, z));
					}
				}
			}
		}
		return M.score();
	}
	
	private void backtrack(int x, int y, int z, Alignment algn) {
		if(x < 0 || y < 0 || z < 0) return;
		
		if(x == 0 && y == 0 && z == 0) {	// we are done, great!
			algnmts.add(algn);
			return;
		}
		
		System.out.println("M["+(x-1)+"]["+(y-1)+"]["+(z-1)+"]="+M.get(x, y, z));
		
		int[] psbl = {0, 0, 0, 0, 0, 0, 0};
		if(M.get(x, y, z) == M.get(x-1, y-1, z-1) 
				+ omega.getScore(seq1.charAt(x), seq2.charAt(y), seq3.charAt(z)))
			psbl[0] = 1;
		if(M.get(x, y, z) == M.get(x-1, y-1, z) 
				+ omega.getScore( seq1.charAt(x), seq2.charAt(y)) + 2*gapCosts)
			psbl[1] = 1;
		if(M.get(x, y, z) == M.get(x-1, y, z-1) 
				+ omega.getScore( seq1.charAt(x), seq2.charAt(z)) + 2*gapCosts)
			psbl[2] = 1;
		if(M.get(x, y, z) == M.get(x, y-1, z-1) 
				+ omega.getScore( seq1.charAt(y), seq2.charAt(z)) + 2*gapCosts)
			psbl[3] = 1;
		if(M.get(x, y, z) == M.get(x-1, y, z) + 2*gapCosts)
			psbl[4] = 1;
		if(M.get(x, y, z) == M.get(x, y-1, z) + 2*gapCosts)
			psbl[5] = 1;
		if(M.get(x, y, z) == M.get(x, y, z-1) + 2*gapCosts)
			psbl[6] = 1;
//
//		int num_psbl = psbl[0]+psbl[1]+psbl[2]+psbl[3]+psbl[4]+psbl[5]+psbl[6];
//		if(randomBackTrace && num_psbl > 1) {
//			int random;
//			if(num_psbl == 7) {
//				psbl[0] = 0; psbl[1] = 0; psbl[2] = 0; psbl[3] = 0; 
//				psbl[4] = 0; psbl[5] = 0; psbl[6] = 0;
//				psbl[new Random().nextInt(7)] = 1;
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
		
		System.out.println(Arrays.toString(psbl));
		
		char[] column = new char[3];
		char[] match = new char[2];
		
		match[0] = ' '; match[1] = '_';
		
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
		seq1 = "*"+(String) params.elementAt(0).data;
		seq2 = "*"+(String) params.elementAt(1).data;
		seq3 = "*"+(String) params.elementAt(2).data;
		usePAM = (Boolean) params.elementAt(3).data;
		randomBackTrace = (Boolean) params.elementAt(4).data;
		gapCosts = (Double) params.elementAt(5).data;
		
		omega = new SubstitutionMatrix(usePAM, gapCosts);
		
		retVal += "\n mhh.. I assume my default values are fine ! ;) \n";
		
		  // ##########  RUN THE PROGRAM  ###########
		
		System.out.println("\nScore: "+score());
		System.out.println(M.toString());
		
		backtrack(seq1.length()-2, seq2.length()-2, seq3.length()-2, new Alignment(3));
		
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
		SumOfPairs myInstance = new SumOfPairs();
		  // run the example the instance with the default parameters
		BioinfAlgorithm.runAlgorithmDefaults( myInstance );
	}

}
