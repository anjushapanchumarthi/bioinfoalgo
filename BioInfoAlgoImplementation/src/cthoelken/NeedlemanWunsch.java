package cthoelken;

import gui.AlgorithmParameter;
import gui.BioinfAlgorithm;

import java.util.LinkedList;
import java.util.Vector;

 /**
  * Example implementation of the BioinfAlgorithm interface.
  * 
  * @author Martin Mann - http://www.bioinf.uni-freiburg.de/~mmann/
  *
  */
public class NeedlemanWunsch extends BioinfAlgorithm {

	private String seq1;
	private String seq2;
	private boolean usePAM;
	private int gapCosts;
	private SubstitutionMatrix omega;
	private CostMatrix M;
	private LinkedList<Alignment> algnmts = new LinkedList<Alignment>();

	 /**
	  * Constructs an example algorithm and initializes all allowed parameters.
	  */
	public NeedlemanWunsch() { 
		
		// create one example parameter for each type of parameter allowed
		
		super.parameters.add(new AlgorithmParameter(	
				"Sequence 1"
				, "A first sequence of the amino acid alphabet with the length " +
						"1-128 characters is required." 
				, String.class
				, "ABCDEF"));
		super.parameters.add(new AlgorithmParameter(	
				"Sequence 2"
				, "A second sequence of the amino acid alphabet with the length " +
						"1-128 characters is required." 
				, String.class
				, "FABBDEB"));
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
	
	@Override
	public String run(Vector<AlgorithmParameter> params) {
		
		String retVal = new String("");
		
		  // ##########  PARSE INPUT PARAMETERS FOR ERRORS  ###########
		
		seq1 = (String) params.elementAt(0).data;
		seq2 = (String) params.elementAt(1).data;
		usePAM = (Boolean) params.elementAt(2).data;
		gapCosts = (Integer) params.elementAt(3).data;
		
		omega = new SubstitutionMatrix(usePAM, gapCosts);
		M = new CostMatrix(seq1.length(), seq2.length(), omega);
		
		retVal += "\n mhh.. I assume my default values are fine ! ;) \n";
		
		  // ##########  RUN THE PROGRAM  ###########

		 // (JUST TO EXEMPLIFY SOME OUTPUT I REPORT THE INPUT PARAMETERS ...)
		for (int i = 0; i < params.size(); i++) {
			retVal	+= "\n Input : "
					+ params.elementAt(i).name
					+ " = "
					+ params.elementAt(i).defVal.toString();
		}
		
		System.out.println(M.toString());
		
		int top = 0;
		int left = 0;
		int diagonal = 0;
		int max;
		
		for(int i=0; i<seq1.length(); i++) {
			for(int j=0; j<seq2.length(); j++) {
				top = M.get(i, j-1) + omega.getScore(seq1.charAt(i), '_');
				left = M.get(i-1, j) + omega.getScore('_', seq2.charAt(j));
				diagonal = M.get(i-1, j-1) + omega.getScore(seq1.charAt(i), seq2.charAt(j));
				max = top;
				if(max<left) max = left;
				if(max<diagonal) max = diagonal;
				if(i+j!=0) M.set(i, j, max);
			}
		}
		
		System.out.println(M.toString());
		
		new Thread(new Backtracer(new Alignment(2), seq1.length()-1, seq2.length()-1)).start();
		
		for(int k=0; k<algnmts.size(); k++) 
			retVal += algnmts.get(k).toString();
		  // returning the algorithm results ...
		return retVal;
	}
	
	private class Backtracer implements Runnable {
		
		private int x;
		private int y;
		private Alignment algn;
		
		Backtracer(Alignment algn, int x, int y) {
			this.x = x; this.y = y;
			this.algn = algn;
		}

		@Override
		public void run() {
			boolean deadEnd = false;
			System.out.println("Bin frisch geschlüpft: "+algn.toString());
			char[] column = new char[algn.sequences.length];
			while(!deadEnd) {
				if(x == 0 && y == 0) {
					algnmts.add(algn);
					deadEnd = true; break;
				}
				System.out.println("Ich renn! "+algn.toString());
				if(M.get(x, y) == M.get(x-1, y) + omega.getScore(seq1.charAt(x), '_')) {
					column[0] = seq1.charAt(x); column[1] = '_';
					new Thread(new Backtracer(new Alignment(algn).addFirst(column), x-1, y)).start();	
				}
				if(M.get(x, y) == M.get(x, y-1) + omega.getScore('_', seq2.charAt(y))) {
					column[1] = seq1.charAt(y); column[0] = '_';
					new Thread(new Backtracer(new Alignment(algn).addFirst(column), x, y-1)).start();	
				}

				if(M.get(x, y) != M.get(x-1, y-1) + omega.getScore(seq1.charAt(x), seq2.charAt(y)))
					deadEnd = true;
			}
			
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
		NeedlemanWunsch myInstance = new NeedlemanWunsch();
		  // run the example the instance with the default parameters
		BioinfAlgorithm.runAlgorithmDefaults( myInstance );
	}

}
