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

	private String seq1;
	private String seq2;
	private boolean usePAM;
	private int gapCosts;
	private SubstitutionMatrix omega;
	private CostMatrix M;
	private LinkedList<Alignment> algnmts = new LinkedList<Alignment>();

	 /**
	  * Constructor which generates an empty vector of parameters of the needed types.
	  */
	public NeedlemanWunsch() { 
		
		// create all needed parameters for the algorithm to work.
		
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
		seq1 = (String) params.elementAt(0).data;
		seq2 = (String) params.elementAt(1).data;
		usePAM = (Boolean) params.elementAt(2).data;
		gapCosts = (Integer) params.elementAt(3).data;
		
		omega = new SubstitutionMatrix(usePAM, gapCosts);
		M = new CostMatrix(seq1.length()+1, seq2.length()+1, omega);
		
		retVal += "\n mhh.. I assume my default values are fine ! ;) \n";
		
		  // ##########  RUN THE PROGRAM  ###########

		 // (JUST TO EXEMPLIFY SOME OUTPUT I REPORT THE INPUT PARAMETERS ...)
		for (int i = 0; i < params.size(); i++) {
			retVal	+= "\n Input : "
					+ params.elementAt(i).name
					+ " = "
					+ params.elementAt(i).defVal.toString();
		}
		
		//TODO print empty cost matrix to console   VERBOSE!!!
		System.out.println(M.toString());
		
		int top = 0;		//score from above
		int left = 0;		//score from left
		int diagonal = 0;	//score diagonally (top-left)
		int max;			//temporary maximum score
		seq1 = "#" + seq1; seq2 = "#" + seq2; //increase sequence length, disregarded afterwards
		
		// fill the cost matrix row-wise
		for(int i=0; i<seq1.length(); i++) {
			for(int j=0; j<seq2.length(); j++) {
				if(!(i==0 && j==0)) {
					// calculate the three previous cells
					top = M.get(i, j-1) + omega.getScore(seq1.charAt(i), '_');
					left = M.get(i-1, j) + omega.getScore('_', seq2.charAt(j));
					diagonal = M.get(i-1, j-1) + omega.getScore(seq1.charAt(i), seq2.charAt(j));
					
					max = top;						// take
					if(max<left) max = left;		// the
					if(max<diagonal) max = diagonal;// maximum
					
					M.set(i, j, max);	//update score
				}
			}
		}
		
		// print final cost matrix to console   VERBOSE!!!
		//TODO add to return string!!!
		System.out.println(M.toString());
		
		Thread bt = new Thread(new Backtracer(new Alignment(2), seq1.length()-1, seq2.length()-1));
		bt.start();
		try {
			bt.join();
		} catch (InterruptedException e) {
			// TODO some kind of synchronization should be handled!!!
			e.printStackTrace();
		}
		
		// ### print alignments to return string ###
		for(int k=0; k<algnmts.size(); k++) 
			retVal += "\n########### Alignment "+k+":\n"+algnmts.get(k).toString();
		
		  // returning the algorithm results
		return retVal;
	}
	
	// auxilliary class for parallel semi-recursive backtracking
	private class Backtracer implements Runnable {
		
		private int x;	// position in sequence 1
		private int y;	// position in sequence 2
		private Alignment algn;
		
		Backtracer(Alignment algn, int x, int y) {
			this.x = x; this.y = y;
			this.algn = algn;
		}

		@Override
		public void run() {
			
			boolean deadEnd = false;	// termination variable
			
			//TODO could also be three single chars instead of arrays
			char[] column = new char[2];	// Symbols in the sequences
			char[] match = new char[1];		// match, mismatch or gap
			
			// to join indetermined number of threads, create a ThreadGroup
			ThreadGroup tg = new ThreadGroup(algn.toString());
			
			// go diagonally, if other routs (insertion or deletion) are possible
			// start a new thread with temporary alignment from there, die if diagonal
			// move is not possible.
			while(!deadEnd) {
				
				deadEnd = true;		// generally we want to die
				
				if(x < 0 || y < 0) break; // we are BEYOND done (off the cliff)
				
				if(x == 0 && y == 0) {	// we are done, great!
					column[0] = seq1.charAt(x); column[1] = seq2.charAt(y);
					match[0] = (seq1.charAt(x) == seq2.charAt(y)) ? '|' : '*';
					algnmts.add(algn.addFirst(column, match));
					break;
				}

				match[0] = ' ';
				
				// insertion, start new thread from there
				if(M.get(x, y) == M.get(x-1, y) + omega.getScore(seq1.charAt(x), '_')) {
					column[0] = seq1.charAt(x); column[1] = '_'; 
					new Thread(tg, new Backtracer(new Alignment(algn).addFirst(column, match), x-1, y)).start();	
				}
				
				// deletion, start new thread from there
				if(M.get(x, y) == M.get(x, y-1) + omega.getScore('_', seq2.charAt(y))) {
					column[0] = '_'; column[1] = seq2.charAt(y);
					new Thread(tg, new Backtracer(new Alignment(algn).addFirst(column, match), x, y-1)).start();	
				}

				// match/mismatch, lets carry on
				if(M.get(x, y) == M.get(x-1, y-1) + omega.getScore(seq1.charAt(x), seq2.charAt(y))) {
					deadEnd = false;	// HAH! not dead yet! I'm a survivor!
					column[0] = seq1.charAt(x); column[1] = seq2.charAt(y);
					match[0] = (seq1.charAt(x) == seq2.charAt(y)) ? '|' : '*';
					algn.addFirst(column, match);
					x--; y--;
				}
			}
			
			// gather up the thread-mess we created and join them
			Thread[] threads = new Thread[tg.activeCount()];
			tg.enumerate(threads);
			for(int i=0; i<tg.activeCount(); i++)
				try {
					threads[i].join();
				} catch (InterruptedException e) {
					// TODO should be fine when access to algnmts is synchronized
					e.printStackTrace();
				}	// UH.. ugly... 6 lines and now bracketed for-loop!
			
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
