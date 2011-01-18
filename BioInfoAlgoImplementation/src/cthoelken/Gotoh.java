package cthoelken;

import gui.AlgorithmParameter;
import gui.BioinfAlgorithm;

import java.util.Random;
import java.util.Vector;

/**
 * Implementation of the Gotoh Algorithm for pairwise alignments with affine gap
 * costs.
 * 
 * @author Clemens Thoelken
 * 
 */
public class Gotoh extends NeedlemanWunsch {

	protected double gapCostsExt;
	protected CostMatrix H;
	protected CostMatrix V;

	/**
	 * Constructor which generates an empty vector of parameters of the needed
	 * types.
	 */
	public Gotoh() {
		super();

		// create additional parameters for the algorithm to work.
		super.parameters.add(new AlgorithmParameter("Extended gap costs",
				"An integer for the constant gap costs used for scoring.",
				Double.class, new Double(-1.0)));

	}

	@Override
	public String getName() {
		return new String("Gotoh");
	}

	@Override
	public String getDescription() {
		return new String("Implemenation of the Gotoh affine gap cost "
				+ "approach to pairwise alignment with the help of dynamic "
				+ "programming.");
	}

	private double calculate() {

		seq1 = "#" + seq1;
		seq2 = "#" + seq2; // increase sequence length, disregarded afterwards

		// fill the cost matrix row-wise
		for (int i = 0; i < seq1.length(); i++) {
			H.set(i, 0, Double.NEGATIVE_INFINITY); // init H fist row
			System.out.println("\n"+H.toString());
			for (int j = 0; j < seq2.length(); j++) {
				if (i == 0)
					V.set(0, j, Double.NEGATIVE_INFINITY); // init V first column
				System.out.println("\n"+V.toString());
				if (!(i == 0 && j == 0)) {
					// calculate auxiliary matrices
					if (j != 0)
						H.set(i, j, maxValue(M.get(i-1, j)
								+ gapCosts + gapCostsExt,
								H.get(i-1, j) + gapCostsExt));
					if (i != 0)
						V.set(i, j, maxValue(M.get(i, j-1)
								+ gapCosts + gapCostsExt,
								V.get(i, j-1) + gapCostsExt));

					M.set(i, j, maxValue(M.get(i-1, j-1)
							+ omega.getScore(seq1.charAt(i), seq2.charAt(j)),
							H.get(i, j), V.get(i, j)));
				}
			}
		}
		// seq1 = seq1.substring(1); seq2 = seq2.substring(1);
		return M.score();
	}

	private void backtrack(int x, int y, Alignment algn) {

		System.out.println("backtrack: x=" + x + " y=" + y);

		if (x < 0 || y < 0)
			return; // we are out of bounds, shouldn't happen

		char[] column = new char[2];
		char[] match = new char[1];

		if (x == 0 && y == 0) { // we are done, great!
		// column[0] = seq1.charAt(x); column[1] = seq2.charAt(y);
		// match[0] = (seq1.charAt(x) == seq2.charAt(y)) ? '|' : '*';
		// algnmts.add(algn.addFirst(column, match));
			System.out.println("All is good!");
			algnmts.add(algn);
			return;
		}

		int[] psbl = { 0, 0, 0 };

		System.out.println("\n M(x,y) = " + M.get(x, y) + " - M(x,y-1) = "
				+ M.get(x, y - 1) + " - SEQ2: " + seq2.charAt(y)); // TODO

		// insertion
		if (M.get(x, y) == H.get(x, y)) psbl[0] = 1;

		// deletion, start new thread from there
		if (M.get(x, y) == V.get(x, y)) psbl[1] = 1;

		// match/mismatch, lets carry on
		if (M.get(x, y) == M.get(x - 1, y - 1)
				+ omega.getScore(seq1.charAt(x), seq2.charAt(y)))
			psbl[2] = 1;

		System.out.println("\n " + psbl[0] + "" + psbl[1] + "" + psbl[2]); // TODO

		int num_psbl = psbl[0] + psbl[1] + psbl[2];
		if (randomBackTrace && num_psbl > 1) {
			int random;
			if (num_psbl == 3) {
				psbl[0] = 0;
				psbl[1] = 0;
				psbl[2] = 0;
				psbl[new Random().nextInt(3)] = 1;
			}
			if (num_psbl == 2) {
				random = new Random().nextInt(2);
				if (psbl[0] == 1) {
					if (psbl[1] == 1) {
						psbl[0] = random;
						psbl[1] = 1 - random;
					} else {
						psbl[0] = random;
						psbl[2] = 1 - random;
					}
				} else {
					psbl[1] = random;
					psbl[2] = 1 - random;
				}
			}
			System.out.println("\n " + psbl[0] + "" + psbl[1] + "" + psbl[2]); // TODO
		}

		match[0] = ' ';

		if (psbl[0] == 1) {
			System.out.println("\n I go left");
			column[0] = seq1.charAt(x);
			column[1] = '_';
			backtrack(x - 1, y, algn.addFirst(column, match));
		}
		if (psbl[1] == 1) {
			System.out.println("\n I go right");
			column[0] = '_';
			column[1] = seq2.charAt(y);
			backtrack(x, y - 1, algn.addFirst(column, match));
		}
		if (psbl[2] == 1) {
			System.out.println("\n I go diagonally");
			column[0] = seq1.charAt(x);
			column[1] = seq2.charAt(y);
			match[0] = (seq1.charAt(x) == seq2.charAt(y)) ? '|' : '*';
			backtrack(x - 1, y - 1, algn.addFirst(column, match));
		}
	}

	/**
	 * Main method of the algorithm.
	 * 
	 * @param params
	 *            The filled out parameters are entered externally.
	 * 
	 * @return Output string containing the used parameters, the results and
	 *         errors.
	 */
	@Override
	public String run(Vector<AlgorithmParameter> params) {

		String retVal = new String("");

		// ########## PARSE INPUT PARAMETERS FOR ERRORS ###########

		// TODO Parse input for errors!!!
		seq1 = (String) params.elementAt(0).data;
		seq2 = (String) params.elementAt(1).data;
		usePAM = (Boolean) params.elementAt(2).data;
		randomBackTrace = (Boolean) params.elementAt(3).data;
		gapCosts = (Double) params.elementAt(4).data;
		gapCostsExt = (Double) params.elementAt(5).data;

		omega = new SubstitutionMatrix(usePAM, gapCosts);
		M = new CostMatrix(seq1.length() + 1, seq2.length() + 1);
		H = new CostMatrix(seq1.length() + 1, seq2.length() + 1);
		V = new CostMatrix(seq1.length() + 1, seq2.length() + 1);

		retVal += "\n mhh.. I assume my default values are fine ! ;) \n";

		// ########## RUN THE PROGRAM ###########

		// (JUST TO EXEMPLIFY SOME OUTPUT I REPORT THE INPUT PARAMETERS ...)
		for (int i = 0; i < params.size(); i++) {
			retVal += "\n Input : " + params.elementAt(i).name + " = "
					+ params.elementAt(i).defVal.toString();
		}

		retVal += "\n Maximal Score: " + calculate();

		// print final cost matrix to console VERBOSE!!!
		// TODO add to return string!!!
		System.out.println(M.toString());
		System.out.println(H.toString());
		System.out.println(V.toString());

		backtrack(seq1.length() - 1, seq2.length() - 1, new Alignment(2));

		// ### print alignments to return string ###
		for (int k = 0; k < algnmts.size(); k++)
			retVal += "\n########### Alignment " + k + ":\n"
					+ algnmts.get(k).toString();

		// returning the algorithm results
		System.out.println("######################\n"
				+ omega.getScore('b', 'b'));
		return retVal;
	}

	/**
	 * Creates an instance of this class and calls the run method using the
	 * default parameters.
	 * 
	 * @param args
	 *            program parameters (completely ignored)
	 */
	public static void main(String[] args) {
		// create an instance of this class
		Gotoh myInstance = new Gotoh();
		// run the example the instance with the default parameters
		BioinfAlgorithm.runAlgorithmDefaults(myInstance);
	}

}
