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
	
	private Alignment calculate(Alignment algn) {
		
		if(algn.getSeq(0).length() == 0) algn.setSeq(0, "_");
		if(algn.getSeq(1).length() == 0) algn.setSeq(1, "_");
		if(algn.getSeq(2).length() == 0) algn.setSeq(2, "_");
		
		if((algn.getSeq(0) == "_" && algn.getSeq(1) == "_") ||
				(algn.getSeq(0) == "_" && algn.getSeq(2) == "_") ||
				(algn.getSeq(1) == "_" && algn.getSeq(2) == "_")) {
//			System.out.println("\n"+algn.toString());
			algn.setScore(0.0);
			return algn;
		}
		
		if(algn.getSeq(0).charAt(0) == '_' && algn.getSeq(1).charAt(0) == '_' && algn.getSeq(2).charAt(0) == '_') return algn;
		
		LinkedList<Alignment> psblAlgn = new LinkedList<Alignment>();
		double max = Double.NEGATIVE_INFINITY;
		Alignment temp = new Alignment(algn);
		temp.delFirst();
		psblAlgn.add(new Alignment(calculate(temp)));					// 000
		if(!temp.getSeq(0).isEmpty()) {
			temp.setSeq(0, "_" + temp.getSeq(0));
			psblAlgn.add(new Alignment(calculate(temp)));				// 100
		}
		if(!temp.getSeq(0).isEmpty() && !temp.getSeq(1).isEmpty()) {
			temp.setSeq(1, "_" + temp.getSeq(1));
			psblAlgn.add(new Alignment(calculate(temp)));				// 110
		}
		if(!temp.getSeq(1).isEmpty() && !temp.getSeq(0).isEmpty()) {
			temp.setSeq(0, temp.getSeq(0).substring(1));
			psblAlgn.add(new Alignment(calculate(temp)));				// 010
		}
		if(!temp.getSeq(1).isEmpty() && !temp.getSeq(2).isEmpty()) {
			temp.setSeq(2, "_" + temp.getSeq(2));
			psblAlgn.add(new Alignment(calculate(temp)));				// 011
		}
		if(!temp.getSeq(2).isEmpty() && !temp.getSeq(1).isEmpty()) {
			temp.setSeq(1, temp.getSeq(1).substring(1));
			psblAlgn.add(new Alignment(calculate(temp)));				// 001
		}
		if(!temp.getSeq(0).isEmpty() && !temp.getSeq(2).isEmpty()) {
			temp.setSeq(0, "_" + temp.getSeq(0));
			psblAlgn.add(new Alignment(calculate(temp)));				// 101
		}
		
		for(int i = 0; i < psblAlgn.size(); i++)
			if(max < psblAlgn.get(i).getScore()) max = psblAlgn.get(i).getScore();
		System.out.println("\n Max: "+ max);
		
//		for(int i = 0; i < 7; i++) {
//			switch (i) {
//			case 1: temp.setSeq(0, "_" + temp.getSeq(0)); break; 		// _00
//			case 2: temp.setSeq(1, "_" + temp.getSeq(1)); break; 		// __0
//			case 3: temp.setSeq(0, temp.getSeq(0).substring(1)); break; // 0_0
//			case 4: temp.setSeq(2, "_" + temp.getSeq(2)); break; 		// 0__
//			case 5: temp.setSeq(1, temp.getSeq(1).substring(1)); break; // 00_
//			case 6: temp.setSeq(0, "_" + temp.getSeq(0)); break; 		// _0_
//			}
//			psblAlgn.add(new Alignment(calculate(temp)));
//			max = NeedlemanWunsch.maxValue(psblAlgn.get(i).getScore(), max);
//		}
		
		double score = 0.0;
		score += omega.getScore(algn.getSeq(0).charAt(0), algn.getSeq(1).charAt(0));
		score += omega.getScore(algn.getSeq(0).charAt(0), algn.getSeq(2).charAt(0));
		score += omega.getScore(algn.getSeq(1).charAt(0), algn.getSeq(2).charAt(0));

		for(int i = 0; i < psblAlgn.size(); i++)
			if(psblAlgn.get(i).getScore() == max) {
				if(algn.getSeq(0).charAt(0) == '*' && algn.getSeq(1).charAt(0) == '*' && algn.getSeq(2).charAt(0) == '*')
					return new Alignment(psblAlgn.get(i));
				psblAlgn.get(i).setScore(max + score);
				char[] seqColumn = {algn.getSeq(0).charAt(0), algn.getSeq(1).charAt(0), algn.getSeq(2).charAt(0)};
				char[] matchColumn = {' ', ' '};
				return new Alignment(psblAlgn.get(i).addFirst(seqColumn, matchColumn));
			}
		return temp;
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
		sequences = new Alignment(3);
		sequences.setSeq(0, "*"+(String) params.elementAt(0).data);
		sequences.setSeq(1, "*"+(String) params.elementAt(1).data);
		sequences.setSeq(2, "*"+(String) params.elementAt(2).data);
		usePAM = (Boolean) params.elementAt(3).data;
		gapCosts = (Double) params.elementAt(4).data;
		
		omega = new SubstitutionMatrix(usePAM, gapCosts);
		
		retVal += "\n mhh.. I assume my default values are fine ! ;) \n";
		
		  // ##########  RUN THE PROGRAM  ###########
		
		alignment = calculate(sequences);
		alignment.makeConnections();

		 // (JUST TO EXEMPLIFY SOME OUTPUT I REPORT THE INPUT PARAMETERS ...)
		for (int i = 0; i < params.size(); i++) {
			retVal	+= "\n Input : "
					+ params.elementAt(i).name
					+ " = "
					+ params.elementAt(i).defVal.toString();
		}
		
		retVal += "\n\n"+ alignment.toString();
				
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
