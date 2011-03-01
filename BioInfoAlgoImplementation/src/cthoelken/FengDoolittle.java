package cthoelken;

import gui.AlgorithmParameter;
import gui.BioinfAlgorithm;
import gui.StringList;
import java.util.Vector;

 /**
  * Implementation of the Feng-Doolittle Algorithm for progressive multiple alignments
  * with linear gap costs.
  * 
  * @author Clemens Thoelken
  *
  */
public class FengDoolittle extends BioinfAlgorithm {

	protected Alignment sequences;
	protected boolean usePAM;
	protected double gapCosts;
	protected Alignment alignment;
	protected boolean useUPGMA;
	protected Cluster tree;

	 /**
	  * Constructor which generates an empty vector of parameters of the needed 
	  * types.
	  */
	public FengDoolittle() { 
		
		// create all needed parameters for the algorithm to work.
		super.parameters.add(new AlgorithmParameter(
				"Sequences"
				, "Enter Sequences in FASTA format." 
				, StringList.class 
				, new StringList("\n\n;Kommentar 1 \n>" +
						"Sequence 1 [Die Schwarmm�ücke] \"Plutonium Maximum\"" +
						"\n;Kommentar2\nABCDEF\n>" +
						"Sequence 2\nabdefdab\nG*\n>Sequnce 3\n" +
						"dfgabab\n>Sequence 4\nccccgabdegfef")));
		super.parameters.add(new AlgorithmParameter(
				"use PAM (BLOSUM otherwise)"
				, "Choose YES to use PAM or NO to use BLOSUM for scoring." 
				, Boolean.class 
				, new Boolean(true)));
		super.parameters.add(new AlgorithmParameter(
				"use WPGMA (otherwise UPGMA) pairing"
				, "Choose YES to use unweighted or NO to use weighted for pairing." 
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
		return new String("Feng-Doolittle");
	}
	

	@Override
	public String getDescription() {
		return new String("Implemenation of the Feng-Doolittle approach " +
				"to multiple alignment with the help of dynamic programming.");
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
		
		try{
			sequences = Util.parseFasta((StringList) params.elementAt(0).data);
		} catch(IllegalArgumentException e) {
			return e.toString();
		} catch(Exception e) { return "The FASTA data is not valid!"; }
		usePAM = (Boolean) params.elementAt(1).data;
		useUPGMA = (Boolean) params.elementAt(2).data;
		try{
			if(params.elementAt(3).data.getClass() == Double.class)
				gapCosts = (Double) params.elementAt(3).data;
			else return "Gap costs are not a valid decimal value!";
		} catch(Exception e) {return "Gap costs are not a valid decimal value!";}
		
		  // ##########  RUN THE PROGRAM  ###########

		tree = new PGMA().calculate(sequences, usePAM, !useUPGMA, gapCosts);
		
		tree.align(usePAM, gapCosts);
		
		retVal += "\n" + tree.toString();
		retVal += "\n" + tree.generateAlignment(sequences).toString();
		retVal += "\n Sum-of-Pairs Score: "  + SumOfPairs.score(tree.generateAlignment(sequences), usePAM, gapCosts);
		
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
		FengDoolittle myInstance = new FengDoolittle();
		  // run the example the instance with the default parameters
		BioinfAlgorithm.runAlgorithmDefaults( myInstance );
	}

}
