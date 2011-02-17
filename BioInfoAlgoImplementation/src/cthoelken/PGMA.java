package cthoelken;

import gui.AlgorithmParameter;
import gui.BioinfAlgorithm;
import gui.StringList;
import java.util.LinkedList;
import java.util.Vector;

 /**
  * Implementation of the Weighted and Unweighted Pair Group Method using the 
  * Arithmetic mean for multiple alignments.
  * 
  * @author Clemens Thoelken
  *
  */
public class PGMA extends BioinfAlgorithm {

	private Alignment sequences;
	private boolean usePAM;
	private boolean weighted;
	private double gapCosts;


	/**
	 * Constructor which generates an empty vector of parameters of the needed 
	 * types.
	 */
	public PGMA() { 
		// create additional parameters for the algorithm to work.
		super.parameters.add(new AlgorithmParameter(
				"Sequences"
				, "Enter Sequences in FASTA format." 
				, StringList.class 
				, new StringList("\n\n;Kommentar 1 \n>" +
						"Sequence 1 [Die Schwarmm�ücke] \"Plutonium Maximum\"" +
						"\n;Kommentar2\nABCDEF\n>" +
						"Sequence 2\nABGEFGGGGGGGGGGGGGG\nG*\n>Sequnce 3\n" +
						"ABGDEF\n>Sequence 4\nABDEFGGGGGGG")));
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
				, "A decimal value for the constant gap costs used for scoring."
				, Double.class 
				, new Double(-1.0)));
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
	
	
	/** Calculates a tree form a given alignment
	 * @param sequences Input alignment with all sequences
	 * @param usePAM Use PAM substitution
	 * @param weighted Use WPGMA
	 * @param gapCosts Gap costs for the calculation
	 * @return Returns the tree as a cluster
	 */
	public Cluster calculate(Alignment sequences, boolean usePAM, boolean weighted, double gapCosts) {
		this.sequences = sequences; this.usePAM = usePAM; 
		this.weighted = weighted; this.gapCosts = gapCosts;
		return calculate();
	}
	
	/** Initialtizes all sequences into Clusters and slusters all untill only one Cluster is left
	 * @return The resulting Cluster
	 */
	private Cluster calculate() {
		CostMatrix D = new CostMatrix(sequences.size(), sequences.size());
		
		LinkedList<Cluster> nodes = new LinkedList<Cluster>();
		
		for(int i = 0; i < sequences.size(); i++)
			nodes.add(new Cluster(i, sequences.getSeq(i)));
		
		while(nodes.size() > 1) {
			double min = Double.POSITIVE_INFINITY;
			int iMin = 0; int jMin = 1;

			for(int i = 0; i < nodes.size(); i++) {
				for(int j = i+1; j < nodes.size(); j++) {
					D.set(i, j, nodes.get(i).distance(nodes.get(j), weighted, usePAM, gapCosts));
					if(min > D.get(i, j)) {
						min = D.get(i, j);
						iMin = i; jMin = j;
					}
				}
			}

			Cluster temp = new Cluster(nodes.get(iMin), nodes.get(jMin), D.get(iMin, jMin));
			nodes.remove(iMin); nodes.remove(jMin-1);
			nodes.add(temp);
		}
		
		return nodes.get(0);
				
	}
	
	/**
	 * Main method of the algorithm.
	 * 
	 * @param params The filled out parameters are entered externally.
	 * 
	 * @return Output string containing the used parameters, the results and errors.
	 */
	/* (non-Javadoc)
	 * @see gui.BioinfAlgorithm#run(java.util.Vector)
	 */
	@Override
	public String run(Vector<AlgorithmParameter> params) {
		
		String retVal = new String("");
		
		  // ##########  PARSE INPUT PARAMETERS FOR ERRORS  ###########
		
		try{
			sequences = Util.parseFasta((StringList) params.elementAt(0).data);
		} catch(IllegalArgumentException e) {
			return e.toString();
		}
		weighted = (Boolean) params.elementAt(1).data;
		usePAM = (Boolean) params.elementAt(2).data;
		try{
			if(params.elementAt(3).data.getClass() == Double.class)
				gapCosts = (Double) params.elementAt(3).data;
			else return "Gap costs are not a valid decimal value!";
		} catch(Exception e) {return "Gap costs are not a valid decimal value!";}
		
		  // ##########  RUN THE PROGRAM  ###########

		retVal += "\n" + calculate().toString();

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
		PGMA myInstance = new PGMA();
		  // run the example the instance with the default parameters
		BioinfAlgorithm.runAlgorithmDefaults( myInstance );
	}

}
