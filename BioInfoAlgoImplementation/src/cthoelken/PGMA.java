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
	
	public static Alignment parseFasta(StringList fasta) {
		String[] lines = fasta.toString().split("\n");
		LinkedList<String> names = new LinkedList<String>();
		LinkedList<String> sequences = new LinkedList<String>();
		String temp = "";

		for(int i = 0; i < lines.length; i++) {
			if(lines[i].length() > 0) {
				if(lines[i].charAt(0) == ';') {	
				} else if(lines[i].charAt(0) == '>' && names.size() == sequences.size()) {
					names.add(lines[i].substring(1));
				} else if(names.size()-1 == sequences.size()) {
					temp += lines[i];
					if(lines[i].charAt(lines[i].length()-1) == '*') {
						sequences.add(temp.substring(0, temp.length()-1));
						temp = "";
					} else if(i+1 == lines.length) {
						sequences.add(temp);
						temp = "";
					} else if(lines[i+1].charAt(0) == '>' || lines[i+1].length() == 0) {
						sequences.add(temp);
						temp = "";
					}
				} else throw new IllegalArgumentException("Line "+i
						+" is not a valid formated according to FASTA!");
			}
		}
		
		Alignment algn = new Alignment(names.size());
		for(int i = 0; i < names.size(); i++) {
			algn.setName(i, names.get(i));
			for(int j = 0; j < sequences.get(i).length(); j++)
				if(new SubstitutionMatrix().getScore(sequences.get(i).charAt(j), sequences.get(i).charAt(j)) < 0)
					throw new IllegalArgumentException("\""
							+ sequences.get(i).charAt(j)
							+ "\" is not a valid amino acid in sequence "+i
							+ " at position "+(j+1)+"!");
			algn.setSeq(i, sequences.get(i));
		}

		return algn;
	}
	
	public Cluster calculate(Alignment sequences, boolean usePAM, boolean weighted, double gapCosts) {
		this.sequences = sequences; this.usePAM = usePAM; 
		this.weighted = weighted; this.gapCosts = gapCosts;
		return calculate();
	}
	
	/**
	 * @return
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
			sequences = PGMA.parseFasta((StringList) params.elementAt(0).data);
		} catch(IllegalArgumentException e) {
			return e.toString();
		}
		weighted = (Boolean) params.elementAt(1).data;
		usePAM = (Boolean) params.elementAt(2).data;
		gapCosts = (Integer) params.elementAt(3).data;
		
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
