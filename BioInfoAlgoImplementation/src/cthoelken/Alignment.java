package cthoelken;

import java.security.InvalidParameterException;

/**
 * Object to keep all needed sequences, symbols for matches and mismatches 
 * between them and print out the alignment.
 * 
 * @author Clemens Thoelken
 *
 */
public class Alignment {
	
	public String[] names;
	public String[] sequences;
	public String[] matches;	// match/mismatch symbols
	public double score;
	
	public void setSequence(int index, String seq) {
		if(index >= sequences.length || index <= 0)
			throw new IllegalArgumentException("Index out of bounds in Alignment");
		sequences[index] = seq;
	}
	
	public void setName(int index, String name) {
		if(index >= sequences.length || index <= 0)
			throw new IllegalArgumentException("Index out of bounds in Alignment");
		names[index] = name;
	}
	
	public String getName(int index) {
		if(index >= sequences.length || index <= 0)
			throw new IllegalArgumentException("Index out of bounds in Alignment");
		return names[index];
	}
	
	public String getSequence(int index) {
		if(index >= sequences.length || index <= 0)
			throw new IllegalArgumentException("Index out of bounds in Alignment");
		return sequences[index];
	}
	
	public double getScore() {
		return score;
	}
	
	public void setScore(Double score) {
		this.score = score;
	}
	
	/**
	 * Constructor
	 * 
	 * @param seqCount Number of sequences that should be handled
	 */
	Alignment(int seqCount) {
		if(seqCount < 2)
			throw new InvalidParameterException("Too few sequences for " +
					"an alignment!");
		names = new String[seqCount];
		sequences = new String[seqCount];
		matches = new String[seqCount-1];
		score = 0;
		sequences[0] = ""; names[0] = "";
		for(int i=1; i<seqCount; i++) {
			names[i] = ""; sequences[i] = ""; matches[i-1] = "";
		}
	}
	
	/**
	 * Copy-constructor
	 * @param algn Alignment that is to be copied
	 */
	Alignment(Alignment algn) {
		sequences = new String[algn.sequences.length];
		matches = new String[algn.sequences.length-1];
		sequences[0] = new String(algn.sequences[0]);
		for(int i=1; i<algn.sequences.length; i++) {
			sequences[i] = new String(algn.sequences[i]);
			matches[i-1] = new String(algn.matches[i-1]);
		}
	}
	
	/**
	 * Attaches a column to the front of the returned alignment
	 * Keeps the current alignment unchanged
	 * @param seqColumn Column of symbols which were aligned
	 * @param matchColumn Column of match/mismatch symbols for the alignment
	 * @return Returns the updated alignment
	 */
	public Alignment addFirst(char[] seqColumn, char[] matchColumn) {
		if(seqColumn.length != sequences.length || 
				matchColumn.length != matches.length) 
			throw new InvalidParameterException("Array length not correct!");
//		int counter = 0;
//		for(int j=0; j<seqColumn.length; j++)
//			if(seqColumn[j]=='#') counter++;
//		if(counter == seqColumn.length) return this;
		Alignment temp = new Alignment(this);
		temp.sequences[0] = seqColumn[0] + sequences[0];
		for(int i=1; i<sequences.length; i++) {
			temp.sequences[i] = seqColumn[i] + sequences[i];
			temp.matches[i-1] = matchColumn[i-1] + matches[i-1];
		}
		return temp;
	}
	
	/**
	 * Attaches a column to the end of the sequences.
	 * @param seqColumn Column of symbols which were aligned
	 * @param matchColumn Column of match/mismatch symbols for the alignment
	 * @return Returns the updated alignment
	 */
	public Alignment addLast(char[] seqColumn, char[] matchColumn) {
		if(seqColumn.length != sequences.length || 
				matchColumn.length != matches.length) 
			throw new InvalidParameterException("Array length not correct!");
		sequences[0] = seqColumn[0] + sequences[0];
		for(int i=1; i<sequences.length; i++) {
			sequences[i] = sequences[i] + seqColumn[i];
			matches[i-1] = matches[i-1] + matchColumn[i-1];
		}
		return this;
	}
	
	/**
	 * Pretty prints the Alignments including the sequences and match/mismatch between them.
	 */
	public String toString() {
		String retVal = "";
		retVal += "Seq1: " + sequences[0] + "\n";
		for(int i=1; i<sequences.length; i++)
			retVal += "      " + matches[i-1] + "\n" + "Seq" + (i+1) +
					": "+ sequences[i] + "\n";
		return retVal;
	}

	public int size() {
		return sequences.length;
	}
	

}
