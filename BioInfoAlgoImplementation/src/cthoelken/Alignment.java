package cthoelken;

import java.security.InvalidParameterException;

/**
 * Object to keep 2 to n sequences, symbols for matches and mismatches 
 * between them and print out the alignment.
 * 
 * @author Clemens Thoelken
 *
 */
public class Alignment {
	
	public String[] names;		// sequence names
	public String[] sequences;	// sequence strings
	public String[] matches;	// match/mismatch symbols
	public double score;
	
	/** Sets a sequence at the given index to the defined string
	 * @param index	The index of the Sequence
	 * @param seq The sequence-string
	 */
	public void setSeq(int index, String seq) {
		if(index >= sequences.length || index < 0)
			throw new IllegalArgumentException("Index out of bounds in Alignment");
		sequences[index] = seq;
	}
	
	
	/** Setter for the name of the sequence by index
	 * @param index Index of the sequence
	 * @param name The name string to be assigned to the sequence
	 */
	public void setName(int index, String name) {
		if(index >= sequences.length || index < 0)
			throw new IllegalArgumentException("Index out of bounds in Alignment");
		names[index] = name;
	}
	
	/** Getter for the name of the sequence by index
	 * @param index Index of the sequence
	 * @return The sequence's name
	 */
	public String getName(int index) {
		if(index >= sequences.length || index < 0)
			throw new IllegalArgumentException("Index out of bounds in Alignment");
		return names[index];
	}
	
	/** Getter for the sequences by index
	 * @param index Index of the sequence
	 * @return The sequence string
	 */
	public String getSeq(int index) {
		if(index >= sequences.length || index < 0)
			throw new IllegalArgumentException("Index out of bounds in Alignment");
		return sequences[index];
	}
	
	/** Returns the score of the alignment. This value is not necessarily computed!
	 * @return The Score of the overall alignment
	 */
	public double getScore() {
		return score;
	}
	
	/** Setter for the score variable
	 * @param score Value of the score to be assigned to the alignment
	 */
	public void setScore(Double score) {
		this.score = score;
	}
	
//	
//	public void makeConnections() {
//		matches = new String[matches.length];
//		int index = 0;
//		int finishedSeq = 0;
//		String temp = " ";
//		while(finishedSeq < matches.length) {
//			for(int i = 0; i < matches.length; i++) {
//				if(index == 0) matches[i] = "";
//				if(index < sequences[i].length() && index < sequences[i+1].length()) {
//					temp = ".";
//					if(sequences[i].charAt(index) == sequences[i+1].charAt(index)) temp = "|";
//					if(sequences[i].charAt(index) == '_' || sequences[i+1].charAt(index) == '_') temp = " ";
//					matches[i] = matches[i] + "" + temp;
//				}
//				if(sequences[i].length() == index+1) finishedSeq++;
//			}
//			index++;
//		}
//	}
	
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
		score = 0.0;
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
			sequences[i] = new String(algn.getSeq(i));
			matches[i-1] = new String(algn.matches[i-1]);
		}
		score = algn.getScore();
		
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
	
	public void delFirst() {
		for(int i = 0; i < sequences.length; i++)
			if(sequences[i].length() > 0) 
				sequences[i] = sequences[i].substring(1);
	}
	
	/* (non-Javadoc)
	 * Pretty prints the Alignment in a 80 chars limited output with symols for
	 * matches and mismatches
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		for(int i = 1; i < sequences.length; i++)
			if(sequences[i].length() != sequences[0].length())
				return "Sequences in the alignment are not propperly aligned!\n";
		String retVal = "";
		String prefix = "0";
		int limit = 0; int matchlimit = 0;
		for(int times = 0; times <= sequences[0].length()/76; times++) {
			limit = (sequences[0].length() > times*76+75) ? (times*76+75) : (sequences[0].length()%76 + times*76);
			retVal += "01: " + sequences[0].substring(times*76,limit) + "\n";
			for(int i=1; i<sequences.length; i++) {
				prefix = (i > 9) ? "" : "0";
				matchlimit = (matches[i-1].length() > times*76+75) ? (times*76+75) : (matches[i-1].length()%76 + times*76);
				retVal += "    " + matches[i-1].substring(times*76,matchlimit) + "\n" + prefix + (i+1) +
					": "+ sequences[i].substring(times*76,limit) + "\n";
			}
			retVal += "\n";
		}
		return retVal;
	}

	
	/** Determines the number of sequences in this alignment
	 * @return Returns the number of Sequences
	 */
	public int size() {
		return sequences.length;
	}
	

}
