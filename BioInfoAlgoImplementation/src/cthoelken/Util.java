package cthoelken;

import gui.StringList;

import java.util.LinkedList;

/** Utility class for static methods used by multiple algorithms for intermediate
 * parsing and computations
 * @author Clemens Thoelken
 *
 */
/**
 * @author n2
 *
 */
public class Util {
	
	/** Returns the maximum value of the double parameters
	 * @param value Input double array
	 * @return maximal input parameter
	 */
	public static double maxValue(double... value) {
		double max = Double.NEGATIVE_INFINITY;
		for(int i = 0; i < value.length; i++)
			if(max < value[i]) max = value[i];
		return max;
	}
	
	/** Returns the index of the maximal double in the parameters
	 * @param value Input double array
	 * @return Index of the maximal input parameter
	 */
	public static int maxIndex(double... value) {
		double max = 0;
		int index = 0;
		for(int i = 0; i < value.length; i++)
			if(max < value[i]) {
				max = value[i];
				index = i;
			}
		return index;
	}
	
	/** Checks a character whether it is a valid amino acid including gaps or not.
	 * @param Letter code for an amino acid or gap '_'.
	 * @return Returns TRUE if valid, FALSE otherwise.
	 */
    private static boolean isAmino(char a) {
    	// check for upper and lowercase characters
    	switch ((String.valueOf(a)).toUpperCase().charAt(0)) {
	    	case 'A': return true;
	    	case 'R': return true;
	    	case 'N': return true;
	    	case 'D': return true;
	    	case 'C': return true;
	    	case 'Q': return true;
	    	case 'E': return true;
	    	case 'G': return true;
	    	case 'H': return true;
	    	case 'I': return true;
	    	case 'L': return true;
	    	case 'K': return true;
	    	case 'M': return true;
	    	case 'F': return true;
	    	case 'P': return true;
	    	case 'S': return true;
	    	case 'T': return true;
	    	case 'W': return true;
	    	case 'Y': return true;
	    	case 'V': return true;
	    	case 'B': return true;
	    	case 'Z': return true;
	    	case 'X': return true;
	    	case '_': return true;
	    	default: return false;		// if not in the alphabet
    	}
    }
    
    
    /** Checks a sequence string for its validity.
     * @param s Sequence String.
     * @return TRUE if the string only consists of valid amino acids.
     */
    public static boolean isValidSequence(String s) {
    	for(int i = 0; i < s.length(); i++)
    		if (!isAmino(s.charAt(i))) return false;
    	return true;
    }

    
    /** Parses an Alignment from a FASTA String with newlines
	 * @param fasta Input string in FASTA format
	 * @return Alignment with all sequences
	 */
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
			if(!isValidSequence(sequences.get(i)))
					throw new IllegalArgumentException("Sequence \'" 
							+ names.get(i) + "\' is not conform to FASTA!");
			algn.setSeq(i, sequences.get(i));
		}

		return algn;
	}
}
