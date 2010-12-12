/**
 * 
 */
package cthoelken;

import java.security.InvalidParameterException;
import java.util.Arrays;

/**
 * @author n2
 *
 */
public class Alignment {
	
	public String[] sequences;
	
	Alignment(int seqCount) {
		sequences = new String[seqCount];
	}
	
	Alignment(Alignment algn) {
		for(int i=0; i<algn.sequences.length; i++)
			this.sequences[i] = new String(algn.sequences[i]);
	}
	
	public Alignment addFirst(char[] column) {
		if(column.length != sequences.length) 
			throw new InvalidParameterException("Array length not correct!");
		for(int i=0; i<sequences.length; i++) {
			sequences[i] = column[i] + sequences[i];
		}
		return this;
	}
	
	public Alignment addLast(char[] column) {
		if(column.length != sequences.length) 
			throw new InvalidParameterException("Array length not correct!");
		for(int i=0; i<sequences.length; i++) {
			sequences[i] = sequences[i] + column[i];
		}
		return this;
	}
	
	public String toString() {
		String retVal = "";
		for(int i=0; i<sequences.length; i++)
			retVal += sequences[i] + "\n";
		return retVal;
	}
	

}
