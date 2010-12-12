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
public class CostMatrix {
	
	private int[][] m;
	
	@SuppressWarnings("unused")
	private SubstitutionMatrix match;

	CostMatrix(int xLength, int yLength) {
		this(xLength, yLength, new SubstitutionMatrix(true));
	}
	
	CostMatrix(int xLength, int yLength, SubstitutionMatrix sub) {
		if(xLength < 1 || yLength < 1) 
			throw new InvalidParameterException("Sequence length not feasable!");
		m = new int[xLength][yLength];
		m[0][0] = 0;
		match = sub;
	}
	
	public int get(int x, int y) {
		if(x < 0 || y < 0 || x >= m.length || y >= m[0].length) 
			return -Integer.MIN_VALUE+1000;
		return m[x][y];
	}
	
	public void set(int x, int y, int value) {
		if(x < 0 || y < 0 || x >= m.length || y >= m[0].length)
			throw new InvalidParameterException("Index not feasable!");
		m[x][y] = value;
	}
	
	public String toString() {
		String retVal = "";
		for(int i=0; i<m.length; i++)
			retVal += Arrays.toString(m[i]) + "\n";
		return retVal;
	}
	
}
