package cthoelken;

import java.security.InvalidParameterException;
import java.util.Arrays;

/**
 * @author Clemens Thoelken
 *
 */
public class CostMatrix {
	
	private int[][] m;

	/**
	 * Constructor
	 * @param xLength Length of sequence 1
	 * @param yLength Length of sequence 2
	 */
	CostMatrix(int xLength, int yLength) {
		if(xLength < 1 || yLength < 1) 
			throw new InvalidParameterException("Sequence length not feasable!");
		m = new int[xLength][yLength];
		m[0][0] = 0;
	}
	
	/**
	 * Get the costs at the point x, y
	 * @param x Position in sequence 1
	 * @param y Position in sequence 2
	 * @return Costs
	 */
	public int get(int x, int y) {
		if(x < 0 || y < 0 || x >= m.length || y >= m[0].length) 
			return Integer.MIN_VALUE+1000;
		return m[x][y];
	}
	
	public int score() {
		return m[m.length-1][m[0].length-1];
	}
	
	/**
	 * Updates the costs at the point x, y
	 * @param x Position in sequence 1
	 * @param y Position in sequence 2
	 * @param value Value that should be inserted
	 */
	public void set(int x, int y, int value) {
		if(x < 0 || y < 0 || x >= m.length || y >= m[0].length)
			throw new InvalidParameterException("Index not feasable!");
		m[x][y] = value;
	}
	
	/**
	 * Pretty prints the cost matrix
	 */
	public String toString() {
		String retVal = "";
		for(int i=0; i<m.length; i++)
			retVal += Arrays.toString(m[i]) + "\n";
		return retVal;
	}
	
}
