package cthoelken;

import java.security.InvalidParameterException;
import java.util.Arrays;

/** Auxilliary class to store the computed values for two or three sequences
 * @author Clemens Thoelken
 *
 */
public class CostMatrix {
	
	private double[][][] m;

	/** Constructor for three sequences
	 * @param xLength Length of sequence 1
	 * @param yLength Length of sequence 2
	 * @param zLength Length of sequence 3
	 */
	CostMatrix(int xLength, int yLength, int zLength) {
		if(xLength < 1 || yLength < 1 || zLength < 1) 
			throw new InvalidParameterException("Sequence length not feasable!");
		m = new double[xLength][yLength][zLength];
		m[0][0][0] = 0.0;
	}
	
	/** Constructor for two sequences
	 * @param xLength Length of sequence 1
	 * @param yLength Length of sequence 2
	 */
	CostMatrix(int xLength, int yLength) {
		this(xLength, yLength, 1);
	}
	
	/**
	 * Get the costs at the point x, y
	 * @param x Position in sequence 1
	 * @param y Position in sequence 2
	 * @return Costs
	 */
	public double get(int x, int y) {
		return get(x, y, 0);
	}
	
	/**
	 * Get the costs at the point x, y, z
	 * @param x Position in sequence 1
	 * @param y Position in sequence 2
	 * @param z Position in sequence 2
	 * @return Costs
	 */
	public double get(int x, int y, int z) {
		if(x < 0 || y < 0 || z < 0 || 
				x >= m.length || y >= m[0].length || z >= m[0][0].length) 
			return Double.NEGATIVE_INFINITY;
		return m[x][y][z];
	}
	
	/** Returns the computed score for the alignment
	 * @return Final Score
	 */
	public double score() {
		return m[m.length-1][m[0].length-1][m[0][0].length-1];
	}
	
	/** Set a Value at x, y
	 * @param x X coordinate
	 * @param y y coordinate
	 * @param value Value that should be inserted
	 */
	public void set(int x, int y, double value) {
		set(x, y, 0, value);
	}
	
	/**
	 * Updates the costs at the point x, y, z
	 * @param x Position in sequence 1
	 * @param y Position in sequence 2
	 * @param z Position in sequence 3
	 * @param value Value that should be inserted
	 */
	public void set(int x, int y, int z, double value) {
		if(x < 0 || y < 0 || z < 0 || 
				x >= m.length || y >= m[0].length || z >= m[0][0].length)
			throw new InvalidParameterException("Index not feasable! x="+x+" y="+y+" z="+z);
		m[x][y][z] = new Double(value);
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		String retVal = "";
		for(int i=0; i<m.length; i++)
			retVal += Arrays.toString(m[i]) + "\n";
		return retVal;
	}
	
}
