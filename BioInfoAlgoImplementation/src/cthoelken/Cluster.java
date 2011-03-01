package cthoelken;

import java.util.LinkedList;

/** Auxilliary class for tree based clusters in UPGMA and WPGMA
 * @author Clemens Thoelken
 *
 */
public class Cluster {
	
	Cluster left, right;	// child nodes
	int c;					// internal Integer of the sequence
	String s;				// sequence String
	double distance;		// distance (not really needed)
	boolean weighted;		// true for WPGMA
	
	
	/** Constructor for a node with two children
	 * @param c1 Left child Cluster
	 * @param c2 Right child Cluster
	 * @param distance Distance between the two children
	 */
	Cluster(Cluster c1, Cluster c2, double distance) {
		this.left = c1; this.right = c2;
		this.distance = distance;
		this.weighted = c1.weighted;
	}
	
	/** Getter for WPGMA boolean
	 * @return Returns whether WPGMA is used or not
	 */
	boolean isWeighted() {
		return weighted;
	}
	
	/** Checks whether the node is a leaf or not
	 * @return Returns true if the node has no children
	 */
	boolean isLeaf() {
		if(left == null) return true;
		return false;
	}
	
	/** Computes the distance between this cluster and another
	 * @param cluster The other cluster
	 * @param weighted Use WPGMA or not
	 * @param usePAM Use PAM for substituionMatrix
	 * @param gapCosts Gap costs for the calculation
	 * @return Returns the inverted similarity score of the two clusters
	 */
	double distance(Cluster cluster, boolean weighted, boolean usePAM, double gapCosts) {
		NeedlemanWunsch nw = new NeedlemanWunsch();
		if(isLeaf() && cluster.isLeaf())	// both nodes are sequences
			return -nw.getScore(s, cluster.s, usePAM, gapCosts);
		if(isLeaf()) {		// this is a sequence but the other one isn't -> recursion
			return (cluster.left.distance(this, weighted, usePAM, gapCosts)+cluster.right.distance(this, weighted, usePAM, gapCosts));
		}
		if(!weighted)		// both are clusters -> recursion according to mode
			return (cluster.distance(left, weighted, usePAM, gapCosts)+cluster.distance(right, weighted, usePAM, gapCosts))/2;
		return (cluster.distance(left, weighted, usePAM, gapCosts)+cluster.distance(right, weighted, usePAM, gapCosts))/size();
	}
	
	/** Recursively aligns the nodes of the cluster according to their hierarchy
	 * @param usePAM Use PAM for substitution
	 * @param gapCosts Gap costs for the calculation
	 */
	void align(boolean usePAM, double gapCosts) {
		if(isLeaf()) return;

		left.align(usePAM, gapCosts);		// RECURSION
		right.align(usePAM, gapCosts);
		
		double dist;
		double min = Double.POSITIVE_INFINITY;
		int xMin = 0; int yMin = 0;
		for(int x = 0; x < left.size(); x++) {
			for(int y = 0; y < right.size(); y++) {
				dist = left.getLeaves().get(x).distance(right.getLeaves().get(y), weighted, usePAM, gapCosts);
				if(min > dist) {
					min = dist; xMin = x; yMin = y;
				}
			}
		}
		Alignment algn = new Alignment(new NeedlemanWunsch().getAlignment(left.getLeaves().get(xMin).s, right.getLeaves().get(yMin).s, usePAM, gapCosts));;
		for(int i = 0; i < algn.getSeq(0).length(); i++)
			for(int j = 0; j < left.size(); j++)
				if(algn.getSeq(0).charAt(i) == '_')
					if(i == 0) left.getLeaves().get(j).s = "-"+left.getLeaves().get(j).s;
					else if(left.getLeaves().get(j).s.length() > i)
						left.getLeaves().get(j).s = left.getLeaves().get(j).s.substring(0, i) + "-" + left.getLeaves().get(j).s.substring(i, left.getLeaves().get(j).s.length());
					else left.getLeaves().get(j).s = left.getLeaves().get(j).s+"-";
		for(int i = 0; i < algn.getSeq(1).length(); i++)
			for(int j = 0; j < right.size(); j++)
				if(algn.getSeq(1).charAt(i) == '_')
					if(i == 0) right.getLeaves().get(j).s = "-"+right.getLeaves().get(j).s;
					else if(right.getLeaves().get(j).s.length() > i)
						right.getLeaves().get(j).s = right.getLeaves().get(j).s.substring(0, i) + "-" + right.getLeaves().get(j).s.substring(i, right.getLeaves().get(j).s.length());
					else right.getLeaves().get(j).s = right.getLeaves().get(j).s+"-";
	}
	
	/** Generates an Alignment according to the trees ranking
	 * @param algn Input alignment without any gaps
	 * @return Output alignment with the sequence strings from the clustered tree
	 */
	public Alignment generateAlignment(Alignment algn) {
		if(size() != algn.size()) throw new IllegalArgumentException("Cannot generate alignment if tree and alignment do not have the same size!");
		for(int i = 0; i < left.size(); i++)
			algn.setSeq(left.getLeaves().get(i).c, left.getLeaves().get(i).s);
		for(int i = 0; i < right.size(); i++)
			algn.setSeq(right.getLeaves().get(i).c, right.getLeaves().get(i).s);
		return algn;
	}

	/** Returns all leaf nodes from the tree
	 * @return A List of the nodes
	 */
	private LinkedList<Cluster> getLeaves() {
		LinkedList<Cluster> leaves = new LinkedList<Cluster>();
		if(isLeaf()) {
			leaves.add(this);
			return leaves;
		}
		
		for(int i = 0; i < left.size(); i++)
			leaves.add(left.getLeaves().get(i));
		for(int i = 0; i < right.size(); i++)
			leaves.add(right.getLeaves().get(i));
		return leaves;
	}

	/** Returns the numer of trees in this tree
	 * @return Number of leaf nodes
	 */
	int size() {
		if(isLeaf()) return 1;
		return left.size() + right.size();
	}

	/** Constructor for a leaf-node without children
	 * @param content Internal Integer code for the Sequence
	 * @param sequence The sequence's String
	 */
	Cluster(int content, String sequence) {
		this.c = content;
		this.s = sequence;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		if(isLeaf()) return ""+c;
		return "("+left.toString()+","+right.toString()+")";
	}

}
