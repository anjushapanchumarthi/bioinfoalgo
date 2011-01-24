package cthoelken;

import java.util.LinkedList;

public class Cluster {
	
	Cluster left, right;
	int c;
	String s;
	double distance;
	boolean weighted;
	
	Cluster(Cluster c1, Cluster c2, double distance) {
		this.left = c1; this.right = c2;
		this.distance = distance;
		this.weighted = c1.weighted;
	}
	
	boolean isWeighted() {
		return weighted;
	}
	
	boolean isLeaf() {
		if(left == null) return true;
		return false;
	}
	
	double distance(Cluster cluster, boolean weighted, boolean usePAM, double gapCosts) {
		NeedlemanWunsch nw = new NeedlemanWunsch();
		if(isLeaf() && cluster.isLeaf())
			return -nw.getScore(s, cluster.s, usePAM, gapCosts);
		if(isLeaf()) {
			return (cluster.left.distance(this, weighted, usePAM, gapCosts)+cluster.right.distance(this, weighted, usePAM, gapCosts));
		}
		if(!weighted)
			return (cluster.distance(left, weighted, usePAM, gapCosts)+cluster.distance(right, weighted, usePAM, gapCosts))/2;
		return (cluster.distance(left, weighted, usePAM, gapCosts)+cluster.distance(right, weighted, usePAM, gapCosts))/size();
	}
	
	void align(boolean usePAM, double gapCosts) {
		if(isLeaf()) return;

		left.align(usePAM, gapCosts);
		right.align(usePAM, gapCosts);
		double dist;
		double min = Double.POSITIVE_INFINITY;
		int xMin = 0; int yMin = 0;
		for(int x = 0; x < left.size(); x++) {
			for(int y = 0; y < right.size(); y++) {
				dist = left.getLeaves().get(x).distance(right.getLeaves().get(y), usePAM, usePAM, gapCosts);
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
	
	public Alignment generateAlignment(Alignment algn) {
		if(size() != algn.size()) throw new IllegalArgumentException("Cannot generate alignment if tree and alignment do not have the same size!");
		for(int i = 0; i < left.size(); i++)
			algn.setSeq(left.getLeaves().get(i).c, left.getLeaves().get(i).s);
		for(int i = 0; i < right.size(); i++)
			algn.setSeq(right.getLeaves().get(i).c, right.getLeaves().get(i).s);
		return algn;
	}

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

	double size() {
		if(isLeaf()) return 1;
		return left.size() + right.size();
	}

	Cluster(int content, String sequence) {
		this.c = content;
		this.s = sequence;
	}

	public String toString() {
		if(isLeaf()) return ""+c;
		return "("+left.toString()+","+right.toString()+")";
	}

}
