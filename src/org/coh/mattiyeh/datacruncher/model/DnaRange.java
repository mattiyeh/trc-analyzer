package org.coh.mattiyeh.datacruncher.model;

public class DnaRange {

	String chr;
	int start;
	int end;

	public DnaRange(String chr, int start, int end) {
		super();
		this.chr = chr;
		this.start = start;
		this.end = end;
	}

	public String getChr() {
		return chr;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}

	public boolean overlaps(int mutStart, int mutEnd) {
		return isInRange(mutStart) || isInRange(mutEnd);
	}

	private boolean isInRange(int snpPosition) {
		return (start <= snpPosition) && (snpPosition <= end);
	}

	@Override
	public String toString() {
		return "DnaRange [chr=" + chr + ", start=" + start + ", end=" + end + "]";
	}

}
